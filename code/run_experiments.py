"""Forking-path inflation audit for molecular property benchmarks.
New code written during the Apart AI x Epistemics sprint (2026-11-14).
CPU only. Data: public MoleculeNet CSVs (ESOL, FreeSolv, Lipophilicity) from the DeepChem S3 bucket.
"""
import json, sys, time, warnings
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
RDLogger.DisableLog('rdApp.*'); warnings.filterwarnings('ignore')
from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[1])
N_SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 20
TEST_FRAC = 0.2
DATASETS = {
    'ESOL': ('delaney-processed.csv', 'smiles', 'measured log solubility in mols per litre'),
    'FreeSolv': ('SAMPL.csv', 'smiles', 'expt'),
    'Lipophilicity': ('Lipophilicity.csv', 'smiles', 'exp'),
}

def load(name):
    f, s, y = DATASETS[name]
    df = pd.read_csv(f'{ROOT}/data/{f}')
    mols, ys, smi = [], [], []
    for a, b in zip(df[s], df[y]):
        m = Chem.MolFromSmiles(str(a).strip())
        if m is not None and np.isfinite(b):
            mols.append(m); ys.append(float(b)); smi.append(Chem.MolToSmiles(m))
    return mols, np.array(ys), smi

def featurize(mols):
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in mols]
    X_fp = np.zeros((len(mols), 2048), dtype=np.uint8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, X_fp[i])
    names = [n for n, _ in Descriptors.descList]
    X_d = np.array([[fn(m) for _, fn in Descriptors.descList] for m in mols], dtype=np.float64)
    X_d = np.nan_to_num(np.clip(X_d, -1e6, 1e6))
    return fps, X_fp, X_d

def random_split(n, seed):
    rng = np.random.RandomState(seed); idx = rng.permutation(n)
    k = int(round(TEST_FRAC * n)); return idx[k:], idx[:k]

def group_split(groups, n, seed):
    """Assign whole groups to test (random group order) until test >= TEST_FRAC*n. Groups larger
    than half the target test size go to train (as in common 'balanced scaffold' splitters)."""
    rng = np.random.RandomState(seed)
    target = int(round(TEST_FRAC * n))
    gids = list(groups.keys()); rng.shuffle(gids)
    test = []
    for g in gids:
        mem = groups[g]
        if len(mem) > target / 2: continue
        if len(test) + len(mem) <= target * 1.05:
            test.extend(mem)
        if len(test) >= target: break
    test = np.array(sorted(test)); mask = np.ones(n, bool); mask[test] = False
    return np.where(mask)[0], test

def scaffold_groups(mols):
    g = {}
    for i, m in enumerate(mols):
        try: s = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
        except Exception: s = ''
        g.setdefault(s if s else f'acyclic_{i}', []).append(i)  # acyclic molecules: own group
    return g

def butina_groups(fps, cutoff=0.6):
    n = len(fps); d = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]); d.extend([1 - x for x in sims])
    cl = Butina.ClusterData(d, n, cutoff, isDistData=True)
    return {k: list(c) for k, c in enumerate(cl)}

def nn_sim(fps, tr, te):
    trf = [fps[i] for i in tr]
    return np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], trf)) for i in te])

def main():
    rows, permol = [], []
    for name in DATASETS:
        t0 = time.time(); mols, y, smi = load(name); n = len(y)
        fps, X_fp, X_d = featurize(mols)
        groups = {'scaffold': scaffold_groups(mols), 'cluster': butina_groups(fps)}
        n_dup = n - len(set(smi))
        meta = dict(dataset=name, n=n, n_scaffolds=len(groups['scaffold']), n_clusters=len(groups['cluster']), dup_canonical_smiles=n_dup)
        print(meta, flush=True)
        json.dump(meta, open(f'{ROOT}/results/meta_{name}.json', 'w'))
        for split in ['random', 'scaffold', 'cluster']:
            for seed in range(N_SEEDS):
                tr, te = random_split(n, seed) if split == 'random' else group_split(groups[split], n, seed)
                sims = nn_sim(fps, tr, te)
                for feat, X in [('Morgan2048', X_fp), ('RDKitDesc', X_d)]:
                    rf = RandomForestRegressor(n_estimators=200, n_jobs=8, random_state=seed, min_samples_leaf=1)
                    rf.fit(X[tr], y[tr]); p = rf.predict(X[te])
                    r2 = r2_score(y[te], p); rmse = mean_squared_error(y[te], p) ** 0.5
                    rows.append(dict(dataset=name, split=split, seed=seed, feat=feat, n_train=len(tr), n_test=len(te),
                                     r2=r2, rmse=rmse, nn_sim_median=float(np.median(sims)),
                                     frac_test_nn_ge_0p8=float(np.mean(sims >= 0.8)), test_y_sd=float(y[te].std())))
                    for i, s_, yy, pp in zip(te, sims, y[te], p):
                        permol.append((name, split, seed, feat, int(i), float(s_), float(yy), float(pp)))
            print(name, split, 'done', round(time.time() - t0, 1), 's', flush=True)
    pd.DataFrame(rows).to_csv(f'{ROOT}/results/runs.csv', index=False)
    pd.DataFrame(permol, columns=['dataset', 'split', 'seed', 'feat', 'idx', 'nn_sim', 'y', 'pred']).to_csv(f'{ROOT}/results/per_molecule.csv.gz', index=False)

if __name__ == '__main__':
    main()
