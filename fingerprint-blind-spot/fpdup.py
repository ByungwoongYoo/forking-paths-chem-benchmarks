import pandas as pd, numpy as np, json
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
res={}
for name in ['ESOL','FreeSolv','Lipophilicity','BACE']:
    d=pd.read_csv(f'parsed_{name}.csv')
    fps=[gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)).astype(np.uint8) for s in d.can]
    X=np.array(fps); y=d.y.values
    keys=[x.tobytes() for x in X]
    d['fpkey']=[hash(k) for k in keys]
    g=d.groupby('fpkey'); sz=g.size(); dupk=sz[sz>1].index
    spread=g.y.agg(lambda v:v.max()-v.min())[sz>1]
    r={'n':len(d),'fp_identical_groups':int(len(dupk)),'rows_in_fp_groups':int(sz[sz>1].sum()),
       'distinct_structures_in_fp_groups':int(d[d.fpkey.isin(dupk)].groupby('fpkey').can.nunique().gt(1).sum()),
       'fp_group_label_spread_median':float(spread.median()) if len(spread) else None,
       'fp_group_label_spread_max':float(spread.max()) if len(spread) else None,
       'fp_groups_spread_ge1':int((spread>=1).sum())}
    seeds=[];
    for seed in range(10):
        rng=np.random.RandomState(seed); p=rng.permutation(len(d)); ntr=int(.8*len(d)); nva=int(.1*len(d))
        tr=p[:ntr]; te=p[ntr+nva:]
        trkeys=set(d.fpkey.values[tr])
        tw=np.array([d.fpkey.values[i] in trkeys for i in te])
        rf=RandomForestRegressor(n_estimators=200,n_jobs=-1,random_state=seed).fit(X[tr],y[tr])
        pr=rf.predict(X[te])
        e=np.abs(pr-y[te])
        seeds.append(dict(n_test=len(te),n_twin=int(tw.sum()),r2_all=r2_score(y[te],pr),
           r2_no_twin=r2_score(y[te][~tw],pr[~tw]),mae_twin=float(e[tw].mean()) if tw.any() else np.nan,mae_rest=float(e[~tw].mean())))
    s=pd.DataFrame(seeds); r['split_80_10_10_10seeds']=s.mean().to_dict(); r['n_twin_range']=[int(s.n_twin.min()),int(s.n_twin.max())]
    res[name]=r; print(name,json.dumps(r))
json.dump(res,open('fp_twin_stats.json','w'),indent=1)
