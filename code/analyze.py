"""Analysis + figures for the forking-path inflation audit. New sprint code (2026-11-14)."""
import json, itertools
import numpy as np, pandas as pd
from math import comb
from scipy.stats import spearmanr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[1])
R = pd.read_csv(f'{ROOT}/results/runs.csv')
P = pd.read_csv(f'{ROOT}/results/per_molecule.csv.gz')
DS = ['ESOL', 'FreeSolv', 'Lipophilicity']; SP = ['random', 'scaffold', 'cluster']; FE = ['Morgan2048', 'RDKitDesc']
out = {}

# A. summary table
g = R.groupby(['dataset', 'feat', 'split'])
T = g.agg(r2_mean=('r2', 'mean'), r2_sd=('r2', 'std'), r2_min=('r2', 'min'), r2_max=('r2', 'max'),
          rmse_mean=('rmse', 'mean'), rmse_sd=('rmse', 'std'), nn_sim_median=('nn_sim_median', 'mean'),
          frac_nn_ge_0p8=('frac_test_nn_ge_0p8', 'mean'), n_seeds=('seed', 'nunique'), n_test=('n_test', 'mean')).reset_index()
T.to_csv(f'{ROOT}/results/table1_summary.csv', index=False)

# B. split-choice inflation
inf = []
for d in DS:
    for f in FE:
        m = {s: R[(R.dataset == d) & (R.feat == f) & (R.split == s)].r2.mean() for s in SP}
        inf.append(dict(dataset=d, feat=f, r2_random=m['random'], r2_scaffold=m['scaffold'], r2_cluster=m['cluster'],
                        gap_random_minus_scaffold=m['random'] - m['scaffold'], gap_random_minus_cluster=m['random'] - m['cluster']))
inf = pd.DataFrame(inf); inf.to_csv(f'{ROOT}/results/table2_split_gap.csv', index=False)

# C. best-of-K (exact expected maximum of K draws without replacement from the n observed seeds)
def exp_max(vals, K):
    v = np.sort(vals); n = len(v)
    w = np.array([comb(i, K - 1) for i in range(n)], float) / comb(n, K)  # P(max is the i-th smallest)
    return float((w * v).sum())
bok = []
for d in DS:
    for f in FE:
        for s in SP:
            v = R[(R.dataset == d) & (R.feat == f) & (R.split == s)].r2.values
            for K in [1, 2, 5, 10, 20]:
                if K <= len(v): bok.append(dict(dataset=d, feat=f, split=s, K=K, exp_max_r2=exp_max(v, K)))
bok = pd.DataFrame(bok); bok.to_csv(f'{ROOT}/results/table3_best_of_k.csv', index=False)

# C2. full garden: best single number over (split x feat x seed) vs pre-registered honest report
garden = []
for d in DS:
    sub = R[R.dataset == d]
    best = sub.loc[sub.r2.idxmax()]
    honest = sub[sub.split == 'cluster'].groupby('feat').r2.mean()
    garden.append(dict(dataset=d, best_r2=best.r2, best_config=f"{best.split}/{best.feat}/seed{best.seed}",
                       prereg_cluster_mean_r2_bestfeat=honest.max(), prereg_feat=honest.idxmax(),
                       prereg_scaffold_mean_r2_bestfeat=sub[sub.split == 'scaffold'].groupby('feat').r2.mean().max(),
                       worst_r2=sub.r2.min(), n_runs=len(sub)))
garden = pd.DataFrame(garden); garden['inflation_vs_cluster'] = garden.best_r2 - garden.prereg_cluster_mean_r2_bestfeat
garden.to_csv(f'{ROOT}/results/table4_garden.csv', index=False)

# D. similarity vs performance across split instances, and per-molecule error by NN-sim bin
corr = []
for d in DS:
    for f in FE:
        s = R[(R.dataset == d) & (R.feat == f)]
        rho, p = spearmanr(s.nn_sim_median, s.r2); corr.append(dict(dataset=d, feat=f, n_split_instances=len(s), spearman_rho=rho, p=p))
corr = pd.DataFrame(corr); corr.to_csv(f'{ROOT}/results/table5_sim_vs_r2.csv', index=False)
BINS = [0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001]
P['bin'] = pd.cut(P.nn_sim, BINS, right=False, labels=False)
P['se'] = (P.y - P.pred) ** 2; P['ae'] = (P.y - P.pred).abs()
eb = P.groupby(['dataset', 'feat', 'bin']).agg(mae=('ae', 'mean'), n=('ae', 'size')).reset_index()
eb.to_csv(f'{ROOT}/results/table6_error_by_sim_bin.csv', index=False)

# E. reweighting audit: estimate scaffold/cluster RMSE from random-split errors reweighted to the target
# split's NN-similarity histogram (computable without training). Pooled over seeds.
rw = []
for d in DS:
    for f in FE:
        src = P[(P.dataset == d) & (P.feat == f) & (P.split == 'random')]
        bin_mse = src.groupby('bin').se.mean()
        for tgt in ['scaffold', 'cluster']:
            for seed in sorted(P.seed.unique()):
                t = P[(P.dataset == d) & (P.feat == f) & (P.split == tgt) & (P.seed == seed)]
                h = t.bin.value_counts(normalize=True)
                est = 0.0
                for b, w in h.items():
                    if b in bin_mse.index: est += w * bin_mse[b]
                    else:  # nearest populated source bin
                        nb = min(bin_mse.index, key=lambda x: abs(x - b)); est += w * bin_mse[nb]
                rw.append(dict(dataset=d, feat=f, target=tgt, seed=seed, rmse_actual=float(np.sqrt(t.se.mean())),
                               rmse_reweighted_estimate=float(np.sqrt(est)), rmse_random_naive=float(np.sqrt(src[src.seed == seed].se.mean()))))
rw = pd.DataFrame(rw); rw.to_csv(f'{ROOT}/results/table7_reweight_audit_per_seed.csv', index=False)
rws = rw.groupby(['dataset', 'feat', 'target']).mean(numeric_only=True).drop(columns='seed').reset_index()
rws['naive_gap'] = rws.rmse_random_naive - rws.rmse_actual
rws['reweighted_gap'] = rws.rmse_reweighted_estimate - rws.rmse_actual
rws['gap_closed_frac'] = 1 - rws.reweighted_gap.abs() / rws.naive_gap.abs()
rws.to_csv(f'{ROOT}/results/table7_reweight_audit_summary.csv', index=False)

# Figures
col = {'random': '#d62728', 'scaffold': '#1f77b4', 'cluster': '#2ca02c'}
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for a, d in zip(ax, DS):
    pos = 0; ticks = []; labs = []
    for f in FE:
        for s in SP:
            v = R[(R.dataset == d) & (R.feat == f) & (R.split == s)].r2.values
            bp = a.boxplot(v, positions=[pos], widths=0.6, patch_artist=True)
            bp['boxes'][0].set_facecolor(col[s]); ticks.append(pos); labs.append(s[:4])
            pos += 1
        pos += 0.7
    a.set_xticks(ticks); a.set_xticklabels(labs, fontsize=8); a.set_title(d); a.set_ylabel('test R²')
    a.text(1, a.get_ylim()[0], 'Morgan', ha='center', va='bottom', fontsize=8); a.text(4.7, a.get_ylim()[0], 'RDKit desc', ha='center', va='bottom', fontsize=8)
fig.suptitle('Figure 1. Same data, same model: test R² under random / scaffold / cluster splits (20 seeds each)', fontsize=10)
fig.tight_layout(); fig.savefig(f'{ROOT}/figures/fig1_r2_by_split.png', dpi=200); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for a, d in zip(ax, DS):
    for f, ls in zip(FE, ['-', '--']):
        for s in SP:
            b = bok[(bok.dataset == d) & (bok.feat == f) & (bok.split == s)]
            a.plot(b.K, b.exp_max_r2, ls, marker='o', ms=3, color=col[s], label=f'{s}/{f}')
    a.set_xscale('log'); a.set_xticks([1, 2, 5, 10, 20]); a.set_xticklabels([1, 2, 5, 10, 20]); a.set_xlabel('K seeds tried, best reported'); a.set_title(d)
ax[0].set_ylabel('expected reported R²'); ax[2].legend(fontsize=6, loc='lower right')
fig.suptitle('Figure 2. Best-of-K seed selection: expected reported R² when the best of K split seeds is reported', fontsize=10)
fig.tight_layout(); fig.savefig(f'{ROOT}/figures/fig2_best_of_k.png', dpi=200); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
mids = [(BINS[i] + min(BINS[i + 1], 1)) / 2 for i in range(len(BINS) - 1)]
for a, d in zip(ax, DS):
    for f, ls in zip(FE, ['-', '--']):
        for s in SP:
            e = P[(P.dataset == d) & (P.feat == f) & (P.split == s)].groupby('bin').ae.mean()
            a.plot([mids[int(i)] for i in e.index], e.values, ls, marker='o', ms=3, color=col[s], label=f'{s}/{f}')
    a.set_xlabel('test molecule: max Tanimoto to train'); a.set_title(d)
ax[0].set_ylabel('mean absolute error'); ax[2].legend(fontsize=6)
fig.suptitle('Figure 3. Error depends on nearest-train similarity, roughly the same curve whatever the split', fontsize=10)
fig.tight_layout(); fig.savefig(f'{ROOT}/figures/fig3_error_vs_similarity.png', dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.6, 4.2))
for _, r in rws.iterrows():
    ax.scatter(r.rmse_actual, r.rmse_random_naive, marker='x', color='grey')
    ax.scatter(r.rmse_actual, r.rmse_reweighted_estimate, marker='o', color=col[r.target])
lim = [0, max(rws.rmse_actual.max(), rws.rmse_reweighted_estimate.max()) * 1.1]; ax.plot(lim, lim, 'k:', lw=1)
ax.set_xlabel('actual RMSE (scaffold / cluster split)'); ax.set_ylabel('estimate from random-split errors')
ax.set_title('Figure 4. x = naive random-split RMSE;\no = similarity-reweighted estimate', fontsize=9)
fig.tight_layout(); fig.savefig(f'{ROOT}/figures/fig4_reweight_audit.png', dpi=200); plt.close(fig)

pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
for name, df in [('T1', T), ('T2', inf), ('T4', garden), ('T5', corr), ('T7', rws)]:
    print('==', name); print(df.round(3).to_string(index=False))
print('== T3 K=1 vs K=20'); print(bok[bok.K.isin([1, 20])].pivot_table(index=['dataset', 'feat', 'split'], columns='K', values='exp_max_r2').round(3))
