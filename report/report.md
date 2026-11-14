# Forking Paths in Molecular Property Benchmarks: How Far a Reported Score Moves Under Defensible Choices, and a Training-Free Audit Signal

Byungwoong Yoo, BigWonder (independent), Seoul, South Korea

**Abstract.** A benchmark score is a claim, and readers, including AI agents that summarise or extend research, usually see the number without the choices behind it. We measure how far one reported test R² can move under choices that are each individually defensible. We used three public MoleculeNet regression sets (ESOL, FreeSolv, Lipophilicity), two standard feature sets and a random forest. Across 360 train/test runs (3 split strategies × 20 split seeds × 2 feature sets × 3 datasets), random splits gave a higher mean R² than scaffold splits, by 0.050 to 0.133, and than Butina-cluster splits, by 0.035 to 0.250. Group splits, usually read as the rigorous option, were 1.2 to 3.3 times more seed-sensitive. Reporting the best of 20 split seeds raised their expected R² by 0.047 to 0.193. In 6 of 12 settings this was enough to match or beat the honest random-split mean. Reporting the best of all 120 runs per dataset overstated the cluster-split mean by 0.059 to 0.205 R². The median nearest-training-neighbour Tanimoto similarity of a split can be computed without training a model, and it tracked R² across splits (Spearman ρ 0.37 to 0.89). Reweighting random-split errors to a target split's similarity profile cut the mean absolute error of predicting group-split RMSE from 0.186 to 0.122, improving 11 of 12 settings, but it overshot on one small dataset. We propose a two-number disclosure, the split similarity profile plus the seed spread, that makes these forks visible to human and automated readers.

## 1. Introduction

A reported benchmark number compresses many choices into one figure: which split, which random seed, which features, which run to report. Each choice can be defended on its own. Taken together they form a "garden of forking paths" (Gelman & Loken, 2013): a researcher, or an automated research agent optimising a score, can move the headline number a long way without fabricating anything. Readers rarely see the forks. That includes the AI systems that increasingly read, summarise and build on published results. This is a direct instance of the sprint's Track 1 question: which exploits make benchmarks gameable, and can one be demonstrated and closed?

We study the question where it can be measured cheaply and exactly: small-molecule property regression. Here split choice is a known, heavily debated degree of freedom (random vs. scaffold vs. cluster splits). Many papers report a single number, and AI-assisted drug-discovery claims are often taken at face value. We ask three questions. (Q1) How far does one reported R² move across defensible splits and seeds? (Q2) Does the "rigorous" split protect readers if the seed is shopped? (Q3) Is there a training-free signal a reader or reviewer agent could compute to tell how solid a reported number is?

Our main contributions are:

1. **A measured "reportable range"** for three public benchmarks: 360 runs showing that the best single run overstates the pre-registered cluster-split mean by 0.059 to 0.205 R².

2. **Evidence that the harder split alone is not a safeguard.** Scaffold and cluster splits are 1.2 to 3.3 times more seed-sensitive than random splits, so best-of-20 seed reporting on a "rigorous" split matched or beat the honest random-split mean in 6 of 12 settings.

3. **A training-free audit signal and a proposed disclosure.** The median nearest-training-neighbour similarity of a split tracks its R² (Spearman ρ 0.37 to 0.89). Reweighting random-split errors by that profile predicted group-split RMSE better than the naive number in 11 of 12 settings. We propose reporting the similarity profile and the seed spread with any score.

## 2. Related Work

**Split choice in molecular ML.** MoleculeNet (Wu et al., 2018) popularised scaffold splits based on Bemis–Murcko frameworks (Bemis & Murcko, 1996). Sheridan (2013) argued that random splits overstate prospective performance and proposed time splits. Sheridan et al. (2004) showed that similarity to the training set predicts per-molecule accuracy. Guo et al. (2024) argue that even scaffold splits overestimate virtual-screening performance compared with cluster-based splits. Our split-gap and error-vs-similarity results (Figures 1 and 3) replicate this known pattern and are not new. They serve as the measured base for Q2 and Q3.

**Researcher degrees of freedom and leakage.** Simmons et al. (2011) and Gelman & Loken (2013) describe how undisclosed analytic flexibility inflates findings. Kapoor & Narayanan (2023) catalogue data leakage across ML-based science. Dodge et al. (2019) proposed reporting expected maximum validation performance as a function of tuning budget, and our best-of-K estimator is the same order-statistic idea applied to split seeds rather than hyperparameters. Bouthillier et al. (2021) show that data-split variance is a major, often ignored, source of benchmark variance.

**What is different here.** We (i) quantify the combined reportable range, rather than one fork at a time, on public data with fully released per-molecule predictions; (ii) show that the usual remedy, a harder split, increases the room for seed shopping; and (iii) test whether a cheap, training-free similarity profile lets a reader or reviewer agent recover most of the gap. That last property matters for AI-for-epistemics: an agent that summarises papers could compute it from released data without retraining anything.

## 3. Methods

**Data.** Three public MoleculeNet regression CSVs from the DeepChem S3 bucket: ESOL aqueous solubility (Delaney, 2004; n = 1,128), FreeSolv hydration free energy (Mobley & Guthrie, 2014; n = 642) and Lipophilicity logD (ChEMBL-derived; n = 4,200). SMILES were parsed with RDKit, and all rows parsed. ESOL contains 11 duplicate canonical SMILES; we kept them as distributed, as a naive user would. SHA-256 hashes of the files are in the README.

**Splits (80/20 train/test, no validation set, no tuning).** Random: a seeded permutation. Scaffold: Bemis–Murcko scaffold groups. Cluster: Butina clustering on Morgan fingerprints at Tanimoto distance 0.6. For group splits, groups are shuffled with the seed and added to the test set until it reaches 20% (5% tolerance). Groups larger than half the target test size always go to training. Acyclic molecules (empty scaffold) each form their own group, a choice that differs from DeepChem's single-empty-scaffold convention (see Limitations). Seeds 0 to 19 for every split type.

**Models and features.** scikit-learn RandomForestRegressor (200 trees, default settings, random_state = split seed) on either a 2048-bit Morgan fingerprint (radius 2) or all RDKit 2D descriptors (Descriptors.descList, NaN set to 0, clipped to ±1e6). There was no hyperparameter search, so the only forks are split type, split seed and feature set.

**Metrics.** Test R² and RMSE. For every test molecule we store its maximum Tanimoto similarity (Morgan, radius 2) to the training set ("NN-sim"). This needs no model. The split-level summary is the median NN-sim.

**Best-of-K.** For each (dataset, features, split) we have 20 R² values. The expected maximum when K seeds are tried and the best is reported is computed exactly from the order statistics of the 20 observed values, by sampling without replacement (Dodge et al., 2019). We report K = 1, 2, 5, 10, 20.

**Reweighting audit.** We take per-molecule squared errors from the random splits (pooled over 20 seeds) and average them within NN-sim bins (edges 0, 0.3, 0.4, …, 0.9, 1.0). For each target scaffold or cluster split we reweight these bin means by the target split's NN-sim histogram, which needs no training on the target. If a target bin is empty in the source, we use the nearest populated bin. We compare the estimated RMSE with the RMSE actually obtained on that split, averaged over seeds.

**Compute.** One shared 8-core CPU Linux machine, no GPU, no paid services. The full experiment ran in about 38 minutes of wall time.

## 4. Results

**Q1: Split choice.** Table 1 and Figure 1 show the 20-seed distributions. With the same data and model, mean R² fell from random to scaffold splits by 0.050 to 0.133, and from random to cluster splits by 0.035 to 0.250. Paired by seed, the group-split RMSE exceeded the random-split RMSE in 65% to 100% of seeds, depending on the setting. Test-label standard deviations were similar across split types within each dataset (for example Lipophilicity 1.196 / 1.216 / 1.198), so the R² drop is not an artefact of a narrower test set. The gap was largest with fingerprints and on Lipophilicity, where the median NN-sim fell from 0.70 under random splits to 0.37 under cluster splits.

| Dataset       | Features   | R² random     | R² scaffold   | R² cluster    | NN-sim r/s/c   |
|:--------------|:-----------|:--------------|:--------------|:--------------|:---------------|
| ESOL          | Morgan2048 | 0.692 ± 0.034 | 0.559 ± 0.085 | 0.487 ± 0.112 | 0.55/0.47/0.42 |
| ESOL          | RDKitDesc  | 0.908 ± 0.018 | 0.859 ± 0.028 | 0.873 ± 0.026 | 0.55/0.47/0.42 |
| FreeSolv      | Morgan2048 | 0.667 ± 0.081 | 0.551 ± 0.097 | 0.574 ± 0.105 | 0.53/0.49/0.42 |
| FreeSolv      | RDKitDesc  | 0.893 ± 0.049 | 0.814 ± 0.084 | 0.848 ± 0.068 | 0.53/0.49/0.42 |
| Lipophilicity | Morgan2048 | 0.525 ± 0.021 | 0.415 ± 0.028 | 0.274 ± 0.043 | 0.70/0.57/0.37 |
| Lipophilicity | RDKitDesc  | 0.664 ± 0.020 | 0.603 ± 0.032 | 0.499 ± 0.053 | 0.70/0.57/0.37 |

Table 1. Test R² (mean ± sd over 20 split seeds) by split strategy, with the mean of the per-split median nearest-train Tanimoto similarity (random / scaffold / cluster). A larger sd means more room to choose a favourable seed.

![](../figures/fig1_r2_by_split.png)
Figure 1. Test R² for the same data and model under random (red), scaffold (blue) and Butina-cluster (green) splits, 20 seeds each. Left group: Morgan fingerprints. Right group: RDKit descriptors.

**Q2: Seed shopping on the "rigorous" split.** Group splits were 1.2 to 3.3 times more variable across seeds than random splits (sd ratio; for example ESOL Morgan cluster sd 0.112 vs. random 0.034). Figure 2 shows the consequence. Reporting the best of K seeds raises the expected reported R² quickly for group splits: from K = 1 to K = 20 by 0.047 to 0.193, against 0.024 to 0.130 for random splits. The best-of-K group-split score matched or exceeded the honest random-split mean in 0 of 12 settings at K = 2, 3 of 12 at K = 5, 5 of 12 at K = 10 and 6 of 12 at K = 20. A reviewer who checks only that a paper "used a scaffold split" therefore cannot tell an honest scaffold number from a shopped one.

![](../figures/fig2_best_of_k.png)
Figure 2. Expected reported R² when the best of K split seeds is reported (exact order-statistic estimate from 20 observed seeds). Solid lines: Morgan fingerprints. Dashed lines: RDKit descriptors.

**Combined garden.** Table 2 takes the single best of all 120 runs per dataset (3 splits × 20 seeds × 2 feature sets). In every dataset the winner was a random split with RDKit descriptors. It overstated the pre-registered cluster-split mean (best feature set) by 0.059 (ESOL), 0.103 (FreeSolv) and 0.205 (Lipophilicity) R². The worst single runs (0.235, 0.280, 0.199) show the downward range is just as wide.

| Dataset       |   Best single run (of 120) | Config of best run      |   Cluster-split mean, best features |   Scaffold-split mean, best features |   Worst single run |   Inflation vs cluster mean |
|:--------------|---------------------------:|:------------------------|------------------------------------:|-------------------------------------:|-------------------:|----------------------------:|
| ESOL          |                      0.932 | random/RDKitDesc/seed7  |                               0.873 |                                0.859 |              0.235 |                       0.059 |
| FreeSolv      |                      0.95  | random/RDKitDesc/seed18 |                               0.848 |                                0.814 |              0.28  |                       0.103 |
| Lipophilicity |                      0.704 | random/RDKitDesc/seed11 |                               0.499 |                                0.603 |              0.199 |                       0.205 |

Table 2. Reportable range per dataset over 120 runs. "Inflation" = best single run minus the 20-seed cluster-split mean of the better feature set.

**Q3: A training-free signal.** Across the 60 split instances per (dataset, features), median NN-sim correlated with test R²: Spearman ρ = 0.645 and 0.557 (ESOL, Morgan / descriptors), 0.369 and 0.388 (FreeSolv) and 0.894 and 0.865 (Lipophilicity), all p < 0.005. Per molecule, mean absolute error fell with NN-sim in all six settings (Figure 3; for example Lipophilicity Morgan MAE 0.964 below 0.3 vs. 0.341 at ≥ 0.9). One exception: on ESOL with Morgan fingerprints, the ≥ 0.9 bin had a higher MAE (1.020) than the 0.8 to 0.9 bin (0.670). We did not check the cause. Duplicate or near-identical fingerprints with differing labels are a plausible explanation.

![](../figures/fig3_error_vs_similarity.png)
Figure 3. Mean absolute error per test molecule as a function of its maximum Tanimoto similarity to the training set, pooled over seeds. Curves from different split strategies largely overlap, which is what makes reweighting possible.

Because the error-vs-similarity curves overlap across split types, a reader with only random-split errors and the target split's NN-sim histogram can estimate the harder split's RMSE (Figure 4; Appendix Table 3). The mean absolute error of the estimate fell from 0.186 (naive random-split RMSE) to 0.122, and the estimate was closer in 11 of 12 settings (median 55% of the gap closed; above 75% in all four Lipophilicity settings). It failed on FreeSolv with Morgan fingerprints and a cluster target: the estimate overshot by 0.471 RMSE, worse than the naive gap of 0.227. FreeSolv is small (about 129 test molecules), and its low-similarity bin is sparse in random splits, so the bin means there are noisy.

![](../figures/fig4_reweight_audit.png)
Figure 4. Predicting group-split RMSE from random-split errors. Grey ×: the naive random-split RMSE. Coloured ○: the similarity-reweighted estimate (blue = scaffold target, green = cluster target). The dotted line is perfect prediction.

## 5. Discussion and Limitations

**Implications.** For AI-for-epistemics, the point is not that random splits are optimistic, which is known. The point is that (a) the honest-looking choice, a harder split, enlarges the room for undisclosed seed selection, and (b) a single number does not carry enough information for a reader, human or agent, to judge how solid it is. Automated research agents that iterate on a score will traverse these forks by default. A summariser agent that repeats the headline number passes the inflation on. We propose a minimal disclosure that both people and agents can consume: (1) the per-split NN-sim profile (median and fraction ≥ 0.8), which costs no training and can be recomputed from released data, and (2) the score's spread over at least 5 split seeds, with the seed list fixed and hashed in advance. In our data the seed spread directly exposes the room for best-of-K reporting (Figure 2), and the NN-sim profile tracks the split effect (ρ 0.37 to 0.89), though we did not test how readers or agents use these signals. The reweighting estimate is a useful but imperfect sanity check, not a correction.

### Limitations

This is one model family (random forest), two feature sets and three small public datasets, all physicochemical regression. We did not test graph neural networks, pretrained models, classification tasks or time splits, and the size of the gaps may differ for them. Group-split construction involves its own choices: the acyclic-molecule handling, the Butina cutoff of 0.6, the 5% size tolerance and excluding large groups from the test set. These choices shift the numbers, and we tested only one setting of each. Our best-of-K figures assume the seeds tried are exchangeable draws, and real seed shopping may also involve features, hyperparameters and metrics, so the true reportable range is probably larger than measured. The reweighting audit uses only one similarity measure and fixed bins, and it failed in one of 12 settings. We did not evaluate it on independent datasets, so its reliability is not established. We did not run any language model: the claims about AI agents are a motivation, not a measured result.

### Future Work

Next steps: (1) run an LLM research agent with a fixed iteration budget on the same task and log which forks it takes and what it reports; (2) test whether giving a summariser LLM the two-number disclosure changes how it describes the result; (3) extend to GNN and pretrained baselines and to the Polaris and TDC benchmark suites; (4) package the NN-sim profile as a small, dependency-light checker that runs on any released split file.

## 6. Conclusion

On three public molecular property benchmarks, defensible choices of split, seed and features moved the best reportable R² 0.06 to 0.21 above a pre-registered cluster-split mean. The "rigorous" split was the most exposed to seed shopping. A split's nearest-neighbour similarity profile, which needs no training, tracked these differences and helped a reader predict the harder-split error. We suggest that benchmark claims carry this profile and a seed spread, so that people and AI systems reading them can judge how solid they are.

## Code and Data

Code repository: to be published before submission (all code, results and figures are in the project folder; see README). Data: MoleculeNet CSVs from https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/ (delaney-processed.csv, SAMPL.csv, Lipophilicity.csv). All per-molecule predictions (per_molecule.csv.gz) are released so every table can be recomputed without retraining.

Prior work disclosure: the author previously released chem-ood-eval (frozen random, scaffold and cluster splits with random-forest baselines for his own chemistry datasets). No code or data from it was used here. All code in this project was written during the sprint (13 to 15 November 2026).

## References

Bemis, G. W., & Murcko, M. A. (1996). The properties of known drugs. 1. Molecular frameworks. Journal of Medicinal Chemistry, 39(15), 2887–2893.

Bouthillier, X., et al. (2021). Accounting for variance in machine learning benchmarks. Proceedings of MLSys 2021.

Butina, D. (1999). Unsupervised data base clustering based on Daylight's fingerprint and Tanimoto similarity. Journal of Chemical Information and Computer Sciences, 39(4), 747–750.

Delaney, J. S. (2004). ESOL: Estimating aqueous solubility directly from molecular structure. Journal of Chemical Information and Computer Sciences, 44(3), 1000–1005.

Dodge, J., Gururangan, S., Card, D., Schwartz, R., & Smith, N. A. (2019). Show your work: Improved reporting of experimental results. EMNLP 2019.

Gelman, A., & Loken, E. (2013). The garden of forking paths. Unpublished manuscript, Columbia University.

Guo, Q., Hernandez-Hernandez, S., & Ballester, P. J. (2024). Scaffold splits overestimate virtual screening performance. arXiv:2406.00873.

Kapoor, S., & Narayanan, A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. Patterns, 4(9), 100804.

Mobley, D. L., & Guthrie, J. P. (2014). FreeSolv: A database of experimental and calculated hydration free energies. Journal of Computer-Aided Molecular Design, 28(7), 711–720.

Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python. Journal of Machine Learning Research, 12, 2825–2830.

RDKit: Open-source cheminformatics. https://www.rdkit.org

Sheridan, R. P. (2013). Time-split cross-validation as a method for estimating the goodness of prospective prediction. Journal of Chemical Information and Modeling, 53(4), 783–790.

Sheridan, R. P., Feuston, B. P., Maiorov, V. N., & Kearsley, S. K. (2004). Similarity to molecules in the training set is a good discriminator for prediction accuracy in QSAR. Journal of Chemical Information and Computer Sciences, 44(6), 1912–1928.

Simmons, J. P., Nelson, L. D., & Simonsohn, U. (2011). False-positive psychology. Psychological Science, 22(11), 1359–1366.

Wu, Z., et al. (2018). MoleculeNet: A benchmark for molecular machine learning. Chemical Science, 9(2), 513–530.

## Appendix A. Limitations and Dual-Use / Ethical Considerations

**What the result does and does not establish.** It establishes, for these three datasets, one model family and the stated split constructions, the measured size of split-, seed- and selection-driven differences in reported R², and the measured accuracy of one training-free reweighting estimate. It does not establish that any published paper inflated its results, and it does not establish the size of the effect for other model classes, tasks or domains. Nor does it measure any language model's behaviour.

**Dual use.** The same measurements show how to make a weak model look strong: choose the split and report the best seed. This is already widely known and requires no tooling. On balance we expect publishing the size of the effect, together with a cheap check that exposes it, to help readers and reviewers more than it helps someone gaming a benchmark. We release no tool designed to optimise a reported score.

**Data and privacy.** Only public chemical-structure datasets were used. No personal data, human subjects or proprietary data were involved.

## Appendix B. Reweighting audit per setting

| Dataset       | Features   | Target split   |   Actual RMSE |   Naive (random) RMSE |   Reweighted estimate |   Share of gap closed |
|:--------------|:-----------|:---------------|--------------:|----------------------:|----------------------:|----------------------:|
| ESOL          | Morgan2048 | cluster        |         1.457 |                 1.159 |                 1.231 |                  0.24 |
| ESOL          | Morgan2048 | scaffold       |         1.407 |                 1.159 |                 1.231 |                  0.29 |
| ESOL          | RDKitDesc  | cluster        |         0.722 |                 0.63  |                 0.714 |                  0.91 |
| ESOL          | RDKitDesc  | scaffold       |         0.788 |                 0.63  |                 0.692 |                  0.39 |
| FreeSolv      | Morgan2048 | cluster        |         2.465 |                 2.239 |                 2.937 |                 -1.08 |
| FreeSolv      | Morgan2048 | scaffold       |         2.429 |                 2.239 |                 2.613 |                  0.03 |
| FreeSolv      | RDKitDesc  | cluster        |         1.46  |                 1.253 |                 1.584 |                  0.4  |
| FreeSolv      | RDKitDesc  | scaffold       |         1.537 |                 1.253 |                 1.449 |                  0.69 |
| Lipophilicity | Morgan2048 | cluster        |         1.019 |                 0.824 |                 1.067 |                  0.76 |
| Lipophilicity | Morgan2048 | scaffold       |         0.929 |                 0.824 |                 0.914 |                  0.86 |
| Lipophilicity | RDKitDesc  | cluster        |         0.846 |                 0.693 |                 0.863 |                  0.89 |
| Lipophilicity | RDKitDesc  | scaffold       |         0.765 |                 0.693 |                 0.758 |                  0.9  |

Table 3. Actual group-split RMSE (mean over 20 seeds) vs. the naive random-split RMSE and the similarity-reweighted estimate. "Share of gap closed" = 1 − |reweighted − actual| / |naive − actual|. A negative value means the estimate was worse than naive.

## LLM Usage Statement

An AI assistant was used to write the experiment and analysis code and to draft the text of this report under the author's direction. Every number in the report is read from the result files produced by the released scripts (results/*.csv). The author is responsible for the final text.
