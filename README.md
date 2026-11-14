# Forking Paths in Molecular Property Benchmarks

Code, figures and result tables for the [Apart AI x Epistemics](https://apartresearch.com/sprints/ai-x-epistemics-research-sprint-2026-11-13-to-2026-11-15) sprint report (13–15 November 2026).

**Report PDF:** [report/report.pdf](report/report.pdf) (also as [report/report.md](report/report.md)).

Author: Byungwoong Yoo (BigWonder, Seoul). Solo project. All code in this repository was written during the sprint.

The report measures how far a reported test R² on three public MoleculeNet regression sets (ESOL, FreeSolv, Lipophilicity) can move under defensible choices of split, seed and features, and tests a training-free nearest-neighbour similarity profile as an audit signal.

Prior work, not used here: [chem-ood-eval](https://github.com/ByungwoongYoo/chem-ml-datasets/tree/main/chem-ood-eval).

## Reproduce (CPU only, about 40 minutes on 8 cores)

```bash
pip install -r requirements.txt
bash data/download_data.sh
python3 code/run_experiments.py 20   # writes results/runs.csv and results/per_molecule.csv.gz
python3 code/analyze.py              # writes results/table*.csv and figures/fig*.png
```

Figures used in the report are in `figures/` as SVG (same pixels as the originals). `analyze.py` also writes PNG files when re-run.

`code/build_report.py` can regenerate the Word document if you place the official Apart sprint template at `report/apart_template_original.docx`. That template is not included here. The finished PDF and Markdown are already in `report/`.

## Data

The raw MoleculeNet CSVs are not in this repository. `data/download_data.sh` downloads them from the public DeepChem S3 bucket and checks these SHA-256 hashes (files retrieved 2026-11-14):

```
8c06a76f0c6487d29ab0f903e6a7a7139f189ab3c1178f159c8be8964602f189  delaney-processed.csv
ab5895d914ee87cb563bd7b9611e869527bba45bec6b014d34dc495a0f9dcb72  SAMPL.csv
aed41590cb30609d51d8e08ad3ff06495a76e80e211358801f596b10da69bacd  Lipophilicity.csv
```

- ESOL: Delaney (2004), file `delaney-processed.csv`
- FreeSolv: Mobley & Guthrie (2014), file `SAMPL.csv`
- Lipophilicity: ChEMBL-derived, file `Lipophilicity.csv`

The datasets remain under their own terms. Experiment outputs in `results/` (including the `results/per_molecule_*.csv` files, which store predictions and nearest-train Tanimoto similarities, not SMILES) are included so the tables can be recomputed without retraining. `code/run_experiments.py` writes the same table as `results/per_molecule.csv.gz`.

## License

The code in this repository is MIT (see LICENSE).
