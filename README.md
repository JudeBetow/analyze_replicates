# Description

'analyze_replicates.py' is a small, focused Python tool to analyse MD-derived protein–ligand (PL) contact data for hit complexes across replicate simulations (mainly for PPI studies).
It computes a delta disruption metric (the difference in contact occupancy between ligand-bound and apo/baseline ensembles), applies user-defined thresholds to classify residue pairs
as disrupted / partial / preserved, and writes both tabular results and a publication-ready heatmap that visualises disruption across the protein.
The script accepts either precomputed PL-contacts CSVs or MD trajectories (via MDAnalysis) and supports replica offsets, multiple normalisation options,
and a range of output files for downstream analysis.

# Features

- Computes per-frame / per-replicate ligand-occupancy and aggregated occupancy tables.
- Calculates several disruption metrics (including disruption_pct_by_lig) based on the delta between ligand-bound and apo occupancies.
- Classifies residue/residue (or residue/atom) pairs as disrupted / partial / preserved using a user threshold and filters low-confidence pairs.
 Produces CSV outputs at multiple levels and a heatmap PNG summarising the disruption metric and classification.
- Can read precomputed PL-contacts CSVs or compute contacts from trajectories using MDAnalysis.
- Handles replica frame offsets and writes a detailed log for reproducibility.

# Requirements

Python 3.8+
Key Python packages (install via pip):
- pandas
- numpy
- MDAnalysis (optional, required only if computing contacts from trajectories),
- matplotlib
- seaborn
- tqdm (optional, for progress bars)

You can install the main dependencies with:
python -m pip install pandas numpy matplotlib seaborn tqdm MDAnalysis
for trajectory support:
python -m pip install MDAnalysis
scipy (optional — for some statistics/normalizations)

# Installation
Clone your GitHub repo with:
'''
git clone https://github.com/<your-org>/<repo>.git
cd <repo>
# run the script with Python 3.8+
python analyze_replicates.py -h
'''
No formal install is required — the script runs as a single Python module.
