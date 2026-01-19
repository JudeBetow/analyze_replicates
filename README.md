# Description

`analyze_replicates.py` is a small, focused Python tool to analyse MD-derived protein–ligand (PL) contact data for hit complexes across replicate simulations (mainly for PPI studies).
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

```python -m pip install pandas numpy matplotlib seaborn tqdm MDAnalysis```

for trajectory support:

```python -m pip install MDAnalysis scipy # optional — for some statistics/normalizations```


# Installation
Clone your GitHub repo with:
```git clone https://github.com/<your-org>/<repo>.git
cd <repo>
# run the script with Python 3.8+
python analyze_replicates.py -h
```
# Input formats & expected files
# 1) Precomputed PL-contacts CSV mode (recommended if you already computed contacts with Desmond, for example)
- One or more CSV files with per-frame contact information (the format expected by the script is: `frame,residue_A,residue_B,contact_flag` or similar — see example below).
- The script will detect replicate names from the filenames (or you can pass an explicit mapping).
```
frame,residue_A,residue_B,contact
0,ALA45,LYS89,1
1,ALA45,LYS89,0
...
```
# 2) Trajectory mode (compute contacts on the fly)
- Coordinate/trajectory: e.g., traj.xtc, traj.dcd.
- Topology: e.g., top.pdb, top.psf (a file MDAnalysis can read). I strongly recommend VMD to convert your CMS file to PDB format. Maestro converts, but there will be a mismatch in atom count between the PDB and XTC files.
- Selection strings/contact definition (e.g., residue sets or atom groups for the protein and ligand).

MDAnalysis must be installed for this mode.

A. Using precomputed PL-contacts CSVs (replace with your right paths):
```
  --root /path/to/root/dir \ # root containing your replica sub-directories
  --outdir /path/to/outdir/analysis_output/ \
  --pp-contacts /path/to/pp-cont/pp-file.csv/ \ 
  --pp-distances /path/to/pp-dst/dist-file.csv/ \
  --pairs "E417-A30,E449-A38,E493-A34,E505-A353,E495-A353" \ # optional - if not specified, all 14 pairs will be computed as default.
  --metric delta \
  --delta-threshold 10 \
  --min-baseline-freq 0.05 \
  --plot-metric disruption_pct_by_lig \
  --log-level INFO
```
`--metric` chooses which disruption metric to compute (default: disruption_pct_by_lig).
`--threshold` (percent points) controls classification into disrupted/partial/preserved.

B. Computing contacts from trajectories (MDAnalysis required
```
python analyze_replicates.py \
  --traj replicate1.xtc replicate2.xtc replicate3.xtc \
  --top topology.pdb \
  --outdir path/to/outdir/ \
  --pairs "E417-A30,E449-A38,E493-A34,E505-A353,E495-A353" \ # optional - if not specified, all 14 pairs will be computed as default.
  --metric delta \
  --delta-threshold 10 \
  --min-baseline-freq 0.05 \
  --plot-metric disruption_pct_by_lig \
  --log-level INFO
```
- `--traj` accepts one or multiple trajectories (for replicate sets)

# Main outputs (file names used by the script)
`<replicate>_occupancy.csv` — per-replicate ligand-occupancy tables (per residue/pair).
`occupancy_aggregated.csv` — aggregated occupancy across replicates (mean, stdev, counts).
`mapped_pp_contacts.csv` — mapped protein–protein contacts (if mapping was supplied/precomputed).
`disruption_metrics.csv` — final per-pair statistics that include apo occupancy, ligand occupancies, delta, normalized scores, and classification (disrupted/partial/preserved).
`disruption_heatmap.png` — heatmap visualising the chosen metric (e.g., disruption_pct_by_lig) for selected residue pairs (or all 14 pairs if not specified); includes a side colour bar showing the classification.
`analysis.log` — execution log with parameters, timestamps and brief progress/error messages.

Each CSV contains human-readable headers that can be opened in Excel, pandas, or other tools for further analysis.

# Troubleshooting & tips
- MDAnalysis errors: check your topology/trajectory formats; MDAnalysis supports many formats but may require additional optional dependencies.
- Large trajectories: subsample frames (e.g., every 10–50 ps) to reduce runtime and memory. Use the --frames-per-snapshot or equivalent option.
- Unexpected classification: verify the apo baseline file — if baseline occupancy is zero for many pairs, delta will equal occupancy_bound. The script filters very-low baseline occupancies but you can adjust the filter.
- Visual appearance: heatmap colormap and clustering options may be configurable; check the script headers or modify the plotting section if you need publication-specific styling.

# Extending and contributing
Contributions, issue reports, or feature requests are welcome.

If you add features, please:
Fork the repo. Create a feature branch.
Submit a pull request with tests/examples.

# License
This project is released under the MIT License — see `LICENSE`.

# Contact / Citation
If you use this script in published work, please cite the GitHub repository and include a short methods note describing:
(a) PL-contact replica used,

(b) frames sampled per replicate,

(c) threshold for disruption classification, and

(d) whether contacts were precomputed or calculated on-the-fly.

For questions or help: open an issue on the GitHub repo or contact the author (betow.jude@ubuea.cm).
