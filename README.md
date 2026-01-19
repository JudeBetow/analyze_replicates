*analyze_replicates.py — README*
analyze_replicates.py is a small, focused Python tool to analyse MD-derived protein–ligand (PL) contact data for hit complexes across replicate simulations.
It computes a delta disruption metric (the difference in contact occupancy between ligand-bound and apo/baseline ensembles), applies user-defined thresholds to classify residue pairs
as disrupted / partial / preserved, and writes both tabular results and a publication-ready heatmap that visualises disruption across the protein.
The script accepts either precomputed PL-contacts CSVs or MD trajectories (via MDAnalysis) and supports replica offsets, multiple normalisation options,
and a range of output files for downstream analysis.
