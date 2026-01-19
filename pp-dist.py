#!/usr/bin/env python3
"""
Compute per-pair minimum atom-atom distance statistics for a list of key pairs
(using MDAnalysis.distance_array with PBC).

Inputs:
 - topology file (e.g., PSF, PDB)
 - one or more trajectory files (e.g., XTC)
 - optional pairs file (otherwise uses built-in default list)

Outputs:
 - outdir/csv/pp_pair_distances_perrep.csv  (rows: replica, pair, mean_min_dist, std_min_dist, mean_min_dist_when_contact, n_contact_frames, frames_examined)
 - outdir/csv/pp_pair_distances_aggregated.csv (rows: pair, mean_of_rep_means, std_of_rep_means, mean_of_rep_contact_means, n_reps)

Example usage:
    python pp-dist.py \
  --topology 6m0j-topology-1.pdb \
  --trajs 6m0j-traj-rep1.xtc,6m0j-traj-rep2.xtc,6m0j-traj-rep3.xtc \
  --pairs pairs.txt \
  --outdir analysis_output \
  --cutoff 4.0 --subsample 1
Requires MDAnalysis, numpy.
"""
import argparse, csv, os
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances as mda_distances

# default pairs if none provided
DEFAULT_KEYPAIRS = [
    ("E", 417, "A", 30),
    ("E", 449, "A", 38),
    ("E", 493, "A", 35),
    ("E", 493, "A", 31),
    ("E", 493, "A", 34),
    ("E", 498, "A", 353),
    ("E", 500, "A", 355),
    ("E", 501, "A", 353),
    ("E", 501, "A", 41),
    ("E", 505, "A", 353),
    ("E", 486, "A", 82),
    ("E", 487, "A", 24),
    ("E", 487, "A", 83),
    ("E", 495, "A", 353),
]

def read_pairs_file(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                pairs.append((parts[0], int(parts[1]), parts[2], int(parts[3])))
    return pairs

def compute_min_dists_for_pair(u, pair, cutoff=None, subsample=1):
    """
    For a loaded Universe u, compute per-frame minimum atom-atom distance for 'pair'.
    pair = (chain1, resid1, chain2, resid2)
    returns: list of per-frame minima (len = frames used), plus list of boolean is_contact for each frame (min <= cutoff) if cutoff provided
    """
    c1, r1, c2, r2 = pair
    sel1 = u.select_atoms(f"chainID {c1} and resid {r1} and not name H*")
    sel2 = u.select_atoms(f"chainID {c2} and resid {r2} and not name H*")
    # fallback if empty
    if len(sel1) == 0:
        sel1 = u.select_atoms(f"resid {r1} and not name H*")
    if len(sel2) == 0:
        sel2 = u.select_atoms(f"resid {r2} and not name H*")
    if len(sel1) == 0 or len(sel2) == 0:
        # return empty lists; caller should detect
        return [], []
    idx1 = sel1.indices
    idx2 = sel2.indices

    min_dists = []
    contact_flags = []
    for ts in u.trajectory[::subsample]:
        pos = u.atoms.positions
        pos1 = pos[idx1]
        pos2 = pos[idx2]
        try:
            dmat = mda_distances.distance_array(pos1, pos2, box=ts.dimensions)
        except Exception:
            dmat = mda_distances.distance_array(pos1, pos2, box=None)
        minval = float(np.min(dmat))
        min_dists.append(minval)
        contact_flags.append(minval <= cutoff if cutoff is not None else False)
    return min_dists, contact_flags

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--topology", "-t", required=True)
    p.add_argument("--trajs", "-x", required=True, help="comma-separated xtc files")
    p.add_argument("--pairs", help="optional pairs file, format: CH1 RES1 CH2 RES2 per line")
    p.add_argument("--outdir", "-o", default="pp_pair_distances_out")
    p.add_argument("--cutoff", type=float, default=4.0)
    p.add_argument("--subsample", type=int, default=1)
    args = p.parse_args()

    topo = args.topology
    trajs = [t.strip() for t in args.trajs.split(",")]
    outdir = Path(args.outdir)
    csvdir = outdir / "csv"
    csvdir.mkdir(parents=True, exist_ok=True)

    pairs = DEFAULT_KEYPAIRS if not args.pairs else read_pairs_file(args.pairs)

    perrep_rows = []
    agg_by_pair = {pair: [] for pair in pairs}
    agg_contact_means_by_pair = {pair: [] for pair in pairs}

    for tr in trajs:
        if not Path(tr).exists():
            raise FileNotFoundError(f"Trajectory not found: {tr}")
        print("Processing", tr)
        u = mda.Universe(topo, tr)
        frames_examined = len(u.trajectory[::args.subsample])
        for pair in pairs:
            min_dists, flags = compute_min_dists_for_pair(u, pair, cutoff=args.cutoff, subsample=args.subsample)
            if len(min_dists) == 0:
                # could not select atoms
                mean_min = np.nan
                std_min = np.nan
                contact_mean = np.nan
                n_contact = 0
            else:
                arr = np.array(min_dists)
                mean_min = float(np.nanmean(arr))
                std_min = float(np.nanstd(arr))
                contact_flags = np.array(flags)
                n_contact = int(np.sum(contact_flags))
                contact_mean = float(n_contact) / float(len(contact_flags)) if len(contact_flags)>0 else 0.0
            perrep_rows.append({
                'replica': tr, 'chain1': pair[0], 'res1': pair[1], 'chain2': pair[2], 'res2': pair[3],
                'mean_min_dist': mean_min, 'std_min_dist': std_min,
                'mean_min_dist_when_contact': (float(np.nanmean(np.array(min_dists)[np.array(flags)])) if any(flags) else np.nan),
                'n_contact_frames': n_contact, 'frames_examined': frames_examined
            })
            if not np.isnan(mean_min):
                agg_by_pair[pair].append(mean_min)
            if not np.isnan(contact_mean):
                agg_contact_means_by_pair[pair].append(contact_mean)

    # write per-rep CSV
    perrep_csv = csvdir / "pp_pair_distances_perrep.csv"
    with perrep_csv.open("w", newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['replica','chain1','res1','chain2','res2','mean_min_dist','std_min_dist','mean_min_dist_when_contact','n_contact_frames','frames_examined'])
        w.writeheader()
        for r in perrep_rows:
            w.writerow(r)
    print("Wrote", perrep_csv)

    # write aggregated CSV: mean/std across replicas of mean_min_dist and of contact means
    agg_csv = csvdir / "pp_pair_distances_aggregated.csv"
    with agg_csv.open("w", newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['chain1','res1','chain2','res2','mean_of_rep_mean_min_dist','std_of_rep_mean_min_dist','mean_contact_freq','std_contact_freq','n_reps'])
        for pair in pairs:
            means = np.array(agg_by_pair.get(pair, []), dtype=float)
            contact_means = np.array(agg_contact_means_by_pair.get(pair, []), dtype=float)
            mean_of_means = float(np.nanmean(means)) if means.size>0 else np.nan
            std_of_means = float(np.nanstd(means)) if means.size>0 else np.nan
            cm_mean = float(np.nanmean(contact_means)) if contact_means.size>0 else np.nan
            cm_std = float(np.nanstd(contact_means)) if contact_means.size>0 else np.nan
            w.writerow([pair[0], pair[1], pair[2], pair[3], mean_of_means, std_of_means, cm_mean, cm_std, len(means)])
    print("Wrote", agg_csv)
    print("Done.")

if __name__ == "__main__":
    main()
