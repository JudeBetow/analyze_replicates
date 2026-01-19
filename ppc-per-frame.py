#!/usr/bin/env python3
"""
Export per-frame PP contacts (atom-atom min distance <= cutoff) as CSV
suitable for analyse_replicates.py.

Writes either:
 - per-rep CSV with header: replica,chain1,res1,chain2,res2,Frame
   (default, safe)
 - merged/global CSV (no replica column) with header: chain1,res1,chain2,res2,Frame
   (use --merged and supply --rep-frame-lens)

Usage examples:

# per-rep mode (default)
python ppc-per-frame.py \
  --topology 6m0j-topology-1.pdb \
  --trajs 6m0j-traj-rep1.xtc,6m0j-traj-rep2.xtc,6m0j-traj-rep3.xtc \
  --out csv/pp_contacts_perframe.csv \
  --pairs-file pairs.txt \
  --cutoff 4.0 --subsample 1

# merged mode (you concatenated rep1,rep2,rep3 into merged_all.xtc and want global frame indices)
python ppc-per-frame.py \
  --topology 6m0j-topology-1.pdb \
  --trajs merged_all.xtc \
  --out csv/pp_contacts_merged.csv \
  --pairs-file pairs.txt --merged \
  --rep-frame-lens 1001,1001,1001 \
  --cutoff 4.0 --subsample 1

Pairs file format (optional): one pair per line, whitespace separated:
E 493 A 35
E 417 A 30
...
If --pairs-file is omitted the script uses a small default list (RBD-ACE2).

Output columns (per-rep mode):
 replica, chain1, res1, chain2, res2, Frame

Output columns (merged/global mode):
 chain1, res1, chain2, res2, Frame

Note: Frame indices are 0-based to match MDAnalysis.
"""
import argparse, csv, sys
from pathlib import Path
from typing import List, Tuple
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances as mda_distances

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

def read_pairs_file(path: Path) -> List[Tuple[str,int,str,int]]:
    pairs = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                pairs.append((parts[0], int(parts[1]), parts[2], int(parts[3])))
    return pairs

def compute_contacts_for_traj(topo: str, traj: str, keypairs, cutoff: float=4.0, subsample: int=1):
    """
    Generator that yields tuples (replica/traj, pair, frame) for frames where pair is in contact.
    pair is (chain1,res1,chain2,res2)
    """
    u = mda.Universe(topo, traj)
    # Build index arrays for each pair (heavy atoms, fallback to resid)
    pair_index_sets = []
    for pair in keypairs:
        c1,r1,c2,r2 = pair
        sel1 = u.select_atoms(f"chainID {c1} and resid {r1} and not name H*")
        sel2 = u.select_atoms(f"chainID {c2} and resid {r2} and not name H*")
        # fallback if empty
        if len(sel1)==0:
            sel1 = u.select_atoms(f"resid {r1} and not name H*")
        if len(sel2)==0:
            sel2 = u.select_atoms(f"resid {r2} and not name H*")
        if len(sel1)==0 or len(sel2)==0:
            # store None to indicate skip
            pair_index_sets.append((pair, None, None))
            continue
        pair_index_sets.append((pair, sel1.indices.copy(), sel2.indices.copy()))
    total_frames = 0
    for ts in u.trajectory[::subsample]:
        frame = int(ts.frame)
        total_frames += 1
        positions = u.atoms.positions
        box = ts.dimensions if hasattr(ts, 'dimensions') else None
        for pair, idx1, idx2 in pair_index_sets:
            if idx1 is None or idx2 is None:
                continue
            pos1 = positions[idx1]
            pos2 = positions[idx2]
            try:
                dmat = mda_distances.distance_array(pos1, pos2, box=box)
            except Exception:
                dmat = mda_distances.distance_array(pos1, pos2, box=None)
            if np.any(dmat <= cutoff):
                yield (pair, frame)
    return

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--topology", "-t", required=True)
    p.add_argument("--trajs", "-x", required=True, help="comma separated list of traj files (one per replica) OR a single merged xtc when --merged is used")
    p.add_argument("--pairs-file", help="optional file with pairs (CH1 RES1 CH2 RES2 per line)")
    p.add_argument("--out", "-o", default="pp_contacts_out.csv", help="output csv path")
    p.add_argument("--cutoff", type=float, default=4.0, help="atom-atom cutoff (Å)")
    p.add_argument("--subsample", type=int, default=1)
    p.add_argument("--merged", action="store_true", help="If set, produce global frames (no replica column). Requires --rep-frame-lens")
    p.add_argument("--rep-frame-lens", help="comma-separated list of frame counts per replica (only used with --merged). e.g. 1001,1001,1001")
    args = p.parse_args()

    topo = args.topology
    trajs = [t.strip() for t in args.trajs.split(",")]
    outpath = Path(args.out)
    cutoff = args.cutoff
    subsample = max(1, int(args.subsample))

    if args.pairs_file:
        keypairs = read_pairs_file(Path(args.pairs_file))
    else:
        keypairs = DEFAULT_KEYPAIRS

    # If merged mode: expect a single trajectory (merged) and rep-frame-lens provided
    if args.merged:
        if len(trajs) != 1:
            print("When --merged is set, --trajs should give the single merged xtc file.", file=sys.stderr)
            sys.exit(2)
        if not args.rep_frame_lens:
            print("When --merged is set, you must provide --rep-frame-lens with per-rep frame counts.", file=sys.stderr)
            sys.exit(2)
        rep_lens = [int(s) for s in args.rep_frame_lens.split(",")]
        # compute offsets cumulative
        offsets = [0]
        for L in rep_lens[:-1]:
            offsets.append(offsets[-1] + L)
        # We'll iterate the single merged trajectory and figure out which replica the frame belongs to -> global frame = running index
        merged_traj = trajs[0]
        u = mda.Universe(topo, merged_traj)
        # build indices per pair once
        pair_index_sets = []
        for pair in keypairs:
            c1,r1,c2,r2 = pair
            sel1 = u.select_atoms(f"chainID {c1} and resid {r1} and not name H*")
            sel2 = u.select_atoms(f"chainID {c2} and resid {r2} and not name H*")
            if len(sel1)==0:
                sel1 = u.select_atoms(f"resid {r1} and not name H*")
            if len(sel2)==0:
                sel2 = u.select_atoms(f"resid {r2} and not name H*")
            if len(sel1)==0 or len(sel2)==0:
                pair_index_sets.append((pair, None, None))
                continue
            pair_index_sets.append((pair, sel1.indices.copy(), sel2.indices.copy()))
        # write header (no replica column)
        with outpath.open("w", newline='') as csvf:
            writer = csv.writer(csvf)
            writer.writerow(["chain1","res1","chain2","res2","Frame"])
            # iterate merged trajectory
            master_frame = -1
            for ts in u.trajectory[::subsample]:
                master_frame += 1
                positions = u.atoms.positions
                box = ts.dimensions if hasattr(ts, 'dimensions') else None
                for pair, idx1, idx2 in pair_index_sets:
                    if idx1 is None or idx2 is None:
                        continue
                    pos1 = positions[idx1]
                    pos2 = positions[idx2]
                    try:
                        dmat = mda_distances.distance_array(pos1, pos2, box=box)
                    except Exception:
                        dmat = mda_distances.distance_array(pos1, pos2, box=None)
                    if np.any(dmat <= cutoff):
                        c1,r1,c2,r2 = pair
                        writer.writerow([c1, r1, c2, r2, master_frame])
        print("Wrote merged per-frame contacts to", outpath)
        return

    # non-merged (per-rep) mode: produce replica column
    with outpath.open("w", newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow(["replica","chain1","res1","chain2","res2","Frame"])
        for traj in trajs:
            print("Processing", traj)
            for pair, frame in compute_contacts_for_traj(topo, traj, keypairs, cutoff=cutoff, subsample=subsample):
                c1,r1,c2,r2 = pair
                writer.writerow([traj, c1, r1, c2, r2, frame])
    print("Wrote per-frame contacts (per-rep) to", outpath)

if __name__ == "__main__":
    main()
