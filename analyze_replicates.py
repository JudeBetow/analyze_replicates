#!/usr/bin/env python3
"""analyze_replicates.py
Analyze PL-Contacts across replicates and compute disruption metrics.

Usage examples (replace paths as needed):

python analyze_replicates.py \
  --root /home/jude/work/mdruns/reaxys20-plc \
  --outdir /home/jude/work/mdruns/Protein-6m0j/analysis_output/delta_metric \
  --pp-contacts /home/jude/work/mdruns/Protein-6m0j/analysis_output/csv/pp_contacts_perframe.csv \
  --pp-distances /home/jude/work/mdruns/Protein-6m0j/analysis_output/csv/pp_pair_distances_perrep.csv \
  --pairs "E417-A30,E449-A38,E493-A34,E505-A353,E495-A353" \
  --metric delta \
  --delta-threshold 10 \
  --min-baseline-freq 0.05 \
  --plot-metric disruption_pct_by_lig \
  --log-level INFO
"""
import sys
import os
from pathlib import Path
import argparse
import logging
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from collections import defaultdict
from typing import Optional, Dict, Tuple, Set

# Try importing MDAnalysis only when needed
try:
    import MDAnalysis as mda
except Exception:
    mda = None

# -------------------------
# Defaults / canonical pairs
# -------------------------
DEFAULT_CONTACT_FILES = [
    "PL-Contacts_HBond.dat",
    "PL-Contacts_Hydrophobic.dat",
    "PL-Contacts_Pi-Pi.dat",
    "PL-Contacts_Pi-Cation.dat",
    "PL-Contacts_Ionic.dat",
    "PL-Contacts_WaterBridge.dat",
    "PL-Contacts_Metal.dat",
]

# KEY PAIR format: (chain1, res1, chain2, res2)
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

# -------------------------
# Helper functions & mapping utilities
# -------------------------
_INT_RE = re.compile(r"(-?\d+)")

def infer_frame_count_from_dat_dir(repdir: Path) -> Optional[int]:
    maxf = -1
    found_any = False
    for f in repdir.glob("**/*.dat"):
        try:
            with f.open("r") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    m = _INT_RE.search(line)
                    if m:
                        found_any = True
                        try:
                            val = int(m.group(1))
                            if val > maxf:
                                maxf = val
                        except Exception:
                            pass
        except Exception:
            pass
    if not found_any:
        return None
    return int(maxf) + 1

def infer_frame_count_from_xtc(repdir: Path) -> Optional[int]:
    try:
        from MDAnalysis.coordinates.XTC import XTCReader
    except Exception:
        return None
    for xtc in repdir.glob("**/*.xtc"):
        try:
            xr = XTCReader(str(xtc))
            n = xr.n_frames
            xr.close()
            return int(n)
        except Exception:
            continue
    return None

def infer_rep_frame_counts(replicas: list, rep_frame_lens_arg: Optional[str], replica_order_arg: Optional[str]) -> Dict[str,int]:
    # optionally re-order by explicit replica_order
    if replica_order_arg:
        requested_order = [s.strip() for s in replica_order_arg.split(",")]
        name2path = {p.name:p for p in replicas}
        ordered = []
        for n in requested_order:
            if n in name2path:
                ordered.append(name2path[n])
            else:
                raise RuntimeError(f"--replica-order specified unknown replica name: {n}")
        replicas = ordered
    rep_names = [p.name for p in replicas]

    rep_frame_counts: Dict[str,int] = {}
    if rep_frame_lens_arg:
        parts = [int(s.strip()) for s in rep_frame_lens_arg.split(",")]
        if len(parts) != len(replicas):
            raise RuntimeError("Length of --rep-frame-lens does not match number of discovered replicas.")
        for name,cnt in zip(rep_names, parts):
            rep_frame_counts[name] = int(cnt)
        return rep_frame_counts

    # infer from .dat
    for p in replicas:
        cnt = infer_frame_count_from_dat_dir(p)
        if cnt is not None:
            rep_frame_counts[p.name] = int(cnt)

    # infer from xtc header if still missing
    for p in replicas:
        if p.name in rep_frame_counts:
            continue
        cnt = infer_frame_count_from_xtc(p)
        if cnt is not None:
            rep_frame_counts[p.name] = int(cnt)

    # unknowns -> set 0 (warn)
    missing = [p.name for p in replicas if p.name not in rep_frame_counts]
    if missing:
        logging.warning(f"Could not infer frame counts for replicas: {missing}. Set to 0; pass --rep-frame-lens to override.")
        for n in missing:
            rep_frame_counts[n] = 0
    return rep_frame_counts

def compute_offsets(rep_frame_counts: Dict[str,int], replica_names_in_order: list) -> Dict[str,int]:
    offsets = {}
    cur = 0
    for n in replica_names_in_order:
        offsets[n] = cur
        cur += int(rep_frame_counts.get(n,0) or 0)
    logging.info(f"Computed replica offsets: {offsets}")
    return offsets

def load_and_map_pp_contacts(pp_contacts_csv: Path, offsets: Dict[str,int]) -> Dict[Tuple[str,int,str,int], Set[int]]:
    """
    Load CSV that may have either:
    - replica,chain1,res1,chain2,res2,Frame  (per-rep)
    - chain1,res1,chain2,res2,Frame          (merged/global)
    Return mapping: pair -> set(global frames)
    """
    df = pd.read_csv(pp_contacts_csv)
    logging.info(f"Loaded pp_contacts CSV: {pp_contacts_csv} rows={len(df)} columns={list(df.columns)}")
    result = {}
    cols = [c.strip() for c in df.columns]
    df.columns = cols
    if 'replica' in df.columns:
        logging.info("Detected 'replica' column in pp_contacts CSV; mapping per-rep frames to global frames using offsets.")
        for idx, row in df.iterrows():
            repname = str(row['replica'])
            repkey = Path(repname).name
            offset = int(offsets.get(repkey, offsets.get(repname,0) or 0))
            try:
                frame = int(row['Frame'])
            except Exception:
                if 'frame' in df.columns:
                    frame = int(row['frame'])
                else:
                    logging.error(f"Cannot parse Frame value for row {idx}: {row.to_dict()}")
                    continue
            global_frame = frame + offset
            try:
                c1 = str(row['chain1']); c2 = str(row['chain2'])
                r1 = int(row['res1']); r2 = int(row['res2'])
            except Exception as e:
                logging.error(f"Malformed pair row {idx}: {row.to_dict()} -> {e}")
                continue
            key = (c1, r1, c2, r2)
            result.setdefault(key, set()).add(int(global_frame))
    else:
        logging.info("No 'replica' column found; reading Frame as global frame indices directly.")
        for idx,row in df.iterrows():
            try:
                frame = int(row['Frame'])
                c1 = str(row['chain1']); c2 = str(row['chain2'])
                r1 = int(row['res1']); r2 = int(row['res2'])
            except Exception as e:
                logging.error(f"Malformed row {idx}: {row.to_dict()} -> {e}")
                continue
            key = (c1, r1, c2, r2)
            result.setdefault(key, set()).add(int(frame))
    logging.info(f"Mapped pp_contacts into {len(result)} pairs with frame sets.")
    return result

def write_mapped_pp_csv(outpath: Path, pp_frames_mapped: Dict[Tuple[str,int,str,int], Set[int]]):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w", newline='') as fh:
        import csv
        w = csv.writer(fh)
        w.writerow(["chain1","res1","chain2","res2","Frames"])
        for (c1,r1,c2,r2), frames in sorted(pp_frames_mapped.items(), key=lambda x: (x[0][0], x[0][1], x[0][2], x[0][3])):
            flist = sorted(int(x) for x in frames)
            w.writerow([c1, r1, c2, r2, ";".join(str(x) for x in flist)])
    logging.info(f"Wrote mapped pp contacts CSV: {outpath}")

# ----------------- distance aggregation helper -----------------
def aggregate_pp_distance_stats(dist_csv_path: Path) -> Dict[Tuple[str,int,str,int], dict]:
    """
    Aggregate per-rep pp distance CSV into per-pair stats.

    Expects CSV with columns:
      replica,chain1,res1,chain2,res2,mean_min_dist,std_min_dist,mean_min_dist_when_contact,n_contact_frames,frames_examined
    Returns mapping:
      (chain1,res1,chain2,res2) -> {
        'mean_min_dist_mean': float,
        'mean_min_dist_std': float,
        'mean_min_dist_when_contact_mean': float or None,
        'n_contact_frames_total': int,
        'frames_examined_total': int,
        'apo_freq': float
      }
    """
    df = pd.read_csv(dist_csv_path)
    # normalize column names for resilience
    df_cols = {c.lower(): c for c in df.columns}
    def col(name):
        return df_cols.get(name.lower())
    required = [col('chain1'), col('res1'), col('chain2'), col('res2')]
    if not all(required):
        logging.error("pp-distances CSV missing required columns (chain1,res1,chain2,res2).")
        return {}

    # ensure numeric
    for num_col in ['mean_min_dist','std_min_dist','mean_min_dist_when_contact','n_contact_frames','frames_examined']:
        if col(num_col) in df.columns:
            df[col(num_col)] = pd.to_numeric(df[col(num_col)], errors='coerce')

    grouped = df.groupby([col('chain1'), col('res1'), col('chain2'), col('res2')])
    out = {}
    for name, g in grouped:
        c1, r1, c2, r2 = name
        mm_mean = float(g[col('mean_min_dist')].dropna().mean()) if col('mean_min_dist') in df.columns else None
        mm_std = float(g[col('std_min_dist')].dropna().mean()) if col('std_min_dist') in df.columns else None
        if col('mean_min_dist_when_contact') in df.columns:
            mw_mean = float(g[col('mean_min_dist_when_contact')].dropna().mean()) if not g[col('mean_min_dist_when_contact')].dropna().empty else None
        else:
            mw_mean = None
        n_total = int(g[col('n_contact_frames')].dropna().sum()) if col('n_contact_frames') in df.columns else 0
        frames_total = int(g[col('frames_examined')].dropna().sum()) if col('frames_examined') in df.columns else 0
        apo_freq = (n_total / frames_total) if frames_total else 0.0
        out[(str(c1), int(r1), str(c2), int(r2))] = {
            'mean_min_dist_mean': mm_mean,
            'mean_min_dist_std': mm_std,
            'mean_min_dist_when_contact_mean': mw_mean,
            'n_contact_frames_total': n_total,
            'frames_examined_total': frames_total,
            'apo_freq': apo_freq
        }
    return out

# ----------------- PDB/PL-Contacts parsing helpers -----------------
def canonical_cols(cols):
    """Map input column names to canonical names used in the script."""
    out = []
    for c in cols:
        s = str(c).strip()
        s_low = s.lower()
        if "frame" in s_low:
            out.append("Frame")
        elif "resid" in s_low or "residue#" in s_low or re.match(r"^residue$", s_low):
            out.append("Residue")
        elif s_low in ("chain", "chainid"):
            out.append("Chain")
        elif "resname" in s_low:
            out.append("ResName")
        elif "atom" in s_low:
            out.append("AtomName")
        else:
            # keep original but normalized
            out.append(re.sub(r"\W+", "_", s))
    return out

def parse_dat_file(path: Path) -> pd.DataFrame:
    """
    Parse a Desmond PL-Contacts .dat file robustly.
    - Finds the first comment header line (starting with '#') and uses it as header.
    - Skips other comment lines.
    - Splits on whitespace or tab or comma.
    - Returns DataFrame with canonical columns where possible.
    """
    txt = path.read_text(encoding="utf8", errors="ignore").splitlines()
    header_line = None
    data_lines = []
    for ln in txt:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            # take the first comment line with column names as header
            if header_line is None:
                header_line = s.lstrip("#").strip()
            continue
        data_lines.append(s)

    if header_line is None:
        logging.warning(f"{path}: no header line starting with '#' found; attempting to infer columns")
        if not data_lines:
            return pd.DataFrame()
        first_cols = re.split(r"[,\t\s]+", data_lines[0].strip())
        header_line = " ".join(f"col{i}" for i in range(len(first_cols)))

    header_cols = re.split(r"[,\t\s]+", header_line)
    header_cols = canonical_cols(header_cols)

    records = []
    for i, ln in enumerate(data_lines, start=1):
        cols = re.split(r"[,\t\s]+", ln.strip())
        if len(cols) != len(header_cols):
            if len(cols) > len(header_cols):
                cols = cols[:len(header_cols)]
            else:
                logging.debug(f"{path.name}: skip line {i} (columns mismatch {len(cols)} vs {len(header_cols)})")
                continue
        records.append(cols)

    if not records:
        return pd.DataFrame(columns=header_cols)

    df = pd.DataFrame(records, columns=header_cols)

    if "Frame" in df.columns:
        df["Frame"] = pd.to_numeric(df["Frame"], errors="coerce").astype("Int64")
    if "Residue" in df.columns:
        def resid_to_int(x):
            try:
                return int(re.sub(r"\D+", "", str(x)))
            except Exception:
                return pd.NA
        df["Residue"] = df["Residue"].apply(resid_to_int).astype("Int64")
    if "Chain" in df.columns:
        df["Chain"] = df["Chain"].astype(str).str.strip()

    return df

def discover_replicas(root: Path, explicit=None):
    """Return list of replicate directories (Path objects)."""
    if explicit:
        return [Path(p) for p in explicit]
    if not root.exists():
        raise FileNotFoundError(f"root {root} does not exist")
    return sorted([p for p in root.iterdir() if p.is_dir()])

def collect_ligand_frames_from_contact_files(rep_dir: Path, contact_files):
    """
    For a replicate dir, parse all CONTACT dat files and return:
    - dict mapping (chain,resid) -> set(frames)
    - max_frame found in files (int or None)
    """
    frames_per_res = {}
    max_frame = -1
    for fname in contact_files:
        p = rep_dir / fname
        if not p.exists():
            logging.debug(f"Missing contact file {p}")
            continue
        df = parse_dat_file(p)
        if df.empty:
            continue
        frame_col = "Frame" if "Frame" in df.columns else None
        res_col = "Residue" if "Residue" in df.columns else None
        chain_col = "Chain" if "Chain" in df.columns else None

        if frame_col is None or res_col is None or chain_col is None:
            logging.debug(f"{p.name}: missing expected columns; columns={df.columns.tolist()}")
            continue

        for _, row in df.iterrows():
            if pd.isna(row[frame_col]) or pd.isna(row[res_col]) or pd.isna(row[chain_col]):
                continue
            f = int(row[frame_col])
            r = int(row[res_col])
            c = str(row[chain_col]).strip()
            key = (c, r)
            frames_per_res.setdefault(key, set()).add(f)
            if f > max_frame:
                max_frame = f
    return frames_per_res, (max_frame if max_frame >= 0 else None)

def compute_ligand_occupancy_from_frames(frames_per_res: dict, total_frames: int):
    recs = []
    for (chain, resid), frameset in sorted(frames_per_res.items(), key=lambda x: (x[0][0], x[0][1])):
        count = len(frameset)
        occ = (count / total_frames) * 100.0 if total_frames and total_frames > 0 else 0.0
        recs.append({"Chain": chain, "Residue": int(resid), "frames": count, "occupancy_pct": occ})
    return pd.DataFrame(recs)

def compute_pp_contacts_from_trajectory(topology: str, traj: str, keypairs, cutoff=4.0, subsample=1):
    if mda is None:
        raise RuntimeError("MDAnalysis not available. Install it or provide precomputed pp_contacts CSV via --pp-contacts.")
    from MDAnalysis.lib import distances as mda_distances

    u = mda.Universe(topology, traj)
    total_frames = len(u.trajectory)
    logging.info(f"Loaded trajectory: {topology} + {traj} with {total_frames} frames (subsample={subsample})")

    pair_index_sets = []
    for pair in keypairs:
        c1, r1, c2, r2 = pair
        sel1 = u.select_atoms(f"chainID {c1} and resid {r1} and name CA")
        sel2 = u.select_atoms(f"chainID {c2} and resid {r2} and name CA")
        if len(sel1) == 0:
            sel1 = u.select_atoms(f"chainID {c1} and resid {r1} and not name H*")
        if len(sel2) == 0:
            sel2 = u.select_atoms(f"chainID {c2} and resid {r2} and not name H*")
        if len(sel1) == 0 or len(sel2) == 0:
            logging.warning(f"Could not find atoms for pair {pair} in topology. Skipping pair.")
            pair_index_sets.append((pair, None, None))
            continue
        idx1 = sel1.indices.copy()
        idx2 = sel2.indices.copy()
        pair_index_sets.append((pair, idx1, idx2))

    pp_frames = {pair: set() for pair in keypairs}

    for ts in u.trajectory[::subsample]:
        frame = int(ts.frame)
        positions = u.atoms.positions
        box = ts.dimensions if hasattr(ts, "dimensions") else None
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
                pp_frames[pair].add(frame)

    return pp_frames, total_frames

# -------------------------
# compute_disruption_metrics (UPDATED with normalized metrics)
# -------------------------
def compute_disruption_metrics(keypairs, pp_frames, ligand_frames_per_pair, total_frames,
                               threshold=20.0, pp_distance_stats: Optional[dict]=None, min_baseline_freq: float=0.05):
    """
    Compute disruption metrics and include normalized columns:
      - disruption_pct_global (old)
      - disruption_pct_by_lig  (n_blocked / n_lig * 100)  <-- recommended
      - disruption_pct_by_testable (n_blocked / |frames_lig ∪ frames_pp| * 100)
      - disruption_pct_by_examined (if frames_examined_total available)
    Also includes overlap_pct_of_pp and blocked_pct_of_lig for convenience.
    """
    recs = []
    for pair in keypairs:
        frames_pp = pp_frames.get(pair, set())
        frames_lig = ligand_frames_per_pair.get((pair[0], pair[1]), set()).union(
            ligand_frames_per_pair.get((pair[2], pair[3]), set())
        )
        n_pp = len(frames_pp)
        n_lig = len(frames_lig)

        pp_occ = len(frames_pp) / total_frames * 100.0 if total_frames else 0.0
        lig_occ = len(frames_lig) / total_frames * 100.0 if total_frames else 0.0
        both = frames_pp & frames_lig
        n_both = len(both)
        simultaneous = n_both / total_frames * 100.0 if total_frames else 0.0

        # blocked frames (lig present & PP absent)
        blocked_by_ligand = frames_lig - frames_pp
        n_blocked = len(blocked_by_ligand)

        # original (global) pct
        disruption_pct_global = (n_blocked / total_frames * 100.0) if total_frames else 0.0

        # normalized metrics
        disruption_pct_by_lig = (n_blocked / n_lig * 100.0) if n_lig else float("nan")
        frames_testable = len(frames_lig.union(frames_pp))
        disruption_pct_by_testable = (n_blocked / frames_testable * 100.0) if frames_testable else float("nan")

        # distance stats if provided
        mm = None; mmstd = None; mmwhen = None; n_total = 0; frames_total = 0; apo_freq = 0.0
        if pp_distance_stats:
            dstats = pp_distance_stats.get(pair) or pp_distance_stats.get((pair[0], pair[1], pair[2], pair[3]))
            if dstats:
                mm = dstats.get('mean_min_dist_mean')
                mmstd = dstats.get('mean_min_dist_std')
                mmwhen = dstats.get('mean_min_dist_when_contact_mean')
                n_total = dstats.get('n_contact_frames_total', 0)
                frames_total = dstats.get('frames_examined_total', 0)
                apo_freq = dstats.get('apo_freq', (n_total / frames_total) if frames_total else 0.0)

        disruption_pct_by_examined = (n_blocked / frames_total * 100.0) if frames_total else float("nan")

        # relative metrics
        overlap_pct_of_pp = (n_both / n_pp * 100.0) if n_pp else float("nan")
        blocked_pct_of_lig = (n_blocked / n_lig * 100.0) if n_lig else float("nan")

        # determine status with baseline filter
        if apo_freq < float(min_baseline_freq):
            status = "insufficient_baseline"
        else:
            if disruption_pct_by_lig is not None and not math.isnan(disruption_pct_by_lig):
                check_val = disruption_pct_by_lig
            else:
                check_val = disruption_pct_global
            if check_val > threshold:
                status = "disrupted"
            elif check_val > (threshold / 2.0):
                status = "partial"
            else:
                status = "preserved"

        recs.append({
            "pair_label": f"{pair[0]}{pair[1]}-{pair[2]}{pair[3]}",
            "chain1": pair[0], "res1": pair[1],
            "chain2": pair[2], "res2": pair[3],
            "pp_occ_pct": pp_occ,
            "lig_occ_pct": lig_occ,
            "simultaneous_pct": simultaneous,
            # original and global
            "disruption_pct": disruption_pct_global,
            "disruption_pct_global": disruption_pct_global,
            # normalized metrics
            "disruption_pct_by_lig": disruption_pct_by_lig,
            "disruption_pct_by_testable": disruption_pct_by_testable,
            "disruption_pct_by_examined": disruption_pct_by_examined,
            # relative metrics
            "overlap_pct_of_pp": overlap_pct_of_pp,
            "blocked_pct_of_lig": blocked_pct_of_lig,
            # counts & distances
            "pp_frames_count": n_pp,
            "lig_frames_count": n_lig,
            "disruption_count": n_blocked,
            "mean_min_dist_mean": mm,
            "mean_min_dist_std": mmstd,
            "mean_min_dist_when_contact_mean": mmwhen,
            "n_contact_frames_total": n_total,
            "frames_examined_total": frames_total,
            "apo_freq": apo_freq,
            "status": status
        })
    return pd.DataFrame(recs)

# -------------------------
# plotting (updated to accept metric_col)
# -------------------------
def plot_heatmap_disruption(df: pd.DataFrame, outfile: Path, metric_col: str = "disruption_pct_by_lig"):
    """
    Create a vertical heatmap (rows = pair_label) showing chosen metric (metric_col),
    with a coloured side-bar indicating disruption 'status'.
    """
    import matplotlib as mpl
    from matplotlib.patches import Patch

    if df is None or df.shape[0] == 0:
        outfile.parent.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(4, 1))
        plt.text(0.5, 0.5, "No data", ha="center", va="center")
        plt.axis("off")
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()
        return

    df_sorted = df.sort_values(["status", metric_col if metric_col in df.columns else "disruption_pct"], ascending=[False, True]).reset_index(drop=True)

    if "pair_label" not in df_sorted.columns:
        df_sorted["pair_label"] = df_sorted.index.astype(str)
    if metric_col not in df_sorted.columns:
        # fallback
        if "disruption_pct" in df_sorted.columns:
            metric_col = "disruption_pct"
        else:
            raise ValueError("No suitable metric column found for plotting.")

    if "status" not in df_sorted.columns:
        df_sorted["status"] = "unknown"

    status_colors = {
        "disrupted": "#66aadadd",
        #"possible disruption": "#c6285629",
        "partial": "#fb8c00",
        "preserved": "#218204eb",
        "insufficient_baseline": "#c62828",
        "unknown": "#bdbdbd"
    }

    statuses = df_sorted["status"].astype(str).tolist()
    row_colors = []
    for s in statuses:
        col = status_colors.get(s) or status_colors.get(s.lower(), status_colors["unknown"])
        row_colors.append(mpl.colors.to_rgb(col)) 
    colors_arr = np.array(row_colors).reshape(len(row_colors), 1, 3) 

    mat = df_sorted.set_index("pair_label")[[metric_col]] 

    n_rows = len(df_sorted)
    fig = plt.figure(figsize=(6.5, max(3.5, 0.25 * n_rows)), constrained_layout=True)
    gs = fig.add_gridspec(nrows=1, ncols=3, width_ratios=[0.12, 0.5, 5.0], wspace=0.02) 

    ax_pad = fig.add_subplot(gs[0, 0])
    ax_pad.axis("off")

    ax_color = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[0, 2])

    ax_color.imshow(colors_arr, aspect='auto', origin='upper', interpolation='nearest')
    ax_color.set_xticks([]) # no x-axis
    ax_color.set_yticks([]) # no y-axis
    ax_color.invert_yaxis()

    cmap = sns.diverging_palette(240, 10, n=9, as_cmap=True) # blue to red
    sns.heatmap(mat, ax=ax_heat, cmap=cmap, annot=True, fmt=".1f",
                cbar_kws={"label": f"{metric_col} (%)"}, linewidths=0.5, linecolor="#eeeeee",
                yticklabels=True)

    ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0, fontsize=9)
    # ax_heat.set_xlabel(metric_col.replace("_", " ").capitalize())
    ax_heat.set_title("Per-pair disruption (metric: {})".format(metric_col.replace("_", " ")))

    unique_statuses = []
    for s in statuses:
        if s not in unique_statuses:
            unique_statuses.append(s)
    legend_items = []
    for s in unique_statuses:
        col = status_colors.get(s) or status_colors.get(s.lower(), status_colors["unknown"])
        legend_items.append(Patch(facecolor=col, edgecolor="#333333", label=str(s)))
    if legend_items:
        ax_color.legend(handles=legend_items, loc="lower center",
                        bbox_to_anchor=(0.7, -0.15), ncol=min(2, len(legend_items)),
                        frameon=False, fontsize=8)

    outfile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

# -------------------------
# Main CLI
# -------------------------
def main():
    p = argparse.ArgumentParser(description="Analyze PL-Contacts across replicates and compute disruption metrics.")
    p.add_argument("--root", type=Path, required=True, help="Top-level directory containing replicate subfolders")
    p.add_argument("--rep-frame-lens", help="Optional comma-separated frame lengths per replica (overrides autodetection)")
    p.add_argument("--replica-order", help="Optional explicit comma-separated replica directory names in the exact concatenation order")
    p.add_argument("--replica-dirs", type=str, default=None,
                help="Optional comma-separated replicate folder names (relative to --root) to analyze")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--pairs", type=str, default=",".join(f"{a}{b}-{c}{d}" for a,b,c,d in DEFAULT_KEYPAIRS),
                help="Comma-separated list of canonical pairs like E493-A35,E505-A353")
    p.add_argument("--contact-files", type=str, default=",".join(DEFAULT_CONTACT_FILES),
                help="Comma-separated list of PL-Contacts files to parse")
    p.add_argument("--traj", type=str, default=None, help="Trajectory file (e.g. merged .xtc) for PP contact calculation")
    p.add_argument("--top", type=str, default=None, help="Topology file (pdb / cms supported by MDAnalysis) for PP contact calculation")
    p.add_argument("--pp-contacts", type=str, default=None, help="Precomputed pp_contacts CSV (optional)")
    p.add_argument("--pp-distances", type=str, default=None, help="Per-rep pp distances CSV to include baseline stats (optional)")
    p.add_argument("--min-baseline-freq", type=float, default=0.05, help="Minimum apo baseline frequency (fraction) required before evaluating disruption")
    p.add_argument("--pp-cutoff", type=float, default=4.0, help="Distance cutoff (Å) for protein-protein contact")
    p.add_argument("--frames", type=int, default=None, help="Total frames (if you want to override detection)")
    p.add_argument("--disruption-threshold", type=float, default=20.0, help="Disruption pct threshold to call 'disrupted'")
    p.add_argument("--metric", type=str, default="disruption", choices=["disruption","delta","both"], help="Which metric to compute/use")
    p.add_argument("--delta-threshold", type=float, default=10.0, help="Threshold (percentage points) for delta metric (freq_bound - freq_apo)")
    p.add_argument("--plot-metric", type=str, default="disruption_pct_by_lig", help="Which metric column to plot on the heatmap")
    p.add_argument("--log-level", type=str, default="INFO")
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=str(args.outdir / "analysis.log"),
                        level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(asctime)s %(levelname)s: %(message)s")
    logging.info("Starting analysis")
    logging.info(f"Root = {args.root}")

    contact_files = [c.strip() for c in args.contact_files.split(",") if c.strip()]
    # parse pairs string into list of tuples
    pair_tokens = [tok.strip() for tok in args.pairs.split(",") if tok.strip()]
    keypairs = []
    for tok in pair_tokens:
        m = re.match(r"^([A-Za-z])(\d+)-([A-Za-z])(\d+)$", tok)
        if not m:
            logging.warning(f"Cannot parse pair token '{tok}' — skip")
            continue
        keypairs.append((m.group(1), int(m.group(2)), m.group(3), int(m.group(4))))
    if not keypairs:
        logging.error("No valid canonical pairs provided. Exiting.")
        sys.exit(1)

    explicit = [str(args.root / d) for d in args.replica_dirs.split(",")] if args.replica_dirs else None
    replicas = discover_replicas(args.root, explicit=explicit)
    logging.info(f"Found replicate directories: {replicas}")

    all_replicate_occ = []
    ligand_frames_per_pair_global = {}
    global_max_frame = -1

    for rep in replicas:
        rep = Path(rep)
        logging.info(f"Processing replicate {rep}")
        frames_per_res, max_frame = collect_ligand_frames_from_contact_files(rep, contact_files)
        if max_frame is not None and max_frame > global_max_frame:
            global_max_frame = max_frame
        if not frames_per_res:
            logging.warning(f"No ligand contact frames found in replicate {rep}")
        total_frames_rep = args.frames or (max_frame + 1 if max_frame is not None else None)
        if total_frames_rep is None:
            logging.info(f"No total_frames known yet for {rep}; will infer later")
        occ_df = compute_ligand_occupancy_from_frames(frames_per_res, total_frames_rep or 0)
        occ_df["replicate"] = rep.name
        out_csv = args.outdir / "csv" / f"{rep.name}_occupancy.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        occ_df.to_csv(out_csv, index=False)
        all_replicate_occ.append(occ_df)
        for k, s in frames_per_res.items():
            ligand_frames_per_pair_global.setdefault(k, set()).update(s)

    if args.frames:
        total_frames = args.frames
    else:
        if args.traj and args.top:
            if mda is None:
                logging.error("MDAnalysis not available to read trajectories; please install it or provide --frames")
                sys.exit(1)
            utest = mda.Universe(args.top, args.traj)
            total_frames = len(utest.trajectory)
        else:
            total_frames = (global_max_frame + 1) if global_max_frame >= 0 else None

    if total_frames is None:
        logging.error("Could not determine total_frames. Provide --frames or supply --top and --traj.")
        sys.exit(1)

    logging.info(f"Using total_frames = {total_frames}")

    if all_replicate_occ:
        agg = pd.concat(all_replicate_occ, ignore_index=True)
        agg_group = agg.groupby(["Chain", "Residue"]).agg(
            occ_mean=("occupancy_pct", "mean"),
            occ_sd=("occupancy_pct", "std"),
            replicates=("replicate", "nunique")
        ).reset_index()
        agg_group.to_csv(args.outdir / "csv" / "occupancy_aggregated.csv", index=False)
    else:
        logging.warning("No per-replicate occupancy data to aggregate")

    # load or compute pp_frames
    if args.pp_contacts:
        logging.info(f"Loading precomputed protein-protein contact CSV: {args.pp_contacts}")
        pp_contacts_path = Path(args.pp_contacts)
        root_dir = Path(args.root) if args.root else Path(".")
        replicas = discover_replicas(root_dir)
        rep_frame_counts = infer_rep_frame_counts(replicas, getattr(args, 'rep_frame_lens', None), getattr(args, 'replica_order', None))
        replica_names_order = [p.name for p in replicas] if not getattr(args,'replica_order', None) else [s.strip() for s in args.replica_order.split(",")]
        offsets = compute_offsets(rep_frame_counts, replica_names_order)
        pp_frames = load_and_map_pp_contacts(pp_contacts_path, offsets)
        out_mapped = Path(args.outdir) / f"{pp_contacts_path.stem}_mapped.csv"
        write_mapped_pp_csv(out_mapped, pp_frames)
        logging.info(f"Loaded and mapped pp contacts; stored mapping csv at: {out_mapped}")
    else:
        logging.info("No precomputed pp_contacts CSV provided; script will compute PP contacts from trajectories as before.")
        if args.top and args.traj:
            pp_frames, _ = compute_pp_contacts_from_trajectory(args.top, args.traj, keypairs, cutoff=args.pp_cutoff)
        else:
            pp_frames = {}

    pp_distance_stats = None
    if args.pp_distances:
        pp_distance_stats = aggregate_pp_distance_stats(Path(args.pp_distances))
        logging.info(f"Loaded pp distance stats for {len(pp_distance_stats)} pairs from {args.pp_distances}")

    # compute disruption metrics and write CSV
    disruption_df = compute_disruption_metrics(keypairs, pp_frames, ligand_frames_per_pair_global,
                                               total_frames, threshold=args.disruption_threshold,
                                               pp_distance_stats=pp_distance_stats,
                                               min_baseline_freq=args.min_baseline_freq)
    out_metrics = args.outdir / "csv" / "disruption_metrics.csv"
    out_metrics.parent.mkdir(parents=True, exist_ok=True)
    disruption_df.to_csv(out_metrics, index=False)

    # plot (choose metric_col from CLI)
    metric_col = args.plot_metric or "disruption_pct_by_lig"
    try:
        plot_heatmap_disruption(disruption_df, Path(args.outdir) / "plots" / "disruption_heatmap.png", metric_col=metric_col)
    except Exception as e:
        logging.exception("Failed to plot heatmap: %s", e)

    logging.info("Analysis complete. Outputs written to %s", args.outdir)

if __name__ == "__main__":
    main()
