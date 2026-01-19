#!/usr/bin/env python3
"""
pp-frameplot.py

Visualize per-frame protein-protein contacts for a list of key residue pairs.

Input:
 - CSV produced by the exporter with columns:
   replica, chain1, res1, chain2, res2, Frame
   (column names are case-insensitive: Frame or frame accepted)

usage:
# default: reads pp_contacts_perframe.csv in current dir, uses default keypairs
python pp-frameplot.py --pp-csv pp_contacts_perframe.csv --outdir analysis_output/plots

# supply custom pairs (comma-separated, same format as before)
python pp-frameplot.py --pp-csv pp_contacts_perframe.csv \
  --pairs "E417-A30,E449-A38,E493-A34,E505-A353,E495-A353" \
  --outdir analysis_output/plots

# show plots interactively as well as saving them (useful when running locally)
python pp-frameplot.py --pp-csv pp_contacts_perframe.csv --outdir analysis_output/plots --show

Output:
 - One PNG per pair in OUTDIR/plots named <pair_label>_contacts.png
 - A CSV summary at OUTDIR/pp_contacts_perpair_summary.csv
"""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap
import csv
import sys

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

def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to lower-case keys we expect."""
    mapping = {}
    cols = df.columns.tolist()
    for c in cols:
        cl = c.strip()
        low = cl.lower()
        if low == "frame":
            mapping[c] = "Frame"
        elif low == "replica":
            mapping[c] = "replica"
        elif low in ("chain1", "chain"):
            mapping[c] = "chain1"
        elif low in ("chain2",):
            mapping[c] = "chain2"
        elif low in ("res1", "residue1", "resid1", "resnum1"):
            mapping[c] = "res1"
        elif low in ("res2", "residue2", "resid2", "resnum2"):
            mapping[c] = "res2"
        elif low in ("res", "residue", "resid", "resnum"):
            # ambiguous; leave it to data
            mapping[c] = "res"
        else:
            mapping[c] = c
    return df.rename(columns=mapping)

def parse_pairs_string(pairs_str: str):
    pairs = []
    for tok in pairs_str.split(","):
        tok = tok.strip()
        if not tok:
            continue
        # expect format "E417-A30" or "E 417 A 30" or "E,417,A,30"
        tok2 = tok.replace(" ", "").replace(",", "-")
        if "-" in tok2:
            left,right = tok2.split("-",1)
            # left may be E417 -> chain letter(s) plus int
            import re
            m1 = re.match(r"^([A-Za-z]+)(\d+)$", left)
            m2 = re.match(r"^([A-Za-z]+)(\d+)$", right)
            if m1 and m2:
                pairs.append((m1.group(1), int(m1.group(2)), m2.group(1), int(m2.group(2))))
                continue
        # fallback: try space split
        parts = tok.split()
        if len(parts) >= 4:
            pairs.append((parts[0], int(parts[1]), parts[2], int(parts[3])))
            continue
        raise ValueError(f"Cannot parse pair token: {tok}")
    return pairs

def lifetimes_from_frame_list(frames_list):
    """Return list of contiguous-run lengths (counts of consecutive frames)."""
    if not frames_list:
        return []
    frames = sorted(frames_list)
    lengths = []
    cur = frames[0]
    length = 1
    for a, b in zip(frames, frames[1:]):
        if b == a + 1:
            length += 1
        else:
            lengths.append(length)
            length = 1
    lengths.append(length)
    return lengths

def plot_pair_raster(df_pair: pd.DataFrame, pair_label: str, outpath: Path, dpi=200,
                     figsize=(8, 2.5), marker='|', markersize=8, show=False):
    """
    df_pair: dataframe filtered for the pair, must contain 'replica' and 'Frame'
    """
    # get unique replicas in sorted order
    rep_names = sorted(df_pair['replica'].unique())
    if len(rep_names) == 0:
        print(f"[WARN] No frames for pair {pair_label} -> skipping plotting.")
        return None

    # build a mapping of replica -> frames list
    frames_per_rep = {rep: sorted(df_pair[df_pair['replica'] == rep]['Frame'].astype(int).tolist()) for rep in rep_names}
    # compute freq per rep
    counts = {rep: len(frames_per_rep[rep]) for rep in rep_names}
    # try to infer frames examined per replica from max frame encountered + 1
    frames_examined = {rep: (max(frames_per_rep[rep]) + 1 if frames_per_rep[rep] else 0) for rep in rep_names}

    # compute overall x-axis limits
    all_frames = [f for frames in frames_per_rep.values() for f in frames]
    xmin = 0 if not all_frames else min(all_frames)
    xmax = 100 if not all_frames else max(all_frames)

    fig, ax = plt.subplots(figsize=figsize)
    yticks = []
    ylabels = []
    for i, rep in enumerate(rep_names):
        frames = frames_per_rep[rep]
        y = [i] * len(frames)
        if frames:
            ax.plot(frames, y, linestyle='', marker=marker, markersize=markersize, color='#025682', alpha=0.7)
        yticks.append(i)
        # show label with count and freq if frames_examined > 0 else show count only
        fe = frames_examined.get(rep, 0)
        freq_text = f"{counts[rep]}/{fe}" if fe else f"{counts[rep]}"
        ylabels.append(f"{rep} ({freq_text})")

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Frame")
    ax.set_title(f"{pair_label} contact frames per replica")
    ax.set_ylim(-1, len(rep_names))
    # expand x-limits slightly
    if all_frames:
        ax.set_xlim(max(0, xmin - 1), xmax + 1)
    else:
        ax.set_xlim(0, 1)

    # annotate overall stats in the center
    total_contacts = sum(counts.values())
    total_frames_examined = sum(frames_examined.values()) if sum(frames_examined.values()) > 0 else None
    if total_frames_examined:
        overall_freq = total_contacts / total_frames_examined
        stats_text = f"Total contacts: {total_contacts}    Overall freq: {overall_freq:.3f}"
    else:
        stats_text = f"Total contacts: {total_contacts}"
    ax.text(0.01, 0.95, stats_text, transform=ax.transAxes, fontsize=8, va='top')

    plt.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=dpi)
    if show:
        plt.show()
    plt.close(fig)

    # return frames_per_rep and counts for summary
    return frames_per_rep, counts

def main():
    p = argparse.ArgumentParser(description="Visualize PP contact frames per key pair.")
    p.add_argument("--pp-csv", required=True, help="CSV with per-frame contacts (replica,chain1,res1,chain2,res2,Frame)")
    p.add_argument("--pairs", help="Comma-separated pairs e.g. E417-A30,E449-A38. Default: built-in keypairs")
    p.add_argument("--outdir", default="pp_contact_frames_output", help="Output directory for plots and summary CSV")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--figwidth", type=float, default=10.0)
    p.add_argument("--figheight", type=float, default=2.5)
    p.add_argument("--show", action="store_true", help="Show plots interactively")
    p.add_argument("--allow-swap", action="store_true", help="Also match pairs reversed (chain2-res2 vs chain1-res1)")
    args = p.parse_args()

    pp_csv = Path(args.pp_csv)
    if not pp_csv.exists():
        print("pp-csv not found:", pp_csv)
        sys.exit(1)

    df = pd.read_csv(pp_csv)
    df = canonicalize_columns(df)

    # normalize expected columns to lowercase keys in dataframe: 'replica','chain1','res1','chain2','res2','Frame'
    # ensure 'Frame' column exists (case normalized to 'Frame' by canonicalize_columns)
    if 'Frame' not in df.columns and 'frame' in df.columns:
        df = df.rename(columns={'frame': 'Frame'})

    # check we have suitable columns
    required = ['replica', 'chain1', 'res1', 'chain2', 'res2', 'Frame']
    if not all(c in df.columns for c in required):
        print("CSV missing required columns. Found columns:", df.columns.tolist())
        print("Expected at least:", required)
        sys.exit(1)

    # coerce numeric types
    df['res1'] = df['res1'].astype(int)
    df['res2'] = df['res2'].astype(int)
    df['Frame'] = df['Frame'].astype(int)
    df['replica'] = df['replica'].astype(str)

    if args.pairs:
        pairs = parse_pairs_string(args.pairs)
    else:
        pairs = DEFAULT_KEYPAIRS

    outdir = Path(args.outdir)
    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for pair in pairs:
        c1, r1, c2, r2 = pair
        pair_label = f"{c1}{r1}-{c2}{r2}"
        # filter for pair OR (if allow_swap) swapped pair
        mask = (df['chain1'] == str(c1)) & (df['res1'] == int(r1)) & (df['chain2'] == str(c2)) & (df['res2'] == int(r2))
        if args.allow_swap:
            mask_rev = (df['chain1'] == str(c2)) & (df['res1'] == int(r2)) & (df['chain2'] == str(c1)) & (df['res2'] == int(r1))
            mask = mask | mask_rev
        df_pair = df[mask].copy()
        out_png = plots_dir / f"{pair_label}_contacts.png"
        res = plot_pair_raster(df_pair, pair_label, out_png, dpi=args.dpi,
                               figsize=(args.figwidth, args.figheight), show=args.show)
        # compute some lifetimes and summary stats
        frames_per_rep = {}
        counts = {}
        median_lifetime = None
        mean_lifetime = None
        if res is not None:
            frames_per_rep, counts = res
            all_lifetimes = []
            for rep, frames in frames_per_rep.items():
                lifs = lifetimes_from_frame_list(frames)
                if lifs:
                    all_lifetimes.extend(lifs)
            if all_lifetimes:
                median_lifetime = float(np.median(all_lifetimes))
                mean_lifetime = float(np.mean(all_lifetimes))
        # total contacts and frequency estimates across replicas
        total_contacts = sum(counts.values()) if counts else 0
        total_frames_examined = 0
        # try infer per-rep frames_ex (best-effort: use max(frame)+1 per replica)
        per_rep_frames_examined = {}
        for rep in df['replica'].unique():
            rep_frames = df[df['replica'] == rep]['Frame']
            per_rep_frames_examined[rep] = int(rep_frames.max() + 1) if not rep_frames.empty else 0
            total_frames_examined += per_rep_frames_examined[rep]
        overall_freq = (total_contacts / total_frames_examined) if total_frames_examined else 0.0

        summary_rows.append({
            'pair': pair_label,
            'c1': c1, 'r1': r1, 'c2': c2, 'r2': r2,
            'total_contacts': total_contacts,
            'total_frames_examined': total_frames_examined,
            'overall_freq': overall_freq,
            'median_lifetime_frames': median_lifetime,
            'mean_lifetime_frames': mean_lifetime
        })
        print(f"Plotted {pair_label} -> {out_png} (contacts={total_contacts}, freq={overall_freq:.4f})")

    # write summary CSV
    summary_csv = outdir / "pp_contacts_perpair_summary.csv"
    with summary_csv.open("w", newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()) if summary_rows else ['pair'])
        w.writeheader()
        for r in summary_rows:
            w.writerow(r)
    print("Wrote per-pair summary to", summary_csv)
    print("Done.")

if __name__ == "__main__":
    main()
