#!/usr/bin/env python3
"""
plot_defence_heatmap.py
Presence/absence heatmap of defence systems (columns) across genomes (rows).

Input is EITHER:
  --consensus-dir  the per-genome profiles (output/05_consensus/), optionally with
                   --resolutions (your curated unresolved_patterns_CURATED.csv), OR
  --matrix         a matrix from create_final_defence_matrix.py (e.g. one you have
                   renamed/edited yourself).

Systems are shown at type level by default (e.g. RM, CBASS); use --level subtype
for finer columns. Genes left *_unresolved are dropped unless --keep-unresolved.
Columns are ordered by prevalence, rows by number of systems.

Needs pandas + matplotlib (the seaborn image has both).

Usage:
    python3 scripts/plot_defence_heatmap.py --consensus-dir output/05_consensus/ \
        --output defence_heatmap
    python3 scripts/plot_defence_heatmap.py --matrix mydata_matrix.csv \
        --output defence_heatmap
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import pandas as pd

# Reuse the exact loading/curation logic of the final-matrix step
sys.path.insert(0, str(Path(__file__).resolve().parent))
import create_final_defence_matrix as cfm  # noqa: E402

ABSENT = "#f0efec"
PRESENT = "#2a78d6"
TEXT = "#0b0b0b"
TEXT_MUTED = "#52514e"
UNRESOLVED = {"type_unresolved", "subtype_unresolved", "UNMAPPED_TYPE", "nan"}
MAX_ROW_LABELS = 80  # hide genome names above this many genomes


def matrix_from_consensus(consensus_dir, resolutions_file, level):
    resolutions = cfm.load_resolutions(Path(resolutions_file)) if resolutions_file else {}
    genes_df, all_genomes = cfm.load_profiles(Path(consensus_dir))
    if not genes_df.empty:
        genes_df = cfm.apply_resolutions(genes_df, resolutions)
    return cfm.build_matrix(genes_df, all_genomes, level=level, binary=True)


def matrix_from_file(matrix_file, level):
    """Collapse <type>#<subtype>#<outcome> columns to the requested level."""
    m = pd.read_csv(matrix_file, index_col=0)
    if all("#" in str(c) for c in m.columns):
        idx = 0 if level == "type" else 1
        m = m.T.groupby(lambda c: str(c).split("#")[idx]).sum().T
    return (m > 0).astype(int)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--consensus-dir", help="directory of *_defenceprofile.csv files")
    src.add_argument("--matrix", help="matrix CSV from create_final_defence_matrix.py")
    parser.add_argument("--resolutions", help="curated CSV (with --consensus-dir only)")
    parser.add_argument("--level", choices=["type", "subtype"], default="type")
    parser.add_argument("--keep-unresolved", action="store_true",
                        help="keep *_unresolved / unmapped columns")
    parser.add_argument("--output", default="defence_heatmap",
                        help="output path without extension (default: defence_heatmap)")
    args = parser.parse_args()

    if args.matrix:
        m = matrix_from_file(args.matrix, args.level)
    else:
        m = matrix_from_consensus(args.consensus_dir, args.resolutions, args.level)

    if not args.keep_unresolved:
        m = m[[c for c in m.columns if str(c) not in UNRESOLVED]]
    m = m.loc[:, m.sum() > 0]
    if m.empty:
        sys.exit("No defence systems to plot.")

    # Most common systems on the left; genomes with most systems at the top
    m = m[m.sum().sort_values(ascending=False).index]
    m = m.loc[m.sum(axis=1).sort_values(ascending=False, kind="stable").index]

    n_genomes, n_systems = m.shape
    show_rows = n_genomes <= MAX_ROW_LABELS
    width = min(max(6, 1.5 + 0.22 * n_systems), 30)
    height = min(max(3, 1.8 + (0.16 if show_rows else 0.04) * n_genomes), 30)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(m.values, aspect="auto", interpolation="none",
              cmap=ListedColormap([ABSENT, PRESENT]), vmin=0, vmax=1)

    # White gridlines between cells
    ax.set_xticks([x - 0.5 for x in range(1, n_systems)], minor=True)
    if show_rows:  # row lines only when rows are tall enough to see them
        ax.set_yticks([y - 0.5 for y in range(1, n_genomes)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="both", length=0, colors=TEXT_MUTED)

    ax.xaxis.tick_top()
    ax.set_xticks(range(n_systems), m.columns, rotation=90)
    if show_rows:
        ax.set_yticks(range(n_genomes), m.index)
    else:
        ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)

    ax.set_ylabel(f"Genomes (n={n_genomes})", color=TEXT_MUTED)
    ax.set_title(f"Defence systems ({args.level} level): presence / absence",
                 loc="left", color=TEXT, fontweight="bold", pad=12)
    ax.legend(handles=[Patch(color=PRESENT, label="Present"),
                       Patch(color=ABSENT, label="Absent")],
              loc="upper left", bbox_to_anchor=(0, -0.01), ncol=2, frameon=False,
              labelcolor=TEXT_MUTED, handlelength=1)

    for ext in ("png", "pdf"):
        fig.savefig(f"{args.output}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    print(f"Heatmap of {n_genomes} genomes x {n_systems} systems -> "
          f"{args.output}.png / .pdf")


if __name__ == "__main__":
    main()
