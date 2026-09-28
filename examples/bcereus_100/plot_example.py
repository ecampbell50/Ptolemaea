#!/usr/bin/env python3
"""
plot_example.py
Two-panel summary figure for the B. cereus group example run.

  A - Prevalence of each defence-system type: % of genomes carrying >= 1 gene
      of that type (automatically resolved calls only; genes left for curation
      are shown in panel B instead).
  B - How each defence gene was called, per species: share of genes by consensus
      status (AGREE / RESOLVED / SINGLE / BLAST resolve automatically; MAPPING /
      CONFLICT are left for curation).

Runs inside the seaborn image (matplotlib + pandas); called by run_example.sh.

Usage:
    python3 plot_example.py --annotations bcereus100_annotations.csv \
        --summary bcereus100_summary.tsv --accessions accessions.tsv \
        --output bcereus100_defence_overview
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Consensus statuses in pipeline order, one fixed colour each
STATUS_ORDER = ["AGREE", "RESOLVED", "SINGLE", "BLAST", "MAPPING", "CONFLICT"]
STATUS_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
BAR_COLOUR = "#2a78d6"
TEXT = "#0b0b0b"
TEXT_MUTED = "#52514e"
GRID = "#e4e3df"

UNRESOLVED_TYPES = {"type_unresolved", "UNMAPPED_TYPE"}
TOP_N_TYPES = 25
MIN_GENOMES_PER_SPECIES = 5


def style_axis(ax):
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=TEXT_MUTED, length=0, labelsize=8)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--accessions", required=True)
    parser.add_argument("--output", required=True, help="output path without extension")
    args = parser.parse_args()

    genes = pd.read_csv(args.annotations)
    summary = pd.read_csv(args.summary, sep="\t", dtype={"genome_id": str})
    species = pd.read_csv(args.accessions, sep="\t", comment="#")
    species = species.set_index("accession")["species"]
    n_genomes = len(summary)

    # --- Panel A: prevalence of each type across genomes ---------------------
    resolved = genes[~genes["final_type"].isin(UNRESOLVED_TYPES)]
    prevalence = (resolved.groupby("final_type")["genome_id"].nunique()
                  .div(n_genomes).mul(100)
                  .sort_values(ascending=False).head(TOP_N_TYPES)
                  .sort_values())

    # --- Panel B: status mix per species -------------------------------------
    genes["species"] = genes["genome_id"].map(species).fillna("Other")
    genomes_per_species = summary["genome_id"].map(species).value_counts()
    keep = genomes_per_species[genomes_per_species >= MIN_GENOMES_PER_SPECIES].index
    genes.loc[~genes["species"].isin(keep), "species"] = "Other species"
    status = (genes.groupby(["species", "status"]).size().unstack(fill_value=0)
              .reindex(columns=STATUS_ORDER, fill_value=0))
    status_pct = status.div(status.sum(axis=1), axis=0).mul(100)
    order = [s for s in status_pct.sort_values("AGREE").index if s != "Other species"]
    if "Other species" in status_pct.index:
        order = ["Other species"] + order
    status_pct = status_pct.loc[order]

    # --- Draw ------------------------------------------------------------------
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "text.color": TEXT, "axes.labelcolor": TEXT_MUTED})
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(12, 7.5),
        gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.55})

    ax_a.barh(prevalence.index, prevalence.values, color=BAR_COLOUR, height=0.7)
    ax_a.set_xlim(0, 100)
    ax_a.set_xlabel("Genomes carrying the system (%)")
    ax_a.set_title("A  Defence-system prevalence", loc="left", fontweight="bold", pad=10)
    style_axis(ax_a)

    left = pd.Series(0.0, index=status_pct.index)
    for name, colour in zip(STATUS_ORDER, STATUS_COLOURS):
        ax_b.barh(status_pct.index, status_pct[name], left=left, color=colour,
                  height=0.7, label=name, edgecolor="white", linewidth=1)
        left += status_pct[name]
    ax_b.set_xlim(0, 100)
    ax_b.set_xlabel("Defence genes (%)")
    ax_b.set_title("B  How each defence gene was called", loc="left",
                   fontweight="bold", pad=10)
    labels = [f"{s.replace('Bacillus ', 'B. ')}  (n={genomes_per_species.get(s, '')})"
              if s != "Other species" else "Other species" for s in status_pct.index]
    ax_b.set_yticks(range(len(labels)), labels)
    style_axis(ax_b)
    ax_b.legend(ncol=6, frameon=False, fontsize=7.5, loc="upper left",
                bbox_to_anchor=(-0.02, -0.08), handlelength=1, columnspacing=1,
                labelcolor=TEXT_MUTED, alignment="left",
                title="Consensus status (MAPPING and CONFLICT need curation)",
                title_fontsize=7.5)

    n_genes = len(genes)
    fig.suptitle(f"Ptolemaea on {n_genomes} Bacillus cereus group genomes "
                 f"({n_genes:,} defence genes)", x=0.06, ha="left",
                 fontsize=12, fontweight="bold")

    for ext in ("png", "pdf"):
        fig.savefig(f"{args.output}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Wrote {args.output}.png and {args.output}.pdf")


if __name__ == "__main__":
    main()
