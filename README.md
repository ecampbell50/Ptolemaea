# Ptolemaea

**Consensus, comprehensive annotation of antiviral defence systems in bacterial genomes**

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.06.26.734901-b31b1b.svg)](https://www.biorxiv.org/content/10.64898/2026.06.26.734901v1)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> *"In Ptolemaea, the final zone of the ninth circle of Hell, traitors to their guests are punished"* — just like how bacteriophages, guests of the bacterial cell, are mistreated by the antiviral defence arsenal

Ptolemaea runs [PADLOC](https://github.com/padlocbio/padloc),
[DefenseFinder](https://github.com/mdmparis/defense-finder) and a bidirectional BLASTp
against the *B. cereus* defence proteins of
[July & Gillis (2025)](https://doi.org/10.1038/s41598-025-86748-8) over one shared protein
set, then reconciles them into a single consensus annotation per genome.

## Abstract

**Motivation:** Bacteria carry a large repertoire of antiviral defence systems, our knowledge
of which is expanding rapidly. Several bioinformatics tools now exist to identify them. Though
powerful, these tools can differ in the models they use and the nomenclature they return, thus
a single tool could both miss an annotation and disagree with its peers.

**Results:** Here we describe Ptolemaea, a pipeline for harmonising phage-defence annotations
across multiple tools by reconciling PADLOC, DefenseFinder, and a bidirectional BLAST. Over a
common predicted set of proteins, Ptolemaea provides a consensus annotation list per genome.
The pipeline is not intended to outperform or replace its component tools; its purpose is to
maximise the number of defence systems recovered from a genome and to make disagreements
between tools explicit and resolvable. We demonstrate the pipeline on 700 complete genomes
spanning the ESKAPE pathogens and *Escherichia coli*, recovering 32,509 defence annotations,
of which 50.6% were supported by more than one annotation source.

Preprint: [bioRxiv 10.64898/2026.06.26.734901](https://www.biorxiv.org/content/10.64898/2026.06.26.734901v1)

## Install

Every tool runs from an Apptainer image, so you only need
[Apptainer](https://apptainer.org/) (Linux or HPC) and internet access for the one-time setup.

```bash
git clone https://github.com/ecampbell50/Ptolemaea.git
cd Ptolemaea
bash singularity_scripts/setup.sh
```

`setup.sh` pulls the tool images into `images/` and the PADLOC and DefenseFinder databases
into `databases_runtime/`, pinned to the versions used in the paper. On an HPC, run it on a
login node.

## Usage

**1. Run.** Put one nucleotide FASTA per genome (`<id>.fna`) in `genomes/`, then:

```bash
bash singularity_scripts/Ptolemaea_singularity.sh .
```

Each genome goes through Pyrodigal (gene calling), PADLOC, DefenseFinder and bidirectional
BLASTp, and every defence protein gets a consensus name and a status. A master key
(`databases/MASTER_ToolKey.tsv`) maps the two tools' nomenclatures onto shared names.

| Status | Meaning |
|---|---|
| `AGREE` | Both tools detect it and map to the same name |
| `RESOLVED` | Both tools detect it but disagree; BLAST breaks the tie |
| `SINGLE` | Only one tool detects it |
| `BLAST` | Neither tool detects it, but forward and reverse BLAST agree |
| `MAPPING` | Detected, but the system is not in the master key: **needs curation** |
| `CONFLICT` | Both tools map it, but the BLAST vote is tied: **needs curation** |

**2. Curate (optional).** Collect the `MAPPING`/`CONFLICT` genes, grouped so each unique
pattern is decided once:

```bash
python3 scripts/extract_unresolved_patterns.py --consensus-dir output/05_consensus/ \
    --output unresolved_patterns.csv
```

Fill in `TYPE`, `SUBTYPE` and `OUTCOME` (e.g. `RM`, `RM_I`, `Non-abi`) and save as
`unresolved_patterns_CURATED.csv`. Skip this step and those genes are kept as `*_unresolved`.

**3. Build the final tables.**

```bash
python3 scripts/create_final_defence_matrix.py --consensus-dir output/05_consensus/ \
    --resolutions unresolved_patterns_CURATED.csv --output-prefix mydata
```

This writes `mydata_matrix.csv` (genome × system gene counts, columns
`<type>#<subtype>#<outcome>`; add `--binary` for presence/absence),
`mydata_annotations.csv` (one row per defence gene) and `mydata_summary.tsv` (per-genome
counts). The Python steps need pandas; on a cluster without it, prefix them with
`apptainer exec images/pandas.sif`.

Paths, threads (`PTOL_CPUS`) and tool versions are set in
`singularity_scripts/ptolemaea.config`. If your genomes are outside the repo on an HPC, bind
that filesystem, e.g. `export PTOL_BIND=/mnt/scratch2`. The original non-container scripts
(Prokka, conda, SLURM) are kept in `scripts/`.

## Example run

Download 100 complete *Bacillus cereus* group genomes (23 species) from NCBI, run the full
pipeline, and plot a summary figure:

```bash
bash examples/bcereus_100/run_example.sh --slurm     # HPC: one SLURM job per genome
bash examples/bcereus_100/run_example.sh --jobs 4    # single machine, 4 genomes at a time
bash examples/bcereus_100/run_example.sh --limit 5   # quick test on 5 genomes
```

Results go to `examples/bcereus_100/run/`, including the matrix (`bcereus100_matrix.csv`) and
the figure (`bcereus100_defence_overview.png`). The genome list is fixed in
`examples/bcereus_100/accessions.tsv`, so the run is reproducible. SLURM options (partition,
time, memory) can be set with `PTOL_SBATCH_ARGS`, e.g.
`export PTOL_SBATCH_ARGS="--partition=k2-hipri --time=00:30:00 --mem=16G --cpus-per-task=8"`.

## Citation

> Campbell E.B.T., Skvortsov T., Creevey C.J. (2026). Ptolemaea: consensus, comprehensive
> annotation of antiviral defence systems in bacterial genomes. *bioRxiv*.
> doi:[10.64898/2026.06.26.734901](https://doi.org/10.64898/2026.06.26.734901)

Please also cite the tools Ptolemaea wraps:
Pyrodigal ([Larralde 2022](https://doi.org/10.21105/joss.04296)),
Prodigal ([Hyatt *et al.* 2010](https://doi.org/10.1186/1471-2105-11-119)),
PADLOC ([Payne *et al.* 2021](https://doi.org/10.1093/nar/gkab883); [2022](https://doi.org/10.1093/nar/gkac400)),
DefenseFinder ([Tesson *et al.* 2022](https://doi.org/10.1038/s41467-022-30269-9); [2024](https://doi.org/10.24072/pcjournal.470)),
MacSyFinder ([Néron *et al.* 2023](https://doi.org/10.24072/pcjournal.250)),
BLAST+ ([Camacho *et al.* 2009](https://doi.org/10.1186/1471-2105-10-421)),
HMMER ([Eddy 2011](https://doi.org/10.1371/journal.pcbi.1002195)) and the
*B. cereus* defence database ([July & Gillis 2025](https://doi.org/10.1038/s41598-025-86748-8)).

## Acknowledgements

July & Gillis for the *B. cereus* defence-protein set and naming conventions; the PADLOC and
DefenseFinder teams; and the Queen's University Belfast HPC team.

## License

MIT, see [LICENSE](LICENSE).
