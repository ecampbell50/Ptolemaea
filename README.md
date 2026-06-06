# Ptolemaea: Consensus Annotation of Antiviral Defence Systems

**Authors:** Emmet B. T. Campbell, Timofey Skvortsov, Sharon A. Huws, Christopher J. Creevey

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> *"In Ptolemaea, the final zone of the ninth circle of Hell, traitors to their guests are punished"* — much as bacteriophages, unwelcome guests in the cell, meet the bacterial antiviral arsenal.

## Overview

Ptolemaea annotates antiviral defence systems in bacterial genomes by running three
complementary methods over a single, common protein set and reconciling them into one
**consensus annotation per genome**:

- **[PADLOC](https://github.com/padlocbio/padloc)** — HMM-based defence-system detection
- **[DefenseFinder](https://github.com/mdmparis/defense-finder)** — model-based detection (MacSyFinder)
- **Bidirectional BLASTp** against curated *Bacillus cereus* defence proteins from [July & Gillis (2025)](https://doi.org/10.1038/s41598-025-86748-8)

It does not introduce new models or claim higher accuracy than its components. Its purpose is
to **maximise recall** and to make inter-tool disagreements **explicit, traceable, and resolvable**,
with an optional human-in-the-loop curation step.

## How it works

Each genome passes through five stages:

1. **Prokka** — predict proteins (`.faa`) and coordinates (`.gff`); contig headers are normalised automatically.
2. **PADLOC** — detect defence systems from the protein/GFF pair.
3. **DefenseFinder** — detect defence systems from the same proteins.
4. **Bidirectional BLASTp** — forward (genome → *B. cereus* DB) and reverse (*B. cereus* → genome) searches.
5. **Consensus** — collate all calls; assign each protein a consensus name and a **status label**.

A **master key** (`databases/MASTER_ToolKey.tsv`) harmonises the two tools' nomenclatures to a
shared consensus name (e.g. PADLOC `AbiD` + DefenseFinder `Abi2` → `Abi2DF`).

### Consensus status labels

| Status | Meaning | Resolution |
|---|---|---|
| **AGREE** | Both tools detect the protein and map to the same consensus name | Automatic |
| **RESOLVED** | Both tools detect but disagree; BLAST voting breaks the tie | Automatic |
| **SINGLE** | Only one tool detects it, and it has a consensus name | Automatic |
| **BLAST** | Neither tool detects it, but forward + reverse BLAST agree and pass size/coverage filters | Automatic |
| **MAPPING** | Detected, but the system is not in the master key (no consensus name) | **Needs curation** |
| **CONFLICT** | Both tools detect and map, but BLAST voting is tied | **Needs curation** |

(BLAST-only hits that fail the filters are discarded and excluded from the output.)

## Installation

### Prerequisites

| Tool | Version |
|---|---|
| [Prokka](https://github.com/tseemann/prokka) | 1.14+ |
| [PADLOC](https://github.com/padlocbio/padloc) | 2.0.0+ |
| [DefenseFinder](https://github.com/mdmparis/defense-finder) | 2.0.1+ |
| HMMER | 3.4+ |
| NCBI BLAST+ | 2.15.0+ |
| Python | 3.8+ (pandas) |

Prokka and PADLOC are expected in their own conda environments; DefenseFinder, HMMER and
BLAST+ must be on `PATH` (e.g. via pip / environment modules).

### Setup

```bash
git clone https://github.com/ecampbell50/Ptolemaea.git
cd Ptolemaea
chmod +x scripts/Ptolemaea.sh
mkdir -p genomes
```

Before the first run, edit the paths for your system:

| File | Set |
|---|---|
| `scripts/pipeline_functions.sh` | Prokka conda env (`run_prokka`), PADLOC conda env (`run_padloc`) |
| `scripts/Ptolemaea.sh` | `module load` lines (HMMER, BLAST+) |
| `scripts/defence_pipeline_array.sh` | `#SBATCH --partition=`, time/mem, `module load` lines *(HPC only)* |
| `scripts/batch_submit_defence_pipeline.sh` | `WORKING_DIR` *(HPC only — exported to the array jobs)* |

**Input:** one genome per file, nucleotide FASTA with a `.fna` extension (e.g. `1005041.3.fna`).
Place them in `genomes/`.

## Usage

The full workflow is: **run → extract unresolved → curate → build final outputs.**

### Quick start (single machine / few genomes)

```bash
mv /path/to/*.fna genomes/

# 1. Run the pipeline (per-genome consensus profiles in output/05_consensus/)
./scripts/Ptolemaea.sh /path/to/Ptolemaea

# 2. Pull out MAPPING/CONFLICT genes for curation
python3 scripts/extract_unresolved_patterns.py \
    --consensus-dir output/05_consensus/ \
    --output unresolved_patterns.csv

# 3. Edit unresolved_patterns.csv (fill TYPE/SUBTYPE/OUTCOME), save as *_CURATED.csv

# 4. Apply curation and build final matrix + tidy table + summary
python3 scripts/create_final_defence_matrix.py \
    --consensus-dir output/05_consensus/ \
    --resolutions unresolved_patterns_CURATED.csv \
    --output-prefix mydata
```

### Large-scale (HPC / SLURM array)

After setting `WORKING_DIR` and SLURM options (see [Setup](#setup)):

```bash
# 1. Submit array jobs (genomes run in parallel, 1000 per batch)
bash scripts/batch_submit_defence_pipeline.sh

# 2. When jobs finish (squeue -u $USER), run steps 2–4 from Quick start above.
```

> **No conflicts?** If step 2 reports a clean run, skip to step 4 and omit `--resolutions`.
> Uncurated MAPPING/CONFLICT genes are labelled `*_unresolved` so nothing is silently dropped.

## Curation

`extract_unresolved_patterns.py` groups all MAPPING/CONFLICT genes by their unique
(PADLOC, DefenseFinder, BLAST) annotation pattern, so you only decide each pattern **once**.
Open the CSV and fill three columns per row:

| Column | Example values |
|---|---|
| `TYPE` | `CBASS`, `BREX`, `RM` |
| `SUBTYPE` | `CBASS_IIs`, `BREX_I`, `RM_I` |
| `OUTCOME` | `Abi`, `Non-abi`, `Unknown` |

When a value genuinely cannot be determined, use `type_unresolved` / `subtype_unresolved` /
`outcome_unresolved` (distinct from `Unknown`, which means a tool detected it but could not classify it).

## Outputs

### Per genome — `output/`
```
01_prokka/<id>/<id>.faa, <id>.gff
02_padloc/<id>/<id>_padloc.csv
03_defensefinder/<id>/<id>_defense_finder_genes.tsv
04_blast/<id>/<id>_vs_bcereus_{forward,reverse}_cleaned.txt
05_consensus/<id>_defenceprofile.csv      # consensus profile (one row per defence gene)
```

### Final aggregated outputs — `create_final_defence_matrix.py`
- **`<prefix>_matrix.csv`** — wide genome × system matrix, **gene counts** (the unit is genes, not
  collapsed multi-gene systems). Columns encode all three levels of annotation as
  `<type>#<subtype>#<outcome>`. Use `--binary` for presence/absence; `--level subtype` or
  `--level type` for single-level columns.
- **`<prefix>_annotations.csv`** — tidy table, one row per defence gene (incl. `status`, `resolution_source`, final TYPE/SUBTYPE/OUTCOME).
- **`<prefix>_summary.tsv`** — per-genome counts (defence genes, unique systems, status breakdown).

Each matrix column header packs all three annotation levels, separated by `#`:

```
            <type> # <subtype> # <outcome>
   e.g.       RM   #   RM_I     #  Non-abi
```

```
# <prefix>_matrix.csv (gene counts; column headers are type#subtype#outcome)
genome_id,RM#RM_I#Non-abi,CBASS#CBASS_IIs#Abi,BREX#BREX_I#Non-abi,DRT#DRT_Type6#Unknown
1005041.3,2,1,1,0
1214195.3,1,0,1,1
```

## Configuration

| Setting | File |
|---|---|
| SLURM time / memory / CPUs | `defence_pipeline_array.sh` (`#SBATCH` lines) |
| Genomes per batch (default 1000) | `batch_submit_defence_pipeline.sh` (`BATCH_SIZE`) |
| BLAST E-value / hit count | `pipeline_functions.sh` (`-evalue 1e-5`, `-max_target_seqs 1`) |

**Resources (per genome, defaults):** ~30 min, 16 GB RAM, 6 cores. Steps are idempotent —
completed stages are skipped on re-run, so failed jobs can be resubmitted safely.

## Troubleshooting

- **No genomes processed** — check files end in `.fna` and `output/genome_list.txt` was generated.
- **PADLOC fails** — GFF must be Prokka-style (Bakta/Prodigal GFFs untested); check the PADLOC database matches the PADLOC version.
- **DefenseFinder errors** — ensure HMMER is available and DefenseFinder models are downloaded.
- **BLAST database not found** — the `databases/bcereus_db/` files must be present; custom databases must use names consistent with `MASTER_ToolKey.tsv`.

## Citation

Ptolemaea is a wrapper: if you use it, please cite **this tool** and **all of the underlying
tools and data** it depends on.

**Ptolemaea**
```
Campbell, E.B.T., Skvortsov, T., Huws, S.A., Creevey, C.J. (2025).
Ptolemaea: consensus annotation of antiviral defence systems in bacterial genomes.
https://github.com/ecampbell50/Ptolemaea
```

**Underlying tools and data**

- **Prokka** — Seemann T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30(14):2068–2069. doi:[10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153)
- **PADLOC** — Payne L.J. *et al.* (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. *Nucleic Acids Research* 49(19):10868–10878. doi:[10.1093/nar/gkab883](https://doi.org/10.1093/nar/gkab883)
- **PADLOC web server** — Payne L.J. *et al.* (2022). PADLOC: a web server for the identification of antiviral defence systems in microbial genomes. *Nucleic Acids Research* 50(W1):W541–W550. doi:[10.1093/nar/gkac400](https://doi.org/10.1093/nar/gkac400)
- **DefenseFinder** — Tesson F. *et al.* (2022). Systematic and quantitative view of the antiviral arsenal of prokaryotes. *Nature Communications* 13:2561. doi:[10.1038/s41467-022-30269-9](https://doi.org/10.1038/s41467-022-30269-9)
- **DefenseFinder models / webservice** — Tesson F. *et al.* (2024). A comprehensive resource for exploring antiphage defense: DefenseFinder webservice, wiki and databases. *Peer Community Journal* 4:e91. doi:[10.24072/pcjournal.470](https://doi.org/10.24072/pcjournal.470)
- **MacSyFinder v2** (DefenseFinder engine) — Néron B. *et al.* (2023). MacSyFinder v2: improved modelling and search engine to identify molecular systems in genomes. *Peer Community Journal* 3:e28. doi:[10.24072/pcjournal.250](https://doi.org/10.24072/pcjournal.250)
- **BLAST+** — Camacho C. *et al.* (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. doi:[10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **HMMER** — Eddy S.R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology* 7(10):e1002195. doi:[10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)
- **\*B. cereus\* defence database** — July E. & Gillis A. (2025). Antiviral defence arsenal across members of the *Bacillus cereus* group. *Scientific Reports* 15:4958. doi:[10.1038/s41598-025-86748-8](https://doi.org/10.1038/s41598-025-86748-8)

## License

MIT License — see [LICENSE](LICENSE).

## Acknowledgements

July & Gillis for the *B. cereus* defence-protein set and naming conventions; the PADLOC and
DefenseFinder teams; and the Queen's University Belfast HPC team.
