# Ptolemaea: Antiviral Defence System Consolidation in Bacteria

**Authors:** Emmet B. T. Campbell, Timofey Skvortsov, Sharon A. Huws, Christopher J. Creevey

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> *"In Ptolemaea, the fourth and final zone of the ninth circle of Hell, traitors to their guests are punished..."*  
> Much like how bacteriophages, as unwelcome guests in the bacterial cell, face the sophisticated arsenal of antiviral defense systems.

## Overview

Ptolemaea is a comprehensive pipeline for detecting and consolidating antiviral defence systems in bacterial genomes. It achieves consensus naming conventions between two popular antiviral defence detection tools:
- **[PADLOC](https://github.com/padlocbio/padloc)** 
- **[DefenseFinder](https://github.com/mdmparis/defense-finder)**

The pipeline integrates multiple detection methods and compares results against a curated database of *Bacillus cereus* defence genes from [July & Gillis (2025)](https://doi.org/10.1038/s41598-025-86748-8).

**📝 Citation:** If you use Ptolemaea, please cite this tool, PADLOC, DefenseFinder, and the July & Gillis paper.

# Getting started - if only running on a few genomes...
To use Ptolemaea, you must have access to prerequisites (see section for this below).

We hope to eventually make all required packages available in one clean download

### Once you have access to prerequisites...
- Clone the repo
```bash
git clone https://github.com/ecampbell50/Ptolemaea.git
```
- cd into Ptolemaea and create a "genomes" directory
```bash
cd Ptolemaea
chmod +x scripts/Ptolemaea.sh

# NB: ensure dir is called 'genomes' unless willing to edit pipeline scripts
mkdir genomes

# Then move any genomes to analyse into it
mv /path/to/genomes/*.fna /path/to/Ptolemaea/genomes/
# NB: script currently does not have extension handling
# Please ensure genome fastas have the extension '.fna' (not .fa/.fasta/etc)
```
You MUST change these lines to suit your system 
```bash
pipeline_functions.sh:
line 39:     conda run -p /path/to/prokka/env prokka \
line 86:     conda run -p /path/to/padloc/env padloc \

Ptolemaea.sh
line 6:      module load apps/hmmer/3.4/gcc-14.1.0
line 7:      module load apps/ncbiblast/2.15.0/gcc-14.1.0
```
- Now you can run the pipeline
```bash
./Ptolemaea.sh /path/to/Ptolemaea   # <- This sets your working directory
```
And you should then have your consensus output in Ptolemaea/output/05_consensus !

If running on lots of genomes I reccommend using the defence_pipeline_array.sh and batch_subit_defence_pipeline.sh scripts.
Mapping issues for all genomes can be rectified using extract_unresolved_patterns.py (see: "Post-Processing: Resolving Conflicts" below)


## Features

- **Automated genome annotation** using Prokka
- **Dual defence system detection** with PADLOC and DefenseFinder
- **Bidirectional BLAST analysis** against curated *B. cereus* defence genes
- **Consensus profile generation** combining all detection methods
- **High-throughput processing** with SLURM array job support
- **Batch submission** for processing thousands of genomes

## Prerequisites

### Software Dependencies

1. **Annotation:**
   - [Prokka](https://github.com/tseemann/prokka) v1.14+

2. **Defence Detection:**
   - [PADLOC](https://github.com/padlocbio/padloc) v2.0.0+ 
   - [DefenseFinder](https://github.com/mdmparis/defense-finder) v2.0.1+
   - HMMER v3.4+

3. **Sequence Analysis:**
   - NCBI BLAST+ v2.15.0+
   
4. **Scripting:**
   - Python 3.8+
   - Bash

### HPC Requirements

The pipeline is designed for SLURM-based HPC systems. Adjust partition names and resource allocations in the scripts according to your cluster.

## Directory Structure
```
Ptolemaea/
├── scripts/
│   ├── batch_submit_defence_pipeline.sh    # Main submission script
│   ├── defence_pipeline_array.slurm        # SLURM array job script
│   ├── pipeline_functions.sh               # Core pipeline functions
│   └── create_defence_profile_direct.py    # Consensus generation
├── databases/
│   ├── MASTER_ToolKey.tsv                  # Defence system naming key
│   └── bcereus_db/
│       ├── Bcereus_ConsensusDefProts_17Sep25     # BLAST database
│       └── Bcereus_DefenceProts_17Sep25.faa      # Defence protein sequences
├── genomes/                                 # Input genome files (.fna)
└── output/
    ├── 01_prokka/                          # Prokka annotations
    ├── 02_padloc/                          # PADLOC results
    ├── 03_defensefinder/                   # DefenseFinder results
    ├── 04_blast/                           # BLAST results
    └── 05_consensus/                       # Final consensus profiles
```

## Usage

### Large-scale analysis

1. **Place your genome files** in the `genomes/` directory
   - Format: `genomeID.fna` (e.g., `1005041.3.fna`)

2. **Configure the pipeline** by editing paths in `batch_submit_defence_pipeline.sh`:
```bash
   WORKING_DIR="/path/to/your/Ptolemaea"
```

3. **Submit the pipeline:**
```bash
   bash scripts/batch_submit_defence_pipeline.sh
```

The script will:
- Count your genomes
- Calculate required batches (1000 genomes per batch)
- Ask for confirmation
- Submit SLURM array jobs for each batch

### Pipeline Stages

#### Stage 1: Genome Annotation (Prokka)
##### NB: contig headers are normalised to ensure compatibility with PADLOC's gene ID matching
```bash
prokka --outdir output/01_prokka/${GENOME_ID} \
       --prefix ${GENOME_ID} \
       --cpus ${CPUS} \
       ${GENOME_FILE}
```

**Output:** 
- `${GENOME_ID}.faa` - protein sequences
- `${GENOME_ID}.gff` - genome annotations

---

#### Stage 2: PADLOC Defence Detection
```bash
padloc --faa ${GENOME_ID}.faa \
       --gff ${GENOME_ID}.gff \
       --outdir output/02_padloc/${GENOME_ID}
```

**Output:** `${GENOME_ID}_padloc.csv`

---

#### Stage 3: DefenseFinder Analysis
```bash
defense-finder run --out-dir output/03_defensefinder/${GENOME_ID} \
                   ${GENOME_ID}.faa
```

**Output:** `${GENOME_ID}_defense_finder_genes.tsv`

---

#### Stage 4: Bidirectional BLAST

**Forward BLAST** (your genome vs *B. cereus* database):
```bash
blastp -query ${GENOME_ID}.faa \
       -db Bcereus_ConsensusDefProts_17Sep25 \
       -out ${GENOME_ID}_vs_bcereus_forward.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
       -evalue 1e-5 \
       -max_target_seqs 1 \
       -num_threads ${CPUS}
```

**Reverse BLAST** (*B. cereus* vs your genome):
```bash
makeblastdb -in ${GENOME_ID}.faa -dbtype prot -out ${GENOME_ID}_PROTdb

blastp -query Bcereus_DefenceProts_17Sep25.faa \
       -db ${GENOME_ID}_PROTdb \
       -out ${GENOME_ID}_vs_bcereus_reverse.txt \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
       -evalue 1e-5 \
       -max_target_seqs 1 \
       -num_threads ${CPUS}
```

**Output:** 
- `${GENOME_ID}_vs_bcereus_forward_cleaned.txt`
- `${GENOME_ID}_vs_bcereus_reverse_cleaned.txt`

---

#### Stage 5: Consensus Profile Generation

Combines PADLOC, DefenseFinder, and BLAST results:
```bash
python3 create_defence_profile_direct.py \
    --padloc ${GENOME_ID}_padloc.csv \
    --defensefinder ${GENOME_ID}_defense_finder_genes.tsv \
    --forward ${GENOME_ID}_vs_bcereus_forward_cleaned.txt \
    --reverse ${GENOME_ID}_vs_bcereus_reverse_cleaned.txt \
    --master-key MASTER_ToolKey.tsv \
    --output ${GENOME_ID}_defenceprofile.csv
```

**Output:** `${GENOME_ID}_defenceprofile.csv` - final consolidated defence profile

---

---

### Post-Processing: Resolving Conflicts

After generating consensus profiles, some genes may have conflicting annotations between tools. Ptolemaea provides a script to identify and resolve these conflicts through manual curation.

#### Step 1: Extract Unresolved Patterns
```bash
python3 scripts/extract_unresolved_patterns.py \
    --consensus-dir output/05_consensus/ \
    --output unresolved_patterns.csv
```

This script:
- Scans all consensus profiles for genes with `MAPPING` or `CONFLICT` status
- Groups them by unique annotation patterns (PADLOC + DefenseFinder + BLAST)
- Creates a template CSV for manual curation
- Shows you the most common problematic patterns

**Example output:**
```
Total problematic proteins: 342
Unique patterns to review: 45

Top patterns:
 1. n= 89  PADLOC:CBASS_other      DF:CBASS_IIs        Fwd:CBASS_II        Rev:CBASS_II
 2. n= 67  PADLOC:Dynamins         DF:Eleos            Fwd:No_hit          Rev:No_hit
 3. n= 45  PADLOC:No_hit           DF:BREX             Fwd:BREX_I          Rev:BREX_I
```

#### Step 2: Manual Curation

Open `unresolved_patterns.csv` in Excel and fill in three columns for each pattern:

| Column | Purpose | Example Values |
|--------|---------|----------------|
| **TYPE** | Defence system type | `CBASS`, `BREX`, `RM` |
| **SUBTYPE** | Defence system subtype | `CBASS_IIs`, `BREX_I`, `RM_I` |
| **OUTCOME** | Biological outcome | `Abi`, `Unknown`, `Non-abi` |

**Resolution Guidelines:**

1. **When tools agree:** Use the consensus annotation
```
   PADLOC=CBASS_other, DF=CBASS_IIs, Fwd=CBASS_II, Rev=CBASS_II
   → TYPE=CBASS, SUBTYPE=CBASS_IIs, OUTCOME=Abi
```

2. **When BLAST provides clarity:** Trust BLAST if tools disagree
```
   PADLOC=No_hit, DF=BREX, Fwd=BREX_I, Rev=BREX_I
   → TYPE=BREX, SUBTYPE=BREX_I, OUTCOME=Non-abi
```

3. **When genuinely unresolvable:** Use special values
```
   PADLOC=PDC-S30, DF=RM_Type_I, Fwd=No_hit, Rev=No_hit
   → TYPE=type_unresolved, SUBTYPE=subtype_unresolved, OUTCOME=outcome_unresolved
```

**Important:** 
- `type_unresolved`, `subtype_unresolved`, `outcome_unresolved` = curator could not resolve (but still a defence gene!)
- `Unknown` = tool detected but classified as unknown
- These are **different** meanings!

#### Step 3: Apply Resolutions

After curating, create the final defence matrix:
```bash
python3 scripts/create_final_defence_matrix.py \
    --consensus-dir output/05_consensus/ \
    --resolutions unresolved_patterns_CURATED.csv \
    --output-prefix species_defence
```

This generates:
- `species_defence_matrix.csv` - Binary presence/absence matrix
- `species_defence_summary.tsv` - Per-genome summary statistics
- `species_defence_annotations.csv` - Full annotated gene list

**Output Example:**
```
genome_id,TYPE#CBASS#SUBTYPE#CBASS_IIs#OUTCOME#Abi,TYPE#BREX#SUBTYPE#BREX_I#OUTCOME#Non-abi...
1005041.3,1,1...
1214195.3,1,10...
```

---

## Output Files

### Final Output
- **`05_consensus/${GENOME_ID}_defenceprofile.csv`** - Consolidated defence system profile with consensus naming

### Intermediate Files
Each genome generates outputs in subdirectories:
```
output/
├── 01_prokka/${GENOME_ID}/
│   ├── ${GENOME_ID}.faa
│   ├── ${GENOME_ID}.gff
│   └── ...
├── 02_padloc/${GENOME_ID}/
│   └── ${GENOME_ID}_padloc.csv
├── 03_defensefinder/${GENOME_ID}/
│   └── ${GENOME_ID}_defense_finder_genes.tsv
└── 04_blast/${GENOME_ID}/
    ├── ${GENOME_ID}_vs_bcereus_forward_cleaned.txt
    └── ${GENOME_ID}_vs_bcereus_reverse_cleaned.txt
```

## Resource Requirements

**Per genome (default settings):**
- **Time:** 30 minutes (safe time for ~2Mb S. suis genome)
- **Memory:** 16 GB
- **CPUs:** 6 cores

**For 2000+ genomes:**
- Submitted in batches of 1000
- Automatic batch management
- ~24-48 hours total runtime (depending on queue)

## Customization

### Adjusting Resources

Edit `defence_pipeline_array.slurm`:
```bash
#SBATCH --time=00:30:00      # Increase for large genomes
#SBATCH --mem=16G            # Increase for complex genomes
#SBATCH --cpus-per-task=6    # Adjust based on availability
```

### Changing Batch Size

Edit `batch_submit_defence_pipeline.sh`:
```bash
BATCH_SIZE=1000  # Change to 500, 2000, etc.
```

### BLAST Parameters

Modify in `pipeline_functions.sh`:
```bash
-evalue 1e-5           # Stringency threshold
-max_target_seqs 1     # Number of hits to report
```

## Troubleshooting

### Common Issues

1. **"No genome ID found"**
   - Check genome files are named correctly: `genomeID.fna`
   - Verify `genome_list.txt` was generated

2. **PADLOC failures**
   - Ensure Prokka GFF format is compatible (not tested with Bakta/Prodigal gffs yet)
   - Check PADLOC database is installed and matches PADLOC version

3. **DefenseFinder errors**
   - Verify HMMER is loaded/installed
   - Check DefenseFinder models are downloaded

4. **BLAST database not found**
   - Ensure *B. cereus* database files are in `databases/bcereus_db/`
   - Have not tested with custom databases, ensure they fit a coherent naming convention that matches the master key file


## Performance Tips

- **Pre-annotate genomes** if you already have Prokka outputs
- **Adjust batch sizes** based on queue limits
- **Use high-priority partitions** when available
- **Monitor failed jobs** and resubmit individually if needed

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

MIT License

## Citation

If you use Ptolemaea in your research, please cite:
```
Campbell, E.B.T., Skvortsov, T., Creevey, C.J. (2025). Ptolemaea: A comprehensive pipeline 
for antiviral defence system detection and consolidation in bacterial genomes. 
GitHub: https://github.com/ecampbell50/Ptolemaea
```

And please cite the underlying tools and databases:
- **PADLOC:** Payne et al. (2021)
- **DefenseFinder:** Tesson et al. (2022)
- **B. cereus database:** July & Gillis (2025)

## Support

Please let me know about any problems you run into or any questions!

## Acknowledgments

- July & Gillis for the *B. cereus* antiviral defence genes and naming conventions
- PADLOC and DefenseFinder development teams
- Queen's University Belfast HPC team
- All contributors to this project
