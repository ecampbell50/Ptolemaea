# Ptolemaea: Antiviral Defence System Consolidation in Bacteria

**Authors:** Emmet B. T. Campbell, Timofey Skvortsov, Christopher J. Creevey

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> *"In Ptolemaea, the fourth and final zone of the ninth circle of Hell, traitors to their guests are punished..."*  
> Much like how bacteriophages, as unwelcome guests in the bacterial cell, face the sophisticated arsenal of antiviral defense systems.

## Overview

Ptolemaea is a tool for achieving consensus naming conventions between two popular antiviral defence detection tools:
- **[PADLOC](https://github.com/padlocbio/padloc)** 
- **[DefenseFinder](https://github.com/mdmparis/defense-finder)**

This pipeline builds upon the work of [July & Gillis (2025)](https://doi.org/10.1038/s41598-025-86748-8), who published a comprehensive list of defence genes across *Bacillus cereus*.

**📝 Citation:** If you use Ptolemaea, please cite this tool, PADLOC, DefenseFinder, and the July & Gillis paper.

## How It Works

Ptolemaea creates a consensus between PADLOC and DefenseFinder outputs by BLASTing proteins against a curated database of antiviral genes. The current version requires you to run PADLOC, DefenseFinder, and BLAST separately, then process outputs through this pipeline.

**🚧 Future Development:** We plan to expand this into a complete tool that accepts FNA/FAA input directly.

## Prerequisites

Before running Ptolemaea, you'll need:

1. **Genome annotation files:**
   - FAA (protein sequences)
   - Matching GFF file
   - **Recommended:** [Prokka](https://github.com/tseemann/prokka) for annotation
   - **Alternative:** Prodigal or Bakta (see PADLOC documentation for FAA/GFF compatibility)

2. **Software dependencies:**
   - [PADLOC](https://github.com/padlocbio/padloc) v2.0.0+ 
   - [DefenseFinder](https://github.com/mdmparis/defense-finder) v2.0.1+ (uses HMMER)
   - NCBI BLAST+ suite
   - Python 3.8+

## Installation

```bash
git clone https://github.com/ecampbell50/Ptolemaea.git
cd Ptolemaea
```

## Pipeline Workflow

### ⚠️ Important: File Naming Convention
All files must follow this naming format: `genomeID.extension` (e.g., `123456.faa`, `123456.gff`)

---

### Step 1: Run PADLOC

```bash
padloc --faa 123456.faa --gff 123456.gff
```

**Fix PADLOC output to include genome IDs:**
```bash
python3 1_genomeID_to_output.py 123456_padloc.csv
```

---

### Step 2: Prepare FAA File

Add genome ID prefix to locus tags:
```bash
python3 2_add_genomeID_to_faa.py 123456.faa
```

**Example transformation:**
```
Before: >locus_001 hypothetical protein
After:  >123456@locus_001 hypothetical protein
```

---

### Step 3: Run DefenseFinder

```bash
defense-finder run 123456.faa
```

---

### Step 4: BLAST Analysis

#### Forward BLAST (*B. cereus* as subject database)

```bash
blastp \
    -query 123456.faa \
    -db Bcereus_data/PROT_ConsensusBcereus_20Mar25 \
    -out 123456_vs_bcereus_forward.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -num_threads 8  # Adjust based on your system
```

#### Reverse BLAST (*B. cereus* as query)

```bash
# Create protein database from your genome
makeblastdb -in 123456.faa -dbtype prot -out 123456_PROTdb

# Run reverse BLAST
blastp \
    -query BcereusDFseqs_wHASH_FRAME1AA.faa \
    -db 123456_PROTdb \
    -out 123456_vs_bcereus_reverse.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -num_threads 8  # Adjust based on your system
```

---

## Expected Output Files

So far, you should have:

- `123456_padloc.csv` (processed)
- `123456_defense_finder_genes.tsv`
- `123456_vs_bcereus_forward.txt`
- `123456_vs_bcereus_reverse.txt`
- Original files with `.original` backups

## Troubleshooting

### Common Issues

1. **File naming errors:** Ensure all files follow the `genomeID.extension` format
2. **PADLOC/GFF compatibility:** Use Prokka for best results
3. **Missing dependencies:** PADLOC installs well with mamba, and DefenseFinder installs with pip, but I had to load HMMER from my HPC's available modules

### Getting Help

If you encounter issues:
1. Check file naming conventions
2. Verify software versions match requirements
3. Ensure input files are properly formatted

## Contributing

Please feel free to submit issues or pull requests for improvements!

## License

MIT License

## Acknowledgments

- July & Gillis for the *B. cereus* antiviral defence genes and names
- PADLOC and DefenseFinder development teams
- All contributors to this project
