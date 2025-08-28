# Ptolemaea: Antiviral Defence System Consolidation in bacteria
Emmet B. T. Campbell - Timofey Skvortsov - Christopher J. Creevey

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> *"In Ptolemaea, the fourth and final zone of the ninth circle of Hell, traitors to their guests are punished..."*  
> Much like how bacteriophages, as unwelcome guests in the bacterial cell, face the sophisticated arsenal of antiviral defense systems.

## Overview

Ptolemaea is a tool for achieving a consensus naming convention between the two popular antiviral defence detection tools PADLOC and DefenseFinder

https://github.com/padlocbio/padloc
https://github.com/mdmparis/defense-finder

It piggy-backs off the work of July & Gillis (https://doi.org/10.1038/s41598-025-86748-8), who published a list of defence genes across _Bacillus cereus_.

If you use Ptolemaea, please also cite these papers!

## How it works
Ptolemaea as of now is a pipeline to get a consensus antiviral defence gene name between PADLOC and DefenseFinder by BLASTing proteins against a database of antiviral genes. In the near future we hope to expand it to a full tool that takes FNA/FAA input. Currently, you will need to run PADLOC, DefenseFinder and BLAST yourself and run the outputs through this pipeline.

## What you need before running
1. An FAA and matching GFF from an annotation tool like prokka. Prodigal and Bakta may also work, but please refer to PADLOC's documentation on FAA/GFF matching with tools other than prokka! Prokka will be used for the rest of this pipeline.
2. A PADLOC ...padloc.csv output (not the .gff), DefenseFinder's ...genes.tsv output file, forward, and reverse BLAST against the _B. cereus_ database provided. see further below on the BLAST portion.

As part of this pipeline, we edit the locus tag of all proteins to include the genome ID. If you are running one genome this may be inconsequential, however for >1 genome (in my case, tens of thousands!), this alleviates a lot of problems when running downstream anlayses.

"Unedited .faa"  = FAA file where sequence headers do not contain the genome ID
                  >locus_tag  gene_information
"Edited .faa"    = FAA file where sequence headers are edited to include the genome ID
                  >genomeID@locus_tag  gene_information

### Running PADLOC
With the **unedited** .faa and .gff file from prokka, run PADLOC with the --faa/--gff option. 
the ...padloc.csv can then have the locus tags edited after (PADLOC requires the sequence header and GFF "ID=" to match)

### Running DefenseFinder
With the **edited** FAA file, run DefenseFinder (the --antidefense flag is not supported in this pipeline at the minute, so it may interfere with the pipeline)

### Running BLAST
BLAST is used for two reasons: to find a consensus name, and to find putative genes under strict parameters.
It is a protein database, and therefore will run a blastp.
#### Forward BLAST
Here, your genome is the query, and the _B. cereus_ database is the subject. This essentially matches proteins across your entire genome to any antiviral defence genes in the database.
#### Reverse BLAST
Here, you must make a protein BLAST database out of your FAA file. The provided BcereusDFseqs_wHASH_FRAME1AA.faa should then be BLASTed as a query against your genome database.
Running a forward and reverse BLAST ensures that not only was a particular protein from your genome matched to an antiviral defence gene in the database, but also an antiviral defence gene matched with that protein out of your whole genome (not just defence genes).

Once you have all of these files, Ptolemaea can work to get you a consensus csv of each protein ID with an antiviral defence gene!

