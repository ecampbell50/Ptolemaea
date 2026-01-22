#!/bin/bash
#SBATCH --job-name=defence_pipeline
#SBATCH --error=logs/defence_%A_%a.err
#SBATCH --output=logs/defence_%A_%a.out
#SBATCH --partition= # Specify partition
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=6
#SBATCH --array=1-1000

# Load modules
# Ensure hmmer and blast are available in path for defense-finder and blast steps
module load apps/hmmer/3.4/gcc-14.1.0
module load apps/ncbiblast/2.15.0/gcc-14.1.0

WORKING_DIR=""
GENOME_DIR="${WORKING_DIR}/genomes"
OUTPUT_BASE="${WORKING_DIR}/output"
MASTER_KEY="${WORKING_DIR}/databases/MASTER_ToolKey.tsv"
BCEREUS_DB="${WORKING_DIR}/databases/bcereus_db/Bcereus_ConsensusDefProts_17Sep25"
BCEREUS_FAA="${WORKING_DIR}/databases/bcereus_db/Bcereus_DefenceProts_17Sep25.faa"
CONSENSUS_SCRIPT="${WORKING_DIR}/scripts/create_defence_profile_direct.py"
GENOME_LIST="${OUTPUT_BASE}/genome_list.txt"

# Source functions
source "${WORKING_DIR}/scripts/pipeline_functions.sh"

# Make output directory and genome list if missing
mkdir -p "${OUTPUT_BASE}"
if [[ ! -f "$GENOME_LIST" ]]; then
    ls "${GENOME_DIR}"/*.fna | xargs -n1 basename | sed 's/\.fna$//' > "$GENOME_LIST"
fi

# Calculate the actual genome index for this task
GENOME_INDEX=$((${BATCH_START:-0} + ${SLURM_ARRAY_TASK_ID}))

# Get genome ID for this index
GENOME_ID=$(sed -n "${GENOME_INDEX}p" "$GENOME_LIST")

if [[ -z "$GENOME_ID" ]]; then
    echo "No genome ID found for index $GENOME_INDEX (batch start: ${BATCH_START:-0}, task: ${SLURM_ARRAY_TASK_ID}), skipping."
    exit 0
fi

echo "Processing genome $GENOME_INDEX: $GENOME_ID"

GENOME_FILE="${GENOME_DIR}/${GENOME_ID}.fna"

## 1 - PHASE ONE - Run PROKKA
run_prokka "${GENOME_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/01_prokka" ${SLURM_CPUS_PER_TASK}

## 2 - PHASE TWO - Run PADLOC
FAA_FILE="${OUTPUT_BASE}/01_prokka/${GENOME_ID}/${GENOME_ID}.faa"
GFF_FILE="${OUTPUT_BASE}/01_prokka/${GENOME_ID}/${GENOME_ID}.gff"
run_padloc "${FAA_FILE}" "${GFF_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/02_padloc" ${SLURM_CPUS_PER_TASK}

## 3 - PHASE THREE - Run DefenseFinder
run_defensefinder "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/03_defensefinder" ${SLURM_CPUS_PER_TASK}

## 4 - PHASE FOUR - Run Bidirectional BLAST
run_bidirectional_blast "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/04_blast" \
    "${BCEREUS_DB}" "${BCEREUS_FAA}" ${SLURM_CPUS_PER_TASK}

clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward.txt" forward > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse.txt" reverse > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"

## 5 - PHASE FIVE - Combine outputs into consensus profile
PADLOC_FILE="${OUTPUT_BASE}/02_padloc/${GENOME_ID}/${GENOME_ID}_padloc.csv"
DF_FILE="${OUTPUT_BASE}/03_defensefinder/${GENOME_ID}/${GENOME_ID}_defense_finder_genes.tsv"
FORWARD_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
REVERSE_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"
CONSENSUS_FILE="${OUTPUT_BASE}/05_consensus/${GENOME_ID}_defenceprofile.csv"

create_consensus_profile "${PADLOC_FILE}" "${DF_FILE}" "${FORWARD_FILE}" "${REVERSE_FILE}" \
    "${MASTER_KEY}" "${CONSENSUS_FILE}" "${CONSENSUS_SCRIPT}"

echo "Pipeline complete for ${GENOME_ID}"
