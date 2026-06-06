#!/bin/bash
# Load modules
# Ensure hmmer, defense-finder, and blast are available in path for defense-finder and blast steps
# My system is a hpc and had hmmer and blast installed as modules
# I downloaded defensefinder with pip, so its in path

module load apps/hmmer/3.4/gcc-14.1.0
module load apps/ncbiblast/2.15.0/gcc-14.1.0

# --- Argument check ---
if [[ -z "$1" ]]; then
    echo "Usage: $0 <working_directory>"
    exit 1
fi

# Set directories
WORKING_DIR=$1
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

# EXAMPLE
#$ ls ../genomes/*.fna | xargs -n1 basename | sed 's/\.fna$//'
#       C67_HBHC_prefixed
#       D32_LBLC_prefixed
#       D68_LBHC_prefixed
# Gives just file name within extension

# -- RUN PIPELINE -------------------------------------------------------------------
for GENOME_ID in $(cat "${GENOME_LIST}");
do
        echo "Processing ${GENOME_ID}"
        GENOME_FILE="${GENOME_DIR}"/"${GENOME_ID}.fna"

        ## 1 - PHASE ONE - Run PROKKA -----------------------------------------------
        run_prokka "${GENOME_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/01_prokka" 6                 # 6 = cpus per task, change as needed

        ## 2 - PHASE TWO - Run PADLOC -----------------------------------------------
        FAA_FILE="${OUTPUT_BASE}/01_prokka/${GENOME_ID}/${GENOME_ID}.faa"
        GFF_FILE="${OUTPUT_BASE}/01_prokka/${GENOME_ID}/${GENOME_ID}.gff"
        run_padloc "${FAA_FILE}" "${GFF_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/02_padloc" 6      # 6 = cpus per task, change as needed

        ## 3 - PHASE THREE - Run DefenseFinder --------------------------------------
        run_defensefinder "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/03_defensefinder" 6      # 6 = cpus per task, change as needed

        ## 4 - PHASE FOUR - Run Bidirectional BLAST ---------------------------------
        run_bidirectional_blast "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/04_blast" \
                                "${BCEREUS_DB}" "${BCEREUS_FAA}" 6                              # 6 = cpus per task, change as needed

        clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward.txt" forward > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
        clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse.txt" reverse > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"

        ## 5 - PHASE FIVE - Combine outputs into consensus profile ------------------
        PADLOC_FILE="${OUTPUT_BASE}/02_padloc/${GENOME_ID}/${GENOME_ID}_padloc.csv"
        DF_FILE="${OUTPUT_BASE}/03_defensefinder/${GENOME_ID}/${GENOME_ID}_defense_finder_genes.tsv"
        FORWARD_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
        REVERSE_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"
        CONSENSUS_FILE="${OUTPUT_BASE}/05_consensus/${GENOME_ID}_defenceprofile.csv"

        create_consensus_profile "${PADLOC_FILE}" "${DF_FILE}" "${FORWARD_FILE}" "${REVERSE_FILE}" \
                                "${MASTER_KEY}" "${CONSENSUS_FILE}" "${CONSENSUS_SCRIPT}"

        echo "Pipeline complete for ${GENOME_ID}"
        echo "Use extract_unresolved_patterns.py to review and edit CONFLICT/MAPPING issues"
done