#!/bin/bash
# Ptolemaea_singularity.sh
# Containerised version of Ptolemaea.sh — every tool runs from an Apptainer .sif
# image (see singularity_scripts/ptolemaea.config). No conda envs, no hmmer/blast
# modules; the only module needed is Apptainer itself.
#
# Usage: bash Ptolemaea_singularity.sh <working_directory>
# The working directory must contain: genomes/ , databases/ , scripts/ .

# --- Apptainer must be on PATH (Kelvin: load the module) ---------------------
# Comment this out / adjust on systems where apptainer is already available.
module load apps/apptainer/1.3.4 2>/dev/null || true

# --- Argument check ----------------------------------------------------------
if [[ -z "$1" ]]; then
    echo "Usage: $0 <working_directory>"
    exit 1
fi

# --- Directories -------------------------------------------------------------
WORKING_DIR=$1
GENOME_DIR="${WORKING_DIR}/genomes"
OUTPUT_BASE="${WORKING_DIR}/output"
MASTER_KEY="${WORKING_DIR}/databases/MASTER_ToolKey.tsv"
BCEREUS_DB="${WORKING_DIR}/databases/bcereus_db/Bcereus_ConsensusDefProts_17Sep25"
BCEREUS_FAA="${WORKING_DIR}/databases/bcereus_db/Bcereus_DefenceProts_17Sep25.faa"
CONSENSUS_SCRIPT="${WORKING_DIR}/scripts/create_defence_profile_direct.py"
GENOME_LIST="${OUTPUT_BASE}/genome_list.txt"

CPUS=8   # threads handed to each tool

# --- Load the containerised functions (this also sources ptolemaea.config) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/pipeline_functions_v2.sh"

# --- Build output dir + genome list if missing -------------------------------
mkdir -p "${OUTPUT_BASE}"
if [[ ! -f "$GENOME_LIST" ]]; then
    ls "${GENOME_DIR}"/*.fna | xargs -n1 basename | sed 's/\.fna$//' > "$GENOME_LIST"
fi

# --- Run the pipeline, one genome at a time ----------------------------------
for GENOME_ID in $(cat "${GENOME_LIST}"); do
    echo "=== Processing ${GENOME_ID} ==="
    GENOME_FILE="${GENOME_DIR}/${GENOME_ID}.fna"

    ## 1 - Gene calling (Pyrodigal)
    run_pyrodigal "${GENOME_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/01_pyrodigal" ${CPUS}

    FAA_FILE="${OUTPUT_BASE}/01_pyrodigal/${GENOME_ID}/${GENOME_ID}.faa"
    GFF_FILE="${OUTPUT_BASE}/01_pyrodigal/${GENOME_ID}/${GENOME_ID}.gff"

    ## 2 - PADLOC
    run_padloc "${FAA_FILE}" "${GFF_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/02_padloc" ${CPUS}

    ## 3 - DefenseFinder
    run_defensefinder "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/03_defensefinder" ${CPUS}

    ## 4 - Bidirectional BLAST + cleaning
    run_bidirectional_blast "${FAA_FILE}" "${GENOME_ID}" "${OUTPUT_BASE}/04_blast" \
                            "${BCEREUS_DB}" "${BCEREUS_FAA}" ${CPUS}

    clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward.txt" forward \
        > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
    clean_blast "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse.txt" reverse \
        > "${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"

    ## 5 - Consensus profile
    PADLOC_FILE="${OUTPUT_BASE}/02_padloc/${GENOME_ID}/${GENOME_ID}_padloc.csv"
    DF_FILE="${OUTPUT_BASE}/03_defensefinder/${GENOME_ID}/${GENOME_ID}_defense_finder_genes.tsv"
    FORWARD_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward_cleaned.txt"
    REVERSE_FILE="${OUTPUT_BASE}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse_cleaned.txt"
    CONSENSUS_FILE="${OUTPUT_BASE}/05_consensus/${GENOME_ID}_defenceprofile.csv"

    create_consensus_profile "${PADLOC_FILE}" "${DF_FILE}" "${FORWARD_FILE}" "${REVERSE_FILE}" \
                             "${MASTER_KEY}" "${CONSENSUS_FILE}" "${CONSENSUS_SCRIPT}"

    echo "=== Pipeline complete for ${GENOME_ID} ==="
done
