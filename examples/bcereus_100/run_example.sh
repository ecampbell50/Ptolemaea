#!/bin/bash
# run_example.sh
# End-to-end worked example: download 100 Bacillus cereus group genomes from NCBI,
# run the full Ptolemaea pipeline on them, build the final matrix and plot a
# summary figure.
#
# Usage (from anywhere):
#   bash examples/bcereus_100/run_example.sh [--slurm] [--jobs N] [--limit N]
#
#   --slurm     run each genome as a SLURM array task, then a dependent job that
#               builds the matrix + figure (recommended on HPC)
#   --jobs N    local mode: genomes processed in parallel (default 1; each uses
#               $PTOL_CPUS threads, default 8)
#   --limit N   only use the first N genomes (e.g. --limit 5 for a quick test)
#
# Everything is written to examples/bcereus_100/run/ (override with
# PTOL_EXAMPLE_DIR). Re-running is safe: finished genomes/steps are skipped.
#
# SLURM resources per genome can be changed with PTOL_SBATCH_ARGS
# (default: "--time=00:30:00 --mem=16G --cpus-per-task=$PTOL_CPUS").
#
# Downloading needs internet: on HPC, run this script on a login node - it
# downloads there and then submits the compute jobs.

set -o pipefail

EX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${EX_DIR}/../../singularity_scripts/ptolemaea.config"
module load apps/apptainer/1.3.4 2>/dev/null || true

WORK="${PTOL_EXAMPLE_DIR:-${EX_DIR}/run}"
ACCESSIONS="${EX_DIR}/accessions.tsv"
PREFIX="bcereus100"

MODE="local"; JOBS=1; LIMIT=100
while [[ $# -gt 0 ]]; do
    case "$1" in
        --slurm)    MODE="slurm" ;;
        --jobs)     JOBS="$2"; shift ;;
        --limit)    LIMIT="$2"; shift ;;
        --task)     MODE="task" ;;        # internal: one SLURM array task
        --finalise) MODE="finalise" ;;    # internal: matrix + figure only
        -h|--help)  sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p "${WORK}/genomes" "${WORK}/output" "${WORK}/logs" "${WORK}/lists"
GENOME_LIST="${WORK}/output/genome_list.txt"

# Run one genome through the pipeline (all five stages)
run_genome() {
    local GID=$1
    echo "${GID}" > "${WORK}/lists/${GID}.txt"
    bash "${PTOL_REPO_DIR}/singularity_scripts/Ptolemaea_singularity.sh" \
        "${WORK}" "${WORK}/lists/${GID}.txt" > "${WORK}/logs/${GID}.log" 2>&1
    if [[ -f "${WORK}/output/05_consensus/${GID}_defenceprofile.csv" ]]; then
        echo "  done: ${GID}"
    else
        echo "  FAILED: ${GID} (see ${WORK}/logs/${GID}.log)"
    fi
}

# Matrix, unresolved patterns and figure from whatever consensus profiles exist
finalise() {
    echo "=== Building final matrix ==="
    cd "${WORK}" || exit 1
    ptol_python "${PTOL_REPO_DIR}/scripts/extract_unresolved_patterns.py" \
        --consensus-dir output/05_consensus/ \
        --output "${PREFIX}_unresolved_patterns.csv" || exit 1
    ptol_python "${PTOL_REPO_DIR}/scripts/create_final_defence_matrix.py" \
        --consensus-dir output/05_consensus/ \
        --output-prefix "${PREFIX}" || exit 1

    echo "=== Plotting ==="
    ${PTOL_APPTAINER} "${PTOL_IMAGE_DIR}/${PTOL_IMG_SEABORN}" python \
        "${EX_DIR}/plot_example.py" \
        --annotations "${PREFIX}_annotations.csv" \
        --summary "${PREFIX}_summary.tsv" \
        --accessions "${ACCESSIONS}" \
        --output "${PREFIX}_defence_overview" || exit 1

    local N_DONE N_ALL
    N_DONE=$(ls output/05_consensus/*_defenceprofile.csv 2>/dev/null | wc -l)
    N_ALL=$(wc -l < "${GENOME_LIST}")
    echo
    echo "Example complete: ${N_DONE}/${N_ALL} genomes profiled. Results in ${WORK}/"
    echo "  ${PREFIX}_defence_overview.png / .pdf   summary figure"
    echo "  ${PREFIX}_matrix.csv                   genome x system matrix"
    echo "  ${PREFIX}_annotations.csv              one row per defence gene"
    echo "  ${PREFIX}_summary.tsv                  per-genome counts"
    echo "  ${PREFIX}_unresolved_patterns.csv      MAPPING/CONFLICT patterns to curate"
}

# --- Internal modes (called by SLURM) ----------------------------------------
if [[ "${MODE}" == "task" ]]; then
    run_genome "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${GENOME_LIST}")"
    exit 0
fi
if [[ "${MODE}" == "finalise" ]]; then
    finalise
    exit $?
fi

# --- 1. Tools and databases ----------------------------------------------------
echo "=== 1. Checking images and databases ==="
bash "${PTOL_REPO_DIR}/singularity_scripts/setup.sh" || exit 1

# --- 2. Download genomes from NCBI ---------------------------------------------
echo "=== 2. Downloading genomes from NCBI ==="
grep -v '^#' "${ACCESSIONS}" | tail -n +2 | cut -f1 | head -n "${LIMIT}" > "${GENOME_LIST}"

TO_GET=()
while read -r ACC; do
    [[ -s "${WORK}/genomes/${ACC}.fna" ]] || TO_GET+=("${ACC}")
done < "${GENOME_LIST}"
echo "$(( $(wc -l < "${GENOME_LIST}") - ${#TO_GET[@]} )) already downloaded, ${#TO_GET[@]} to fetch"

API="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
BATCH=10
for (( i=0; i<${#TO_GET[@]}; i+=BATCH )); do
    CHUNK=$(IFS=,; echo "${TO_GET[*]:i:BATCH}")
    ZIP="${WORK}/genomes/_batch.zip"
    echo "  fetching ${CHUNK//,/ }"
    curl -sSfL --retry 3 ${NCBI_API_KEY:+-H "api-key: ${NCBI_API_KEY}"} \
        -o "${ZIP}" "${API}/${CHUNK}/download?include_annotation_type=GENOME_FASTA" \
        || { echo "ERROR: download failed"; exit 1; }
    # Each assembly is ncbi_dataset/data/<acc>/<acc>_<asm>_genomic.fna -> genomes/<acc>.fna
    unzip -q -o -j "${ZIP}" '*_genomic.fna' -d "${WORK}/genomes/_tmp" || exit 1
    for F in "${WORK}"/genomes/_tmp/*_genomic.fna; do
        ACC=$(basename "$F" | cut -d_ -f1,2)
        mv "$F" "${WORK}/genomes/${ACC}.fna"
    done
    rm -rf "${ZIP}" "${WORK}/genomes/_tmp"
done

MISSING=0
while read -r ACC; do
    [[ -s "${WORK}/genomes/${ACC}.fna" ]] || { echo "  missing: ${ACC}"; MISSING=1; }
done < "${GENOME_LIST}"
[[ ${MISSING} -eq 0 ]] || { echo "ERROR: some genomes did not download"; exit 1; }

# --- 3. Run the pipeline ---------------------------------------------------------
N=$(wc -l < "${GENOME_LIST}")
echo "=== 3. Running Ptolemaea on ${N} genomes (${MODE}) ==="

if [[ "${MODE}" == "slurm" ]]; then
    SBATCH_ARGS="${PTOL_SBATCH_ARGS:---time=00:30:00 --mem=16G --cpus-per-task=${PTOL_CPUS}}"
    ARRAY_ID=$(sbatch --parsable ${SBATCH_ARGS} --job-name=ptol_example \
        --array=1-${N} --output="${WORK}/logs/slurm_%A_%a.out" \
        --wrap="bash '${EX_DIR}/run_example.sh' --task") || exit 1
    FINAL_ID=$(sbatch --parsable --time=00:30:00 --mem=8G --job-name=ptol_example_final \
        --dependency=afterany:${ARRAY_ID} --output="${WORK}/logs/slurm_final_%j.out" \
        --wrap="bash '${EX_DIR}/run_example.sh' --finalise") || exit 1
    echo "Submitted array job ${ARRAY_ID} (${N} genomes) and final job ${FINAL_ID}."
    echo "Monitor with: squeue -u \$USER"
    echo "When done, the figure is ${WORK}/${PREFIX}_defence_overview.png"
    exit 0
fi

export -f run_genome
export WORK PTOL_REPO_DIR
xargs -P "${JOBS}" -I{} bash -c 'run_genome {}' < "${GENOME_LIST}"

# --- 4. Final matrix + figure ------------------------------------------------
finalise
