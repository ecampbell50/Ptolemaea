#!/bin/bash
# submit_slurm.sh
# Run the containerised pipeline as a SLURM array: one task per genome in
# <working_directory>/genomes/. Results land in <working_directory>/output/,
# exactly as with Ptolemaea_singularity.sh.
#
# Usage: bash singularity_scripts/submit_slurm.sh <working_directory>
#
# SLURM options (partition, time, memory) via PTOL_SBATCH_ARGS, e.g.
#   export PTOL_SBATCH_ARGS="--partition=k2-hipri --time=00:30:00 --mem=16G --cpus-per-task=8"
# Re-submitting is safe: finished genomes/steps are skipped.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPT_DIR}/ptolemaea.config"

if [[ -z "$1" ]]; then
    echo "Usage: $0 <working_directory>"
    exit 1
fi
WORK="$(cd "$1" && pwd -P)"

mkdir -p "${WORK}/output/lists" "${WORK}/logs"
GENOME_LIST="${WORK}/output/genome_list.txt"
ls "${WORK}"/genomes/*.fna 2>/dev/null | xargs -n1 basename | sed 's/\.fna$//' > "${GENOME_LIST}"
N=$(wc -l < "${GENOME_LIST}")
if [[ ${N} -eq 0 ]]; then
    echo "ERROR: no .fna files in ${WORK}/genomes/"
    exit 1
fi

# Each array task writes a one-genome list and runs the normal pipeline on it
TASK='GID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "'"${GENOME_LIST}"'"); '\
'echo "${GID}" > "'"${WORK}"'/output/lists/${GID}.txt"; '\
'bash "'"${SCRIPT_DIR}"'/Ptolemaea_singularity.sh" "'"${WORK}"'" "'"${WORK}"'/output/lists/${GID}.txt"'

SBATCH_ARGS="${PTOL_SBATCH_ARGS:---time=00:30:00 --mem=16G --cpus-per-task=${PTOL_CPUS}}"
JOB_ID=$(sbatch --parsable ${SBATCH_ARGS} --job-name=ptolemaea --array=1-${N} \
    --output="${WORK}/logs/ptolemaea_%A_%a.out" --wrap="${TASK}") || exit 1

echo "Submitted array job ${JOB_ID}: ${N} genomes. Monitor with: squeue -u \$USER"
echo "Consensus profiles will appear in ${WORK}/output/05_consensus/"
