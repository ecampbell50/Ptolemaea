#!/bin/bash
WORKING_DIR=""   # <-- SET THIS to your Ptolemaea working directory (absolute path)
BATCH_SIZE=1000  # Set number of genomes to process per job

# --- Safety check: WORKING_DIR must be set ---
if [[ -z "$WORKING_DIR" ]]; then
    echo "ERROR: WORKING_DIR is empty. Edit this script and set it before running."
    exit 1
fi

GENOME_DIR="${WORKING_DIR}/genomes"
OUTPUT_BASE="${WORKING_DIR}/output"
GENOME_LIST="${OUTPUT_BASE}/genome_list.txt"
ARRAY_SCRIPT="${WORKING_DIR}/scripts/defence_pipeline_array.sh"

# --- Safety check: array script must exist ---
if [[ ! -f "$ARRAY_SCRIPT" ]]; then
    echo "ERROR: array job script not found: $ARRAY_SCRIPT"
    exit 1
fi

mkdir -p "${WORKING_DIR}/logs"
mkdir -p "${OUTPUT_BASE}"

# Ensure genome list exists
if [[ ! -f "$GENOME_LIST" ]]; then
    ls "${GENOME_DIR}"/*.fna | xargs -n1 basename | sed 's/\.fna$//' > "$GENOME_LIST"
fi

TOTAL_GENOMES=$(wc -l < "$GENOME_LIST")

echo "Total genomes: $TOTAL_GENOMES"
NUM_BATCHES=$(( (TOTAL_GENOMES + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Number of batches needed: $NUM_BATCHES"

read -p "Submit all $NUM_BATCHES batches? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

for ((batch=0; batch<NUM_BATCHES; batch++)); do
    BATCH_START=$((batch * BATCH_SIZE))

    REMAINING=$((TOTAL_GENOMES - BATCH_START))
    GENOMES_IN_BATCH=$(( REMAINING < BATCH_SIZE ? REMAINING : BATCH_SIZE ))

    echo "Submitting batch $((batch+1))/$NUM_BATCHES, genomes $((BATCH_START+1)) to $((BATCH_START+GENOMES_IN_BATCH))"

    # --export=ALL keeps the current environment (conda init, defense-finder on PATH,
    # loaded modules) and adds the two pipeline variables. WORKING_DIR is passed through
    # so it only needs to be set here, not also in the array script.
    JOB_ID=$(sbatch --chdir="${WORKING_DIR}" \
                    --export=ALL,BATCH_START=$BATCH_START,WORKING_DIR="${WORKING_DIR}" \
                    --array=1-$GENOMES_IN_BATCH \
                    "${ARRAY_SCRIPT}" | awk '{print $4}')
    echo "  Job ID: $JOB_ID"

    sleep 3
done

echo "All batches submitted."
