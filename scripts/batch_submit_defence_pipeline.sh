#!/bin/bash
WORKING_DIR=""
GENOME_DIR="${WORKING_DIR}/genomes"
OUTPUT_BASE="${WORKING_DIR}/output"
GENOME_LIST="${OUTPUT_BASE}/genome_list.txt"
BATCH_SIZE=1000 # Set number of genomes to process per job

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

    JOB_ID=$(sbatch --export=BATCH_START=$BATCH_START --array=1-$GENOMES_IN_BATCH "${WORKING_DIR}/scripts/defence_pipeline_array.slurm" | awk '{print $4}')
    echo "  Job ID: $JOB_ID"

    sleep 3
done

echo "All batches submitted."
