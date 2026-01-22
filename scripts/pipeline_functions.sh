#!/bin/bash
# pipeline_functions.sh
# Reusable functions for defence system pipeline
# Source this file in your SLURM scripts with: source pipeline_functions.sh
##### SET CONDA ENVS SPECIFIC TO YOUR SYSTEM (lines 27, 70,  


# Function to run Prokka annotation
run_prokka() {
    local GENOME_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local CPUS=${4:-6}

    local PROKKA_OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Check if already complete
    if [[ -f "${PROKKA_OUT}/${GENOME_ID}.gff" ]] && [[ -f "${PROKKA_OUT}/${GENOME_ID}.faa" ]]; then
        echo "Prokka already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running Prokka for ${GENOME_ID}"
    mkdir -p "${PROKKA_OUT}"

    # Run exactly as in original script - using conda run
    conda run -p /path/to/prokka/env/prokka prokka \
        --outdir "${PROKKA_OUT}" \
        --prefix "${GENOME_ID}" \
        --kingdom Bacteria \
        --gcode 11 \
        --cpus ${CPUS} \
        --force \
        "${GENOME_FILE}"

    # Check success
    if [[ -f "${PROKKA_OUT}/${GENOME_ID}.gff" ]]; then
        local GENE_COUNT=$(grep -c "CDS" "${PROKKA_OUT}/${GENOME_ID}.gff" || echo "0")
        echo "Prokka complete: ${GENOME_ID} with ${GENE_COUNT} genes"
        return 0
    else
        echo "ERROR: Prokka failed for ${GENOME_ID}"
        return 1
    fi
}

# Function to run PADLOC
run_padloc() {
    local FAA_FILE=$1
    local GFF_FILE=$2
    local GENOME_ID=$3
    local OUTPUT_DIR=$4
    local CPUS=${5:-6}

    local PADLOC_OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Check if already complete with genome ID added
    if [[ -f "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" ]]; then
        if grep -q "@" "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" 2>/dev/null; then
            echo "PADLOC already complete for ${GENOME_ID}"
            return 0
        fi
    fi

    echo "Running PADLOC for ${GENOME_ID}"
    mkdir -p "${PADLOC_OUT}"
    cd "${PADLOC_OUT}"

    # Run exactly as in original script - using conda run
    conda run -p /path/to/padloc/env/padloc padloc \
        --faa "${FAA_FILE}" \
        --gff "${GFF_FILE}" \
        --cpu ${CPUS} \
        --outdir .

    cd - > /dev/null

    # Check and fix output
    if [[ -f "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" ]]; then
        fix_padloc_output "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" "${GENOME_ID}"
        echo "PADLOC complete: ${GENOME_ID}"
        return 0
    else
        echo "PADLOC produced no output for ${GENOME_ID} (no defence systems found)"
        touch "${PADLOC_OUT}/${GENOME_ID}_padloc.csv"
        return 0
    fi
}

# Function to fix PADLOC output (add genome ID)
fix_padloc_output() {
    local PADLOC_CSV=$1
    local GENOME_ID=$2

    # Check if genome ID already added
    if grep -q "@" "${PADLOC_CSV}" 2>/dev/null; then
        return 0
    fi

    echo "Adding genome ID to PADLOC output"

    python3 -c "
import pandas as pd
df = pd.read_csv('${PADLOC_CSV}')
if 'target.name' in df.columns:
    df['target.name'] = '${GENOME_ID}@' + df['target.name'].astype(str)
    df.to_csv('${PADLOC_CSV}', index=False)
    print('  Added genome ID prefix to ${PADLOC_CSV}')
"
}

# Function to add genome ID to FAA headers
add_genome_id_to_faa() {
    local INPUT_FAA=$1
    local OUTPUT_FAA=$2
    local GENOME_ID=$3

    python3 -c "
with open('${INPUT_FAA}', 'r') as infile:
    with open('${OUTPUT_FAA}', 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line[1:].strip()
                if not header.startswith('${GENOME_ID}@'):
                    outfile.write(f'>${GENOME_ID}@{header}\n')
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
"
}

# Function to run DefenseFinder
run_defensefinder() {
    local FAA_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local CPUS=${4:-6}

    local DF_OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Check if already complete
    if [[ -f "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" ]]; then
        echo "DefenseFinder already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running DefenseFinder for ${GENOME_ID}"
    mkdir -p "${DF_OUT}"

    # Create modified FAA with genome ID
    local MODIFIED_FAA="${DF_OUT}/${GENOME_ID}_modified.faa"
    add_genome_id_to_faa "${FAA_FILE}" "${MODIFIED_FAA}" "${GENOME_ID}"

    # Run DefenseFinder (pip installed, in PATH)
    defense-finder run \
        --workers ${CPUS} \
        --out-dir "${DF_OUT}" \
        "${MODIFIED_FAA}"

    # DefenseFinder creates files based on input filename, so we need to rename them
    # From: 1307.762_modified_defense_finder_genes.tsv
    # To:   1307.762_defense_finder_genes.tsv

    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_genes.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_genes.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv"
        echo "Renamed DefenseFinder output file"
    fi

    # Also rename other DefenseFinder output files if they exist (for consistency)
    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_systems.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_systems.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_systems.tsv"
    fi

    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_hmmer.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_hmmer.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_hmmer.tsv"
    fi

    # Clean up modified FAA
    rm -f "${MODIFIED_FAA}"

    # Check success
    if [[ -f "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" ]]; then
        local DF_COUNT=$(tail -n +2 "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" 2>/dev/null | wc -l || echo "0")
        echo "DefenseFinder complete: ${GENOME_ID} with ${DF_COUNT} hits"
        return 0
    else
        echo "DefenseFinder produced no output for ${GENOME_ID}"
        touch "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv"
        return 0
    fi
}

# Function to run bidirectional BLAST
run_bidirectional_blast() {
    local FAA_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local BCEREUS_DB=$4
    local BCEREUS_FAA=$5
    local CPUS=${6:-6}

    local BLAST_OUT="${OUTPUT_DIR}/${GENOME_ID}"
    local FORWARD_BLAST="${BLAST_OUT}/${GENOME_ID}_vs_bcereus_forward.txt"
    local REVERSE_BLAST="${BLAST_OUT}/${GENOME_ID}_vs_bcereus_reverse.txt"

    # Check if already complete
    if [[ -f "${FORWARD_BLAST}" ]] && [[ -f "${REVERSE_BLAST}" ]]; then
        echo "BLAST already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running bidirectional BLAST for ${GENOME_ID}"
    mkdir -p "${BLAST_OUT}"

    # Create modified FAA with genome ID
    local MODIFIED_FAA="${BLAST_OUT}/${GENOME_ID}_modified.faa"
    add_genome_id_to_faa "${FAA_FILE}" "${MODIFIED_FAA}" "${GENOME_ID}"

    # Forward BLAST
    if [[ ! -f "${FORWARD_BLAST}" ]]; then
        echo "  Running forward BLAST..."
        blastp \
            -query "${MODIFIED_FAA}" \
            -db "${BCEREUS_DB}" \
            -out "${FORWARD_BLAST}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads ${CPUS}

        local FORWARD_COUNT=$(wc -l < "${FORWARD_BLAST}")
        echo "  Forward BLAST complete with ${FORWARD_COUNT} hits"
    fi

    # Reverse BLAST
    if [[ ! -f "${REVERSE_BLAST}" ]]; then
        echo "  Creating temporary BLAST database..."
        local TEMP_DB="${BLAST_OUT}/${GENOME_ID}_temp_db"
        makeblastdb -in "${MODIFIED_FAA}" -dbtype prot -out "${TEMP_DB}"

        echo "  Running reverse BLAST..."
        blastp \
            -query "${BCEREUS_FAA}" \
            -db "${TEMP_DB}" \
            -out "${REVERSE_BLAST}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads ${CPUS}

        local REVERSE_COUNT=$(wc -l < "${REVERSE_BLAST}")
        echo "  Reverse BLAST complete with ${REVERSE_COUNT} hits"

        # Clean up
        rm -f ${TEMP_DB}.*
    fi

    # Clean up modified FAA
    rm -f "${MODIFIED_FAA}"

    echo "BLAST complete: ${GENOME_ID}"
    return 0
}


# --------------------------------------------
# Function: clean_blast
# Usage: clean_blast <input_file> <mode>
#   input_file: path to BLAST tabular output file (outfmt 6)
#   mode: "forward" or "reverse"
#
# Output:
#   Writes cleaned (and filtered for reverse mode) BLAST results to stdout
#
# Description:
#   - Cleans trailing '_1' from specific columns based on mode
#   - For forward mode: cleans subject IDs (column 2)
#   - For reverse mode: cleans query IDs (column 1) and keeps only
#     best hit per genome protein (subject ID, column 2) by bit score (column 12)
# --------------------------------------------
clean_blast() {
    local input_file="$1"
    local mode="$2"

    # Check file exists
    if [[ ! -f "$input_file" ]]; then
        echo "Error: Input file not found: $input_file" >&2
        return 1
    fi

    # Validate mode
    if [[ "$mode" != "forward" && "$mode" != "reverse" ]]; then
        echo "Error: Mode must be 'forward' or 'reverse'" >&2
        return 1
    fi

    # Exception names: these legitimate names keep their trailing '_1'
    # e.g. DISARM_1 should NOT lose '_1'
    local exceptions=("DISARM" "PD-T7-5" "GAO_19")

    # The clean function replicates your Python logic:
    # Remove trailing '_1' unless base is in exceptions.
    # Also removes only one trailing '_1', so DISARM_1_1 → DISARM_1
    clean_name_func='
    function clean(name, exceptions_arr, n_exceptions) {
        # Check if name ends with "_1"
        if (name ~ /_1$/) {
            base = substr(name, 1, length(name)-2)  # Remove last 2 chars "_1"
            # Check if base is in exceptions
            for (i=1; i<=n_exceptions; i++) {
                if (base == exceptions_arr[i]) {
                    return name  # Keep original if exception
                }
            }
            return base  # Remove trailing "_1" if not exception
        }
        return name  # Return as is if no trailing "_1"
    }
    '

    if [[ "$mode" == "forward" ]]; then
        # Forward BLAST: Clean subject IDs (column 2)
        awk -v OFS="\t" -v n_exceptions="${#exceptions[@]}" \
            -v exceptions_arr="${exceptions[*]}" \
            -v clean_func="$clean_name_func" \
            '
            BEGIN {
                split(exceptions_arr, exceptions_list, " ")
                n = n_exceptions
                # Define the clean function from variable clean_func (awk -f trick)
                # Instead, we define inline below (see below)
            }

            # Define the cleaning function inline:
            function clean(name) {
                if (name ~ /_1$/) {
                    base = substr(name, 1, length(name)-2)
                    for (i = 1; i <= n; i++) {
                        if (base == exceptions_list[i]) {
                            return name
                        }
                    }
                    return base
                }
                return name
            }

            {
                $2 = clean($2)   # Clean subject ID (column 2)
                print
            }
            ' "$input_file"

    elif [[ "$mode" == "reverse" ]]; then
        # Reverse BLAST:
        # 1. Clean query IDs (column 1)
        # 2. Keep only best hit per genome protein (subject ID in col 2) by bit score (col 12)

        awk -v OFS="\t" -v n_exceptions="${#exceptions[@]}" \
            -v exceptions_arr="${exceptions[*]}" \
            '
            BEGIN {
                split(exceptions_arr, exceptions_list, " ")
                n = n_exceptions
            }

            function clean(name) {
                if (name ~ /_1$/) {
                    base = substr(name, 1, length(name)-2)
                    for (i = 1; i <= n; i++) {
                        if (base == exceptions_list[i]) {
                            return name
                        }
                    }
                    return base
                }
                return name
            }

            {
                $1 = clean($1)                  # Clean query ID (col 1)

                key = $2                       # Use subject ID (genome protein) as key
                score = $12 + 0               # Bit score (field 12), numeric

                # Keep best hit per genome protein by bit score
                if (!(key in best_score) || score > best_score[key]) {
                    best_score[key] = score
                    best_line[key] = $0
                }
            }

            END {
                # Print best hit line for each genome protein
                for (k in best_line) {
                    print best_line[k]
                }
            }
            ' "$input_file"
    fi
}


# Function to create consensus defence profile
create_consensus_profile() {
    local PADLOC_FILE=$1
    local DF_FILE=$2
    local FORWARD_BLAST=$3
    local REVERSE_BLAST=$4
    local MASTER_KEY=$5
    local OUTPUT_FILE=$6
    local CONSENSUS_SCRIPT=$7

    # Check if already complete
    if [[ -f "${OUTPUT_FILE}" ]]; then
        echo "Consensus profile already exists"
        return 0
    fi

    echo "Creating consensus defence profile"

    mkdir -p "$(dirname "${OUTPUT_FILE}")"

    # Ensure all input files exist (create empty if missing)
    [[ ! -f "${PADLOC_FILE}" ]] && touch "${PADLOC_FILE}"
    [[ ! -f "${DF_FILE}" ]] && touch "${DF_FILE}"
    [[ ! -f "${FORWARD_BLAST}" ]] && touch "${FORWARD_BLAST}"
    [[ ! -f "${REVERSE_BLAST}" ]] && touch "${REVERSE_BLAST}"

    # Run consensus script
    python3 "${CONSENSUS_SCRIPT}" \
        --padloc "${PADLOC_FILE}" \
        --defensefinder "${DF_FILE}" \
        --forward-blast "${FORWARD_BLAST}" \
        --reverse-blast "${REVERSE_BLAST}" \
        --master-key "${MASTER_KEY}" \
        --output "${OUTPUT_FILE}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Consensus script execution failed"
        return 1
    fi

    if [[ -f "${OUTPUT_FILE}" ]]; then
        local PROTEIN_COUNT=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
        echo "Consensus profile created with ${PROTEIN_COUNT} defence proteins"
        return 0
    else
        echo "ERROR: Consensus profile creation failed"
        return 1
    fi
}

# Function to check step completion
check_step_complete() {
    local STEP=$1
    local GENOME_ID=$2
    local BASE_DIR=$3

    case ${STEP} in
        prokka)
            [[ -f "${BASE_DIR}/01_prokka/${GENOME_ID}/${GENOME_ID}.gff" ]] && \
            [[ -f "${BASE_DIR}/01_prokka/${GENOME_ID}/${GENOME_ID}.faa" ]]
            ;;
        padloc)
            [[ -f "${BASE_DIR}/02_padloc/${GENOME_ID}/${GENOME_ID}_padloc.csv" ]] && \
            grep -q "@" "${BASE_DIR}/02_padloc/${GENOME_ID}/${GENOME_ID}_padloc.csv" 2>/dev/null
            ;;
        defensefinder)
            [[ -f "${BASE_DIR}/03_defensefinder/${GENOME_ID}/${GENOME_ID}_defense_finder_genes.tsv" ]]
            ;;
        blast)
            [[ -f "${BASE_DIR}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_forward.txt" ]] && \
            [[ -f "${BASE_DIR}/04_blast/${GENOME_ID}/${GENOME_ID}_vs_bcereus_reverse.txt" ]]
            ;;
        consensus)
            [[ -f "${BASE_DIR}/05_consensus/${GENOME_ID}_defenceprofile.csv" ]]
            ;;
        *)
            return 1
            ;;
    esac
}
