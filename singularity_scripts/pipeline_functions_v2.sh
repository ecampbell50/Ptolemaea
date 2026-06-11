#!/bin/bash
# pipeline_functions.sh
# functions defined for the containerised Ptolemaea pipeline
#
# All tools run from Apptainer .sif images
# Image locations, host bind path, and databases are configured
# ----------------------------------------- in ptolemaea.config

# - Load configuration ----------------------------------------

# Determine directory script resides, cd to it, and assign to _PF_DIR
_PF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_PF_DIR}/ptolemaea.config"

# Function to run gene calling with Pyrodgial (replaced Prokka)
# Produces ${GENOME_ID}.faa (proteins) and ${GENOME_ID}.gff
run_pyrodigal() {
    # Using 'local' ensures global pipeline variables are not changed
    # within the function
    local GENOME_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local CPUS=${4:-8}
    local OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Check to see if pyrodigal output already exists
    if [[ -f "${OUT}/${GENOME_ID}.gff" ]] && [[ -f "${OUT}/${GENOME_ID}.faa" ]]; then
        echo "Pyrodigal already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running Pyrodgial for ${GENOME_ID}"
    mkdir -p "${OUT}"

    # Run pyrodigal
    ${PTOL_APPTAINER} "${PTOL_IMAGE_DIR}/${PTOL_IMG_PYRODIGAL}" pyrodigal \
        -i "${GENOME_FILE}" \
        -a "${OUT}/${GENOME_ID}.faa" \
        -o "${OUT}/${GENOME_ID}.gff" -f gff \
        -p single -g 11

    # Check if pyrodigal ran
    if [[ -f "${OUT}/${GENOME_ID}.gff" ]]; then
        local GENE_COUNT
        # if grep fails, throw away noise and just report '0'
        GENE_COUNT=$(grep -c $'\tCDS\t' "${OUT}/${GENOME_ID}.gff" 2>/dev/null || echo "0")
        echo "Pyrodigal complete: ${GENOME_ID} with ${GENE_COUNT} CDSs"
        return 0
    else
        echo "ERROR: Pyrodigal failed for ${GENOME_ID}"
        return 1
    fi
}

# Function to run PADLOC on Pyrodigal output
# Produces ${GENOME_ID}_padloc.csv (defence systems found by PADLOC's HMMs)
run_padloc() {
    local FAA_FILE=$1
    local GFF_FILE=$2
    local GENOME_ID=$3
    local OUTPUT_DIR=$4
    local CPUS=${5:-8}

    local PADLOC_OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Skip if already done. Two "done" states:
    #  - a CSV containing "@" (the genome ID we prepend below = real results), OR
    #  - an empty CSV (PADLOC ran but found no defence systems)
    if [[ -f "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" ]]; then
        if grep -q "@" "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" 2>/dev/null; then
            echo "PADLOC already complete for ${GENOME_ID}"
            return 0
        elif [[ ! -s "${PADLOC_OUT}/${GENOME_ID}_padloc.csv" ]]; then
            echo "PADLOC already run for ${GENOME_ID} (no defence systems found)"
            return 0
        fi
    fi

    echo "Running PADLOC for ${GENOME_ID}"
    mkdir -p "${PADLOC_OUT}"

    # PADLOC wants absolute paths. realpath turns whatever we were given into a
    # full /mnt/scratch2/... path, which is visible inside the container thanks
    # to the --bind in $PTOL_APPTAINER.
    local FAA_ABS GFF_ABS PADLOC_OUT_ABS
    FAA_ABS=$(realpath "${FAA_FILE}")
    GFF_ABS=$(realpath "${GFF_FILE}")
    PADLOC_OUT_ABS=$(realpath "${PADLOC_OUT}")

    # The image's database is at /usr/local/data, which is READ-ONLY inside the
    # .sif. Overlay the writable host DB dir on top of that path so PADLOC finds
    # the database we downloaded once with `padloc --db-update`.
    ${PTOL_APPTAINER} --bind "${PTOL_PADLOC_DATA}:/usr/local/data" \
        "${PTOL_IMAGE_DIR}/${PTOL_IMG_PADLOC}" padloc \
        --faa "${FAA_ABS}" \
        --gff "${GFF_ABS}" \
        --cpu ${CPUS} \
        --force \
        --outdir "${PADLOC_OUT_ABS}"

    # PADLOC writes the CSV only if it found something. Tag it with the genome ID,
    # or leave an empty CSV as a "ran, found nothing" marker.
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

# Helper: prepend the genome ID to PADLOC's target.name column, so proteins stay
# traceable to their genome once everything is merged later (e.g. NZ_..._60 ->
# GCF_046031115.2@NZ_..._60). Uses pandas (via the pandas image) to edit the CSV.
fix_padloc_output() {
    local PADLOC_CSV=$1
    local GENOME_ID=$2

    # "@" present => already tagged this file, nothing to do
    if grep -q "@" "${PADLOC_CSV}" 2>/dev/null; then
        return 0
    fi

    echo "Adding genome ID to PADLOC output"
    # Only tag if PADLOC actually found something (non-empty file)
    if [[ -s "${PADLOC_CSV}" ]]; then
        ptol_python -c "
import pandas as pd
df = pd.read_csv('${PADLOC_CSV}')
if 'target.name' in df.columns:
    df['target.name'] = '${GENOME_ID}@' + df['target.name'].astype(str)
    df.to_csv('${PADLOC_CSV}', index=False)
    print('  Added genome ID prefix to ${PADLOC_CSV}')
"
    fi
}

# Helper: prepend "${GENOME_ID}@" to every FASTA header in a .faa, so DefenseFinder
# and BLAST results carry the genome ID too. Done in pure awk (no container)
add_genome_id_to_faa() {
    local INPUT_FAA=$1
    local OUTPUT_FAA=$2
    local GENOME_ID=$3

    awk -v gid="${GENOME_ID}" '
        /^>/ {                                  # header lines start with ">"
            h = substr($0, 2)                   # drop the ">"
            if (index(h, gid "@") == 1)         # already tagged? leave it
                print $0
            else
                print ">" gid "@" h             # otherwise add the prefix
            next
        }
        { print }                               # sequence lines pass through
    ' "${INPUT_FAA}" > "${OUTPUT_FAA}"
}

# Function to run DefenseFinder on Pyrodigal output
# Produces ${GENOME_ID}_defense_finder_genes.tsv (and systems/hmmer tsvs)
run_defensefinder() {
    local FAA_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local CPUS=${4:-8}

    local DF_OUT="${OUTPUT_DIR}/${GENOME_ID}"

    # Skip if already done
    if [[ -f "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" ]]; then
        echo "DefenseFinder already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running DefenseFinder for ${GENOME_ID}"
    mkdir -p "${DF_OUT}"

    # Give the proteins genome-tagged headers before searching
    local MODIFIED_FAA="${DF_OUT}/${GENOME_ID}_modified.faa"
    add_genome_id_to_faa "${FAA_FILE}" "${MODIFIED_FAA}" "${GENOME_ID}"

    # Run DefenseFinder from its image. --models-dir points at the host models
    # installed once with `defense-finder update` (must match the DF version — see
    # the coupling note in ptolemaea.config).
    ${PTOL_APPTAINER} "${PTOL_IMAGE_DIR}/${PTOL_IMG_DEFENSEFINDER}" defense-finder run \
        --models-dir "${PTOL_DF_MODELS}" \
        --workers ${CPUS} \
        --out-dir "${DF_OUT}" \
        "${MODIFIED_FAA}"

    # DefenseFinder names its outputs after the input file (..._modified_...), so
    # rename them back to plain ${GENOME_ID}_... for the rest of the pipeline.
    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_genes.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_genes.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv"
        echo "Renamed DefenseFinder output file"
    fi

    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_systems.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_systems.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_systems.tsv"
    fi

    if [[ -f "${DF_OUT}/${GENOME_ID}_modified_defense_finder_hmmer.tsv" ]]; then
        mv "${DF_OUT}/${GENOME_ID}_modified_defense_finder_hmmer.tsv" "${DF_OUT}/${GENOME_ID}_defense_finder_hmmer.tsv"
    fi

    # Tidy up the temporary tagged FAA
    rm -f "${MODIFIED_FAA}"

    # Report (an empty marker file means "ran, found nothing")
    if [[ -f "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" ]]; then
        local DF_COUNT
        DF_COUNT=$(tail -n +2 "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv" 2>/dev/null | wc -l || echo "0")
        echo "DefenseFinder complete: ${GENOME_ID} with ${DF_COUNT} hits"
        return 0
    else
        echo "DefenseFinder produced no output for ${GENOME_ID}"
        touch "${DF_OUT}/${GENOME_ID}_defense_finder_genes.tsv"
        return 0
    fi
}

# Function to run bidirectional BLASTp against the B. cereus defence DB
# Forward = genome proteins vs B. cereus ; Reverse = B. cereus vs genome proteins
run_bidirectional_blast() {
    local FAA_FILE=$1
    local GENOME_ID=$2
    local OUTPUT_DIR=$3
    local BCEREUS_DB=$4
    local BCEREUS_FAA=$5
    local CPUS=${6:-8}

    local BLAST_OUT="${OUTPUT_DIR}/${GENOME_ID}"
    local FORWARD_BLAST="${BLAST_OUT}/${GENOME_ID}_vs_bcereus_forward.txt"
    local REVERSE_BLAST="${BLAST_OUT}/${GENOME_ID}_vs_bcereus_reverse.txt"

    # Skip if both directions already done
    if [[ -f "${FORWARD_BLAST}" ]] && [[ -f "${REVERSE_BLAST}" ]]; then
        echo "BLAST already complete for ${GENOME_ID}"
        return 0
    fi

    echo "Running bidirectional BLAST for ${GENOME_ID}"
    mkdir -p "${BLAST_OUT}"

    # Genome-tagged proteins (so hit IDs stay traceable)
    local MODIFIED_FAA="${BLAST_OUT}/${GENOME_ID}_modified.faa"
    add_genome_id_to_faa "${FAA_FILE}" "${MODIFIED_FAA}" "${GENOME_ID}"

    # Reuse the same blast image for every call below
    local BLAST_IMG="${PTOL_IMAGE_DIR}/${PTOL_IMG_BLAST}"

    # Forward: each genome protein searched against the prebuilt B. cereus DB
    if [[ ! -f "${FORWARD_BLAST}" ]]; then
        echo "  Running forward BLAST..."
        ${PTOL_APPTAINER} "${BLAST_IMG}" blastp \
            -query "${MODIFIED_FAA}" \
            -db "${BCEREUS_DB}" \
            -out "${FORWARD_BLAST}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads ${CPUS}

        local FORWARD_COUNT
        FORWARD_COUNT=$(wc -l < "${FORWARD_BLAST}")
        echo "  Forward BLAST complete with ${FORWARD_COUNT} hits"
    fi

    # Reverse: build a temporary DB from the genome proteins, then search the
    # B. cereus proteins against it.
    if [[ ! -f "${REVERSE_BLAST}" ]]; then
        echo "  Creating temporary BLAST database..."
        local TEMP_DB="${BLAST_OUT}/${GENOME_ID}_temp_db"
        ${PTOL_APPTAINER} "${BLAST_IMG}" makeblastdb -in "${MODIFIED_FAA}" -dbtype prot -out "${TEMP_DB}"

        echo "  Running reverse BLAST..."
        ${PTOL_APPTAINER} "${BLAST_IMG}" blastp \
            -query "${BCEREUS_FAA}" \
            -db "${TEMP_DB}" \
            -out "${REVERSE_BLAST}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen" \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -num_threads ${CPUS}

        local REVERSE_COUNT
        REVERSE_COUNT=$(wc -l < "${REVERSE_BLAST}")
        echo "  Reverse BLAST complete with ${REVERSE_COUNT} hits"

        # Remove the temp DB files (.phr/.pin/.psq etc.)
        rm -f ${TEMP_DB}.*
    fi

    # Tidy up the tagged FAA
    rm -f "${MODIFIED_FAA}"

    echo "BLAST complete: ${GENOME_ID}"
    return 0
}


# --------------------------------------------------------------------------
# Function: clean_blast <input_file> <mode>
#   mode = "forward" or "reverse"
# Strips a trailing "_1" that the B. cereus protein names pick up (an artefact
# of how they were built), EXCEPT for a few names where "_1" is real (DISARM_1
# etc.). For reverse mode it also keeps only the single best hit (by bit score,
# column 12) per genome protein. Pure awk — writes the cleaned table to stdout.
# --------------------------------------------------------------------------
clean_blast() {
    local input_file="$1"
    local mode="$2"

    if [[ ! -f "$input_file" ]]; then
        echo "Error: Input file not found: $input_file" >&2
        return 1
    fi

    if [[ "$mode" != "forward" && "$mode" != "reverse" ]]; then
        echo "Error: Mode must be 'forward' or 'reverse'" >&2
        return 1
    fi

    # Names that legitimately END in "_1" and must NOT be trimmed
    local exceptions=("DISARM" "PD-T7-5" "GAO_19")

    if [[ "$mode" == "forward" ]]; then
        # Forward: clean the subject ID (column 2 = the B. cereus protein)
        awk -v OFS="\t" -v n_exceptions="${#exceptions[@]}" \
            -v exceptions_arr="${exceptions[*]}" \
            '
            BEGIN {
                split(exceptions_arr, exceptions_list, " ")
                n = n_exceptions
            }

            function clean(name) {
                if (name ~ /_1$/) {                       # ends in "_1"?
                    base = substr(name, 1, length(name)-2)
                    for (i = 1; i <= n; i++) {            # is it an exception?
                        if (base == exceptions_list[i]) return name
                    }
                    return base                           # trim the "_1"
                }
                return name
            }

            { $2 = clean($2); print }
            ' "$input_file"

    elif [[ "$mode" == "reverse" ]]; then
        # Reverse: clean the query ID (column 1) AND keep only the best hit per
        # genome protein (subject ID, column 2) by bit score (column 12).
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
                        if (base == exceptions_list[i]) return name
                    }
                    return base
                }
                return name
            }

            {
                $1 = clean($1)              # clean query ID
                key = $2                    # genome protein
                score = $12 + 0             # bit score, forced numeric

                # remember the highest-scoring line for each genome protein
                if (!(key in best_score) || score > best_score[key]) {
                    best_score[key] = score
                    best_line[key]  = $0
                }
            }

            END {
                for (k in best_line) print best_line[k]
            }
            ' "$input_file"
    fi
}


# Function to merge all four tool outputs into one consensus profile per genome
# Calls create_defence_profile_direct.py (pandas) inside the pandas image.
create_consensus_profile() {
    local PADLOC_FILE=$1
    local DF_FILE=$2
    local FORWARD_BLAST=$3
    local REVERSE_BLAST=$4
    local MASTER_KEY=$5
    local OUTPUT_FILE=$6
    local CONSENSUS_SCRIPT=$7

    # Skip if already done
    if [[ -f "${OUTPUT_FILE}" ]]; then
        echo "Consensus profile already exists"
        return 0
    fi

    echo "Creating consensus defence profile"
    mkdir -p "$(dirname "${OUTPUT_FILE}")"

    # The python script expects all four inputs to exist; create empty
    # placeholders for any tool that found nothing, so it doesn't error.
    [[ ! -f "${PADLOC_FILE}" ]]   && touch "${PADLOC_FILE}"
    [[ ! -f "${DF_FILE}" ]]       && touch "${DF_FILE}"
    [[ ! -f "${FORWARD_BLAST}" ]] && touch "${FORWARD_BLAST}"
    [[ ! -f "${REVERSE_BLAST}" ]] && touch "${REVERSE_BLAST}"

    ptol_python "${CONSENSUS_SCRIPT}" \
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
        local PROTEIN_COUNT
        PROTEIN_COUNT=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
        echo "Consensus profile created with ${PROTEIN_COUNT} defence proteins"
        return 0
    else
        echo "ERROR: Consensus profile creation failed"
        return 1
    fi
}

# Function to test whether a given stage already finished for a genome.
# Used by the SLURM scripts so re-submitted jobs skip completed work.
check_step_complete() {
    local STEP=$1
    local GENOME_ID=$2
    local BASE_DIR=$3

    case ${STEP} in
        pyrodigal)
            [[ -f "${BASE_DIR}/01_pyrodigal/${GENOME_ID}/${GENOME_ID}.gff" ]] && \
            [[ -f "${BASE_DIR}/01_pyrodigal/${GENOME_ID}/${GENOME_ID}.faa" ]]
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