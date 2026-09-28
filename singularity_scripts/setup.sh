#!/bin/bash
# setup.sh
# One-time setup for the containerised pipeline: pulls the tool images into
# $PTOL_IMAGE_DIR and downloads the PADLOC + DefenseFinder databases into
# $PTOL_PADLOC_DATA / $PTOL_DF_MODELS (all set in ptolemaea.config).
#
# Needs internet, so on an HPC run it on a login / data-mover node, NOT in a
# compute job. Safe to re-run: anything already present is skipped.
#
# Usage: bash singularity_scripts/setup.sh

module load apps/apptainer/1.3.4 2>/dev/null || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/ptolemaea.config"

if ! command -v apptainer >/dev/null 2>&1; then
    echo "ERROR: apptainer not found on PATH. Install or load Apptainer first."
    exit 1
fi

# Keep Apptainer's (large) cache out of $HOME unless the user already chose a location
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${PTOL_REPO_DIR}/.apptainer_cache}"
mkdir -p "${PTOL_IMAGE_DIR}" "${PTOL_PADLOC_DATA}" "${PTOL_DF_MODELS}" "${APPTAINER_CACHEDIR}"

# --- 1. Images (paper/thesis versions) ---------------------------------------
declare -A IMAGES=(
    ["${PTOL_IMG_PYRODIGAL}"]="docker://quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1"
    ["${PTOL_IMG_PADLOC}"]="docker://quay.io/biocontainers/padloc:2.0.0--hdfd78af_1"
    ["${PTOL_IMG_DEFENSEFINDER}"]="docker://quay.io/biocontainers/defense-finder:2.0.1--pyhdfd78af_0"
    ["${PTOL_IMG_BLAST}"]="docker://quay.io/biocontainers/blast:2.16.0--h66d330f_5"
    ["${PTOL_IMG_PANDAS}"]="docker://quay.io/biocontainers/pandas:2.2.1"
    ["${PTOL_IMG_SEABORN}"]="docker://quay.io/biocontainers/seaborn:0.13.2"
)

for IMG in "${!IMAGES[@]}"; do
    if [[ -f "${PTOL_IMAGE_DIR}/${IMG}" ]]; then
        echo "Image present: ${IMG}"
    else
        echo "Pulling ${IMG} <- ${IMAGES[$IMG]}"
        apptainer pull "${PTOL_IMAGE_DIR}/${IMG}" "${IMAGES[$IMG]}" || exit 1
    fi
done

# --- 2. PADLOC database ------------------------------------------------------
if [[ -n "$(ls -A "${PTOL_PADLOC_DATA}" 2>/dev/null)" ]]; then
    echo "PADLOC database present: ${PTOL_PADLOC_DATA}"
else
    echo "Downloading PADLOC database -> ${PTOL_PADLOC_DATA}"
    ${PTOL_APPTAINER} --bind "${PTOL_PADLOC_DATA}:/usr/local/data" \
        "${PTOL_IMAGE_DIR}/${PTOL_IMG_PADLOC}" padloc --db-update || exit 1
fi

# --- 3. DefenseFinder models (pinned versions, see ptolemaea.config) ---------
# `defense-finder update` cannot pin a version and always fetches the latest
# models, so install the exact versions with macsydata (ships in the DF image).
if [[ -n "$(ls -A "${PTOL_DF_MODELS}" 2>/dev/null)" ]]; then
    echo "DefenseFinder models present: ${PTOL_DF_MODELS}"
else
    echo "Installing defense-finder-models ${PTOL_DF_MODELS_VERSION} + CasFinder ${PTOL_CASFINDER_VERSION} -> ${PTOL_DF_MODELS}"
    DF_IMG="${PTOL_IMAGE_DIR}/${PTOL_IMG_DEFENSEFINDER}"
    ${PTOL_APPTAINER} "${DF_IMG}" macsydata install --org mdmparis \
        --target "${PTOL_DF_MODELS}" "defense-finder-models==${PTOL_DF_MODELS_VERSION}" || exit 1
    ${PTOL_APPTAINER} "${DF_IMG}" macsydata install \
        --target "${PTOL_DF_MODELS}" "CasFinder==${PTOL_CASFINDER_VERSION}" || exit 1
fi

echo "Setup complete."
