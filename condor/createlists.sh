#!/bin/bash

# Check if run number is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <runnumber>"
    exit 1
fi

RUNNUMBER=$1
mkdir -p lists

BASE_DIR_CALO="/sphenix/lustre01/sphnxpro/production2/run3pp/physics/calofitting/new_newcdbtag_v008/"
BASE_DIR_CLUSTER="/sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_2025p009_v001/DST_TRKR_CLUSTER/"

OUTPUT_FILE_CALO="lists/DST_CALOFITTING_${RUNNUMBER}.list"
OUTPUT_FILE_CLUSTER="lists/DST_CLUSTER_${RUNNUMBER}.list"
OUTPUT_FILE_MAPPED="lists/CALOFITTING_TRACK_${RUNNUMBER}.list"

# Determine the appropriate folder range
RUN_BASE=$((RUNNUMBER / 100 * 100))
RUN_FOLDER="run_000${RUN_BASE}_000$((RUN_BASE + 100))"

DIRECTORY_CALO="${BASE_DIR_CALO}/${RUN_FOLDER}"
DIRECTORY_CLUSTER="${BASE_DIR_CLUSTER}/${RUN_FOLDER}"

# Check directories
if [ ! -d "$DIRECTORY_CALO" ]; then
    echo "Error: Directory $DIRECTORY_CALO does not exist."
    exit 1
fi

if [ ! -d "$DIRECTORY_CLUSTER" ]; then
    echo "Error: Directory $DIRECTORY_CLUSTER does not exist."
    exit 1
fi

# Create file lists
ls ${DIRECTORY_CALO}/*${RUNNUMBER}* > ${OUTPUT_FILE_CALO}
echo "File list saved to $OUTPUT_FILE_CALO"

#find "$DIRECTORY_CALO" -type f -name "*${RUNNUMBER}*" > "$OUTPUT_FILE_CALO"
#echo "File list saved to $OUTPUT_FILE_CALO"

#ls ${DIRECTORY_CLUSTER}/DST_TRKR_CLUSTER*${RUNNUMBER}*.root > ${OUTPUT_FILE_CLUSTER}
#echo "File list saved to $OUTPUT_FILE_CLUSTER"

max_rows=1000

#ls ${DIRECTORY_CLUSTER}/DST_TRKR_CLUSTER*${RUNNUMBER}*.root | head -n ${max_rows} > ${OUTPUT_FILE_CLUSTER}

ls ${DIRECTORY_CLUSTER}/DST_TRKR_CLUSTER*${RUNNUMBER}*.root > ${OUTPUT_FILE_CLUSTER}

#find "$DIRECTORY_CLUSTER" -type f -name "DST_TRKR_CLUSTER*${RUNNUMBER}*.root" \
#  | head -n ${max_rows} > "$OUTPUT_FILE_CLUSTER"

echo "File list saved to $OUTPUT_FILE_CLUSTER"

# Mapping
echo "Creating mapping file: $OUTPUT_FILE_MAPPED"

mapfile -t cluster_files < "$OUTPUT_FILE_CLUSTER"

> "$OUTPUT_FILE_MAPPED"

# Output directory for ROOT files
OUTDIR="output/run_${RUNNUMBER}"
#mkdir -p "${OUTDIR}"

# Directory passed to Fun4All_TrackSeeding
CLUSTERDIR="${DIRECTORY_CLUSTER}/"

count=0

for cluster in "${cluster_files[@]}"; do

#    ((count++))
#    [[ $count -gt $max_rows ]] && break

    cluster_base=$(basename "$cluster")

    echo "${OUTPUT_FILE_CALO} ${cluster_base} ${CLUSTERDIR} ${OUTDIR}" \
        >> "$OUTPUT_FILE_MAPPED"
done

echo "Mapping file saved to $OUTPUT_FILE_MAPPED"

