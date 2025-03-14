#!/bin/bash

BILLING_PROJECT="my_billing_project"
BATCH_REGIONS="us-central1"
REMOTE_TMPDIR="gs://MY_BUCKET/test-glimpse/batch/"

REFERENCE_DIR="gs://broad-dsde-methods-bge-resources-public/GlimpseImputation/ReferencePanels/1000G_HGDP_no_singletons/chunks/cm_4/boost_1_78_0/bin/"
BINARY_REFERENCE_FILE_REGEX="reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"

# use this to update existing chunk files with the actual number of common and rare variants in the chunk
# we needed to use this option because we filtered the reference panel after splitting the reference and
# generating the original chunk files
CHUNK_INFO_DIR="--chunk-info-dir gs://MY_BUCKET/reference_panel/chunks/"

OUTPUT_DIR="gs://MY_BUCKET/regenerate-chunks-1000G_HGDP_no_singletons/"

python3 -m glimpse_hail_batch.regenerate_chunk_metadata \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions $BATCH_REGIONS \
    --reference-dir $REFERENCE_DIR \
    $CHUNK_INFO_DIR \
    --output-directory $OUTPUT_DIR \
    --gcs-requester-pays-configuration $BILLING_PROJECT \
    --binary-reference-file-regex $BINARY_REFERENCE_FILE_REGEX
