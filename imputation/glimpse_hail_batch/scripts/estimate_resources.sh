#!/bin/bash

BILLING_PROJECT="my_billing_project"

CHUNK_INFO_DIR="gs://MY_BUCKET/glimpse_reference/chunks/"
BINARY_REFERENCE_DIR="gs://MY_BUCKET/glimpse_reference/binary_reference/"
FILE_REGEX="ref_chunk_(?P<contig>.+)_(?P<chunk_index>\d+).bin"
CHUNK_FILE_REGEX="chunks_(?P<contig>.+)\.txt"

N_SAMPLES=30000

python3 -m glimpse_hail_batch.resource_estimator \
    --binary-reference-file-regex $FILE_REGEX \
    --reference-dir $BINARY_REFERENCE_DIR \
    --chunk-file-regex $CHUNK_FILE_REGEX \
    --chunk-info-dir $CHUNK_INFO_DIR \
    --n-samples $N_SAMPLES \
    --gcs-requester-pays-configuration $BILLING_PROJECT
