#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"
REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/test-glimpse/batch/"
SPLIT_REFERENCE_PANEL_OUTPUT_DIR="gs://bge-dragen-imputation/hgdp1kg/"
OUTPUT_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks-chrX/"

python3 -m glimpse_hail_batch.regenerate_chunk_metadata \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --split-reference-dir $SPLIT_REFERENCE_PANEL_OUTPUT_DIR \
    --output-directory $OUTPUT_DIR
