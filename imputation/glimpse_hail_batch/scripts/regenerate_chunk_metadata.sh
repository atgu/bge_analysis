#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"
REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/test-glimpse/batch/"
#SPLIT_REFERENCE_PANEL_OUTPUT_DIR="gs://bge-dragen-imputation/hgdp1kg/"
SPLIT_REFERENCE_PANEL_OUTPUT_DIR="gs://broad-dsde-methods-bge-resources-public/GlimpseImputation/ReferencePanels/1000G_HGDP_no_singletons/chunks/cm_4/boost_1_78_0/bin/"
OUTPUT_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks-1000G_HGDP_no_singletons/"

python3 -m glimpse_hail_batch.regenerate_chunk_metadata \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --reference-dir $SPLIT_REFERENCE_PANEL_OUTPUT_DIR \
    --output-directory $OUTPUT_DIR \
    --gcs-requester-pays-configuration $BILLING_PROJECT \
    --binary-reference-file-regex "reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"
