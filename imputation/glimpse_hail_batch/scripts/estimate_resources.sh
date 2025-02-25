#!/bin/bash

SPLIT_REFERENCE_PANEL_DIR="gs://jigold-batch-tmp-ezxyx/split-hgdp-reference/"

python3 -m glimpse_hail_batch.resource_estimator \
    --split-reference-dir $SPLIT_REFERENCE_PANEL_DIR \
    --max-runtime-mins 60 \
    --target-runtime-mins 30 \
    --n-samples 30000
