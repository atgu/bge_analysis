#!/bin/bash

INPUT_PATH="gs://MY_BUCKET/glimpse-outputs/PROJECT_NAME/PROJECT_NAME.chr22.mt"
OUTPUT_PATH="gs://MY_BUCKET/glimpse-outputs/PROJECT_NAME/PROJECT_NAME.chr22.annotated.mt"
SAMPLE_ANNOTATIONS="sample_manifest.tsv"
SAMPLE_ID_COL_ANN="entity:sample_id"
GROUP_BY_COLS="reported_sex research_project"

python3 -m glimpse_hail_batch.add_impute_info_scores \
    --input-path $INPUT_PATH \
    --output-path $OUTPUT_PATH \
    --sample-annotations $SAMPLE_ANNOTATIONS \
    --sample-id-col-input "s" \
    --sample-id-col-ann $SAMPLE_ID_COL_ANN \
    --group-by-col $GROUP_BY_COLS \
    --hl-init-kwarg "backend=batch" \
    --overwrite
