#!/bin/bash

BILLING_PROJECT="my_billing_project"
REMOTE_TMPDIR="gs://MY_BUCKET/batch/"
REGIONS_FILE="data/split_reference_panel_genetic_map_inputs.tsv"
OUTPUT_DIR="gs://MY_BUCKET/split-hgdp-reference/"
DOCKER="us-central1-docker.pkg.dev/MY_PROJECT/glimpse2/glimpse2-gcloud:odelaneau_bd93ade"

python3 -m glimpse_hail_batch.split_reference \
    --billing-project $BILLING_PROJECT \
    --remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --regions-file $REGIONS_FILE \
    --region "chr22" \
    --seed 34243253 \
    --docker $DOCKER \
    --cpu 1 \
    --memory "standard" \
    --output-dir $OUTPUT_DIR \
    --batch-name glimpse_split_reference
