#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"
REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/test-glimpse/batch/"
SPLIT_REFERENCE_PANEL_OUTPUT_DIR="gs://jigold-batch-tmp-ezxyx/split-hgdp-reference/"
SAMPLE_MANIFEST="/Users/jigold/Downloads/sample_700.tsv"
STAGING_REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/glimpse-staging-cohort1/"
FASTA="gs://jigold-batch-tmp-ezxyx/fasta/Homo_sapiens_assembly38.fasta"

OUTPUT_FILE="gs://jigold-batch-tmp-ezxyx/glimpse-outputs/final_merged_cohort1.mt"

python3 -m glimpse_hail_batch.imputation \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --batch-name "glimpse" \
    --docker-hail "hailgenetics/hail:0.2.133" \
    --docker-glimpse "us-central1-docker.pkg.dev/neale-pumas-bge/glimpse2/glimpse2-gcloud:odelaneau_bd93ade" \
    --fasta $FASTA \
    --staging-remote-tmpdir $STAGING_REMOTE_TMPDIR \
    --split-reference-dir $SPLIT_REFERENCE_PANEL_OUTPUT_DIR \
    --samples-per-copy-group 50 \
    --sample-group-size 175 \
    --phase-cpu 4 \
    --phase-memory "standard" \
    --ligate-cpu 8 \
    --ligate-memory "standard" \
    --merge-vcf-cpu 8 \
    --merge-vcf-storage "20Gi" \
    --sample-manifest $SAMPLE_MANIFEST \
    --sample-id-col "collaborator_sample_id" \
    --cram-path-col "genome_cram_path" \
    --cram-index-path-col "genome_crai_path" \
    --output-file $OUTPUT_FILE \
    --save-checkpoints
