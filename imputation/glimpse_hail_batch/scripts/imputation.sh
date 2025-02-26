#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"
REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/test-glimpse/batch/"

CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks-1000G_HGDP_no_singletons/"
BINARY_REFERENCE_DIR="gs://broad-dsde-methods-bge-resources-public/GlimpseImputation/ReferencePanels/1000G_HGDP_no_singletons/chunks/cm_4/boost_1_78_0/bin/"
FILE_REGEX="reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"

#SAMPLE_MANIFEST="/Users/jigold/Downloads/Vietnam_SCZ_BGE_Cases.sample_manifest.tsv"
SAMPLE_MANIFEST="/Users/jigold/Downloads/sample_700.tsv"

#STAGING_REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/staging/Vietnam_SCZ_BGE_Cases-test/"
STAGING_REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/staging/test/"

FASTA="gs://jigold-batch-tmp-ezxyx/fasta/Homo_sapiens_assembly38.fasta"

#OUTPUT_FILE="gs://jigold-batch-tmp-ezxyx/glimpse-outputs/Vietnam_SCZ_BGE_Cases.mt"
OUTPUT_FILE="gs://jigold-batch-tmp-ezxyx/glimpse-outputs/test.mt"

python3 -m glimpse_hail_batch.imputation \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --batch-name "glimpse" \
    --docker-hail "hailgenetics/hail:0.2.133" \
    --docker-glimpse "us-central1-docker.pkg.dev/neale-pumas-bge/glimpse2/glimpse2-gcloud:odelaneau_bd93ade" \
    --fasta $FASTA \
    --staging-remote-tmpdir $STAGING_REMOTE_TMPDIR \
    --samples-per-copy-group 50 \
    --sample-group-size 350 \
    --phase-cpu 16 \
    --phase-memory "standard" \
    --ligate-cpu 4 \
    --ligate-memory "standard" \
    --merge-vcf-cpu 16 \
    --merge-vcf-storage "20Gi" \
    --sample-manifest $SAMPLE_MANIFEST \
    --sample-id-col "entity:sample_id" \
    --cram-path-col "genome_cram_path" \
    --cram-index-path-col "genome_crai_path" \
    --output-file $OUTPUT_FILE \
    --binary-reference-file-regex $FILE_REGEX \
    --reference-dir $BINARY_REFERENCE_DIR \
    --chunk-info-dir $CHUNK_INFO_DIR \
    --gcs-requester-pays-configuration $BILLING_PROJECT \
    --save-checkpoints \
    --contig "chr22" \
    --n-samples 350
