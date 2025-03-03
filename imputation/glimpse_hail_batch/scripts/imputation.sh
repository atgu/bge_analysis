#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"
REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/test-glimpse/batch/"

#CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks-1000G_HGDP_no_singletons/"
#BINARY_REFERENCE_DIR="gs://broad-dsde-methods-bge-resources-public/GlimpseImputation/ReferencePanels/1000G_HGDP_no_singletons/chunks/cm_4/boost_1_78_0/bin/"
#FILE_REGEX="reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"

CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks/"
BINARY_REFERENCE_DIR="gs://bge-dragen-imputation/hgdp1kg/binary_reference/"
FILE_REGEX="ref_chunk_(?P<contig>.*)_(?P<chunk_index>\d+).bin"
SAMPLE_GROUP_SIZE=300
PHASE_CPU=16

#SAMPLE_MANIFEST="/Users/jigold/Downloads/Vietnam_SCZ_BGE_Cases.sample_manifest.tsv"
SAMPLE_MANIFEST="/Users/jigold/Downloads/gcp_pumas_scz.sample_manifest.tsv"

#STAGING_REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/staging/Vietnam_SCZ_BGE_Cases-test/"
STAGING_REMOTE_TMPDIR="gs://jigold-batch-tmp-ezxyx/staging/test_50/"

FASTA="gs://jigold-batch-tmp-ezxyx/fasta/Homo_sapiens_assembly38.fasta"

#OUTPUT_FILE="gs://jigold-batch-tmp-ezxyx/glimpse-outputs/Vietnam_SCZ_BGE_Cases.mt"
OUTPUT_FILE="gs://jigold-batch-tmp-ezxyx/glimpse-outputs/test_50.mt"


python3 -m glimpse_hail_batch.imputation \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions "us-central1" \
    --batch-name "glimpse-test-50" \
    --docker-hail "hailgenetics/hail:0.2.133" \
    --docker-glimpse "us-central1-docker.pkg.dev/neale-pumas-bge/glimpse2/glimpse2-gcloud:odelaneau_bd93ade" \
    --fasta $FASTA \
    --staging-remote-tmpdir $STAGING_REMOTE_TMPDIR \
    --samples-per-copy-group 100 \
    --sample-group-size $SAMPLE_GROUP_SIZE \
    --phase-cpu $PHASE_CPU \
    --phase-memory "standard" \
    --ligate-cpu 4 \
    --ligate-memory "standard" \
    --merge-vcf-cpu 4 \
    --merge-vcf-memory "standard" \
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
    --use-checkpoints \
    --contig "chr22" \
    --n-samples 50
