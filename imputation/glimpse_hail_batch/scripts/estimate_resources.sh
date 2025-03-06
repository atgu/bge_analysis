#!/bin/bash

BILLING_PROJECT="neale-pumas-bge"

#CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks-1000G_HGDP_no_singletons/"
#BINARY_REFERENCE_DIR="gs://broad-dsde-methods-bge-resources-public/GlimpseImputation/ReferencePanels/1000G_HGDP_no_singletons/chunks/cm_4/boost_1_78_0/bin/"
#FILE_REGEX="reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"

CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/test-regenerate-chunks/"
BINARY_REFERENCE_DIR="gs://bge-dragen-imputation/hgdp1kg/binary_reference/"
FILE_REGEX="ref_chunk_(?P<contig>.+)_(?P<chunk_index>\d+).bin"
CHUNK_FILE_REGEX="chunks_(?P<contig>.+)\.txt"

#CHUNK_INFO_DIR="gs://jigold-batch-tmp-ezxyx/chunk_info/"
#BINARY_REFERENCE_DIR="gs://jigold-batch-tmp-ezxyx/reference_data/"
#FILE_REGEX="reference_panel_contigindex_(?P<contig_index>\d+)_chunkindex_(?P<chunk_index>\d+)_(?P<contig>.*)_(?P<start>\d+)_(?P<end>\d+).bin"


python3 -m glimpse_hail_batch.resource_estimator \
    --binary-reference-file-regex $FILE_REGEX \
    --reference-dir $BINARY_REFERENCE_DIR \
    --chunk-file-regex $CHUNK_FILE_REGEX \
    --chunk-info-dir $CHUNK_INFO_DIR \
    --n-samples 30000 \
    --gcs-requester-pays-configuration $BILLING_PROJECT
