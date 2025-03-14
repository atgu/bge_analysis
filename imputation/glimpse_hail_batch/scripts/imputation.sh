#!/bin/bash

BILLING_PROJECT="my_billing_project"

# Make this your own temporary bucket and a path you like; your batch service account needs read/write access
REMOTE_TMPDIR="gs://MY_BUCKET/test-glimpse/batch/"

BATCH_REGIONS="us-central1"

# you need to make sure these images exist
DOCKER_HAIL="us-central1-docker.pkg.dev/MY_PROJECT/glimpse2/hail-with-job-groups:0.0.4"
DOCKER_GLIMPSE="us-central1-docker.pkg.dev/MY_PROJECT/glimpse2/glimpse2-gcloud:odelaneau_bd93ade"

# change these only with a different reference panel
FASTA="gs://MY_FASTA_BUCKET/fasta/Homo_sapiens_assembly38.fasta"
CHUNK_INFO_DIR="gs://MY_BUCKET/reference_panel/chunks/"
BINARY_REFERENCE_DIR="gs://MY_BUCKET/reference_panel/binary_reference/"

BINARY_REFERENCE_FILE_REGEX="ref_chunk_(?P<contig>.+)_(?P<chunk_index>\d+).bin"
CHUNK_FILE_REGEX="chunks_(?P<contig>.+)\.txt"

NON_PAR_CONTIGS="chrX_nonPAR"

# Make sure all of these are correct for your sample manifest (just a tsv file) and the female code
SAMPLE_MANIFEST="sample_manifest.tsv"
SAMPLE_ID_COL="entity:sample_id"
CRAM_PATH_COL="genome_cram_path"
CRAM_INDEX_PATH_COL="genome_crai_path"
SEX_COL="reported_sex"
FEMALE_CODE="Female"

# change this to a name that makes sense to you for what you are running
BATCH_NAME="my_project_name"

# Make this your own temporary bucket and a path you like; your batch service account needs read/write access
STAGING_REMOTE_TMPDIR="gs://MY_BUCKET/staging/${BATCH_NAME}/"

# This should probably be in bge-dragen-imputation unless you want to copy the finished matrix table over later
# Notice the jinja templating with the {{ }}. My code fills in that file path with the contig name.
OUTPUT_FILE="gs://MY_BUCKET/glimpse-outputs/my_project_name/my_project_name.{{contig}}.mt"

SAMPLE_GROUP_SIZE=750
SAMPLES_PER_COPY_GROUP=300

PHASE_CPU=16
LIGATE_CPU=4
MERGE_VCF_CPU=2

PHASE_MEMORY="standard"
LIGATE_MEMORY="standard"
MERGE_VCF_MEMORY="standard"

# I'd always have --save-checkpoints on, but you will need to remember to delete the temporary files when you are done
# Use_checkpoints checks whether each file exists; will be faster to submit the first time you run the pipeline without this turned on
SAVE_CHECKPOINTS="--save-checkpoints"
USE_CHECKPOINTS="--use-checkpoints"

# Leave these as empty strings for no filters
# My suggestion is to start with this small test set to make sure everything runs before running the full thing
CONTIG="--contig chr22"
CHUNK_INDEX=""
N_SAMPLES="--n-samples 30"

python3 -m glimpse_hail_batch.imputation.submit \
    --gcs-requester-pays-configuration $BILLING_PROJECT \
    --billing-project $BILLING_PROJECT \
    --batch-remote-tmpdir $REMOTE_TMPDIR \
    --batch-regions $BATCH_REGIONS \
    --batch-name $BATCH_NAME \
    --docker-hail $DOCKER_HAIL \
    --docker-glimpse $DOCKER_GLIMPSE \
    --fasta $FASTA \
    --staging-remote-tmpdir $STAGING_REMOTE_TMPDIR \
    --samples-per-copy-group $SAMPLES_PER_COPY_GROUP \
    --sample-group-size $SAMPLE_GROUP_SIZE \
    --phase-cpu $PHASE_CPU \
    --phase-memory $PHASE_MEMORY \
    --ligate-cpu $LIGATE_CPU \
    --ligate-memory $LIGATE_MEMORY \
    --merge-vcf-cpu $MERGE_VCF_CPU \
    --merge-vcf-memory $MERGE_VCF_MEMORY \
    --sample-manifest $SAMPLE_MANIFEST \
    --sample-id-col $SAMPLE_ID_COL \
    --cram-path-col $CRAM_PATH_COL \
    --cram-index-path-col $CRAM_INDEX_PATH_COL \
    --sex-col $SEX_COL \
    --female-code $FEMALE_CODE \
    --output-file $OUTPUT_FILE \
    --binary-reference-file-regex $BINARY_REFERENCE_FILE_REGEX \
    --chunk-file-regex $CHUNK_FILE_REGEX \
    --reference-dir $BINARY_REFERENCE_DIR \
    --chunk-info-dir $CHUNK_INFO_DIR \
    --non-par-contigs $NON_PAR_CONTIGS \
    $SAVE_CHECKPOINTS \
    $USE_CHECKPOINTS \
    $CONTIG \
    $CHUNK_INDEX \
    $N_SAMPLES
