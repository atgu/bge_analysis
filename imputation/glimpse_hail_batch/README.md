# GLIMPSE Implementation in Hail Batch

## Overview

This repo contains a set of Python scripts to run GLIMPSE2 in Hail Batch designed for data stored in the Terra Data Repository. It is intended to be the same implementation as the WDL pipeline developed by Michael Gatzen and Chris Kachulis available [here](https://github.com/broadinstitute/palantir-workflows/tree/main/GlimpseImputationPipeline). The workflow is as follows:

1. Build Docker Images
2. Register Hail Batch Service Account in Terra
3. Chunk and Split Reference Panel (if not already done)
4. Estimate the Optimal Machine Type and Sample Group Size
5. Run Imputation
    - Run a first pass through all the data
    - Rerun any failed jobs
6. Cleanup Temporary Files

## Differences Between Hail Batch and WDL Implementation

The Hail Batch implemented in this repo is almost identical to that written in the WDL implementation. However, there are some notable differences:

- CRAM files are copied from their original source to a temporary bucket. These files are then mounted to a Batch job using `gcsfuse`. This is to avoid persistent 429 errors in GLIMPSE2 when there are too many concurrent jobs. The htslib reader in GLIMPSE2 seems to not retry 429 errors. Instead, we copy the CRAM files to a temporary bucket using `gcloud storage cp` and then delete the files when a group of samples have completed.
- We use Hail to locally merge all the ligated VCFs per sample group into a temporary MatrixTable rather than implementing a way to chunk and aggregate the VCFs and merge them using a combination of GATK and bcftools. Because this step is run on a single VM, the resources for this job must be large to get enough parallelization and have enough memory to do the merging.
- We use Hail Query on Batch to merge all the merged MatrixTables per sample group into a final merged MatrixTable or VCF file.
- As of the time of writing, this pipeline doesn’t automatically merge coverage statistics.
- Sample QC is not run. It’s up to the user to run this Hail method on the final MatrixTable or VCF file.
- There is no automatic retry of jobs that run out of memory or fail for transient reasons in a single pipeline run. Instead, we use checkpointing and rerun the script with higher memory and CPU requirements for the jobs that failed.

## Considerations

All CRAM files being read either in the Terra Data Repository or a bucket managed outside the Terra Data Repository must be in a **Regional** bucket. The temporary buckets used for both Batch temporary files and the staging remote temporary directory should be located in a **Regional** bucket in the same region as the CRAMs bucket(s). Lastly, the region specified to use with Batch `--batch-regions` must also be in the same region as the storage buckets. For example, to be clear, the following six items must be in the same region such as `us-central1`:

1. CRAM and CRAM index files bucket
2. Batch remote temporary directory `--batch-remote-tmpdir`
3. Staging remote temporary directory `--staging-remote-tmpdir`
4. Regions to use with Batch `--batch-regions`
5. Bucket with the split reference panel `--split-reference-dir`
6. Bucket containing the FASTA file `--fasta`

## Installation

```
git clone https://github.com/atgu/bge_analysis.git
cd bge_analysis/imputation/
pip3 install glimpse_hail_batch/
```

## Docker

We need to build a modified GLIMPSE image that contains `gcloud` and is based on the GLIMPSE image from Chris and Michael with checkpointing implemented `us.gcr.io/broad-dsde-methods/glimpse:odelaneau_bd93ade`. We need to modify the Docker image because in Hail Batch, when a job is preempted, it does not get to reuse the underlying storage. Therefore, we need to periodically upload the checkpoint file to the expected remote temporary file location so it can be downloaded when the job restarts.

If you have not created a "glimpse2" repository in your [Google Artifact Registry](https://cloud.google.com/artifact-registry/docs/repositories/create-repos), then you must do so first.

```
cd glimpse_hail_batch/
docker build -t us-central1-docker.pkg.dev/<MY_PROJECT>/glimpse2/glimpse2-gcloud:odelaneau_bd93ade .
docker push us-central1-docker.pkg.dev/<MY_PROJECT>/glimpse2/glimpse2-gcloud:odelaneau_bd93ade
```

We also need to build a modified version of Hail with unpublished features including job groups and being able to submit QoB jobs into an existing batch.

```
docker build -t us-central1-docker.pkg.dev/<MY_PROJECT>/glimpse2/hail-with-job-groups:0.0.1 -f Dockerfile.hail-job-groups
docker push us-central1-docker.pkg.dev/<MY_PROJECT>/glimpse2/hail-with-job-groups:0.0.1
```

## Register Hail Batch Service Account with Terra

If your data is in the Terra Data Repository, then you need to register your Hail Batch Service Account with Terra in order to access the data. You can do that by running the following script (only needs to be done once):

```
% python3 -m glimpse_hail_batch.add_batch_service_account_to_terra --help
usage: add_batch_service_account_to_terra.py [-h] --email EMAIL --billing-project BILLING_PROJECT --regions REGIONS --remote-tmpdir REMOTE_TMPDIR --first-name FIRST_NAME --last-name LAST_NAME

options:
  -h, --help            show this help message and exit
  --email EMAIL
  --billing-project BILLING_PROJECT
  --regions REGIONS
  --remote-tmpdir REMOTE_TMPDIR
  --first-name FIRST_NAME
  --last-name LAST_NAME
```


## Split Reference Panel

The reference panel must be chunked and split. This step must be done once per reference panel. The reference panel that was used in testing was located here (`gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/`). The code uses a TSV file (`--regions-file`) with three fields to determine the contigs, the reference file for the contig, and the location of the corresponding genetic map file. See `data/split_reference_panel_genetic_map_inputs.tsv` for an example.
`--region` can be used to run the GLIMPSE split_reference for only a specific contig. `gcs_bucket_allow_list` is necessary 

```
% python3 -m glimpse_hail_batch.split_reference --help
usage: split_reference.py [-h] [--billing-project BILLING_PROJECT] [--remote-tmpdir REMOTE_TMPDIR] [--batch-regions BATCH_REGIONS] --regions-file REGIONS_FILE [--region REGION] [--seed SEED]
                          [--min-window-cm MIN_WINDOW_CM] [--uniform-number-of-variants] [--keep-monomorphic-ref-sites] --docker DOCKER [--cpu CPU] [--memory MEMORY] --output-dir OUTPUT_DIR
                          [--batch-name BATCH_NAME]

options:
  -h, --help            show this help message and exit
  --billing-project BILLING_PROJECT
  --remote-tmpdir REMOTE_TMPDIR
  --batch-regions BATCH_REGIONS
  --regions-file REGIONS_FILE
  --region REGION
  --seed SEED
  --min-window-cm MIN_WINDOW_CM
  --uniform-number-of-variants
  --keep-monomorphic-ref-sites
  --docker DOCKER
  --cpu CPU
  --memory MEMORY
  --output-dir OUTPUT_DIR
  --batch-name BATCH_NAME
```

An example script for invoking this command with example input parameters is at `scripts/split_reference.sh`.

## Estimate Resources

We provide a utility script to estimate the optimal sample group size given a specific machine time and then the corresponding runtimes and per sample cost. This script uses the chunk compositions of the reference panel to estimate the runtime and memory (same calculation from Michael and Chris) given the sample size and chunk compositions. For example, for chromosome 22 in the HGDP reference panel, the optimal sample group size for an 8 core machine is ~350. The script also outputs the average memory requirement and the expected cost per sample. Depending on the memory requirement (GB/core), you want “standard” if this number is less than 4 or “highmem” if this number is less than 8. If it’s over 8 GB, then you will want to use a smaller batch size with a higher CPU to make sure there is enough memory. Note, the maximum machine size in this script is 16 cores. If you run into problems, we can add more functionality to put the jobs on a private machine with more cores.

```
% python3 -m glimpse_hail_batch.resource_estimator --help
usage: resource_estimator.py [-h] --reference-dir REFERENCE_DIR --chunk-info-dir CHUNK_INFO_DIR --binary-reference-file-regex BINARY_REFERENCE_FILE_REGEX --chunk-file-regex CHUNK_FILE_REGEX --n-samples N_SAMPLES
                             [--gcs-requester-pays-configuration GCS_REQUESTER_PAYS_CONFIGURATION]

options:
  -h, --help            show this help message and exit
  --reference-dir REFERENCE_DIR
  --chunk-info-dir CHUNK_INFO_DIR
  --binary-reference-file-regex BINARY_REFERENCE_FILE_REGEX
  --chunk-file-regex CHUNK_FILE_REGEX
  --n-samples N_SAMPLES
  --gcs-requester-pays-configuration GCS_REQUESTER_PAYS_CONFIGURATION
```

An example of running this script is at `scripts/estimate_resources.sh`

*Example Output for Chromosome 22:*

```
  cores    n_samples    goal_runtime    max_runtime    optimal_batch_size    max_batch_size    actual_batch_size    max_memory    mean_memory    rough_cost_estimate    rough_cost_estimate_per_sample
-------  -----------  --------------  -------------  --------------------  ----------------  -------------------  ------------  -------------  ---------------------  --------------------------------
      1        30000              30             60               43.9098           58.3334              43.9098        4               3                     164.16                          0.005472
      2        30000              30             60               87.8197          116.667               87.8197        2.5             2                     164.16                          0.005472
      4        30000              30             60              175.639           233.334              175.639         1.75            1.5                   164.16                          0.005472
      8        30000              30             60              351.279           466.667              351.279         1.625           1.375                 165.12                          0.005504
     16        30000              30             60              702.557           933.335              702.557         1.4375          1.25                  165.12                          0.005504
```


## Imputation

Imputation in GLIMPSE runs in three steps:

1. Phase samples for a specific chunk of the genome per sample group.
2. Ligate all the chunks for a specific chromosome per sample group.
3. Union all the ligated chromosomes for all sample groups and use Query on Batch to merge them together into one MatrixTable or VCF file.

Samples are split into "Sample Groups" in order to increase parallelism and reduce memory requirements. All sample groups are merged in a final step that unions all the MatrixTable(s) in Hail Query on Batch into either a MatrixTable or compressed VCF file specified with `--output-file` by chromosome.

The Python interface for running imputation is as follows:

```
% python3 -m glimpse_hail_batch.imputation.submit --help
usage: submit.py [-h] [--billing-project BILLING_PROJECT] [--batch-remote-tmpdir BATCH_REMOTE_TMPDIR] [--batch-regions BATCH_REGIONS] [--batch-name BATCH_NAME] [--batch-id BATCH_ID] [--seed SEED] --docker-glimpse DOCKER_GLIMPSE --docker-hail
                 DOCKER_HAIL --reference-dir REFERENCE_DIR --chunk-info-dir CHUNK_INFO_DIR --binary-reference-file-regex BINARY_REFERENCE_FILE_REGEX --chunk-file-regex CHUNK_FILE_REGEX --sample-manifest SAMPLE_MANIFEST --sample-id-col SAMPLE_ID_COL
                 --cram-path-col CRAM_PATH_COL --cram-index-path-col CRAM_INDEX_PATH_COL [--sex-col SEX_COL] [--female-code FEMALE_CODE] [--n-samples N_SAMPLES] --sample-group-size SAMPLE_GROUP_SIZE [--samples-per-copy-group SAMPLES_PER_COPY_GROUP]
                 --ligate-cpu LIGATE_CPU [--ligate-memory LIGATE_MEMORY] [--ligate-storage LIGATE_STORAGE] [--ligate-ref-dict LIGATE_REF_DICT] [--contig CONTIG] [--chunk-index CHUNK_INDEX] [--sample-group-index SAMPLE_GROUP_INDEX]
                 --staging-remote-tmpdir STAGING_REMOTE_TMPDIR --output-file OUTPUT_FILE [--fasta FASTA] [--use-checkpoints] [--save-checkpoints] [--always-delete-temp-files] --phase-cpu PHASE_CPU --phase-memory PHASE_MEMORY
                 [--phase-impute-reference-only-variants] [--phase-call-indels] [--phase-n-burn-in PHASE_N_BURN_IN] [--phase-n-main PHASE_N_MAIN] [--phase-effective-population-size PHASE_EFFECTIVE_POPULATION_SIZE] --merge-vcf-cpu MERGE_VCF_CPU
                 [--merge-vcf-memory MERGE_VCF_MEMORY] [--merge-vcf-storage MERGE_VCF_STORAGE] [--gcs-requester-pays-configuration GCS_REQUESTER_PAYS_CONFIGURATION] [--non-par-contigs NON_PAR_CONTIGS]

options:
  -h, --help            show this help message and exit
  --billing-project BILLING_PROJECT
  --batch-remote-tmpdir BATCH_REMOTE_TMPDIR
  --batch-regions BATCH_REGIONS
  --batch-name BATCH_NAME
  --batch-id BATCH_ID
  --seed SEED
  --docker-glimpse DOCKER_GLIMPSE
  --docker-hail DOCKER_HAIL
  --reference-dir REFERENCE_DIR
  --chunk-info-dir CHUNK_INFO_DIR
  --binary-reference-file-regex BINARY_REFERENCE_FILE_REGEX
  --chunk-file-regex CHUNK_FILE_REGEX
  --sample-manifest SAMPLE_MANIFEST
  --sample-id-col SAMPLE_ID_COL
  --cram-path-col CRAM_PATH_COL
  --cram-index-path-col CRAM_INDEX_PATH_COL
  --sex-col SEX_COL
  --female-code FEMALE_CODE
  --n-samples N_SAMPLES
  --sample-group-size SAMPLE_GROUP_SIZE
  --samples-per-copy-group SAMPLES_PER_COPY_GROUP
  --ligate-cpu LIGATE_CPU
  --ligate-memory LIGATE_MEMORY
  --ligate-storage LIGATE_STORAGE
  --ligate-ref-dict LIGATE_REF_DICT
  --contig CONTIG
  --chunk-index CHUNK_INDEX
  --sample-group-index SAMPLE_GROUP_INDEX
  --staging-remote-tmpdir STAGING_REMOTE_TMPDIR
  --output-file OUTPUT_FILE
  --fasta FASTA
  --use-checkpoints
  --save-checkpoints
  --always-delete-temp-files
  --phase-cpu PHASE_CPU
  --phase-memory PHASE_MEMORY
  --phase-impute-reference-only-variants
  --phase-call-indels
  --phase-n-burn-in PHASE_N_BURN_IN
  --phase-n-main PHASE_N_MAIN
  --phase-effective-population-size PHASE_EFFECTIVE_POPULATION_SIZE
  --merge-vcf-cpu MERGE_VCF_CPU
  --merge-vcf-memory MERGE_VCF_MEMORY
  --merge-vcf-storage MERGE_VCF_STORAGE
  --gcs-requester-pays-configuration GCS_REQUESTER_PAYS_CONFIGURATION
  --non-par-contigs NON_PAR_CONTIGS
```

An example script for invoking this command with example input parameters is at `scripts/imputation.sh`.

### Submission

When you run the imputation script, it first creates a batch on your local computer and then submits one job to that batch.
That job then will submit jobs after determining which checkpoints exist if requested (see below). This setup was chosen in
order to use an unpublished Hail feature (Job Groups) and accommodate large analyses that might take a long time to submit.

### Docker Images

There are two Docker images required:

1. The image built above with `GLIMPSE` and `gcloud` specified by `--docker-glimpse`
2. An image built above containing the custom version of Hail specified by `--docker-hail`.

### Selecting a Subset of the Data

You can use the `--n-samples` option to limit the number of samples. This is useful for testing on a subset of the dataset or retrying a specific part of the dataset with different compute resources. You can use the `--contig` and optionally a `--chunk-index` argument to limit the contig or even to a specific chunk. You can use the `--sample-group-index` to subset to a specific group of samples.

### Controlling Copy Parallelism

If you increase the sample group size (`--sample-group-size`), you'll want to increase the `--samples-per-copy-group` parameter. This is because we copy the CRAM files into the staging bucket in parallel. The more parallelism here, the higher the chance of transient errors.

### Checkpointing

- For large workflows, there are bound to be jobs that fail especially due to out of memory errors. To be able to repair the imputation process, you'll want to use the `--save-checkpoints` option when running the pipeline (including the initial run!).
- When running the same pipeline again (repair failed jobs), use the `--use-checkpoints` option to avoid rerunning phasing, ligate jobs, and merge VCF jobs that have already succeeded for a given sample group.

### Clean Up

CRAM files are always deleted every time a pipeline is run (not possible to checkpoint these files due to cost concerns). **Make sure you check they were deleted by making sure the *delete* job for each sample group succeeded and the size of the staging bucket is 0 for the CRAM files (see below).

If you want to delete all temporary files (phase, ligate, merged VCF outputs), use the `--always-delete-temp-files` option. This will clean up all sample groups even if they have previously succeeded.

Make sure when you are done and while the pipeline is running that you check the size of your `--staging-remote-tmpdir` to make sure files are being deleted:

```
gcloud storage du -s --total --readable-sizes "gs://<MY_BUCKET>/staging/*"
```

To specifically look at the **CRAM** file sizes:

```
gcloud storage du -s --total --readable-sizes "gs://<MY_BUCKET>/staging/sg-*/crams/*"
```

To specifically look at the **phase** file sizes:

```
gcloud storage du -s --total --readable-sizes "gs://<MY_BUCKET>/staging/sg-*/phase/*/chunk-*/*"
```

To specifically look at the **ligate** file sizes:

```
gcloud storage du -s --total --readable-sizes "gs://<MY_BUCKET>/staging/sg-*/ligate/*/*"
```

## Monitoring Progress

You can use this script to show a live view of imputation progress by supplying the batch ID of interest.

```
% python3 -m glimpse_hail_batch.monitor --help
usage: monitor.py [-h] --batch-id BATCH_ID

options:
  -h, --help           show this help message and exit
  --batch-id BATCH_ID
```

## Adding IMPUTE Info Scores

You can run a post-processing step to add annotations to an existing MatrixTable or VCF with GT and GP defined.
Given a sample annotations TSV file and a list of columns to group by, this script will add annotations for the IMPUTE INFO
score and allele frequency for each category of the grouped column. For example, if grouping by "reported_sex", you would
get the INFO score and AF for both Males and Females.

```
% python3 -m glimpse_hail_batch.add_impute_info_scores --help
usage: add_impute_info_scores.py [-h] --input-path INPUT_PATH --output-path OUTPUT_PATH --sample-annotations SAMPLE_ANNOTATIONS [--sample-id-col-input SAMPLE_ID_COL_INPUT] --sample-id-col-ann SAMPLE_ID_COL_ANN --group-by-col GROUP_BY_COL
                                 [GROUP_BY_COL ...] [--overwrite] [--hl-init-kwarg [INIT_KWARGS ...]]

options:
  -h, --help            show this help message and exit
  --input-path INPUT_PATH
  --output-path OUTPUT_PATH
  --sample-annotations SAMPLE_ANNOTATIONS
  --sample-id-col-input SAMPLE_ID_COL_INPUT
  --sample-id-col-ann SAMPLE_ID_COL_ANN
  --group-by-col GROUP_BY_COL [GROUP_BY_COL ...]
  --overwrite
  --hl-init-kwarg [INIT_KWARGS ...]
```

An example version of running this script is at `scripts/add_impute_info_scores.sh`.
