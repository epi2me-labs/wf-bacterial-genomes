# Bacterial SNP calling Workflow

This repository contains a nextflow workflow and associated Docker
container build for performing haploid variant calling with medaka
from basecalls and a reference file.

> The pipeline is currently functional but contains little
> configuration of minimap2, racon, and medaka beyond setting the
> number of compute threads to use.

## Quickstart

### Building the container

> This step is not necessary if you intend to run the workflow using
> conda environments.

The Docker container image can be built with the following command:
```bash
# build the container
CONTAINER_TAG=template-workflow
docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=epi2melabs/base-workflow-image:latest \
    .
```

The `BASEIMAGE` argument here can be changed to use an alternative image.

### Running the workflow

The template includes a simple workflow that outputs a file with the lengths
of sequences contained in a .fastq.gz file.

**Running the workflow with Docker containers**

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

```
# run the pipeline with the test data
OUTPUT=snp_calling
nextflow run workflow.nf \
    -w ${OUTPUT}/workspace 
    -profile standard
    --reads test_data/subset.fa.gz --reference test_data/reference.subseq.fa.gz 
    --threads 4 --out_dir ${OUTPUT}/snp_calling
```

The output of the pipeline will be found in `./snp_calling` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline: a `medaka_consensus.fasta` file and a
`medaka_consensus.vcf` file.

**Using conda environments**

To run the workflow backed by conda environments, simply provide the
`-profile conda` argument to `nextflow run`.

```
# run the pipeline with the test data
OUTPUT=snp_calling
nextflow run workflow.nf \
    -w ${OUTPUT}/workspace 
    -profile conda
    --reads test_data/subset.fa.gz --reference test_data/reference.subseq.fa.gz 
    --threads 4 --out_dir ${OUTPUT}
```

This will create a conda environment with all required software within the
workspace directory. When running multiple analyses on distinct datasets
it may not be desirable to have Nextflow create a conda environment for each
analysis. To avoid the situation editing the file `nextflow.config` will
be necessary. Search for the term `cacheDir` and set this to a directory
where you wish the conda environment to be placed.

## Useful links

* [medaka](https://www.github.com/nanoporetech/medaka)
