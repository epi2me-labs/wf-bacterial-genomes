# Bacterial SNP calling Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
performing haploid variant calling of whole genome data with
[medaka](https://www.github.com/nanoporetech/medaka) from basecalls and a
reference file. The workflow can optionally run
[prokka](https://github.com/tseemann/prokka) to annotate the resulting
consensus sequence.

> The pipeline is currently functional but contains little
> configuration of minimap2, racon, and medaka beyond setting the
> number of compute threads to use.

## Quickstart

### Running the workflow

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-hap-snps --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a consensus sequence scaffolded from a provided reference sequence,
* a [VCF]() file containing variants in the sample compared to the provided reference,
* an HTML report document detailing QC metrics and the primary findings of the workflow,
* (optionally) an annotation of the consensus sequence using prokka.

**Running the workflow with Docker containers**

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

```
# run the pipeline with the test data
OUTPUT=snp_calling
nextflow run epi2me-labs/wf-hap-snps \
    -w ${OUTPUT}/workspace 
    -profile standard
    --fastq test_data/subset.fa.gz --reference test_data/reference.subseq.fa.gz 
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
nextflow run epi2me-labs/wf-hap-snps \
    -w ${OUTPUT}/workspace 
    -profile conda
    --fastq test_data/subset.fa.gz --reference test_data/reference.subseq.fa.gz 
    --threads 4 --out_dir ${OUTPUT}
```

This will create a conda environment with all required software within the
workspace directory. When running multiple analyses on distinct datasets
it may not be desirable to have Nextflow create a conda environment for each
analysis. To avoid the situation editing the file `nextflow.config` will
be necessary. Search for the term `cacheDir` and set this to a directory
where you wish the conda environment to be placed.


## Development

The following notes are for developers only.

### Building the container

The Docker container image can be built with the following command:
```bash
# build the container
CONTAINER_TAG=ontresearch/bacterial-snps:latest
docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:latest \
    .
```

The `BASEIMAGE` argument here can be changed to use an alternative image.


## Useful links

* [medaka](https://www.github.com/nanoporetech/medaka)
* [prokka](https://github.com/tseemann/prokka)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
