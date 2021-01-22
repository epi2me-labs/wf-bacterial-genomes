# Bacterial SNP calling Workflow

This repository contains a nextflow workflow and associated Docker
container build for performing haploid variant calling with medaka
from basecalls and a reference file.

> The pipeline is currently functional but contains little
> configuration of minimap2, racon, and medaka beyond setting the
> number of compute threads to use.

## Quickstart

```bash
# build the container
CONTAINER_TAG=bacterial-snps
docker build -t ${CONTAINER_TAG} -f Dockerfile  .

# run the pipeline with the test data
nextflow run workflow.nf \
    -w snp_calling_docker/workspace 
    -profile withdocker
    --reads test_data/subset.fa.gz --reference test_data/reference.subseq.fa.gz 
    --threads 4 --out_dir snp_calling
```

The output of the pipeline will be found in `./snp_calling` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline: a `medaka_consensus.fasta` file and a
`medaka_consensus.vcf` file.


## Useful links

* [medaka](https://www.github.com/nanoporetech/medaka)
