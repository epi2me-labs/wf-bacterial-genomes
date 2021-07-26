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

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-hap-snps --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a [FASTA](https://en.wikipedia.org/wiki/FASTA) consensus sequence scaffolded from a provided reference sequence,
* a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file containing variants in the sample compared to the provided reference,
* an HTML report document detailing QC metrics and the primary findings of the workflow,
* (optionally) an annotation of the consensus sequence using prokka.


## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
