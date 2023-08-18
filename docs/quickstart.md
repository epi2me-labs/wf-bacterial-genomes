## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) (default) or
[Singularity](https://sylabs.io/singularity/) (`-profile singularity`) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-bacterial-genomes --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a [FASTA](https://en.wikipedia.org/wiki/FASTA) consensus sequence scaffolded from a provided reference sequence,
* a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file containing variants in the sample compared to the reference (if provided),
* an HTML report document detailing QC metrics and the primary findings of the workflow,
* (optionally) an annotation of the consensus sequence using prokka.
* (optionally) a per-sample ResFinder output directory with various results.