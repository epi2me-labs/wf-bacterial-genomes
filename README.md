# Isoform Workflow

This repository contains a Dockerfile to create a container that will run the
[pipeline-nanopore-ref-isoforms](https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms)
workflow. Example data is included in `test_data/` to demonstrate running of
the container.


## Requirements

* Docker


## Quickstart

```bash
# build the container
docker build -t epi2melabs-ref-isoforms -f Dockerfile  .

# run the pipeline with the test data
docker run \
    --user $(id -u):$(id -g) --group-add 100 \
    -v $(pwd)/test_data/:/input/ -v $(pwd):/output/ \
    epi2melabs-ref-isoforms:latest \
    /input/sample_data.fastq.gz \
    /input/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz \
    /input/Drosophila_melanogaster.BDGP6.95.gtf.gz
    --output_label isoform_analysis
```

The output of the pipeline will be found in `./isoform_analysis` for the above
example. More generically the output will be in `<host_path>/<output_label>`
where `host_path` is the path corresponding to the `/output/` mount, and
`output_label` is that given as the `--output_label` argument as above.


## Useful links

* Original github repo for
  [pipeline-nanopore-ref-isoforms](https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms)

For information concerning the outputs of the workflow, please see the original
repository.
