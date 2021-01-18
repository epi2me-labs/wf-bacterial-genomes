# Containered Workflow template

This repository comprises a template for constructing a workflow container
using micromamba to install most software. It handles user permissions
and establishes a basic IO contract.

CI scripts are included for building and testing. Example data 
should be included in `test_data/` to demonstrate running of
the container.

## Workflow Developer Instructions

1. Edit the Docker file

   For simple cases the only thing to edit in the Dockerfile will be the dependency
   installation step and copying the necessary workflow code.

2. The Dockerfile expects an command-line program to be installed to

       /home/epi2melabs/conda/bin/run_workflow

   with a `--help` option. The example `run_workflow` script demonstrates a way
   to run a `snakemake` workflow which conforms to the (horrible) `snakemake`
   standard of maintaining workflow code alongside workflow outputs. The entrypoint
   script should save its outputs to `/output/<output_label>` for consistency
   across workflows.

3. The CI file will need editing to show usage of the container image with
   the sample data.

4. The ***Quickstart*** section below should be edited to make sense for the
   implemented workflow. This should amount to a copy/paste from the CI file.


## Quickstart

```bash
# build the container
CONTAINER_TAG=epi2melabs-<name>
docker build -t ${CONTAINER_TAG} -f Dockerfile  .

# run the pipeline with the test data
docker run \
    --user $(id -u):$(id -g) --group-add 100 \
    -v $(pwd)/test_data/:/input/ -v $(pwd):/output/ \
    ${CONTAINER_TAG}:latest \
    /input/<input_file1.fa.gz> \
    /input/<input_file2.fa.gz> \
    --output_label <my_analysis>
```

The output of the pipeline will be found in `./<my_analysis>` for the above
example. More generically the output will be in `<host_path>/<output_label>`
where `host_path` is the path corresponding to the `/output/` mount, and
`output_label` is that given as the `--output_label` argument as above.


## Useful links

* Include useful links to e.g. the source workflow
