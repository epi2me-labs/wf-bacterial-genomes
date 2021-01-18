# Containered Workflow template

This repository comprises a template for constructing a workflow container
using micromamba to install most software. It handles user permissions
and establishes a basic IO contract.

CI scripts are included for building and testing. Example data 
should be included in `test_data/` to demonstrate running of
the container.


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
