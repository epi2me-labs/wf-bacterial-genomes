# Bacterial SNP calling Workflow

TODO: Fill me in

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

docker run \
    --user $(id -u):$(id -g) --group-add 100 
    -v $(pwd)/test_data/:/input/ -v $(pwd):/output/ 
    ${CONTAINER_TAG}:latest
    -w /output/snp_calling/workspace \
    --reads /input/subset.fa.gz \
    --reference /input/reference.subseq.fa.gz \
    --threads 4 \
    --out_dir /output/snp_calling
```

The output of the pipeline will be found in `./snp_calling` for the above
example. More generically the output will be in `<host_path>/<output_dir>`
where `host_path` is the path corresponding to the `/output/` mount, and
`output_dir` is that given as the `--output_dir` argument as above.


## Useful links

* TODO: Link to medaka
