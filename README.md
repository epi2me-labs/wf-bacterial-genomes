# Bacterial assembly and annotation workflow

Assembly, variant calling, and annotation of bacterial genomes.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

This workflow is primarily used to assemble genomes from bacterial reads and provide information on features of interest within those assemblies through annotations.

The workflow can provide additional information about the assembly, such as antimicrobial resistance (AMR) analysis and sequence typing through an optional `--isolates` mode. 

In brief, this workflow will perform the following: 

+ De novo (or reference-based) assembly of bacterial genomes 
+ Annotation of regions of interest within the assembly
+ Species identification and sequence typing (`--isolates` mode only)
+ Identify genes and SNVs associated with AMR (`--isolates` mode only)



## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 32GB

Minimum requirements:

+ CPUs = 4
+ Memory = 16GB

Approximate run time: 20-40 minutes per sample with ~50x coverage using minimum requirements

ARM processor support: True




## Install and run

<!---Nextflow text remains the same across workflows, update example cmd and demo data sections.--->

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-bacterial-genomes –-help
```
A demo dataset is provided for testing of the workflow. It can be downloaded using:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-bacterial-genomes/wf-bacterial-genomes-demo.tar.gz
tar -xzvf wf-bacterial-genomes-demo.tar.gz
```
The workflow can be run with the demo data using:
```
nextflow run epi2me-labs/wf-bacterial-genomes \
    --fastq wf-bacterial-genomes-demo/isolates_fastq \
    --isolates \
    --reference_based_assembly \
    --reference wf-bacterial-genomes-demo/ref/ref.fasta.gz \
    --sample_sheet wf-bacterial-genomes-demo/isolates_sample_sheet.csv \
    -profile standard
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/



## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Optimal DNA extraction will be dependent on the gram status of the organism. Some useful protocols are provided below:
+ [Gram-positive bacteria](https://community.nanoporetech.com/extraction_method_groups/gram-positive-bacterial-gnda)
+ [Gram-negative bacteria](https://community.nanoporetech.com/extraction_methods/gram-ve-dna)


Find more related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| reference_based_assembly | boolean | Enable reference guided assembly instead of de novo assembly. | By default de novo assembly will be performed with Flye. Enable this to instead perform a reference-based consensus. A reference must be provided. | False |
| reference | string | Reference sequence FASTA file. | The reference sequence is used when performing reference-based assembly. |  |
| basecaller_cfg | string | Name of the model that was used to basecall signal data, used to select an appropriate Medaka model. | The basecaller configuration is used to automatically select the appropriate Medaka model. The automatic selection can be overridden with the `medaka_variant_model` and `medaka_consensus_model` parameters. The model list only shows models that are compatible with this workflow. | dna_r10.4.1_e8.2_400bps_sup@v4.2.0 |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Isolate options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| isolates | boolean | Run the Isolates pipeline on the assembly results if set to True. | Isolates mode adds further analysis options to the workflow such as multi-locus sequence typing and antimicrobial resistance calling, as well as producing single reports for each sample in the run. | False |
| resfinder_version | string | ResFinder version to use. | ResFinder is the tool used to check for antimicrobial resistance genes in isolates of bacteria. | 4.3.2 |
| resfinder_threshold | string | Threshold of required identity to report a match between a gene in the ResFinder database and the assembly. Valid interval: 0.00-1.00 | Identity refers to the ratio of base pairs that match between the sequence in your assembly and that of the sequence in the ResFinder database. Increasing the threshold will results in fewer, but more accurate hits against the database. | 0.8 |
| resfinder_coverage | string | Minimum coverage (breadth-of) threshold required to report a match between a gene in the ResFinder database and the assembly. Valid interval: 0.00-1.00 | The amount of an AMR gene that has to be present within the assembly as compared to the reference in the ResFinder database. | 0.6 |
| mlst_version | string | MLST version to use. |  | 2.23.0 |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| run_prokka | boolean | Run prokka on consensus sequence | Will provide an output file with a list of annotations for your sequence. Optional because it can take some time. | True |
| prokka_opts | string | Command-line arguments for prokka | [Command line arguments](https://github.com/tseemann/prokka#command-line-options) which can be used to alter prokka output annotation files. |  |
| flye_genome_size | integer | Estimated genome size for de novo assembly in non-SI prefix format (e.g 5000000 for 5 Mb genome) | This setting is used in conjunction with `flye_asm_coverage` to subsample the reads used in the initial disjointig step only; all reads are used in subsequent steps. The values in the two parameters are used to calculate a target yield used to subsample the longest reads in the dataset, see [Flye docs](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) for more information. *Note* For runs with mixed genome sizes, preference the larger genome size. |  |
| flye_asm_coverage | integer | Target coverage to use for subsampling in de novo assembly | This setting is used in conjunction with `flye_genome_size` to subsample the reads used in the initial disjointig step only; all reads are used in subsequent steps. The values in the two parameters are used to calculate a target yield used to subsample the longest reads in the dataset, see [Flye docs](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) for more information. |  |
| flye_opts | string | Command-line arguments for flye | [Command line arguments](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) which can be used to alter the de novo assembly process. Enter the command as quoted string (e.g '--meta --iterations 2'). Flye's `--genome-size` and `--asm-coverage` parameters can be set directly in the workflow with `--flye_genome_size` and `--flye_asm_coverage`, respectively. |  |
| medaka_consensus_model | string | The name of a Medaka consensus model to use. This name will override the model automatically chosen based on the provided basecaller configuration. | The workflow will attempt to map the basecalling model used to a suitable Medaka consensus model. You can override this by providing a model with this option instead. |  |
| medaka_variant_model | string | The name of a Medaka variant model to use. This name will override the model automatically chosen based on the provided basecaller configuration. | The workflow will attempt to map the basecalling model used to a suitable Medaka variant model. You can override this by providing a model with this option instead. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads. | Provided to alignment, flye assembly and prokka steps to improve performance. | 3 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | wf-bacterial-genomes-report.html | Report for all samples | aggregated |
| Workflow results | results.json | Structured workflow results for internal/onward use. | aggregated |
| workflow checkpoints | checkpoints.json | Structured workflow checkpoints for internal/onward use. | aggregated |
| Concatenated sequence data | {{ alias }}/seqs.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |
| Per file read stats | {{ alias }}/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, per sample. | per-sample |
| Per read stats | {{ alias }}/fastcat_stats/per-read-stats.tsv.gz | A TSV with per read stats, per sample | per-sample |
| alignment | {{ alias }}.bam | Aligned reads for the sample in BAM format. | per-sample |
| alignment index | {{ alias }}.bam.bai | An index file for the alignment in BAI format. | per-sample |
| Draft assembly FASTA file | {{ alias }}.medaka.fasta.gz | Consensus file generated from either de-novo assembly or reference variant calling pipeline. | per-sample |
| Variants VCF file | {{ alias }}.medaka.vcf.gz | VCF file of variants detected against the provided reference (Reference mode only). | per-sample |
| Variants summary | {{ alias }}.variants.stats | TSV file of summary statistics for variants in sample (Reference mode only). | per-sample |
| Annotations files | {{ alias }}.prokka.gff | Annotations of regions of interest in assembly in GBK and GFF format. | per-sample |
| Sequence typing results | {{ alias }}.mlst.json | Sequence typing results in JSON format (isolates mode only). | per-sample |
| AMR calling results | {{ alias }}_resfinder_results | Resfinder results for AMR calling (isolates mode only). | per-sample |
| isolates per sample report | {{ alias }}-isolates-report.html | Per sample report isolates mode | per-sample |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2a. De-novo assembly

#### i. Assembly

[Flye](https://github.com/fenderglass/Flye) is used to create a draft assembly from the FASTQ reads. This will run by default on the `--nano-raw` paramter for flye. Additional configuration can be performed using `--flye_opts` parameter. 

#### ii. Polishing

The draft assembly from flye is then polished using [Medaka](https://github.com/nanoporetech/medaka). This step will attempt to correct any errors that were introduced during the de-novo assembly process. 

The workflow selects the appropriate [Medaka models](https://github.com/nanoporetech/medaka#models) based on the basecaller configuration that was used to process the signal data.
You can use the parameter `--basecaller_cfg` to provide this information (e.g. `dna_r10.4.1_e8.2_400bps_hac`).
Alternatively, you can choose the [Medaka](https://github.com/nanoporetech/medaka) model directly with `--medaka_consensus_model`.


### 2b. Variant calling mode

#### i. Align reads

Reads are aligned against the provided reference with [mini_align](https://github.com/nanoporetech/pomoxis/).

#### ii. Call variants

After alignment, haploid variants are called with [Medaka](https://github.com/nanoporetech/medaka).

The workflow selects the appropriate [Medaka models](https://github.com/nanoporetech/medaka#models) based on the basecaller configuration that was used to process the signal data.
You can use the parameter `--basecaller_cfg` to provide this information (e.g. `dna_r10.4.1_e8.2_400bps_hac`).
Alternatively, you can choose the [Medaka](https://github.com/nanoporetech/medaka) model directly with `--medaka_model`.

#### iii. Use the variants to generate a consensus

The variants passing the depth filter are then incorporated in the reference to create the consensus sequence. Variant stats are also created at this point.

### 3. Annotations

Regions of interest within your assembly are identified and annotated using [Prokka](https://github.com/tseemann/prokka). By default, prokka will run with it's default databases, but users can refine the annotation using the `--prokka_opts` command. **NOTE** The workflow does not current accept any additional files sent to prokka such as GBK or GFF files.

### 4. Isolates mode (optional)

#### i. Multi-locus sequence typing (MLST)

MLST is a common technique used to help characterise your bacterial isolate, by using allelic variation from internal DNA fragments of 6-7 house keeping genes. Typing schemes for specific species and genera are found on [PubMLST](https://pubmlst.org/) and are pre-loaded into this workflow. [MLST](https://github.com/tseemann/mlst) will try to infer the correct typing scheme to use by scanning the assembly and subsequently identify the allele variant found.

#### ii. Antimicrobial resistance (AMR) calling

[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) is used to identify genes/SNVs associated with AMR in your assembly. Assemblies of any species will be searched for the detection of acquired resistance genes, however SNVs conferring resistance are only available to a few well characterised species/genera. These are:
* Campylobacter spp.
* Enterococcus faecalis
* Enterococcus faecium
* Escherichia coli
* Helicobacter pylori
* Klebsiella spp. 
* Mycobacterium tuberculosis
* Neisseria gonorrhoeae
* Salmonella spp.
* Staphylococcus aureus

The species/genera of your assembly will be detected from the results of the MLST step and SNV will be selected automatically if applicable.








## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

### No results for multi-locus sequence typing?
This usually occurs if the assembly is incomplete and does not have sufficient coverage to identify the house keeping genes of the typing scheme. Another, rarer scenario is if the assembly is from an organism with no typing scheme. A list of the available typing schemes can be found [here](https://github.com/tseemann/mlst/tree/master/db/pubmlst). In both scenarios, AMR calling will still be performed but only for acquired resistance genes.



If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-bacterial-genomes/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

<!---Any other sections that are relevant specifically to this workflow and may be useful to users eg. ## Related blog posts. ## Learning center links.--->

## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



