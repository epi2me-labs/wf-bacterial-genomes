### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| reference_based_assembly | boolean | Enable reference guided assembly instead of de novo assembly. | By default de novo assembly will be performed with Flye. Enable this to instead perform a reference-based consensus. A reference must be provided. | False |
| reference | string | Reference sequence FASTA file. | The reference sequence is used when performing reference-based assembly. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Isolate options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| isolates | boolean | Run the Isolates pipeline on the assembly results if set to True. | Isolates mode adds further analysis options to the workflow such as multi-locus sequence typing and antimicrobial resistance calling, as well as producing single reports for each sample in the run. | True |
| sourmash_db_exclude_list | string | Path to file containing assembly IDs (GCA_/GCF_) to exclude from Sourmash species identification, one per line. | A text file containing assembly identifiers (starting with GCA_ or GCF_) that should be excluded from the sourmash database during species identification. Each line should contain one assembly ID. Invalid entries will be skipped. |  |
| resfinder_threshold | string | Threshold of required identity to report a match between a gene in the ResFinder database and the assembly. Valid interval: 0.00-1.00 | Identity refers to the ratio of base pairs that match between the sequence in your assembly and that of the sequence in the ResFinder database. Increasing the threshold will result in fewer, but more accurate hits against the database. | 0.8 |
| resfinder_coverage | string | Minimum coverage (breadth-of) threshold required to report a match between a gene in the ResFinder database and the assembly. Valid interval: 0.00-1.00 | The amount of an AMR gene that has to be present within the assembly as compared to the reference in the ResFinder database. | 0.6 |
| pointfinder_ignore_indels | boolean | Ignore indels when calling point mutations in PointFinder analysis. | When enabled, PointFinder will skip detection of insertion and deletion mutations that may confer antimicrobial resistance. This is useful when sequencing errors or assembly artifacts might cause false positive indel calls, but may miss genuine resistance-conferring indels. | False |
| pointfinder_ignore_stop_codons | boolean | Ignore stop codons when calling point mutations in PointFinder analysis. | When enabled, PointFinder will not report premature stop codons as resistance mutations. This can be useful to reduce noise from sequencing errors that introduce false stop codons, but may miss genuine resistance mutations that create truncated proteins. | False |
| resfinder_db | string | Path to custom ResFinder database directory. | Override the default ResFinder database with a custom one. This should point to a directory containing the ResFinder database files. If not specified, ResFinder will use its bundled ResFinder database. |  |
| pointfinder_db | string | Path to custom PointFinder database directory. | Override the default PointFinder database with a custom one. This should point to a directory containing the PointFinder database files. If not specified, ResFinder will use its bundled PointFinder database. |  |
| sourmash_pointfinder_mapping | string | CSV file mapping Sourmash species names to PointFinder species names. Must contain 'Sourmash species' and 'PointFinder species' columns. | Optional CSV file to customize the mapping between Sourmash taxonomic classifications and PointFinder species for antimicrobial resistance analysis. If not provided, default mappings will be used. The CSV should have two columns: 'Sourmash species' and 'PointFinder species'. |  |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |
| run_bakta | boolean | Run Bakta on consensus sequence | Will provide an output file with a list of annotations for your sequence. Optional because it can take some time. | True |
| bakta_opts | string | Command-line arguments for Bakta | [Command line arguments](https://github.com/oschwengers/bakta#usage) which can be used to alter Bakta output annotation files. |  |
| bakta_db | string | Path to existing Bakta database directory. If not provided, the workflow will download the database. | Provide the path to a Bakta database directory containing both Bakta and AMRFinderPlus databases. If not specified, the workflow will automatically download the database type specified by bakta_db_type. |  |
| bakta_db_type | string | Type of Bakta database to download if bakta_db is not provided. | Choose 'light' for the standard database (~1.3 GB) or 'full' for the complete database (~30 GB). Only used when bakta_db is not specified. | light |
| flye_genome_size | integer | Estimated genome size for de novo assembly in non-SI prefix format (e.g 5000000 for 5 Mb genome) | This setting is used in conjunction with `flye_asm_coverage` to subsample the reads used in the initial disjointig step only; all reads are used in subsequent steps. The values in the two parameters are used to calculate a target yield used to subsample the longest reads in the dataset, see [Flye docs](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) for more information. *Note* For runs with mixed genome sizes, preference the larger genome size. |  |
| flye_asm_coverage | integer | Target coverage to use for subsampling in de novo assembly | This setting is used in conjunction with `flye_genome_size` to subsample the reads used in the initial disjointig step only; all reads are used in subsequent steps. The values in the two parameters are used to calculate a target yield used to subsample the longest reads in the dataset, see [Flye docs](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) for more information. |  |
| flye_opts | string | Command-line arguments for flye | [Command line arguments](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage) which can be used to alter the de novo assembly process. Enter the command as quoted string (e.g '--meta --iterations 2'). Flye's `--genome-size` and `--asm-coverage` parameters can be set directly in the workflow with `--flye_genome_size` and `--flye_asm_coverage`, respectively. |  |
| min_read_length | integer | Reads below this threshold will be removed from analysis. |  | 1000 |
| max_coverage_plots | integer | Maximum number of largest contigs for which depth coverage plots will be displayed in the summary report. | When specified, only the N largest contigs (by total length) will be shown in the genome coverage depth plots. Set to 'null' to display all contigs. | 20 |
| run_plasmid_id | boolean | Enable plasmid identification and characterization using MOB-suite. | MOB-suite will identify plasmid sequences in the assembled genomes, reconstruct plasmids, and predict plasmid mobility based on the presence of mobilization and replication genes. Results include plasmid sequences, mobility predictions, and incompatibility group typing. | True |
| plasmid_id_opts | string | Command-line arguments for mob_recon | Additional command-line arguments that can be passed to mob_recon. For example, `--filter_db` to specify custom databases. See [MOB-suite](https://github.com/phac-nml/mob-suite) documentation for available options. |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Number of CPU threads. | Provided to alignment, Flye assembly, Dnaapler, aligment, Bakta, and MOB-suite workflow steps to improve performance. | 3 |


