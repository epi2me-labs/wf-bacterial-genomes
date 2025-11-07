# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Fixed
- Changed Nextflow data flow logic to avoid mismatch between sample labels and read statistics in the report.
- Fixed a defect that lead to medakaInference process failing upon retry.

## [v1.4.4]
This patch release of wf-bacterial-genomes correctly labels the workflow as unsupported on ARM architectures to prevent users installing the workflow on unsupported devices via EPI2ME Desktop.

## [v1.4.3]
This patch release of wf-bacterial-genomes updates the workflow title for display in EPI2ME Desktop. This patch does not affect any workflow outputs. Users do not need to adopt this release.

### Changed
- Updated to wf-template v5.6.2 to maintain compliance with our latest wf-template standard: this does not impact the workflow.
    - Reduce verbosity of debug logging from fastcat which can occasionally occlude errors found in FASTQ files during ingress.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.
    - pre-commit configuration to resolve an internal dependency problem with flake8.

### Fixed
- Updated to wf-template v5.6.2 to maintain compliance with our latest wf-template standard:
    - dacite.exceptions.WrongTypeError during report generation when barcode is null.
    - Sequence summary read length N50 incorrectly displayed minimum read length, it now correctly shows the N50.
    - Sequence summary component alignment and coverage plots failed to plot under some conditions.

## [v1.4.2]
### Changed
- Reconciled workflow with wf-template v5.5.0.
- Increasing memory retries for flye deNovo process
### Fixed
- `makeReport` process used to fail with unaligned (u)BAMs as input because unaligned bams were excluded by default by `xam_ingress`. We fixed this by setting `xam_ingress` to keep the unaligned files as well.
- `bam` and `bai` files from mapping the reads to the reference or draft assembly, depending on mode used, were previously missing from the output directory. They are now included.
- Single sample reports incorrectly displaying 'de novo assembly failed' on successful samples.

## [v1.4.1]
### Added
- `Flye` stats output and report
### Changed
- Reconciled workflow with wf-template v5.3.0
- Filtering samples with no data on meta.n_seqs instead of presence of fastq file post ingress

## [v1.4.0]
### Changed
- Reconciled workflow with wf-template v5.2.5.
- Report wording.
- Updated Medaka to v2.0.0.

## [v1.3.0]
### Added
- Workflow now accepts BAM as well as FASTQ files as input (using the `--bam` or `--fastq` parameters, respectively)
- `--override_basecaller_cfg` parameter for cases where automatic basecall model detection fails or users wish to override the automatic choice.
### Removed
- The `--basecaller_cfg`, `--medaka_consensus_model`, and `--medaka_variant_model` parameters as the appropriate Medaka model is now automatically determined from the input data.
### Changed
- New docker images for [resfinder](https://hub.docker.com/r/ontresearch/resfinder) and [mlst](https://hub.docker.com/r/ontresearch/mlst).

## [v1.2.0]
### Added
- `client_fields` parameter to allow input of a JSON file of key value pairs to display on output reports.
- `min_read_length` parameter to remove reads below specified length (default 1000bp) from downstream analysis, to improve de novo assembly process
- Salmonella serotyping with `SeqSero2`

### Fixed
- Duplicate entries in `Pointfinder` processing

## [v1.1.1]
### Fixed
- Report generation when `Resfinder` fails

## [v1.1.0]
### Added
- Sample results aggregated into `results.json`
- `flye_genome_size` and `flye_asm_coverage` parameters for controlling the initial downsampling step before the de novo assembly

### Changed
- De novo assembly mode uses `--nano-hq` rather than `--nano-raw`
- Some formatting in github issue template.
- Retry and memory bump if de novo assembly fails first time

### Fixed
- Workflow now runs to completion when a sample fails the de novo assembly
- The workflow requesting too little memory for some processes.

## [v1.0.2]
### Added
- FASTA now includes basecaller model in headers

### Changed
- Updated to most recent version of Medaka container.

## [v1.0.1]
### Changed
- Minimum compute requirements

## [v1.0.0]
### Added
- Cloud support for the workflow within the EPI2ME Application.

### Changed
- Documentation

## [v0.4.0]
### Added
- MacOS ARM64 support
- New parameter `--flye_opts` for passing additional arguments to `flye`.

### Changed
- Clarify docker is default in README

### Fixed
- De novo assembly failing due to low coverage in some situations.

## [v0.3.3]
### Fixed
- Overwrites in Nextflow config implemented incorrectly

### Changed
- Updated Medaka to 1.9.1.

## [v0.3.2]
### Fixed
- Edge case where medaka variant output is unsorted and causes medaka annotate to exit

### Changed
- Bumped minimum required Nextflow version to 23.04.2.
- Now uses Medaka v1.8.2 with updated models.
- Options for the `--basecaller_cfg` parameter were updated. The default is now `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`.

## [v0.3.1]
### Changed
- GitHub issue templates
- Output GFF and GBK files from Prokka
- Updated resfinder version to 4.3.2
- Removed mutation of unknown effect in SNP-mediated AMR genes output

## [v0.3.0]	
### Added
- Isolate single sample reports
- Include disinfectant resistance results in the report.
- MLST core gene analysis added to `--isolates` parameter.

### Changed
- `species` parameter is removed, valid pointfinder species will be inferred from MSLT results.
- In case `flye` fails due to low coverage, the workflow will continue and this will be indicated in the report.
- Bumped minimum required Nextflow version to 22.10.8
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`

## [v0.2.14]
### Added
- `--isolates` parameter that will run the ResFinder tool on the final assembly to output antimicrobial resistance genes.
- Configuration for running demo data in AWS

### Changed
- Report is now created with `ezcharts`.

## [v0.2.13]
### Fixed
- Rows with too few / too many columns in `medaka_models.tsv`.
- Check sample sheet script.

### Changed
- Now uses new `fastq_ingress` implementation.

## [v0.2.12]
### Fixed
- Medaka models added to container

### Removed
- QUAST

## [v0.2.11]
### Changed
- `--basecall_cfg` is now used to determine suitable Medaka model, alternatively provide the name of a model with `--medaka_consensus_model` and `--medaka_variant_model` to override automatic selection.

## [v0.2.10]
### Fixed
- sample_sheet format in schema to expect a file

## [v0.2.9]
### Changed
- Updated description in manifest

## [v0.2.8]
### Changed
- Output QUAST stats for reference and denovo based assembly

## [v0.2.7]
### Changes
- Replace QUAST with MetaQUAST
- Add species ID to run summary table
- For reference based assembly `--reference_based_assembly` parameter should now be provided with a `--reference`. The default is to use denovo assembly.
- Tidy up presentation in report
- `-profile conda` is no longer supported, users should use `-profile standard` (Docker) or `-profile singularity` instead
- Docs update

### Added
- `nextflow run epi2me-labs/wf-bacterial-genomes --version` will now print the workflow version number and exit

### Fixes
- Prokka only runs in denovo assembly mode
- Tidy up report code

## [v0.2.6]
### Changes
- Added QUAST for assembly stats
- Remove sanitize option

### Fixes
- Update syntax to fix reference error

## [v0.2.5]
### Changes
- Better help text on cli
- Fastqingress metadata map
- Use groovy script to ping after workflow has run

### Fixes
- Output medaka vcf
- Remove reliance on simpleName

## [v0.2.4]
### Fixed
- Amend report name.

## [v0.2.3]
### Added
- Add read me docs.

## [v0.2.2]
### Added
- Option to add suffix to HTML report name.
- Visualisation of prokka output.
- Choice of de novo assembly or alignment.
- Supports multibarcodes

### Updated
- Medaka version.
- Depth coverage graphs.
- Use mosdepth and fastcat.

## [v0.2.1]
### Changed
- Update project with latest practices from wf-template.
- Use `mamba` by default when using conda profile.

### Fixed
- Incorrect specification of conda environment file location.

## [v0.2.0]
### Changed
- Rework workflow to use new medaka methodology for improved robustness of results.

## [v0.1.3]
### Changed
- Standardised report name.

## [v0.1.2]
### Fixed
- Resolved confusion between documentation and workflow: workflow now requires a directory as input.

## [v0.1.1]
### Fixed
- aplanat import error

## [v0.1.0]
### Added
- Prokka can be optionally run to annotate consensus sequence.

### Changed
- Variant call summary produced using aplanat report component.

## [v0.0.1]
* Initial release

### Added
- Basic running of medaka variant calling and report.

