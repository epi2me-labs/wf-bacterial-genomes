# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

