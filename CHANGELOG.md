# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## unreleased
### Fixes
- Prokka only runs in denovo assembly mode

## [v0.2.6]
### Changes
- Added quast for assembly stats
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

Initial release

### Added
- Basic running of medaka variant calling and report.
