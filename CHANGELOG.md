<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.1] - 2020-06-27

### Added

- Add MultiplexScheme parameters to debug log output
- Add version number to debug log and run report

### Changed

- Improve --amplicon-size help text
- Move reference info logging into cli

### Fixed

- Fix outpath option in README
- Fix --amplicon-size % variation in help text

### Removed

- Remove logger name from log file output

## [1.2.0] - 2020-06-26

### Added

- Add `--pinned, -p` option: reverts to designing primers exclusively from first reference
- Add short option names
- Add option to pass single target amplicon size and auto-determine a suitable range
- Add web interface link to README

### Changed

- Switch from Argparse to Click
- Rename options
- Allow `--amplicon-size, -a` to be passed once to set a target size, or twice to set an exact size range
- Refactor test_cli.py to use CliRunner throughout

### Removed

- Remove argparse
- Remove non-functional --step-distance option
- Remove automatic reference sort on length (no longer serves a purpose)

### Fixed

- Fix edge case where last-but-one region steps right, putting last region out of bounds
- Fix (usually) inconsequential out-by-one error to make Window.slice_end coord inclusive

## [1.1.2] - 2020-06-22

### Fixed

- Fix number of reported candidates considered
- Fix missing exception message where no suitable pair found
- Fix incorrect logging of --no-sort

## [1.1.1] - 2020-06-15

### Added

- Define ProgressTracker interface as abstract base class (to facilitate web interface progress)

### Changed

- Reorganise Region and Window into region.py

## [1.1.0] - 2020-06-14

### Added

- Drop secondary references on repeated flank alignment failure
- Add primary_ref, secondary_refs and excluded_refs to json run report

### Changed

- Reduce complexity of Region find_primers() method
- Use progress library
- Add region and considered primers to progress message

### Removed

- DIY progress bar

## [1.0.0] - Nighthawk - 2020-06-12

### Added

- Use all reference sequences for primer generation
- 3' distance-weighted mismatch scoring
- Avoid stable heterodimers in same pool
- Sort inputs by length, to make longest the primary by default (--no-sort reverts)
- Inserts bed file output
- Run report json output

### Changed

- Infer mismatches from aligned references and primer coordinates
- Merge mismatch and base primer scoring
- Generate primers by digesting all references into K-mers at region flanks
- Use own implementation of Primer3 base penalty algorithm
- Update parameters for thermodyamic calculations
- Rationalise the config file
- Remove any remaining constants from the code base
- Reduce code complexity
- Consolidate exception handling
- Improve candidate sorting to ensure deterministic output.

### Removed

- primer3 FindPrimers function
- Redundant classes
- Redundant exceptions

## [0.3.0] - Blackbird - 2020-05-21

### Added

- Windows compatibility
- Tests
- Valid BED file strand column

### Changed

- Use parasail for alignment
- Improve stdout and logging
- Improve code readability
- Improved CLI
- Candidate ranking improvements
- Various bug fixes
- Speed improvements
