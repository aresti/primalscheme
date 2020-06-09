# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - Nighthawk - 2020-06-09
### Added
- Use all reference sequences for primer generation
- 3' distance weighted mismatch scoring
- Avoid stable heteodimers in same pool

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

### Removed
- primer3 FindPrimers function
- Redundant classes

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









