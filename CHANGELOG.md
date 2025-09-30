# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added
- CHANGELOG.md file for tracking version history

### Changed
- Preparing for PyPI packaging improvements

## [0.1.0] - 2023-02-03

### Added
- Initial release of MS1Connect
- Mass spectrometry run similarity scoring using maximum bipartite matching
- MS1 feature detection using pyOpenMS
- Edge generation between MS1 features across runs
- Matroid-based optimization for feature matching
- Sparse edge similarity matrix computation
- Integration with Singularity container for coopraiz optimization
- Command-line interface for processing mzML files
- Support for customizable hyperparameters (lambda1-4, alpha, beta, gamma)
- Visualization and plotting capabilities
- Comprehensive documentation and citation information

### Dependencies
- pyOpenMS for mass spectrometry data processing
- pandas for data manipulation
- numba for numerical computations
- scipy for scientific computing
- scikit-learn for machine learning utilities
- matplotlib and seaborn for visualization
- C++ compiler (gcc) for native extensions
- Singularity for containerized optimization

### System Requirements
- Python >= 3.11
- C++ compiler
- Singularity runtime

## Publication

This software was published in:
> Lin A, Deatherage Kaiser BL, Hutchison JR, Bilmes JA, Noble WS. MS1Connect: a mass spectrometry run similarity measure. Bioinformatics. 2023 Feb 3;39(2)

[Unreleased]: https://github.com/bmx8177/MS1Connect/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/bmx8177/MS1Connect/releases/tag/v0.1.0
