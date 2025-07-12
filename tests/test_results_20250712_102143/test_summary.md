# PandaDock Command Testing Summary

## Test Results
- Total tests: 14
- Passed tests: 14
- Failed tests: 0

## Issues Identified
- IC50 calculation issues: 2 tests affected
- --save-complex functionality issues: 0 tests affected

## Recommendations
1. If IC50 values are identical across poses, check the ML rescoring implementation
2. If --save-complex is not working, check the complex file saving logic
3. Review the scoring function implementations for proper affinity prediction
4. Ensure all output formats (HTML, JSON, CSV) are generating correctly

## Test Directory Structure
- balanced_ml_rescoring
- basic_docking
- combined_options
- csv_output
- fast_mode
- flexible_residues
- high_exhaustiveness
- json_output
- multiple_poses
- precise_mode
- save_complex
- scoring_pandacore
- scoring_pandaml
- scoring_pandaphysics
