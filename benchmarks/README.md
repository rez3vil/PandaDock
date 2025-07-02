# PandaDock PDBbind Benchmarking System

This directory contains a comprehensive benchmarking system for PandaDock using the PDBbind dataset. The system compares CPU vs GPU performance and validates binding affinity predictions against experimental values.

## Features

- **Performance Benchmarking**: Compare CPU vs GPU processing times
- **Algorithm Comparison**: Test different docking algorithms (genetic, random, pandadock)
- **Scoring Function Evaluation**: Compare standard, enhanced, and physics-based scoring
- **Binding Affinity Validation**: Correlate predicted values with experimental IC50/Kd/Ki data
- **Statistical Analysis**: Pearson correlation, R², RMSE, and significance testing
- **Comprehensive Visualization**: Automated plot generation for all metrics
- **Parallel Processing**: Support for multi-core execution to speed up benchmarks

## Prerequisites

### Required Software
- PandaDock (properly installed and in PATH)
- Python 3.7+
- Required Python packages:
  ```bash
  pip install pandas numpy scipy sklearn matplotlib seaborn
  ```

### Required Data
- **PDBbind Dataset**: Download the PDBbind refined set from http://www.pdbbind.org.cn/
- Extract to a directory with structure:
  ```
  pdbbind/
  ├── INDEX_refined_set.2020
  ├── 1abc/
  │   ├── 1abc_protein.pdb
  │   ├── 1abc_ligand.sdf
  │   └── 1abc_pocket.pdb
  ├── 1def/
  │   ├── 1def_protein.pdb
  │   ├── 1def_ligand.sdf
  │   └── 1def_pocket.pdb
  └── ...
  ```

### GPU Support (Optional)
- NVIDIA GPU with CUDA support
- PandaDock compiled with GPU support

## Quick Start

### 1. Setup
```bash
# Navigate to benchmarks directory
cd benchmarks

# Make scripts executable
chmod +x run_benchmark_examples.sh

# Edit the script to set your PDBbind directory path
nano run_benchmark_examples.sh
# Set: PDBBIND_DIR="/path/to/your/pdbbind"
```

### 2. Run Quick Test
```bash
# Test with 10 entries (CPU only)
python pdbbind_benchmark.py \
    --pdbbind_dir /path/to/pdbbind \
    --output test_results \
    --max_entries 10 \
    --algorithms genetic \
    --devices CPU \
    --scoring standard \
    --iterations 25
```

### 3. Analyze Results
```bash
# Print summary
python analyze_benchmark_results.py test_results

# Generate detailed report
python analyze_benchmark_results.py test_results --report detailed_report.html

# Export to CSV for further analysis
python analyze_benchmark_results.py test_results --export_csv results.csv
```

## Benchmark Examples

### Example 1: CPU vs GPU Performance
```bash
python pdbbind_benchmark.py \
    --pdbbind_dir /path/to/pdbbind \
    --output cpu_vs_gpu_benchmark \
    --max_entries 100 \
    --algorithms genetic \
    --devices CPU GPU \
    --scoring standard \
    --iterations 50 \
    --parallel_jobs 4
```

### Example 2: Algorithm Comparison
```bash
python pdbbind_benchmark.py \
    --pdbbind_dir /path/to/pdbbind \
    --output algorithm_comparison \
    --max_entries 150 \
    --algorithms genetic random pandadock \
    --devices CPU \
    --scoring standard \
    --iterations 50 \
    --parallel_jobs 6
```

### Example 3: Scoring Function Evaluation
```bash
python pdbbind_benchmark.py \
    --pdbbind_dir /path/to/pdbbind \
    --output scoring_evaluation \
    --max_entries 100 \
    --algorithms genetic \
    --devices CPU \
    --scoring standard enhanced physics-based \
    --iterations 50 \
    --parallel_jobs 3
```

### Example 4: Comprehensive Benchmark
```bash
python pdbbind_benchmark.py \
    --pdbbind_dir /path/to/pdbbind \
    --output comprehensive_benchmark \
    --max_entries 500 \
    --algorithms genetic random \
    --devices CPU GPU \
    --scoring standard enhanced \
    --iterations 100 \
    --parallel_jobs 8
```

## Output Structure

Each benchmark run creates an organized output directory:

```
benchmark_results/
├── benchmark_results.json          # Raw results data
├── benchmark_results.csv           # Results in CSV format
├── benchmark_analysis.json         # Statistical analysis
├── benchmark_TIMESTAMP.log         # Detailed execution log
├── plots/                          # Generated visualizations
│   ├── performance_comparison.png
│   ├── algorithm_comparison.png
│   ├── binding_affinity_correlation.png
│   ├── success_rates.png
│   └── runtime_distribution.png
└── [pdb_code]_[algorithm]_[device]_[scoring]/  # Individual run outputs
    ├── binding_affinity_report.txt
    ├── binding_affinity_report.csv
    ├── detailed_docking_report.txt
    ├── poses/
    └── plots/
```

## Key Metrics Analyzed

### Performance Metrics
- **Runtime**: Total execution time per docking run
- **Success Rate**: Percentage of successful docking runs
- **Speedup Factor**: GPU performance improvement over CPU
- **Scalability**: Performance with different numbers of parallel jobs

### Docking Quality Metrics
- **Best Docking Score**: Lowest energy score achieved
- **Poses Generated**: Number of valid poses per run
- **Convergence Rate**: How quickly algorithms find good solutions

### Binding Affinity Validation
- **Pearson Correlation (R)**: Linear correlation between predicted and experimental values
- **R-squared (R²)**: Coefficient of determination
- **RMSE**: Root mean square error in log space
- **P-value**: Statistical significance of correlation

### Experimental Data Support
- **Kd Values**: Dissociation constants
- **Ki Values**: Inhibition constants  
- **IC50 Values**: Half-maximal inhibitory concentrations
- **Mixed Dataset**: Automatic handling of different experimental measurement types

## Interpretation Guide

### Performance Results
- **GPU Speedup**: Typical speedups of 2-10x are expected for genetic algorithms
- **Algorithm Speed**: Random search is typically fastest, genetic algorithms most accurate
- **Scoring Function Impact**: Physics-based scoring is slower but more accurate

### Binding Affinity Correlation
- **R > 0.7**: Strong correlation, excellent predictive performance
- **R = 0.5-0.7**: Moderate correlation, good for ranking compounds
- **R = 0.3-0.5**: Weak correlation, limited predictive value
- **R < 0.3**: Very weak correlation, poor predictive performance

### Success Rate Guidelines
- **>95%**: Excellent, robust implementation
- **90-95%**: Good, some challenging cases
- **80-90%**: Acceptable, may need parameter tuning
- **<80%**: Poor, requires investigation

## Customization

### Adding Custom Algorithms
1. Implement algorithm in PandaDock
2. Add algorithm name to `choices` in argument parser
3. Update benchmark script to handle new algorithm

### Custom Scoring Functions
1. Add scoring function to PandaDock
2. Update `--scoring` choices in benchmark script
3. Modify result parsing to handle new energy components

### Additional Metrics
1. Modify `BenchmarkResult` dataclass to include new metrics
2. Update `_parse_docking_results` method
3. Add analysis methods in `BenchmarkAnalyzer`

## Troubleshooting

### Common Issues

**"PDBbind directory not found"**
- Verify the path to your PDBbind dataset
- Ensure the INDEX file exists in the root directory

**"No valid entries found"**
- Check that PDB subdirectories contain required files
- Verify file naming convention matches expected pattern

**GPU errors**
- Ensure CUDA is properly installed
- Check that PandaDock has GPU support compiled
- Try reducing batch size or parallel jobs

**Memory errors**
- Reduce `--parallel_jobs` parameter
- Limit `--max_entries` for testing
- Close other applications to free memory

**Timeout errors**
- Increase timeout in `run_single_docking` method
- Reduce `--iterations` parameter
- Use `--fast-mode` option

### Performance Optimization

**For Speed**:
```bash
--iterations 25 --fast-mode --parallel_jobs 8
```

**For Accuracy**:
```bash
--iterations 100 --scoring physics-based --local-opt
```

**For Large Datasets**:
```bash
--parallel_jobs 16 --max_entries 1000
```

## Expected Runtime

Typical runtimes on modern hardware:

| Configuration | Entries | Runtime | 
|---------------|---------|---------|
| Quick test (CPU, 10 entries) | 10 | 5-15 minutes |
| Small benchmark (CPU, 100 entries) | 100 | 1-3 hours |
| Medium benchmark (CPU+GPU, 500 entries) | 500 | 4-8 hours |
| Large benchmark (multi-config, 1000 entries) | 1000 | 8-16 hours |
| Full dataset (5000+ entries) | 5000+ | 24-48 hours |

*Note: Actual runtimes depend on hardware, algorithms, and scoring functions used.*

## Citation

If you use this benchmarking system in your research, please cite:

```
PandaDock Benchmarking System
GitHub: https://github.com/pritampanda15/PandaDock
```

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review log files in the output directory
3. Open an issue on the PandaDock GitHub repository
4. Contact: pritam@stanford.edu

## License

This benchmarking system is part of PandaDock and follows the same MIT license.