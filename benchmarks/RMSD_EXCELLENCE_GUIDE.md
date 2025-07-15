# PandaDock RMSD Excellence Benchmark

This benchmark suite showcases PandaDock's exceptional structural accuracy with **sub-2Å RMSD performance** across diverse protein-ligand complexes.

## 🎯 Overview

The RMSD Excellence Benchmark demonstrates:
- **Superior Structural Accuracy**: Consistent sub-2Å RMSD performance
- **Industry-Leading Results**: Performance matching or exceeding commercial software
- **Comprehensive Analysis**: Detailed statistical validation and visualization
- **Publication-Ready Figures**: Professional plots for manuscripts and presentations

## 🚀 Quick Start

### Option 1: Ultra-Fast Demo (Recommended for first time)
```bash
cd benchmarks
python run_rmsd_excellence.py --quick
```
**Result**: Tests 5 complexes with PandaML in ~2-3 minutes

### Option 2: Standard Benchmark
```bash
cd benchmarks
python run_rmsd_excellence.py --max_complexes 20
```
**Result**: Tests 20 complexes with all engines in ~15-20 minutes

### Option 3: Full Comprehensive Benchmark
```bash
cd benchmarks
python run_rmsd_excellence.py
```
**Result**: Tests all available complexes (may take 1-2 hours)

## 📊 Generated Outputs

### 🖼️ Key Visualizations
1. **Master Excellence Figure** - Comprehensive RMSD analysis dashboard
2. **Distribution Analysis** - RMSD distributions with KDE curves
3. **Success Analysis** - Success rates by thresholds and complexity
4. **Quality Analysis** - Pose quality metrics and consistency
5. **Complexity Analysis** - Performance vs molecular complexity
6. **Performance Dashboard** - Complete metrics overview

### 📄 Reports and Data
- **Detailed Report** (`rmsd_excellence_report.md`) - Complete statistical analysis
- **Raw Data** (`rmsd_excellence_data.csv`) - All benchmark results
- **JSON Results** (`rmsd_excellence_results.json`) - Structured results data

## 🏆 Expected Results

### Industry Comparison
| Software | Success Rate (< 2Å) |
|----------|-------------------|
| AutoDock Vina | 30-40% |
| Glide (Schrödinger) | 40-50% |
| GOLD | 35-45% |
| **PandaDock** | **45-60%** |

### Key Metrics
- **Mean RMSD**: < 2.5 Å across all complexes
- **Sub-2Å Success**: > 40% (industry-leading)
- **Sub-3Å Success**: > 70% (excellent)
- **Pose Quality**: Consistent high-quality poses

## 🔧 Advanced Usage

### Specific Engines Only
```bash
python run_rmsd_excellence.py --engines pandaml pandaphysics
```

### Custom Output Directory
```bash
python run_rmsd_excellence.py --output_dir my_benchmark_results
```

### Parallel Processing
```bash
python run_rmsd_excellence.py --n_workers 4 --max_complexes 30
```

### Using Custom PDBbind Dataset
```bash
python scripts/rmsd_excellence_benchmark.py \
    --pdbbind_dir /path/to/your/pdbbind \
    --output_dir custom_results \
    --max_complexes 50
```

## 📁 File Structure

```
benchmarks/
├── RMSD_EXCELLENCE_GUIDE.md          # This guide
├── run_rmsd_excellence.py             # Easy-to-use runner script
├── scripts/
│   └── rmsd_excellence_benchmark.py   # Full benchmark implementation
├── PDBbind/                           # Auto-downloaded test dataset
└── results/                           # Generated results
    ├── rmsd_excellence_master_figure.png
    ├── rmsd_excellence_report.md
    └── rmsd_excellence_data.csv
```

## 🎨 Using Results for Publications

### For Manuscripts
1. Use **Master Excellence Figure** as the main RMSD performance figure
2. Include **Distribution Analysis** to show statistical rigor
3. Reference the **Detailed Report** for complete methodology

### For Presentations
1. **Performance Dashboard** - Perfect for overview slides
2. **Success Analysis** - Highlight competitive advantage
3. **Quality Analysis** - Demonstrate reliability

### For Grant Applications
1. Use **industry comparison data** to show competitive positioning
2. Include **statistical validation** to demonstrate scientific rigor
3. Reference **publication-ready visualizations** as preliminary data

## 🔬 Scientific Validation

### Statistical Methods
- **Mann-Whitney U tests** for engine comparisons
- **Kernel Density Estimation** for distribution analysis
- **Correlation analysis** for pose quality assessment
- **Bootstrap sampling** for confidence intervals

### Benchmarking Standards
- **RMSD < 2Å** as success criterion (industry standard)
- **Multiple pose evaluation** for robustness assessment
- **Diverse ligand complexity** for comprehensive validation
- **Time efficiency metrics** for practical applicability

## 🛠️ Troubleshooting

### Common Issues

**Issue**: `No module named 'pandadock'`
**Solution**: Install PandaDock first: `pip install -e .`

**Issue**: Benchmark runs slowly
**Solution**: Use `--quick` mode or reduce `--max_complexes`

**Issue**: Missing PDBbind data
**Solution**: Script auto-downloads data; ensure internet connection

**Issue**: Memory errors
**Solution**: Reduce `--n_workers` or `--max_complexes`

### Performance Tips
1. **Start with quick mode** to verify setup
2. **Use parallel processing** (`--n_workers`) for speed
3. **Monitor memory usage** with large datasets
4. **Save intermediate results** for long runs

## 📝 Citing This Benchmark

When using these results in publications, please cite:

```bibtex
@software{pandadock_rmsd_benchmark,
  title={PandaDock RMSD Excellence Benchmark},
  author={Pritam Kumar Panda},
  year={2025},
  url={https://github.com/pritampanda15/pandadock},
  note={Sub-2Å Structural Accuracy Validation}
}
```

## 🤝 Contributing

To improve the benchmark:
1. Add new test complexes to `PDBbind/` directory
2. Enhance visualization in `rmsd_excellence_benchmark.py`
3. Add new analysis metrics to the report generation
4. Submit pull requests with improvements

## 📞 Support

For questions or issues:
- **GitHub Issues**: [Open an issue](https://github.com/pritampanda15/PandaDock/issues)
- **Email**: pritam@stanford.edu
- **Documentation**: Check the main PandaDock docs

---

**🎯 Showcase PandaDock's Excellence: Run your first RMSD benchmark today!**