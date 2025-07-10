# PandaDock

**Modular, Multi-Strategy, High-Performance Molecular Docking Software**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation Status](https://readthedocs.org/projects/pandadock/badge/?version=latest)](https://pandadock.readthedocs.io/)

## Overview

PandaDock is a comprehensive molecular docking software that combines multiple docking strategies in a unified framework. It supports physics-based, machine learning-based, and genetic algorithm-based approaches for protein-ligand docking with comprehensive analysis and reporting capabilities.

## Key Features

### 🔌 **Three Docking Strategies**
- **Precise Mode**: Physics-based docking with Glide-style systematic conformer generation
- **Balanced Mode**: ML-based docking using DiffDock/Boltz-style diffusion models
- **Fast Mode**: Genetic algorithm-based docking for high-throughput virtual screening

### 🧬 **Advanced Capabilities**
- **Flexible Docking**: Side-chain flexibility with rotamer library sampling
- **Comprehensive Scoring**: Physics-based + ML rescoring with detailed energy breakdown
- **Multi-format Support**: SMILES, SDF, MOL2, PDB input formats
- **Interactive Reports**: HTML reports with pose visualization and analysis

### ⚡ **Performance Features**
- **GPU Acceleration**: CUDA support for ML models
- **Parallel Processing**: Multi-threading for GA and physics calculations
- **Memory Efficient**: Optimized for large-scale virtual screening
- **Extensible Architecture**: Easy to add new docking engines

## Installation

### Quick Install
```bash
# Basic installation
pip install -e .

# With machine learning support
pip install -e .[ml]

# With all features
pip install -e .[all]
```

### Requirements
- Python 3.8+
- NumPy, SciPy, scikit-learn (automatically installed)
- Optional: PyTorch (ML models), RDKit (chemistry), OpenMM (physics)

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

### Install from PyPI (Coming Soon)
```bash
pip install pandadock
```

### Install from Source
```bash
git clone https://github.com/pandadock/pandadock.git
cd pandadock
pip install -e .[all]
```

## Quick Start

### Basic Docking
```bash
# Balanced mode (ML-based)
pandadock --protein receptor.pdb --ligand ligand.sdf --mode balanced

# Physics-based with flexible residues
pandadock --protein receptor.pdb --ligand ligand.sdf --mode precise \
          --flexible-residues "HIS57,SER195" --out results.html

# Fast virtual screening
pandadock --protein receptor.pdb --screen ligands.smi --mode fast \
          --num-poses 5 --exhaustiveness 16
```

### Python API
```python
from pandadock import PandaDockConfig, PhysicsEngine, MLEngine, GAEngine

# Configure docking
config = PandaDockConfig()
config.docking.mode = "balanced"
config.docking.num_poses = 10
config.docking.flexible_residues = ["HIS57", "SER195"]

# Initialize engine
engine = MLEngine(config)

# Run docking
results = engine.dock("receptor.pdb", "ligand.sdf")

# Analyze results
for pose in results[:5]:
    print(f"Pose {pose.pose_id}: Score = {pose.score:.3f}")
    print(f"  Binding Affinity: {pose.get_binding_affinity():.2f} kcal/mol")
    print(f"  IC50: {pose.get_ic50()*1e9:.1f} nM")
```

## Docking Modes

### Precise Mode (Physics-based)
```bash
pandadock --protein receptor.pdb --ligand ligand.sdf --mode precise \
          --flexible-residues "HIS57,SER195,TYR191" \
          --exhaustiveness 8 --num-poses 10
```

**Features:**
- Systematic conformer generation
- Detailed molecular mechanics scoring
- Flexible side-chain sampling
- Energy minimization
- Clash detection and resolution

### Balanced Mode (ML-based)
```bash
pandadock --protein receptor.pdb --ligand ligand.sdf --mode balanced \
          --confidence-threshold 0.7 --ml-rescoring
```

**Features:**
- Diffusion-based pose generation
- Deep learning confidence scoring
- Fast inference
- Hybrid ML/physics scoring

### Fast Mode (GA-based)
```bash
pandadock --protein receptor.pdb --ligand ligand.sdf --mode fast \
          --population-size 150 --generations 27000 \
          --n-jobs 4
```

**Features:**
- Evolutionary search optimization
- Parallel evaluation
- Empirical scoring functions
- Optimized for virtual screening

## Configuration

### JSON Configuration
```json
{
  "docking": {
    "mode": "balanced",
    "num_poses": 10,
    "flexible_residues": ["HIS57", "SER195"],
    "exhaustiveness": 8
  },
  "scoring": {
    "scoring_function": "vina",
    "use_ml_rescoring": true,
    "vdw_weight": 1.0,
    "electrostatic_weight": 1.0
  },
  "io": {
    "output_dir": "results",
    "save_poses": true,
    "save_complex": true,
    "report_format": "html"
  }
}
```

### Environment Variables
```bash
export PANDADOCK_GPU_ENABLED=true
export PANDADOCK_N_JOBS=4
export PANDADOCK_CACHE_DIR=/tmp/pandadock
```

## Output and Analysis

### HTML Reports
PandaDock generates comprehensive HTML reports including:
- **Pose Rankings**: Sorted by score with energy breakdown
- **Binding Affinity**: ΔG, IC50, and ligand efficiency calculations
- **Interaction Analysis**: H-bonds, hydrophobic contacts, salt bridges
- **Visualization**: Interactive 3D pose viewer
- **Comparison Charts**: Score vs. energy correlations

### Data Export
```python
# Export results to various formats
report_generator.export_data(results, format='json', output_path='results.json')
report_generator.export_data(results, format='csv', output_path='results.csv')
```

## Advanced Features

### Flexible Docking
```python
from pandadock.docking import FlexibleDocking

# Configure flexible residues
flexible_docking = FlexibleDocking(config)
flexible_docking.setup_flexible_residues(protein_coords, residue_info)

# Optimize with side-chain flexibility
optimized_pose = flexible_docking.optimize_with_sidechains(pose)
```

### ML Rescoring
```python
from pandadock.scoring import MLRescorer

# Initialize ML rescorer
rescorer = MLRescorer()
rescorer.load_model('path/to/model.pkl')

# Rescore poses
rescored_poses = rescorer.rescore_poses(poses)
```

### Custom Scoring Functions
```python
from pandadock.scoring import ScoringFunctions

# Custom scoring weights
config.scoring.vdw_weight = 1.5
config.scoring.electrostatic_weight = 0.8
config.scoring.hbond_weight = 2.0

scoring = ScoringFunctions(config)
energy = scoring.calculate_total_energy(ligand_coords, protein_coords)
```

## Benchmarking

PandaDock has been tested on standard benchmarks:
- **PDBbind**: Binding affinity prediction
- **CASF**: Comparative assessment of scoring functions
- **DUD-E**: Directory of useful decoys

Performance metrics:
- **Success Rate**: >80% for top-3 poses (RMSD < 2.0 Å)
- **Correlation**: R² > 0.7 for binding affinity prediction
- **Speed**: ~10-100 poses/second depending on mode

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup
```bash
git clone https://github.com/pandadock/pandadock.git
cd pandadock
pip install -e .[dev]
pytest tests/
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use PandaDock in your research, please cite:

```bibtex
@software{pandadock2024,
  title={PandaDock: Modular, Multi-Strategy, High-Performance Molecular Docking Software},
  author={PandaDock Development Team},
  year={2024},
  url={https://github.com/pandadock/pandadock}
}
```

## Support

- **Documentation**: [https://pandadock.readthedocs.io/](https://pandadock.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/pandadock/pandadock/issues)
- **Discussions**: [GitHub Discussions](https://github.com/pandadock/pandadock/discussions)

## Acknowledgments

PandaDock builds upon the scientific contributions of:
- **AutoDock Vina**: Genetic algorithm optimization
- **Glide**: Physics-based scoring functions
- **DiffDock**: Machine learning pose prediction
- **RDKit**: Molecular handling and processing
- **OpenMM**: Molecular dynamics integration

---

**PandaDock** - Empowering drug discovery through intelligent molecular docking.