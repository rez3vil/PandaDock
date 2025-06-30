# 🐼 PandaDock

**A Python-based GPU/CPU-accelerated molecular docking platform for computational drug discovery and bioinformatics.**

<p align="center">
  <a href="https://github.com/pritampanda15/PandaDock">
    <img src="https://github.com/pritampanda15/PandaDock/blob/main/logo/Pandadock_logo.png" width="600" alt="PandaDock Logo"/>
  </a>
</p>

<p align="center">
  <a href="https://pypi.org/project/pandadock/">
    <img src="https://img.shields.io/pypi/v/pandadock.svg" alt="PyPI Version">
  </a>
  <a href="https://github.com/pritampanda15/PandaDock/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/pritampanda15/PandaDock" alt="License">
  </a>
  <a href="https://github.com/pritampanda15/PandaDock/stargazers">
    <img src="https://img.shields.io/github/stars/pritampanda15/PandaDock?style=social" alt="GitHub Stars">
  </a>
  <a href="https://github.com/pritampanda15/PandaDock/issues">
    <img src="https://img.shields.io/github/issues/pritampanda15/PandaDock" alt="GitHub Issues">
  </a>
  <a href="https://pepy.tech/project/pandadock">
    <img src="https://static.pepy.tech/badge/pandadock" alt="Downloads">
  </a>
</p>

---

## 🚀 Overview

**PandaDock** is a modern, modular molecular docking toolkit that combines cutting-edge algorithms with high-performance computing capabilities. Built for drug discovery, computational chemistry, and bioinformatics workflows.

### ✨ Key Features

- 🎯 **Multiple Docking Algorithms**: Genetic Algorithm, Monte Carlo, PANDADOCK, Random Search, Metal-based docking
- ⚡ **Hardware Acceleration**: Native GPU (PyTorch/CUDA) and multi-core CPU parallelization
- 📊 **Advanced Scoring**: Physics-based scoring with MM-GBSA inspired energy decomposition
- 🧬 **Flexible Residues**: Automatic and manual flexible residue detection
- 📈 **Virtual Screening**: High-throughput screening with batch processing
- 🔬 **Comprehensive Analysis**: Binding affinity calculations, pose clustering, interaction analysis
- 🎨 **Rich Visualization**: Interactive HTML reports with energy breakdowns and plots
- 🛠️ **Modular Architecture**: Clean, extensible Python API for custom workflows

---

## 📦 Installation

### Quick Install

```bash
# Install PandaDock
pip install pandadock

# For GPU acceleration (optional)
pip install pandadock[gpu]

# For RDKit integration (recommended)
pip install pandadock[rdkit]

# Install all optional dependencies
pip install pandadock[gpu,rdkit,dev]
```

### Development Install

```bash
# Clone the repository
git clone https://github.com/pritampanda15/PandaDock.git
cd PandaDock

# Install in development mode
pip install -e .
```

### System Requirements

- **Python**: 3.8+
- **OS**: Linux, macOS, Windows
- **Memory**: 4GB+ RAM (8GB+ recommended)
- **GPU**: CUDA-compatible GPU (optional, for acceleration)

---

## 🔥 Quick Start

### Basic Docking

```bash
# Simple protein-ligand docking
pandadock -p protein.pdb -l ligand.sdf -o results

# With binding site specification
pandadock -p protein.pdb -l ligand.sdf -o results -s -15.7 -17.7 8.1 -r 10.0

# Enhanced docking with GPU acceleration
pandadock -p protein.pdb -l ligand.sdf -o results --enhanced --use-gpu
```

### Physics-Based Docking

```bash
# High-accuracy physics-based scoring
pandadock -p protein.pdb -l ligand.sdf -o results --physics-based

# PANDADOCK algorithm with simulated annealing
pandadock -p protein.pdb -l ligand.sdf -o results -a pandadock --physics-based

# With flexible residues
pandadock -p protein.pdb -l ligand.sdf -o results --physics-based --auto-flex
```

### Virtual Screening

```bash
# High-throughput screening
pandadock -p protein.pdb -l ligand_library/ -o screening_results --enhanced

# Batch processing with parallel execution
pandadock -p protein.pdb -l compounds.sdf -o batch_results --cpu-workers 8
```

---

## 🐍 Python API

### Basic Usage

```python
from pandadock.core import DockingEngine
from pandadock.molecules import LigandHandler, ProteinHandler

# Load molecules
protein_handler = ProteinHandler()
ligand_handler = LigandHandler()

protein = protein_handler.load_protein("protein.pdb")
ligand = ligand_handler.load_ligand("ligand.sdf")

# Configure docking
config = {
    'algorithm': 'genetic',
    'scoring': 'enhanced',
    'site_center': [10.0, 20.0, 30.0],
    'site_radius': 12.0,
    'iterations': 100
}

# Run docking
engine = DockingEngine(config)
results = engine.dock(protein, ligand)

# Analyze results
for pose in results['poses']:
    print(f"Pose {pose['rank']}: Score = {pose['score']:.2f}")
```

### Advanced Workflow

```python
from pandadock.scoring import ScoringFunctionFactory
from pandadock.algorithms import AlgorithmFactory
from pandadock.analysis import BindingAffinityCalculator

# Create custom scoring function
scoring_factory = ScoringFunctionFactory()
scorer = scoring_factory.create_scoring_function(
    type='physics-based',
    use_gpu=True
)

# Use specialized algorithm
algorithm_factory = AlgorithmFactory()
algorithm = algorithm_factory.create_algorithm(
    type='pandadock',
    population_size=200,
    temperature=500
)

# Calculate binding affinities
affinity_calc = BindingAffinityCalculator()
affinities = affinity_calc.calculate_batch_affinities(results['poses'])

for affinity in affinities:
    print(f"ΔG: {affinity.delta_g:.2f} kcal/mol, Kd: {affinity.kd:.2e} M")
```

---

## 📊 Output & Analysis

Each docking run generates comprehensive results:

```
results/
├── poses/                          # Individual pose PDB files
│   ├── pose_1_score_-14.25.pdb
│   └── pose_2_score_-14.07.pdb
├── complexes/                      # Protein-ligand complexes
│   └── complex_pose_1_score_-14.25.pdb
├── reports/                        # Analysis reports
│   ├── docking_report.html         # Interactive HTML report
│   ├── binding_affinity_report.csv
│   └── detailed_analysis.txt
├── plots/                          # Visualization plots
│   ├── score_distribution.png
│   ├── energy_breakdown.png
│   └── binding_modes.png
└── docking_results.json           # Machine-readable results
```

### Sample Analysis Report

```
🔬 DOCKING ANALYSIS REPORT
═══════════════════════════════════════

ALGORITHM: Genetic Algorithm (Enhanced)
SCORING: Physics-based with MM-GBSA
POSES GENERATED: 20
BEST SCORE: -14.25 kcal/mol

📊 BINDING AFFINITY ANALYSIS
─────────────────────────────
Pose 1: ΔG = -14.25 kcal/mol, Kd = 3.60e-11 M, IC50 = 7.20e-11 M
Pose 2: ΔG = -14.07 kcal/mol, Kd = 4.87e-11 M, IC50 = 9.75e-11 M
Pose 3: ΔG = -13.94 kcal/mol, Kd = 6.02e-11 M, IC50 = 1.20e-10 M

🧬 INTERACTION ANALYSIS
─────────────────────────
Hydrogen Bonds: 3
Hydrophobic Contacts: 12
π-π Stacking: 1
Salt Bridges: 0
```

---

## 🧠 Algorithms

### Available Algorithms

| Algorithm | Description | Best For |
|-----------|-------------|----------|
| **Genetic** | Evolutionary optimization | General docking, large ligands |
| **Monte Carlo** | Metropolis sampling | Conformational sampling |
| **PANDADOCK** | Simulated annealing + MD | High accuracy, flexible ligands |
| **Random Search** | Stochastic exploration | Quick screening, benchmarking |
| **Metal Docking** | Coordination constraints | Metalloproteins, metal cofactors |

### Scoring Functions

- **Basic**: VDW + Hydrogen bonding
- **Enhanced**: + Electrostatics + Desolvation + Hydrophobic
- **Physics-based**: Full MM-GBSA with entropy estimation
- **Metal-aware**: Coordination geometry constraints

---

## ⚙️ Configuration

### Command Line Options

```bash
# Essential options
-p, --protein PROTEIN       Protein PDB file
-l, --ligand LIGAND         Ligand MOL/SDF file
-o, --output OUTPUT         Output directory
-s, --site X Y Z            Binding site center
-r, --radius RADIUS         Binding site radius (Å)

# Algorithm selection
-a, --algorithm {genetic,monte-carlo,pandadock,random}
-i, --iterations N          Number of generations/steps
--exhaustiveness N          Independent runs

# Hardware acceleration
--use-gpu                   Enable GPU acceleration
--cpu-workers N             Number of CPU cores

# Scoring options
--enhanced-scoring          Enhanced scoring function
--physics-based            Physics-based scoring (MM-GBSA)
--local-opt                 Local optimization

# Flexibility
--auto-flex                 Automatic flexible residues
--flex-residues RES1 RES2   Manual flexible residues

# Output formats
--report-format {html,csv,json,all}
--detailed-energy           Energy component breakdown
```

---

## 🏗️ Architecture

PandaDock features a clean, modular architecture:

```
pandadock/
├── core/                   # Core docking engine
├── algorithms/             # Docking algorithms
├── scoring/               # Scoring functions
├── molecules/             # Molecule handling
├── analysis/              # Post-processing & analysis
├── hardware/              # GPU/CPU abstraction
├── screening/             # Virtual screening
├── cli/                   # Command-line interface
└── utils/                 # Utilities & helpers
```

---

## 🧪 Testing

```bash
# Run tests
pytest

# Run with coverage
pytest --cov=pandadock

# Run specific test categories
pytest tests/test_basic.py
pytest tests/test_scoring.py
```

---

## 📚 Documentation

- **GitHub Wiki**: Comprehensive guides and tutorials
- **API Reference**: Full Python API documentation
- **Examples**: Sample scripts and workflows
- **Paper**: [Cite PandaDock] (if published)

---

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md).

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📞 Support

- **GitHub Issues**: [Report bugs or request features](https://github.com/pritampanda15/PandaDock/issues)
- **Discussions**: [Community Q&A](https://github.com/pritampanda15/PandaDock/discussions)
- **Email**: pritam@stanford.edu

---

## 📖 Citation

If you use PandaDock in your research, please cite:

```bibtex
@software{pandadock2025,
  title={PandaDock: A Python-based Molecular Docking Platform},
  author={Panda, Pritam Kumar},
  year={2025},
  url={https://github.com/pritampanda15/PandaDock},
  version={3.0.0}
}
```

---

## 🙏 Acknowledgments

- RDKit community for cheminformatics tools
- PyTorch team for GPU acceleration framework
- AutoDock Vina for inspiration
- Scientific Python ecosystem

---

**⚗️ Dock Smarter. Discover Faster. Build Better.**