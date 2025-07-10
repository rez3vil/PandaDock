# PandaDock Installation Guide

## Quick Installation

### Option 1: Basic Installation (Recommended)
Install PandaDock with core dependencies only:

```bash
pip install -e .
```

### Option 2: Install with Specific Features
Install PandaDock with optional feature sets:

```bash
# For machine learning capabilities
pip install -e .[ml]

# For visualization support  
pip install -e .[viz]

# For chemistry libraries
pip install -e .[chem]

# For GPU acceleration (Linux/Windows only)
pip install -e .[gpu]

# For development
pip install -e .[dev]

# Install all optional dependencies
pip install -e .[all]
```

### Option 3: Manual Dependency Installation
If you encounter issues with the automatic installation, install dependencies manually:

```bash
# Install core dependencies first
pip install -r requirements-core.txt

# Then install PandaDock
pip install -e . --no-deps
```

## Installation Options by Use Case

### For Basic Molecular Docking
```bash
pip install -e .[chem]
```

### For Machine Learning-based Docking
```bash
pip install -e .[ml,chem]
```

### For Visualization and Analysis
```bash
pip install -e .[viz,chem]
```

### For Full-featured Installation
```bash
pip install -e .[all]
```

### For Development
```bash
pip install -e .[dev,all]
```

## Troubleshooting

### Common Issues

#### 1. PyMOL Installation Fails
PyMOL can be problematic via pip. Use conda instead:
```bash
conda install -c conda-forge pymol-open-source
pip install -e .[ml,chem] # Install without viz group
```

#### 2. RDKit Installation Issues
For RDKit problems, use conda:
```bash
conda install -c conda-forge rdkit
pip install -e . --no-deps
pip install -r requirements-core.txt
```

#### 3. GPU Installation on macOS
GPU acceleration is not available on macOS. Use:
```bash
pip install -e .[ml,chem,viz] # Skip GPU dependencies
```

#### 4. Commercial Software Dependencies
Commercial packages (OpenEye, Schrödinger) require separate licenses:
```bash
pip install -e .[commercial] # Only if you have licenses
```

#### 5. CI/CD and GitHub Actions
For automated builds, use the CI-safe requirements:
```bash
pip install -r requirements-ci.txt
pip install -e . --no-deps
```

## Verification

Test your installation:

```bash
python -c "import pandadock; print('PandaDock installed successfully!')"
pandadock --help
```

## Dependencies

### Core Dependencies (Always Required)
- numpy>=1.21.0
- scipy>=1.7.0
- pandas>=1.3.0
- scikit-learn>=1.0.0
- matplotlib>=3.5.0
- pyyaml>=5.4.0
- click>=8.0.0
- tqdm>=4.62.0
- psutil>=5.8.0
- requests>=2.25.0

### Optional Dependencies

#### Machine Learning (`[ml]`)
- torch>=1.9.0
- torch-geometric>=2.0.0
- transformers>=4.0.0
- xgboost>=1.5.0

#### Chemistry (`[chem]`)
- rdkit>=2022.03.1
- biopython>=1.79
- openmm>=7.7.0

#### Visualization (`[viz]`)
- pymol-open-source>=2.5.0 (install via conda if pip fails)
- py3dmol>=1.8.0
- nglview>=3.0.0

#### GPU Acceleration (`[gpu]`)
- cupy-cuda11x>=9.0.0 (Linux/Windows only)

#### Web Interface (`[web]`)
- flask>=2.0.0
- fastapi>=0.70.0
- uvicorn>=0.15.0

## System Requirements

- Python 3.8 or higher
- 4GB RAM minimum (8GB+ recommended)
- 2GB disk space
- CUDA-compatible GPU (optional, for GPU acceleration)

## Platform-Specific Notes

### Linux
All features supported. For GPU acceleration:
```bash
pip install -e .[gpu,all]
```

### macOS
GPU acceleration not available. Use:
```bash
pip install -e .[ml,chem,viz,web]
```

### Windows
All features supported. May require Visual Studio Build Tools for some packages.

## Virtual Environment (Recommended)

Create an isolated environment:

```bash
python -m venv pandadock-env
source pandadock-env/bin/activate  # Linux/macOS
# or
pandadock-env\Scripts\activate     # Windows

pip install -e .[all]
```