# -*- coding: utf-8 -*-
"""
Setup script for PandaDock
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read core requirements
def read_requirements(filename='requirements-core.txt'):
    requirements = []
    req_file = os.path.join(this_directory, filename)
    
    if not os.path.exists(req_file):
        # Fallback to basic requirements if file doesn't exist
        return [
            'numpy>=1.21.0',
            'scipy>=1.7.0',
            'pandas>=1.3.0',
            'scikit-learn>=1.0.0',
            'matplotlib>=3.5.0',
            'pyyaml>=5.4.0',
            'click>=8.0.0',
            'tqdm>=4.62.0',
            'psutil>=5.8.0',
            'requests>=2.25.0'
        ]
    
    with open(req_file, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if line and not line.startswith('#'):
                # Remove inline comments and handle environment markers properly
                if ';' in line:
                    # Keep the environment marker for proper handling
                    req_line = line
                else:
                    req_line = line
                
                if req_line:
                    requirements.append(req_line)
    return requirements

requirements = read_requirements()

setup(
    name="pandadock",
    version="2.4.0",
    author="Pritam Kumar Panda",
    author_email="pritam@stanford.edu",
    description="Modular, Multi-Strategy, High-Performance Molecular Docking Software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pritampanda15/pandadock",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "gpu": ["cupy-cuda11x>=9.0.0; platform_system != 'Darwin'"],
        "ml": [
            "torch>=1.9.0", 
            "torch-geometric>=2.0.0", 
            "transformers>=4.0.0",
            "xgboost>=1.5.0"
        ],
        "viz": [
            "pymol-open-source>=2.5.0", 
            "rdkit>=2022.03.1",
            "py3dmol>=1.8.0",
            "nglview>=3.0.0"
        ],
        "chem": [
            "rdkit>=2022.03.1",
            "biopython>=1.79",
            "openmm>=7.7.0"
        ],
        "web": [
            "flask>=2.0.0",
            "fastapi>=0.70.0",
            "uvicorn>=0.15.0"
        ],
        "dev": [
            "pytest>=6.0.0", 
            "pytest-cov>=3.0.0",
            "black>=22.0.0", 
            "flake8>=4.0.0", 
            "mypy>=0.910",
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0"
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "sphinx-autodoc-typehints>=1.12.0",
            "myst-parser>=0.18.0"
        ],
        "commercial": [
            "openeye-toolkits>=2021.10.0",
            "schrodinger-suite>=2021.4"
        ],
        "all": [
            "cupy-cuda11x>=9.0.0; platform_system != 'Darwin'",
            "torch>=1.9.0",
            "torch-geometric>=2.0.0", 
            "transformers>=4.0.0",
            "xgboost>=1.5.0",
            "rdkit>=2022.03.1",
            "biopython>=1.79",
            "openmm>=7.7.0",
            "pymol-open-source>=2.5.0",
            "py3dmol>=1.8.0",
            "nglview>=3.0.0",
            "flask>=2.0.0",
            "fastapi>=0.70.0"
        ]
    },
    entry_points={
        "console_scripts": [
            "pandadock=pandadock.__main__:main",
        ],
    },
    include_package_data=True,
    package_data={
        "pandadock": [
            "data/*.json",
            "data/*.csv",
            "templates/*.html",
            "examples/*.py",
            "examples/*.pdb",
            "examples/*.sdf",
        ],
    },
    project_urls={
        "Documentation": "https://pandadock.readthedocs.io/en/latest/",
        "Source": "https://github.com/pritampkp15/pandadock",
        "Tracker": "https://github.com/pritampkp15/pandadock/issues",
    },
    keywords=[
        "molecular docking",
        "drug discovery",
        "computational chemistry",
        "protein-ligand binding",
        "virtual screening",
        "machine learning",
        "genetic algorithm",
        "physics-based modeling",
    ],
    zip_safe=False,
)