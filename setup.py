"""
Setup script for PandaDock
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open(os.path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name="pandadock",
    version="1.0.0",
    author="PandaDock Development Team",
    author_email="info@pandadock.org",
    description="Modular, Multi-Strategy, High-Performance Molecular Docking Software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pandadock/pandadock",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
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
        "gpu": ["cupy-cuda11x>=9.0.0"],
        "ml": ["torch>=1.9.0", "torch-geometric>=2.0.0", "transformers>=4.0.0"],
        "viz": ["pymol-open-source>=2.5.0", "rdkit>=2022.03.1"],
        "dev": ["pytest>=6.0.0", "black>=22.0.0", "flake8>=4.0.0", "mypy>=0.910"],
    },
    entry_points={
        "console_scripts": [
            "pandadock=pandadock.main:main",
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
        "Documentation": "https://pandadock.readthedocs.io/",
        "Source": "https://github.com/pandadock/pandadock",
        "Tracker": "https://github.com/pandadock/pandadock/issues",
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