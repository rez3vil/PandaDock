"""
PandaDock: Modular, Multi-Strategy, High-Performance Molecular Docking Software

PandaDock is a Python-based molecular docking software that provides multiple
docking strategies including physics-based, machine learning-based, and genetic
algorithm-based approaches.

Author: Pritam Kumar Panda
Email: pritam@stanford.edu
Version: 3.0.0
"""

__version__ = "3.0.0"
__author__ = "Pritam Kumar Panda"
__email__ = "pritam@stanford.edu"

from .config import PandaDockConfig

__all__ = [
    'PandaDockConfig',
    '__version__',
    '__author__',
    '__email__'
]