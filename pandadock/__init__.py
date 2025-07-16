"""
PandaDock: Modular, Multi-Strategy, High-Performance Molecular Docking Software

PandaDock is a Python-based molecular docking software that provides multiple
docking strategies including physics-based, machine learning-based, and genetic
algorithm-based approaches.

Author: Pritam Kumar Panda
Email: pritam@stanford.edu
Version: 3.0.0
"""

__version__ = "2.2.0"
__author__ = "Pritam Kumar Panda"
__email__ = "pritam@stanford.edu"

try:
    from .config import PandaDockConfig
    __all__ = [
        'PandaDockConfig',
        '__version__',
        '__author__',
        '__email__'
    ]
except ImportError:
    # Allow basic package import even if config fails
    __all__ = [
        '__version__',
        '__author__',
        '__email__'
    ]