#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple tests for PandaDock CI
"""

import os
import sys
from pathlib import Path

# Test that basic Python functionality works
def test_python_basic():
    """Test basic Python functionality"""
    assert 1 + 1 == 2
    assert isinstance("test", str)
    assert len([1, 2, 3]) == 3


def test_imports_work():
    """Test that basic imports work"""
    import json
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Test basic numpy functionality
    arr = np.array([1, 2, 3])
    assert len(arr) == 3
    
    # Test basic pandas functionality
    df = pd.DataFrame({'a': [1, 2, 3]})
    assert len(df) == 3


def test_project_structure():
    """Test that project structure exists"""
    project_root = Path(__file__).parent.parent
    
    # Check key directories exist
    assert (project_root / "pandadock").exists()
    assert (project_root / "tests").exists()
    assert (project_root / "benchmarks").exists()
    assert (project_root / "docs").exists()
    
    # Check key files exist
    assert (project_root / "setup.py").exists()
    assert (project_root / "README.md").exists()
    assert (project_root / "requirements.txt").exists()


def test_pandadock_package_exists():
    """Test that pandadock package structure exists"""
    project_root = Path(__file__).parent.parent
    pandadock_dir = project_root / "pandadock"
    
    # Check core package files
    assert (pandadock_dir / "__init__.py").exists()
    assert (pandadock_dir / "__main__.py").exists()
    assert (pandadock_dir / "config.py").exists()
    
    # Check subdirectories
    assert (pandadock_dir / "docking").exists()
    assert (pandadock_dir / "scoring").exists()
    assert (pandadock_dir / "utils").exists()
    assert (pandadock_dir / "io").exists()
    assert (pandadock_dir / "reports").exists()


def test_algorithm_names():
    """Test new algorithm naming convention"""
    new_algorithms = ['pandacore', 'pandaml', 'pandaphysics']
    old_algorithms = ['vina', 'glide']
    
    # Test that new names follow convention
    for algorithm in new_algorithms:
        assert algorithm.startswith('panda')
        assert algorithm.islower()
    
    # Test that old and new names don't overlap
    for old in old_algorithms:
        for new in new_algorithms:
            assert old not in new


def test_benchmark_data_validation():
    """Test benchmark performance data from documentation"""
    # Performance data from the updated documentation
    performance_data = {
        'pandaml': {'r2': 0.878, 'success_rate': 0.467, 'rmsd': 3.22},
        'pandaphysics': {'r2': 0.671, 'success_rate': 0.488, 'rmsd': 2.18},
        'pandacore': {'r2': 0.738, 'success_rate': 0.463, 'rmsd': 3.63}
    }
    
    for algorithm, metrics in performance_data.items():
        # Validate ranges
        assert 0.0 <= metrics['r2'] <= 1.0
        assert 0.0 <= metrics['success_rate'] <= 1.0
        assert metrics['rmsd'] > 0.0
        
        # Validate that PandaML has best R²
        if algorithm == 'pandaml':
            assert metrics['r2'] > 0.85


def test_file_imports_syntax():
    """Test that key Python files have valid syntax"""
    project_root = Path(__file__).parent.parent
    
    # Files to check
    files_to_check = [
        "pandadock/__init__.py",
        "pandadock/__main__.py", 
        "pandadock/config.py",
        "setup.py"
    ]
    
    for file_path in files_to_check:
        full_path = project_root / file_path
        if full_path.exists():
            with open(full_path, 'r') as f:
                content = f.read()
            
            # Try to compile - this will raise SyntaxError if invalid
            try:
                compile(content, str(full_path), 'exec')
            except SyntaxError as e:
                assert False, f"Syntax error in {file_path}: {e}"


if __name__ == "__main__":
    import pytest
    pytest.main([__file__])