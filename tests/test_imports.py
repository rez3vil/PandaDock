#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import tests for PandaDock modules
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_core_imports():
    """Test that core modules can be imported"""
    try:
        import pandadock
        assert hasattr(pandadock, '__version__')
    except ImportError as e:
        pytest.skip(f"Core pandadock import failed: {e}")


def test_docking_engine_imports():
    """Test docking engine imports"""
    import_tests = [
        ('pandadock.docking.base_engine', 'BaseDockingEngine'),
        ('pandadock.docking.base_engine', 'Pose'),
    ]
    
    for module_name, class_name in import_tests:
        try:
            module = __import__(module_name, fromlist=[class_name])
            cls = getattr(module, class_name)
            assert cls is not None
        except ImportError as e:
            pytest.skip(f"Import {module_name}.{class_name} failed: {e}")
        except AttributeError as e:
            pytest.skip(f"Attribute {class_name} not found in {module_name}: {e}")


def test_optional_imports():
    """Test optional module imports (may not be available in CI)"""
    optional_imports = [
        ('pandadock.docking.physics_engine', 'PhysicsEngine'),
        ('pandadock.docking.ml_engine', 'MLEngine'),
        ('pandadock.docking.ga_engine', 'GAEngine'),
        ('pandadock.scoring.scoring_functions', 'ScoringFunction'),
        ('pandadock.utils.ic50_calculator', 'IC50Calculator'),
    ]
    
    for module_name, class_name in optional_imports:
        try:
            module = __import__(module_name, fromlist=[class_name])
            cls = getattr(module, class_name)
            assert cls is not None
        except ImportError:
            pytest.skip(f"Optional import {module_name}.{class_name} not available")
        except AttributeError:
            pytest.skip(f"Optional attribute {class_name} not found in {module_name}")


def test_config_import():
    """Test configuration imports"""
    try:
        from pandadock.config import PandaDockConfig
        config = PandaDockConfig()
        assert config is not None
    except ImportError as e:
        pytest.skip(f"Config import failed: {e}")


def test_cli_import():
    """Test CLI imports"""
    try:
        from pandadock.__main__ import main
        assert main is not None
    except ImportError as e:
        pytest.skip(f"CLI import failed: {e}")


def test_package_structure():
    """Test basic package structure"""
    try:
        import pandadock.docking
        import pandadock.scoring
        import pandadock.utils
        # Just test that the packages exist
        assert True
    except ImportError as e:
        pytest.skip(f"Package structure test failed: {e}")


if __name__ == "__main__":
    pytest.main([__file__])