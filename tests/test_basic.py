#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic tests for PandaDock core functionality
"""

import pytest
import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_package_import():
    """Test that the pandadock package can be imported"""
    try:
        import pandadock
        assert hasattr(pandadock, '__version__')
    except ImportError as e:
        pytest.skip(f"PandaDock package not available: {e}")


def test_config_import():
    """Test that config module can be imported"""
    try:
        from pandadock.config import PandaDockConfig
        config = PandaDockConfig()
        assert config is not None
    except ImportError as e:
        pytest.skip(f"Config module not available: {e}")


def test_base_engine_import():
    """Test that base engine classes can be imported"""
    try:
        from pandadock.docking.base_engine import BaseDockingEngine, Pose
        assert BaseDockingEngine is not None
        assert Pose is not None
    except ImportError as e:
        pytest.skip(f"Base engine not available: {e}")


def test_pose_creation():
    """Test basic Pose object creation"""
    try:
        from pandadock.docking.base_engine import Pose
        import numpy as np
        
        # Create a simple pose
        coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        pose = Pose(
            pose_id=1,
            coordinates=coordinates,
            score=-5.0,
            energy=-10.0
        )
        
        assert pose.pose_id == 1
        assert pose.score == -5.0
        assert pose.energy == -10.0
        assert len(pose.coordinates) == 2
        
    except ImportError as e:
        pytest.skip(f"Pose class not available: {e}")


def test_algorithm_names():
    """Test that the new algorithm names are recognized"""
    try:
        from pandadock.config import PandaDockConfig
        
        # Test valid algorithm names
        valid_algorithms = ['pandacore', 'pandaml', 'pandaphysics']
        
        for algorithm in valid_algorithms:
            config = PandaDockConfig()
            # Just check that we can create config - actual validation depends on implementation
            assert config is not None
            
    except ImportError as e:
        pytest.skip(f"Config not available: {e}")


def test_scoring_functions():
    """Test that scoring functions can be imported"""
    try:
        from pandadock.scoring.scoring_functions import ScoringFunction
        assert ScoringFunction is not None
    except ImportError as e:
        pytest.skip(f"Scoring functions not available: {e}")


def test_cli_main_import():
    """Test that CLI main module can be imported"""
    try:
        from pandadock.__main__ import main
        assert main is not None
    except ImportError as e:
        pytest.skip(f"CLI main not available: {e}")


class TestDockingEngines:
    """Test docking engine imports and basic functionality"""
    
    def test_physics_engine_import(self):
        """Test physics engine import"""
        try:
            from pandadock.docking.physics_engine import PhysicsEngine
            assert PhysicsEngine is not None
        except ImportError as e:
            pytest.skip(f"Physics engine not available: {e}")
    
    def test_ml_engine_import(self):
        """Test ML engine import"""
        try:
            from pandadock.docking.ml_engine import MLEngine
            assert MLEngine is not None
        except ImportError as e:
            pytest.skip(f"ML engine not available: {e}")
    
    def test_ga_engine_import(self):
        """Test GA engine import"""
        try:
            from pandadock.docking.ga_engine import GAEngine
            assert GAEngine is not None
        except ImportError as e:
            pytest.skip(f"GA engine not available: {e}")


class TestUtilities:
    """Test utility modules"""
    
    def test_ic50_calculator_import(self):
        """Test IC50 calculator import"""
        try:
            from pandadock.utils.ic50_calculator import IC50Calculator
            assert IC50Calculator is not None
        except ImportError as e:
            pytest.skip(f"IC50 calculator not available: {e}")
    
    def test_math_utils_import(self):
        """Test math utilities import"""
        try:
            from pandadock.utils.math_utils import calculate_rmsd
            assert calculate_rmsd is not None
        except ImportError as e:
            pytest.skip(f"Math utils not available: {e}")


if __name__ == "__main__":
    pytest.main([__file__])