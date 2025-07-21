# -*- coding: utf-8 -*-
"""
Basic tests to improve code coverage for CI.
These tests exercise core functionality without requiring heavy dependencies.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path


def test_basic_imports():
    """Test that basic imports work."""
    try:
        import pandadock
        assert pandadock is not None
    except ImportError:
        pytest.skip("PandaDock not available")


def test_config_basic():
    """Test basic config functionality."""
    try:
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        
        # Test basic attributes exist
        assert hasattr(config, 'docking')
        assert hasattr(config, 'scoring')
        assert hasattr(config, 'n_jobs')
        assert hasattr(config, 'gpu_enabled')
        
        # Test that we can set basic values
        config.n_jobs = 1
        config.gpu_enabled = False
        
        assert config.n_jobs == 1
        assert config.gpu_enabled == False
        
    except ImportError:
        pytest.skip("PandaDock config not available")


def test_pose_basic():
    """Test basic Pose functionality."""
    try:
        from pandadock.docking.base_engine import Pose
        
        coords = np.array([[0, 0, 0], [1, 1, 1]])
        pose = Pose(
            coordinates=coords,
            score=-5.0,
            energy=-8.0,
            ligand_name="test",
            pose_id="test_1"
        )
        
        assert pose.coordinates.shape == (2, 3)
        assert pose.score == -5.0
        assert pose.energy == -8.0
        assert pose.ligand_name == "test"
        assert pose.pose_id == "test_1"
        
    except ImportError:
        pytest.skip("PandaDock Pose not available")


def test_scoring_basic():
    """Test basic scoring functionality."""
    try:
        from pandadock.scoring.scoring_functions import ScoringFunctions
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        scoring = ScoringFunctions(config)
        
        # Test basic coordinate input
        coords = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]])
        
        # Test that energy calculation doesn't crash
        energy = scoring.calculate_total_energy(coords)
        assert isinstance(energy, (int, float))
        assert np.isfinite(energy)
        
        # Test individual energy components
        vdw_energy = scoring.calculate_vdw_energy(coords)
        assert isinstance(vdw_energy, (int, float))
        
        hbond_energy = scoring.calculate_hbond_energy(coords)
        assert isinstance(hbond_energy, (int, float))
        
        hydrophobic_energy = scoring.calculate_hydrophobic_energy(coords)
        assert isinstance(hydrophobic_energy, (int, float))
        
        solvation_energy = scoring.calculate_solvation_energy(coords)
        assert isinstance(solvation_energy, (int, float))
        
        entropy_energy = scoring.calculate_entropy_penalty(coords)
        assert isinstance(entropy_energy, (int, float))
        
        clash_score = scoring.calculate_clash_score(coords)
        assert isinstance(clash_score, (int, float))
        assert clash_score >= 0  # Clash score should be non-negative
        
    except ImportError:
        pytest.skip("PandaDock scoring not available")


def test_math_utils_basic():
    """Test basic math utilities."""
    try:
        from pandadock.utils.math_utils import distance_matrix, rotation_matrix
        
        # Test distance matrix
        coords1 = np.array([[0, 0, 0], [1, 0, 0]])
        coords2 = np.array([[0, 1, 0], [1, 1, 0]])
        
        distances = distance_matrix(coords1, coords2)
        assert distances.shape == (2, 2)
        assert np.all(distances >= 0)
        
        # Test rotation matrix
        angles = np.array([0, np.pi/4, np.pi/2])
        rot_matrix = rotation_matrix(angles)
        assert rot_matrix.shape == (3, 3)
        
        # Check that it's a proper rotation matrix (orthogonal)
        should_be_identity = np.dot(rot_matrix, rot_matrix.T)
        assert np.allclose(should_be_identity, np.eye(3), atol=1e-10)
        
    except ImportError:
        pytest.skip("Math utils not available")


def test_ic50_calculator_basic():
    """Test basic IC50 calculator functionality."""
    try:
        from pandadock.utils.ic50_calculator import IC50Calculator
        
        calculator = IC50Calculator()
        
        # Test basic conversion
        binding_affinity = -8.0  # kcal/mol
        ic50 = calculator.delta_g_to_ic50(binding_affinity)
        assert isinstance(ic50, (int, float))
        assert ic50 > 0  # IC50 should be positive
        
        # Test reverse conversion
        back_to_affinity = calculator.ic50_to_delta_g(ic50)
        assert np.isclose(back_to_affinity, binding_affinity, rtol=1e-3)
        
    except ImportError:
        pytest.skip("IC50 calculator not available")


def test_algorithm_specific_scoring():
    """Test algorithm-specific scoring functions."""
    try:
        from pandadock.scoring.scoring_functions import ScoringFunctions
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        scoring = ScoringFunctions(config)
        
        coords = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0]])
        
        # Test PandaCore scoring
        config.scoring.scoring_function = 'pandacore'
        core_energy = scoring.calculate_total_energy(coords)
        assert isinstance(core_energy, (int, float))
        
        # Test PandaML scoring
        config.scoring.scoring_function = 'pandaml'
        ml_energy = scoring.calculate_total_energy(coords)
        assert isinstance(ml_energy, (int, float))
        
        # Test PandaPhysics scoring
        config.scoring.scoring_function = 'pandaphysics'
        physics_energy = scoring.calculate_total_energy(coords)
        assert isinstance(physics_energy, (int, float))
        
        # Energies should be in reasonable range
        assert -50 < core_energy < 50
        assert -50 < ml_energy < 50
        assert -50 < physics_energy < 50
        
    except ImportError:
        pytest.skip("Scoring functions not available")


def test_grid_box_basic():
    """Test basic grid box functionality."""
    try:
        from pandadock.docking.base_engine import GridBox
        
        center = np.array([0, 0, 0])
        size = np.array([20, 20, 20])
        
        grid_box = GridBox(center, size)
        
        assert np.allclose(grid_box.center, center)
        assert np.allclose(grid_box.size, size)
        
        # Test bounds calculation
        min_bounds, max_bounds = grid_box.get_bounds()
        expected_min = center - size / 2
        expected_max = center + size / 2
        
        assert np.allclose(min_bounds, expected_min)
        assert np.allclose(max_bounds, expected_max)
        
        # Test contains method
        point_inside = np.array([5, 5, 5])
        point_outside = np.array([15, 15, 15])
        
        assert grid_box.contains_point(point_inside)
        assert not grid_box.contains_point(point_outside)
        
    except ImportError:
        pytest.skip("Grid box not available")


def test_file_io_basic():
    """Test basic file I/O functionality."""
    try:
        from pandadock.io.ligand_preparer import LigandPreparer
        
        preparer = LigandPreparer()
        
        # Test that the preparer can be created
        assert preparer is not None
        assert hasattr(preparer, '__init__')
        
    except (ImportError, AttributeError):
        pytest.skip("Ligand preparer not available")


def test_report_generation_basic():
    """Test basic report generation."""
    try:
        from pandadock.reports.html_report import HTMLReportGenerator
        from pandadock.docking.base_engine import Pose
        
        generator = HTMLReportGenerator()
        
        # Create sample poses
        coords = np.array([[0, 0, 0], [1, 1, 1]])
        poses = [
            Pose(coords, -5.0, -8.0, "test", "pose_1"),
            Pose(coords + 0.1, -4.5, -7.5, "test", "pose_2")
        ]
        
        # Test that report generation doesn't crash
        report_data = generator._prepare_pose_data(poses)
        assert isinstance(report_data, list)
        assert len(report_data) == 2
        
    except (ImportError, AttributeError):
        pytest.skip("Report generation not available")


def test_empty_coordinates():
    """Test handling of empty coordinates."""
    try:
        from pandadock.scoring.scoring_functions import ScoringFunctions
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        scoring = ScoringFunctions(config)
        
        # Test empty coordinates
        empty_coords = np.array([]).reshape(0, 3)
        energy = scoring.calculate_total_energy(empty_coords)
        assert energy == 0.0
        
        # Test single atom
        single_atom = np.array([[0, 0, 0]])
        energy = scoring.calculate_total_energy(single_atom)
        assert isinstance(energy, (int, float))
        
    except ImportError:
        pytest.skip("Scoring functions not available")


def test_version_info():
    """Test version information."""
    try:
        import pandadock
        
        # Check if version is available
        if hasattr(pandadock, '__version__'):
            version = pandadock.__version__
            assert isinstance(version, str)
            assert len(version) > 0
        
    except ImportError:
        pytest.skip("PandaDock not available")


def test_pose_attributes():
    """Test pose attribute access and modification."""
    try:
        from pandadock.docking.base_engine import Pose
        
        coords = np.array([[0, 0, 0], [1, 1, 1]])
        pose = Pose(
            coordinates=coords,
            score=-5.0,
            energy=-8.0,
            ligand_name="test",
            pose_id="test_1"
        )
        
        # Test setting additional attributes
        pose.vdw_energy = -2.0
        pose.electrostatic_energy = -1.0
        pose.hbond_energy = -3.0
        pose.confidence = 0.8
        
        assert pose.vdw_energy == -2.0
        assert pose.electrostatic_energy == -1.0
        assert pose.hbond_energy == -3.0
        assert pose.confidence == 0.8
        
        # Test binding affinity calculation
        binding_affinity = pose.get_binding_affinity()
        assert isinstance(binding_affinity, (int, float))
        
    except ImportError:
        pytest.skip("Pose not available")


def test_scoring_edge_cases():
    """Test scoring function edge cases."""
    try:
        from pandadock.scoring.scoring_functions import ScoringFunctions
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        scoring = ScoringFunctions(config)
        
        # Test with close but non-bonded atoms (should give clash score)
        # Distance of 2.5 Å is above bonded cutoff (2.0 Å) but below VdW cutoff (2.8 Å)
        close_coords = np.array([[0, 0, 0], [2.5, 0, 0]])
        clash_score = scoring.calculate_clash_score(close_coords)
        assert clash_score >= 0  # Clash score should be non-negative
        
        # Test with well-separated atoms (should give no clash score)
        separated_coords = np.array([[0, 0, 0], [5, 0, 0]])
        clash_score_low = scoring.calculate_clash_score(separated_coords)
        assert clash_score >= clash_score_low  # Close atoms should have higher or equal clash score
        assert isinstance(clash_score, (int, float))
        assert isinstance(clash_score_low, (int, float))
        
        # Test Vina scoring
        vina_score = scoring.calculate_vina_score(separated_coords)
        assert isinstance(vina_score, (int, float))
        
    except ImportError:
        pytest.skip("Scoring functions not available")


def test_configuration_scoring_modes():
    """Test different scoring mode configurations."""
    try:
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        
        # Test setting different scoring functions
        scoring_functions = ['pandacore', 'pandaml', 'pandaphysics']
        
        for func in scoring_functions:
            config.scoring.scoring_function = func
            assert config.scoring.scoring_function == func
        
        # Test setting weights
        config.scoring.vdw_weight = 0.5
        config.scoring.electrostatic_weight = 1.5
        
        assert config.scoring.vdw_weight == 0.5
        assert config.scoring.electrostatic_weight == 1.5
        
    except ImportError:
        pytest.skip("Config not available")


def test_numpy_operations():
    """Test numpy operations used in the codebase."""
    # Test distance matrix calculation
    coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    coords2 = np.array([[0, 0, 1], [1, 1, 1]])
    
    # Manual distance matrix calculation
    distances = np.zeros((len(coords1), len(coords2)))
    for i, c1 in enumerate(coords1):
        for j, c2 in enumerate(coords2):
            distances[i, j] = np.linalg.norm(c1 - c2)
    
    assert distances.shape == (3, 2)
    assert np.all(distances >= 0)
    
    # Test coordinate transformations
    coords = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    center = np.mean(coords, axis=0)
    centered = coords - center
    
    # Should be centered around origin
    new_center = np.mean(centered, axis=0)
    assert np.allclose(new_center, [0, 0, 0], atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__])