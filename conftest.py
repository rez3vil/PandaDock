# -*- coding: utf-8 -*-
"""
Pytest configuration file for PandaDock test suite.

This file contains shared fixtures, configuration settings,
and test utilities used across all test files.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from typing import Dict, Any, List
import shutil

# Import PandaDock modules for testing
try:
    from pandadock.config import PandaDockConfig
    from pandadock.docking.base_engine import Pose
    from pandadock.docking.ga_engine import GAEngine
    from pandadock.docking.ml_engine import MLEngine
    from pandadock.docking.physics_engine import PhysicsEngine
    from pandadock.scoring.scoring_functions import ScoringFunctions
except ImportError:
    # Handle missing imports gracefully for CI environments
    PandaDockConfig = None
    Pose = None
    GAEngine = None
    MLEngine = None
    PhysicsEngine = None
    ScoringFunctions = None


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "gpu: marks tests requiring GPU support"
    )
    config.addinivalue_line(
        "markers", "ml: marks tests requiring ML dependencies"
    )
    config.addinivalue_line(
        "markers", "chem: marks tests requiring chemistry libraries"
    )
    config.addinivalue_line(
        "markers", "viz: marks tests requiring visualization libraries"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test names."""
    for item in items:
        # Mark slow tests
        if "slow" in item.name or "benchmark" in item.name:
            item.add_marker(pytest.mark.slow)
        
        # Mark integration tests
        if "integration" in item.name or "test_integration" in str(item.fspath):
            item.add_marker(pytest.mark.integration)
        
        # Mark GPU tests
        if "gpu" in item.name or "cuda" in item.name:
            item.add_marker(pytest.mark.gpu)
        
        # Mark ML tests
        if any(keyword in item.name for keyword in ["ml", "diffusion", "transformer"]):
            item.add_marker(pytest.mark.ml)


# Fixtures for test configuration
@pytest.fixture(scope="session")
def test_config():
    """Create a test configuration for PandaDock."""
    if PandaDockConfig is None:
        pytest.skip("PandaDock not available")
    
    config = PandaDockConfig()
    
    # Override settings for testing
    config.docking.num_poses = 5
    config.docking.exhaustiveness = 2
    config.docking.population_size = 20
    config.docking.generations = 10
    config.docking.minimization_steps = 50
    config.scoring.vdw_weight = 1.0
    config.scoring.electrostatic_weight = 1.0
    config.scoring.hbond_weight = 1.0
    config.gpu_enabled = False
    config.n_jobs = 1
    
    return config


@pytest.fixture(scope="session")
def temp_dir():
    """Create a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp(prefix="pandadock_test_")
    yield Path(temp_dir)
    # Cleanup
    shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.fixture
def sample_ligand_coords():
    """Generate sample ligand coordinates for testing."""
    # Create a small molecule with realistic coordinates
    coords = np.array([
        [0.0, 0.0, 0.0],      # Central atom
        [1.5, 0.0, 0.0],      # Bonded atom
        [-1.5, 0.0, 0.0],     # Bonded atom
        [0.0, 1.5, 0.0],      # Bonded atom
        [0.0, 0.0, 1.5],      # Bonded atom
        [2.5, 0.5, 0.5],      # Extended atom
        [-2.5, -0.5, -0.5],   # Extended atom
    ])
    return coords


@pytest.fixture
def sample_protein_coords():
    """Generate sample protein coordinates for testing."""
    # Create a small protein-like structure
    coords = np.array([
        [10.0, 10.0, 10.0],   # Backbone atom
        [11.5, 10.0, 10.0],   # Side chain
        [10.0, 11.5, 10.0],   # Backbone atom
        [10.0, 10.0, 11.5],   # Side chain
        [8.5, 10.0, 10.0],    # Backbone atom
        [10.0, 8.5, 10.0],    # Side chain
        [12.0, 12.0, 12.0],   # Distant atom
        [8.0, 8.0, 8.0],      # Distant atom
    ])
    return coords


@pytest.fixture
def sample_pose(sample_ligand_coords):
    """Create a sample pose for testing."""
    if Pose is None:
        pytest.skip("PandaDock not available")
    
    pose = Pose(
        coordinates=sample_ligand_coords,
        score=-5.2,
        energy=-8.7,
        ligand_name="test_ligand",
        pose_id="test_pose_1"
    )
    
    # Add some realistic energy terms
    pose.vdw_energy = -3.2
    pose.electrostatic_energy = -1.5
    pose.hbond_energy = -2.8
    pose.hydrophobic_energy = -1.2
    pose.solvation_energy = 1.8
    pose.entropy_energy = 2.3
    pose.confidence = 0.85
    
    return pose


@pytest.fixture
def sample_poses(sample_ligand_coords):
    """Create multiple sample poses for testing."""
    if Pose is None:
        pytest.skip("PandaDock not available")
    
    poses = []
    for i in range(5):
        # Slightly perturb coordinates for each pose
        perturbed_coords = sample_ligand_coords + np.random.normal(0, 0.1, sample_ligand_coords.shape)
        
        pose = Pose(
            coordinates=perturbed_coords,
            score=-5.0 + i * 0.5,
            energy=-8.0 + i * 0.3,
            ligand_name="test_ligand",
            pose_id=f"test_pose_{i+1}"
        )
        
        pose.confidence = 0.9 - i * 0.1
        poses.append(pose)
    
    return poses


@pytest.fixture
def sample_pdb_content():
    """Sample PDB file content for testing."""
    return """HEADER    TEST PROTEIN                            01-JAN-25   TEST            
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1      11.500  10.000  10.000  1.00 20.00           C  
ATOM      3  C   ALA A   1      10.000  11.500  10.000  1.00 20.00           C  
ATOM      4  O   ALA A   1      10.000  10.000  11.500  1.00 20.00           O  
ATOM      5  CB  ALA A   1       8.500  10.000  10.000  1.00 20.00           C  
TER       6      ALA A   1                                                      
END                                                                             
"""


@pytest.fixture
def sample_sdf_content():
    """Sample SDF file content for testing."""
    return """
  Test Molecule
  
  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.5000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.5000    0.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5000   -0.5000   -0.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  3  7  1  0  0  0  0
M  END
$$$$
"""


@pytest.fixture
def pdb_file(temp_dir, sample_pdb_content):
    """Create a temporary PDB file for testing."""
    pdb_file = temp_dir / "test_protein.pdb"
    pdb_file.write_text(sample_pdb_content)
    return str(pdb_file)


@pytest.fixture
def sdf_file(temp_dir, sample_sdf_content):
    """Create a temporary SDF file for testing."""
    sdf_file = temp_dir / "test_ligand.sdf"
    sdf_file.write_text(sample_sdf_content)
    return str(sdf_file)


@pytest.fixture
def scoring_functions(test_config):
    """Create a ScoringFunctions instance for testing."""
    if ScoringFunctions is None:
        pytest.skip("PandaDock not available")
    
    return ScoringFunctions(test_config)


@pytest.fixture
def ga_engine(test_config):
    """Create a GAEngine instance for testing."""
    if GAEngine is None:
        pytest.skip("PandaDock not available")
    
    return GAEngine(test_config)


@pytest.fixture
def ml_engine(test_config):
    """Create an MLEngine instance for testing."""
    if MLEngine is None:
        pytest.skip("PandaDock not available")
    
    return MLEngine(test_config)


@pytest.fixture
def physics_engine(test_config):
    """Create a PhysicsEngine instance for testing."""
    if PhysicsEngine is None:
        pytest.skip("PandaDock not available")
    
    return PhysicsEngine(test_config)


# Simple mock fixtures for testing without dependencies
@pytest.fixture
def mock_torch_model():
    """Mock PyTorch model for ML testing."""
    class MockModel:
        def eval(self):
            return self
        def forward(self, *args, **kwargs):
            return None
    return MockModel()


@pytest.fixture
def mock_rdkit_mol():
    """Mock RDKit molecule for chemistry testing."""
    class MockMol:
        def GetNumAtoms(self):
            return 7
        def GetConformers(self):
            return []
    return MockMol()


# Performance testing fixtures
@pytest.fixture
def benchmark_data():
    """Data for benchmark testing."""
    return {
        "small_ligand": np.random.randn(10, 3),
        "medium_ligand": np.random.randn(50, 3),
        "large_ligand": np.random.randn(100, 3),
        "protein": np.random.randn(500, 3)
    }


# Skip conditions for optional dependencies
def requires_torch():
    """Skip test if PyTorch is not available."""
    try:
        import torch
        return False
    except ImportError:
        return True


def requires_rdkit():
    """Skip test if RDKit is not available."""
    try:
        import rdkit
        return False
    except ImportError:
        return True


def requires_gpu():
    """Skip test if GPU is not available."""
    try:
        import torch
        return not torch.cuda.is_available()
    except ImportError:
        return True


# Pytest markers for conditional skipping
skip_if_no_torch = pytest.mark.skipif(requires_torch(), reason="PyTorch not available")
skip_if_no_rdkit = pytest.mark.skipif(requires_rdkit(), reason="RDKit not available")
skip_if_no_gpu = pytest.mark.skipif(requires_gpu(), reason="GPU not available")


# Utility functions for testing
def assert_pose_valid(pose):
    """Assert that a pose object is valid."""
    assert pose is not None
    assert pose.coordinates is not None
    assert pose.coordinates.shape[1] == 3
    assert len(pose.coordinates) > 0
    assert isinstance(pose.score, (int, float))
    assert isinstance(pose.energy, (int, float))
    assert isinstance(pose.ligand_name, str)
    assert isinstance(pose.pose_id, str)


def assert_coordinates_realistic(coords):
    """Assert that coordinates are realistic."""
    assert coords is not None
    assert coords.shape[1] == 3
    assert len(coords) > 0
    assert np.all(np.isfinite(coords))
    assert np.max(np.abs(coords)) < 1000  # Reasonable coordinate range


def assert_energy_realistic(energy):
    """Assert that energy values are realistic."""
    assert isinstance(energy, (int, float))
    assert np.isfinite(energy)
    assert -50.0 < energy < 50.0  # Reasonable energy range for molecular systems


# Custom assertion helpers
class PandaDockAssertions:
    """Custom assertions for PandaDock testing."""
    
    @staticmethod
    def assert_algorithm_consistency(poses1, poses2, tolerance=1e-6):
        """Assert that two sets of poses are consistent."""
        assert len(poses1) == len(poses2)
        
        for p1, p2 in zip(poses1, poses2):
            assert np.allclose(p1.coordinates, p2.coordinates, atol=tolerance)
            assert abs(p1.score - p2.score) < tolerance
            assert abs(p1.energy - p2.energy) < tolerance
    
    @staticmethod
    def assert_scoring_function_consistency(scores1, scores2, tolerance=0.1):
        """Assert that scoring functions produce consistent results."""
        assert len(scores1) == len(scores2)
        
        for s1, s2 in zip(scores1, scores2):
            assert abs(s1 - s2) < tolerance


@pytest.fixture
def pandadock_assertions():
    """Provide custom PandaDock assertions."""
    return PandaDockAssertions


# Test data validation
def validate_test_environment():
    """Validate that the test environment is properly set up."""
    required_modules = ["numpy", "scipy", "pandas", "matplotlib"]
    missing_modules = []
    
    for module in required_modules:
        try:
            __import__(module)
        except ImportError:
            missing_modules.append(module)
    
    if missing_modules:
        pytest.skip(f"Missing required modules: {', '.join(missing_modules)}")


# Session-level setup
@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up the test environment."""
    validate_test_environment()
    
    # Set random seeds for reproducible tests
    np.random.seed(42)
    
    # Configure logging for tests
    import logging
    logging.basicConfig(level=logging.WARNING)
    
    # Suppress warnings during testing
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=PendingDeprecationWarning)


# Cleanup fixture
@pytest.fixture(autouse=True)
def cleanup_after_test():
    """Clean up after each test."""
    yield
    # Any necessary cleanup code here
    pass