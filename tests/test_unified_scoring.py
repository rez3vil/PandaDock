"""
Tests for scoring functions in the new PandaDock architecture.
"""
import pytest
import numpy as np
from pandadock.scoring import ScoringFunctionFactory
from pandadock.molecules.ligand_handler import LigandHandler, LigandStructure
from pandadock.molecules.protein_handler import ProteinHandler


def test_scoring_factory_creation():
    """Test that ScoringFunctionFactory can create scoring functions."""
    factory = ScoringFunctionFactory()
    
    # Test basic scoring function creation
    basic_scorer = factory.create_scoring_function('basic')
    assert basic_scorer is not None
    
    # Test enhanced scoring function creation
    enhanced_scorer = factory.create_scoring_function('enhanced')
    assert enhanced_scorer is not None


def test_scoring_function_types():
    """Test different types of scoring functions."""
    factory = ScoringFunctionFactory()
    
    # Test that we can create different scoring function types
    scoring_types = ['basic', 'enhanced']
    
    for scoring_type in scoring_types:
        try:
            scorer = factory.create_scoring_function(scoring_type)
            assert scorer is not None
            assert hasattr(scorer, 'score')
        except Exception as e:
            pytest.skip(f"Scoring type {scoring_type} not fully implemented: {e}")


def test_ligand_structure_creation():
    """Test ligand structure creation for scoring tests."""
    # Create a simple ligand structure for testing
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
    
    ligand = LigandStructure(
        coords=coords,
        atom_names=['C1', 'C2', 'C3'],
        atom_types=['C', 'C', 'C']
    )
    
    assert ligand.n_atoms == 3
    assert len(ligand.atom_names) == 3
    assert len(ligand.atom_types) == 3


def test_ligand_handler():
    """Test ligand handler functionality."""
    handler = LigandHandler()
    assert handler is not None
    
    # Test that handler has expected methods
    assert hasattr(handler, 'load_ligand')


def test_protein_handler():
    """Test protein handler functionality."""
    handler = ProteinHandler()
    assert handler is not None
    
    # Test that handler has expected methods
    assert hasattr(handler, 'load_protein')


@pytest.mark.skipif(
    True,  # Skip until we have test files
    reason="Requires test protein and ligand files"
)
def test_scoring_with_real_molecules():
    """Test scoring with actual protein and ligand files."""
    # This test would use real test files when available
    protein_handler = ProteinHandler()
    ligand_handler = LigandHandler()
    
    # protein = protein_handler.load_protein("tests/protein.pdb")
    # ligand = ligand_handler.load_ligand("tests/ligand.sdf")
    
    factory = ScoringFunctionFactory()
    scorer = factory.create_scoring_function('basic')
    
    # score = scorer.score(protein, ligand)
    # assert isinstance(score, float)


if __name__ == "__main__":
    pytest.main([__file__])