"""
Basic tests for PandaDock modular architecture.
"""
import unittest
import pytest


class TestBasicImports(unittest.TestCase):
    """Test that the new modular architecture imports work correctly."""

    def test_core_imports(self):
        """Test core module imports."""
        import pandadock
        from pandadock.core import DockingEngine, PoseGenerator, ResultProcessor
        self.assertTrue(hasattr(pandadock, '__version__'))

    def test_algorithm_imports(self):
        """Test algorithm module imports."""
        from pandadock.algorithms import AlgorithmFactory
        from pandadock.algorithms.genetic_algorithm_enhanced import GeneticAlgorithm
        from pandadock.algorithms.random_search import RandomSearchAlgorithm
        self.assertIsNotNone(AlgorithmFactory)

    def test_scoring_imports(self):
        """Test scoring module imports."""
        from pandadock.scoring import ScoringFunctionFactory
        from pandadock.scoring.base_scoring import BaseScoringFunction
        self.assertIsNotNone(ScoringFunctionFactory)

    def test_molecule_imports(self):
        """Test molecule handling imports."""
        from pandadock.molecules.ligand_handler import LigandHandler, LigandStructure
        from pandadock.molecules.protein_handler import ProteinHandler
        self.assertIsNotNone(LigandHandler)

    def test_analysis_imports(self):
        """Test analysis module imports."""
        from pandadock.analysis.binding_affinity import BindingAffinityCalculator
        from pandadock.analysis.pose_clusterer import PoseClusterer
        self.assertIsNotNone(BindingAffinityCalculator)

    def test_hardware_imports(self):
        """Test hardware abstraction imports."""
        from pandadock.hardware import DeviceManager, ComputeBackend
        self.assertIsNotNone(DeviceManager)

    def test_cli_imports(self):
        """Test CLI module imports."""
        from pandadock.cli import create_argument_parser
        self.assertIsNotNone(create_argument_parser)


class TestVersion(unittest.TestCase):
    """Test version information."""

    def test_version_exists(self):
        """Test that version is defined."""
        import pandadock
        self.assertTrue(hasattr(pandadock, '__version__'))
        self.assertIsInstance(pandadock.__version__, str)

    def test_version_format(self):
        """Test that version follows semantic versioning."""
        import pandadock
        version_parts = pandadock.__version__.split('.')
        self.assertGreaterEqual(len(version_parts), 2)


if __name__ == '__main__':
    unittest.main()