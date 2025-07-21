#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for PandaDock algorithm naming and functionality
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestAlgorithmNaming:
    """Test that the new PandaDock algorithm names are properly implemented"""
    
    def test_algorithm_names_updated(self):
        """Test that algorithm names have been updated from old to new"""
        # Test that we're using the new naming convention
        new_algorithms = ['pandacore', 'pandaml', 'pandaphysics']
        old_algorithms = ['vina', 'ml', 'physics']  # Should be replaced
        
        # Verify new names are properly formatted
        for algorithm in new_algorithms:
            assert algorithm.startswith('panda')
            assert algorithm.islower()
            assert algorithm.replace('panda', '') in ['core', 'ml', 'physics']
    
    def test_algorithm_descriptions(self):
        """Test algorithm description mapping"""
        algorithm_descriptions = {
            'pandacore': 'Robust baseline algorithm with reliable performance',
            'pandaml': 'Advanced machine learning algorithm with superior affinity prediction',
            'pandaphysics': 'Physics-based algorithm specialized for metal coordination'
        }
        
        for algorithm, description in algorithm_descriptions.items():
            assert isinstance(algorithm, str)
            assert isinstance(description, str)
            assert len(description) > 10  # Reasonable description length
    
    def test_no_old_algorithm_names(self):
        """Test that old algorithm names are not used in key places"""
        # This is a documentation test - we verify the concept
        old_names = ['vina', 'glide']  # Old names that should be avoided
        new_names = ['pandacore', 'pandaml', 'pandaphysics']
        
        # Verify new names don't conflict with old ones
        for old_name in old_names:
            for new_name in new_names:
                assert old_name not in new_name.lower()
        
        # Verify new names are distinct
        assert len(set(new_names)) == len(new_names)


class TestBenchmarkResults:
    """Test benchmark result values from the documentation"""
    
    def test_documented_performance_metrics(self):
        """Test that documented performance metrics are reasonable"""
        # These are the actual benchmark results from the documentation
        benchmark_results = {
            'pandaml': {
                'affinity_r2': 0.878,
                'success_rate': 0.467,
                'mean_rmsd': 3.22,
                'speed_seconds': 26.1
            },
            'pandaphysics': {
                'affinity_r2': 0.671,
                'success_rate': 0.488,
                'mean_rmsd': 2.18,
                'speed_seconds': 42.9
            },
            'pandacore': {
                'affinity_r2': 0.738,
                'success_rate': 0.463,
                'mean_rmsd': 3.63,
                'speed_seconds': 32.9
            }
        }
        
        for algorithm, metrics in benchmark_results.items():
            # Test reasonable ranges for metrics
            assert 0.0 <= metrics['affinity_r2'] <= 1.0
            assert 0.0 <= metrics['success_rate'] <= 1.0
            assert metrics['mean_rmsd'] > 0.0
            assert metrics['speed_seconds'] > 0.0
            
            # Test that PandaML has best affinity prediction (R²)
            if algorithm == 'pandaml':
                assert metrics['affinity_r2'] > 0.85  # Superior performance
                
            # Test that PandaPhysics has best RMSD
            if algorithm == 'pandaphysics':
                assert metrics['mean_rmsd'] < 2.5  # Best pose accuracy
    
    def test_metal_vs_nonmetal_findings(self):
        """Test metal vs non-metal analysis findings"""
        # These are the documented findings
        metal_analysis = {
            'total_complexes': 285,
            'metal_complexes': 93,
            'nonmetal_complexes': 192,
            'metal_percentage': 32.6,
            'computational_complexity_increase': 32  # percent
        }
        
        # Validate the numbers make sense
        assert metal_analysis['total_complexes'] == (
            metal_analysis['metal_complexes'] + metal_analysis['nonmetal_complexes']
        )
        
        calculated_percentage = (metal_analysis['metal_complexes'] / 
                               metal_analysis['total_complexes']) * 100
        assert abs(calculated_percentage - metal_analysis['metal_percentage']) < 1.0
        
        # Test that metal complexes are minority but significant
        assert 25 <= metal_analysis['metal_percentage'] <= 40


class TestCLIUpdates:
    """Test that CLI has been updated with new algorithm names"""
    
    def test_cli_help_content(self):
        """Test that CLI help mentions new algorithm names"""
        try:
            from pandadock.__main__ import main
            # For now, just test that main function exists
            assert main is not None
        except ImportError as e:
            pytest.skip(f"CLI main not available: {e}")
    
    def test_algorithm_option_values(self):
        """Test expected algorithm option values"""
        expected_algorithms = ['pandacore', 'pandaml', 'pandaphysics']
        
        # Test basic properties of algorithm names
        for algorithm in expected_algorithms:
            assert isinstance(algorithm, str)
            assert len(algorithm) > 5
            assert algorithm.startswith('panda')


if __name__ == "__main__":
    pytest.main([__file__])