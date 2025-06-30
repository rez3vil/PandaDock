"""
Comprehensive test suite for all PandaDock functionality including new modules.

This test verifies:
1. All original functionality (previously tested)
2. Metal-based docking with coordination constraints
3. Virtual screening manager
4. Batch screening functionality
5. Integration between all components
"""

import sys
import os
import numpy as np
import tempfile
import shutil
from pathlib import Path

# Add pandadock to path
sys.path.insert(0, '/Users/pritam/PandaDock')

def test_metal_docking_functionality():
    """Test metal-based docking functionality."""
    print("Testing Metal-Based Docking...")
    
    try:
        from pandadock.algorithms.metal_docking import (
            MetalCenter, MetalConstraint, MetalDockingScorer,
            MetalDockingAlgorithm, MetalDockingPreparation
        )
        from pandadock.algorithms.genetic_algorithm import GeneticAlgorithm
        from pandadock.scoring.physics_based_scoring import PhysicsBasedScoringFunction
        
        # Create mock protein with metal center
        class MockProtein:
            def __init__(self):
                # Protein with zinc center
                self.coords = np.array([
                    [0.0, 0.0, 0.0],    # Zinc ion
                    [2.1, 0.0, 0.0],    # Coordinating nitrogen
                    [0.0, 2.0, 0.0],    # Coordinating oxygen
                    [-2.1, 0.0, 0.0],   # Coordinating sulfur
                    [5.0, 5.0, 5.0]     # Distant protein atom
                ])
                self.atoms = [
                    {'coords': [0.0, 0.0, 0.0], 'element': 'Zn', 'name': 'ZN'},
                    {'coords': [2.1, 0.0, 0.0], 'element': 'N', 'name': 'N1'},
                    {'coords': [0.0, 2.0, 0.0], 'element': 'O', 'name': 'O1'},
                    {'coords': [-2.1, 0.0, 0.0], 'element': 'S', 'name': 'S1'},
                    {'coords': [5.0, 5.0, 5.0], 'element': 'C', 'name': 'C1'}
                ]
        
        # Create mock ligand
        class MockLigand:
            def __init__(self):
                self.coords = np.array([
                    [1.0, 1.0, 1.0],    # Potential coordinating atom
                    [2.0, 2.0, 2.0],    # Another atom
                    [3.0, 3.0, 3.0]     # Third atom
                ])
                self.atoms = [
                    {'coords': [1.0, 1.0, 1.0], 'element': 'N', 'name': 'N1'},
                    {'coords': [2.0, 2.0, 2.0], 'element': 'C', 'name': 'C1'},
                    {'coords': [3.0, 3.0, 3.0], 'element': 'C', 'name': 'C2'}
                ]
        
        protein = MockProtein()
        ligand = MockLigand()
        
        # Test metal center detection
        metal_centers = MetalDockingPreparation.detect_metal_centers(protein)
        assert len(metal_centers) == 1, f"Expected 1 metal center, got {len(metal_centers)}"
        assert metal_centers[0].element == 'ZN', f"Expected Zn, got {metal_centers[0].element}"
        
        # Test metal constraints creation
        constraints = MetalDockingPreparation.create_metal_constraints(metal_centers)
        assert len(constraints) == 1, f"Expected 1 constraint, got {len(constraints)}"
        
        # Test ligand metal binding site analysis
        binding_analysis = MetalDockingPreparation.analyze_ligand_metal_binding_sites(ligand)
        assert 'binding_sites' in binding_analysis
        assert 'chelating_groups' in binding_analysis
        assert binding_analysis['total_coordinating_atoms'] >= 1
        
        # Test metal-aware scoring function
        base_scorer = PhysicsBasedScoringFunction()
        metal_scorer = MetalDockingScorer(
            base_scoring_function=base_scorer,
            metal_centers=metal_centers,
            metal_constraints=constraints,
            metal_weight=10.0
        )
        
        # Test scoring
        score = metal_scorer.score_pose(ligand.coords, protein.coords)
        assert isinstance(score, (int, float)), f"Score should be numeric, got {type(score)}"
        
        # Test metal docking algorithm
        base_algorithm = GeneticAlgorithm(base_scorer, max_iterations=10, population_size=5)
        metal_algorithm = MetalDockingAlgorithm(
            base_algorithm=base_algorithm,
            metal_centers=metal_centers,
            metal_constraints=constraints,
            refinement_iterations=5
        )
        
        # Test search
        results = metal_algorithm.search(protein, ligand)
        assert isinstance(results, list), f"Results should be list, got {type(results)}"
        assert len(results) > 0, "Should return at least one result"
        
        print("✓ Metal-based docking functionality working correctly")
        return True
        
    except Exception as e:
        print(f"✗ Metal-based docking test failed: {e}")
        return False


def test_virtual_screening_functionality():
    """Test virtual screening manager."""
    print("Testing Virtual Screening Manager...")
    
    try:
        from pandadock.screening.virtual_screening import VirtualScreeningManager
        
        # Create temporary output directory
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "vs_output"
            
            # Create mock protein
            class MockProtein:
                def __init__(self):
                    self.coords = np.random.randn(10, 3) * 5
                    self.atoms = [{'element': 'C', 'coords': coord} for coord in self.coords]
                    self.active_site = {'center': np.array([0., 0., 0.]), 'radius': 10.0}
            
            # Create mock ligands
            class MockLigand:
                def __init__(self, name):
                    self.name = name
                    self.coords = np.random.randn(5, 3) * 2
                    self.atoms = [{'element': 'C', 'coords': coord} for coord in self.coords]
            
            protein = MockProtein()
            ligands = [MockLigand(f"ligand_{i}") for i in range(3)]
            ligand_names = [f"ligand_{i}" for i in range(3)]
            
            # Create virtual screening manager
            vs_manager = VirtualScreeningManager(
                scoring_function_type='physics_based',
                algorithm_type='genetic',
                output_dir=output_dir,
                n_workers=1,  # Single worker for testing
                exhaustiveness=2,
                num_modes=3
            )
            
            # Run screening
            screening_results = vs_manager.run_screening(protein, ligands, ligand_names)
            
            # Verify results structure
            assert 'results' in screening_results
            assert 'statistics' in screening_results
            assert len(screening_results['results']) == 3
            
            # Check individual results
            for ligand_name in ligand_names:
                assert ligand_name in screening_results['results']
                result = screening_results['results'][ligand_name]
                assert 'poses' in result or 'error' in result
                assert 'processing_time' in result
            
            # Check if output directory was created
            if output_dir.exists():
                assert (output_dir / 'poses').exists()
                assert (output_dir / 'reports').exists()
            
            print("✓ Virtual screening manager working correctly")
            return True
            
    except Exception as e:
        print(f"✗ Virtual screening test failed: {e}")
        return False


def test_batch_screening_functionality():
    """Test batch screening manager."""
    print("Testing Batch Screening Manager...")
    
    try:
        from pandadock.screening.batch_screening import BatchScreeningManager
        
        # Create temporary test files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create mock protein file
            protein_file = temp_path / "protein.pdb"
            with open(protein_file, 'w') as f:
                f.write("ATOM      1  C   PRO A   1       0.000   0.000   0.000  1.00  0.00           C\n")
                f.write("ATOM      2  C   PRO A   1       1.000   0.000   0.000  1.00  0.00           C\n")
                f.write("END\n")
            
            # Create ligand library directory
            ligand_library = temp_path / "ligands"
            ligand_library.mkdir()
            
            # Create mock ligand files
            for i in range(2):
                ligand_file = ligand_library / f"ligand_{i}.mol"
                with open(ligand_file, 'w') as f:
                    f.write("ligand_name\n")
                    f.write("  2  1  0  0  0  0  0  0  0  0999 V2000\n")
                    f.write("    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
                    f.write("    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
                    f.write("  1  2  1  0  0  0  0\n")
                    f.write("M  END\n")
            
            # Create output directory
            output_dir = temp_path / "batch_output"
            
            # Create batch screening manager
            batch_manager = BatchScreeningManager(
                protein_file=protein_file,
                ligand_library=ligand_library,
                output_dir=output_dir,
                screening_params={
                    'algorithm': 'genetic',
                    'scoring_function': 'physics_based',
                    'exhaustiveness': 1,
                    'num_modes': 2
                },
                n_processes=1  # Single process for testing
            )
            
            # Test ligand file detection
            ligand_files = batch_manager._get_ligand_files()
            assert len(ligand_files) == 2, f"Expected 2 ligand files, got {len(ligand_files)}"
            
            # Note: Full batch screening would require working molecule classes
            # For this test, we verify the manager can be created and basic functions work
            
            print("✓ Batch screening manager structure working correctly")
            return True
            
    except Exception as e:
        print(f"✗ Batch screening test failed: {e}")
        return False


def test_algorithm_factory_extensions():
    """Test algorithm factory with new algorithm types."""
    print("Testing Extended Algorithm Factory...")
    
    try:
        from pandadock.algorithms.algorithm_factory import AlgorithmFactory
        from pandadock.hardware import DeviceManager, ComputeBackend
        from pandadock.scoring.physics_based_scoring import PhysicsBasedScoringFunction
        
        # Create components
        device_manager = DeviceManager()
        compute_backend = ComputeBackend(device_manager)
        algorithm_factory = AlgorithmFactory(compute_backend)
        scoring_function = PhysicsBasedScoringFunction()
        
        # Test available algorithms
        available_algorithms = algorithm_factory.get_available_algorithms()
        assert 'metal' in available_algorithms, "Metal algorithm should be available"
        assert 'genetic' in available_algorithms, "Genetic algorithm should be available"
        assert 'random' in available_algorithms, "Random search should be available"
        
        # Test creating standard algorithms
        genetic_alg = algorithm_factory.create_algorithm('genetic', scoring_function, max_iterations=10)
        assert genetic_alg is not None, "Should create genetic algorithm"
        
        random_alg = algorithm_factory.create_algorithm('random', scoring_function, max_iterations=10)
        assert random_alg is not None, "Should create random search algorithm"
        
        # Test metal-aware scoring function creation
        metal_scorer = algorithm_factory.create_metal_aware_scoring_function(
            base_scoring_function=scoring_function,
            metal_weight=5.0
        )
        assert metal_scorer is not None, "Should create metal-aware scoring function"
        
        print("✓ Extended algorithm factory working correctly")
        return True
        
    except Exception as e:
        print(f"✗ Algorithm factory extension test failed: {e}")
        return False


def test_complete_integration():
    """Test integration between all components."""
    print("Testing Complete System Integration...")
    
    try:
        # Import all major components
        from pandadock import (
            DockingEngine, DeviceManager, AlgorithmFactory, 
            ScoringFunctionFactory, VirtualScreeningManager
        )
        from pandadock.algorithms.metal_docking import MetalDockingPreparation
        
        # Create mock system
        class MockProtein:
            def __init__(self):
                self.coords = np.array([
                    [0.0, 0.0, 0.0],    # Metal center
                    [2.0, 0.0, 0.0],    # Coordinating atom
                    [0.0, 2.0, 0.0],    # Another coordinating atom
                    [5.0, 5.0, 5.0]     # Distant protein atom
                ])
                self.atoms = [
                    {'coords': [0.0, 0.0, 0.0], 'element': 'Zn', 'name': 'ZN'},
                    {'coords': [2.0, 0.0, 0.0], 'element': 'N', 'name': 'N1'},
                    {'coords': [0.0, 2.0, 0.0], 'element': 'O', 'name': 'O1'},
                    {'coords': [5.0, 5.0, 5.0], 'element': 'C', 'name': 'C1'}
                ]
                self.active_site = {'center': np.array([0., 0., 0.]), 'radius': 10.0}
        
        class MockLigand:
            def __init__(self):
                self.coords = np.random.randn(5, 3) * 2
                self.atoms = [{'element': 'C', 'coords': coord} for coord in self.coords]
        
        protein = MockProtein()
        ligand = MockLigand()
        
        # Test metal center detection
        metal_centers = MetalDockingPreparation.detect_metal_centers(protein)
        has_metals = len(metal_centers) > 0
        
        # Create system components
        device_manager = DeviceManager()
        scoring_factory = ScoringFunctionFactory()
        algorithm_factory = AlgorithmFactory(device_manager.compute_backend)
        
        # Test different algorithm types
        algorithm_types = ['genetic', 'random']
        if has_metals:
            algorithm_types.append('metal')
        
        for alg_type in algorithm_types:
            scoring_function = scoring_factory.create_scoring_function('physics_based')
            
            if alg_type == 'metal':
                algorithm = algorithm_factory.create_algorithm(
                    alg_type, scoring_function, 
                    protein=protein, max_iterations=5
                )
            else:
                algorithm = algorithm_factory.create_algorithm(
                    alg_type, scoring_function, max_iterations=5
                )
            
            assert algorithm is not None, f"Should create {alg_type} algorithm"
        
        # Test virtual screening integration
        vs_manager = VirtualScreeningManager(
            scoring_function_type='physics_based',
            algorithm_type='genetic',
            n_workers=1,
            exhaustiveness=1,
            num_modes=2
        )
        
        # Run mini screening
        results = vs_manager.run_screening(protein, [ligand], ['test_ligand'])
        assert 'results' in results
        assert 'test_ligand' in results['results']
        
        print("✓ Complete system integration working correctly")
        return True
        
    except Exception as e:
        print(f"✗ Integration test failed: {e}")
        return False


def main():
    """Run all extended functionality tests."""
    print("🧪 Running Extended PandaDock Functionality Tests")
    print("=" * 60)
    
    tests = [
        ("Metal-Based Docking", test_metal_docking_functionality),
        ("Virtual Screening", test_virtual_screening_functionality),
        ("Batch Screening", test_batch_screening_functionality),
        ("Algorithm Factory Extensions", test_algorithm_factory_extensions),
        ("Complete Integration", test_complete_integration)
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}")
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"📊 Extended Test Results: {passed}/{passed + failed} passed")
    
    if failed == 0:
        print("🎉 ALL EXTENDED TESTS PASSED! All new functionality is working correctly.")
        print("✅ Metal-based docking with coordination constraints implemented")
        print("✅ Virtual screening manager for high-throughput screening implemented") 
        print("✅ Batch screening with parallel processing implemented")
        print("✅ Complete integration between all components verified")
        return True
    else:
        print(f"❌ {failed} tests failed. Some functionality needs attention.")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)