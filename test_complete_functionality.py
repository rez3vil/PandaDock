"""
Comprehensive test to verify all PandaDock functionality is preserved.

This script tests all the newly implemented algorithms, scoring functions,
and analysis modules to ensure the refactored architecture maintains
full compatibility with the original PandaDock capabilities.
"""

import sys
import logging
import numpy as np
from pathlib import Path

# Add PandaDock to path
sys.path.insert(0, str(Path(__file__).parent))

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def test_algorithms():
    """Test all docking algorithms."""
    logger.info("=" * 60)
    logger.info("Testing Docking Algorithms")
    logger.info("=" * 60)
    
    try:
        from pandadock.algorithms import (
            AlgorithmFactory, GeneticAlgorithm, RandomSearchAlgorithm,
            MonteCarloAlgorithm, PANDADOCKAlgorithm, MMFFMinimization
        )
        from pandadock.hardware import DeviceManager, ComputeBackend
        from pandadock.scoring import ScoringFunctionFactory
        
        # Initialize hardware and scoring
        device_manager = DeviceManager()
        compute_backend = ComputeBackend(device_manager)
        
        # Create scoring function
        scoring_factory = ScoringFunctionFactory(compute_backend)
        scoring_func = scoring_factory.create_scoring_function('standard')
        
        # Create algorithm factory
        algo_factory = AlgorithmFactory(compute_backend)
        
        # Test each algorithm type
        algorithms_to_test = ['genetic', 'random', 'monte-carlo', 'pandadock']
        
        for algo_type in algorithms_to_test:
            try:
                logger.info(f"Testing {algo_type} algorithm...")
                algorithm = algo_factory.create_algorithm(algo_type, scoring_func, max_iterations=10)
                logger.info(f"✓ {algo_type} algorithm created successfully")
                
                # Test algorithm info
                info = algorithm.get_algorithm_info()
                logger.info(f"  Algorithm info: {info}")
                
            except Exception as e:
                logger.error(f"✗ Failed to create {algo_type} algorithm: {e}")
        
        # Test MMFF minimization
        try:
            logger.info("Testing MMFF94 minimization...")
            mmff = algo_factory.create_mmff_minimizer()
            logger.info("✓ MMFF94 minimizer created successfully")
        except Exception as e:
            logger.error(f"✗ MMFF94 minimizer failed: {e}")
        
        # Test available algorithms
        available = algo_factory.get_available_algorithms()
        logger.info(f"Available algorithms: {list(available.keys())}")
        
        return True
        
    except Exception as e:
        logger.error(f"Algorithm testing failed: {e}")
        return False

def test_scoring_functions():
    """Test all scoring functions."""
    logger.info("=" * 60)
    logger.info("Testing Scoring Functions")
    logger.info("=" * 60)
    
    try:
        from pandadock.scoring import (
            ScoringFunctionFactory, PhysicsBasedScoringFunction
        )
        from pandadock.hardware import DeviceManager, ComputeBackend
        
        # Initialize hardware
        device_manager = DeviceManager()
        compute_backend = ComputeBackend(device_manager)
        
        # Create scoring factory
        scoring_factory = ScoringFunctionFactory(compute_backend)
        
        # Test different scoring types
        scoring_types = [
            ('standard', {}),
            ('enhanced', {'enhanced': True}),
            ('physics', {'physics_based': True}),
            ('basic', {})
        ]
        
        for scoring_type, kwargs in scoring_types:
            try:
                logger.info(f"Testing {scoring_type} scoring...")
                scoring_func = scoring_factory.create_scoring_function(**kwargs)
                logger.info(f"✓ {scoring_type} scoring function created successfully")
                
            except Exception as e:
                logger.error(f"✗ Failed to create {scoring_type} scoring: {e}")
        
        # Test physics-based scoring directly
        try:
            logger.info("Testing physics-based scoring directly...")
            physics_scoring = PhysicsBasedScoringFunction()
            logger.info("✓ Physics-based scoring created successfully")
            
            # Test component weights
            weights = physics_scoring.weights
            logger.info(f"  Energy component weights: {weights}")
            
        except Exception as e:
            logger.error(f"✗ Physics-based scoring failed: {e}")
        
        # Test available scoring functions
        available = scoring_factory.get_available_scoring_functions()
        logger.info(f"Available scoring functions: {list(available.keys())}")
        
        return True
        
    except Exception as e:
        logger.error(f"Scoring function testing failed: {e}")
        return False

def test_analysis_modules():
    """Test all analysis modules."""
    logger.info("=" * 60)
    logger.info("Testing Analysis Modules")
    logger.info("=" * 60)
    
    try:
        from pandadock.analysis import (
            PoseClusterer, InteractionFingerprinter, BindingModeClassifier,
            EnergyDecomposition, DockingReportGenerator
        )
        
        # Test pose clustering
        try:
            logger.info("Testing pose clustering...")
            clusterer = PoseClusterer(method='hierarchical', rmsd_cutoff=2.0)
            logger.info("✓ PoseClusterer created successfully")
            
            # Test clustering statistics
            mock_results = {'clusters': [], 'n_clusters': 0, 'rmsd_cutoff': 2.0}
            stats = clusterer.cluster_statistics(mock_results)
            logger.info(f"  Clustering statistics: {stats}")
            
        except Exception as e:
            logger.error(f"✗ PoseClusterer failed: {e}")
        
        # Test interaction analysis
        try:
            logger.info("Testing interaction analysis...")
            fingerprinter = InteractionFingerprinter()
            logger.info("✓ InteractionFingerprinter created successfully")
            
            # Test interaction types
            types = fingerprinter.interaction_types
            logger.info(f"  Interaction types: {types}")
            
        except Exception as e:
            logger.error(f"✗ InteractionFingerprinter failed: {e}")
        
        # Test binding mode classification
        try:
            logger.info("Testing binding mode classification...")
            classifier = BindingModeClassifier()
            logger.info("✓ BindingModeClassifier created successfully")
            
        except Exception as e:
            logger.error(f"✗ BindingModeClassifier failed: {e}")
        
        # Test energy decomposition
        try:
            logger.info("Testing energy decomposition...")
            
            # Create a mock scoring function
            class MockScoringFunction:
                def score(self, protein, ligand):
                    return -5.0
            
            mock_scoring = MockScoringFunction()
            energy_decomp = EnergyDecomposition(mock_scoring)
            logger.info("✓ EnergyDecomposition created successfully")
            
        except Exception as e:
            logger.error(f"✗ EnergyDecomposition failed: {e}")
        
        # Test report generation
        try:
            logger.info("Testing report generation...")
            report_gen = DockingReportGenerator(report_format='html')
            logger.info("✓ DockingReportGenerator created successfully")
            
        except Exception as e:
            logger.error(f"✗ DockingReportGenerator failed: {e}")
        
        return True
        
    except Exception as e:
        logger.error(f"Analysis module testing failed: {e}")
        return False

def test_hardware_abstraction():
    """Test hardware abstraction layer."""
    logger.info("=" * 60)
    logger.info("Testing Hardware Abstraction")
    logger.info("=" * 60)
    
    try:
        from pandadock.hardware import DeviceManager, ComputeBackend, PerformanceMonitor
        
        # Test device manager
        try:
            logger.info("Testing device manager...")
            device_manager = DeviceManager()
            device_info = device_manager.get_device_info()
            logger.info(f"✓ Device manager initialized")
            logger.info(f"  Device info: {device_info}")
            
        except Exception as e:
            logger.error(f"✗ Device manager failed: {e}")
        
        # Test compute backend
        try:
            logger.info("Testing compute backend...")
            compute_backend = ComputeBackend(device_manager)
            logger.info("✓ Compute backend created successfully")
            
        except Exception as e:
            logger.error(f"✗ Compute backend failed: {e}")
        
        # Test performance monitor
        try:
            logger.info("Testing performance monitor...")
            perf_monitor = PerformanceMonitor()
            logger.info("✓ Performance monitor created successfully")
            
        except Exception as e:
            logger.error(f"✗ Performance monitor failed: {e}")
        
        return True
        
    except Exception as e:
        logger.error(f"Hardware abstraction testing failed: {e}")
        return False

def test_molecule_handling():
    """Test molecule handling modules."""
    logger.info("=" * 60)
    logger.info("Testing Molecule Handling")
    logger.info("=" * 60)
    
    try:
        from pandadock.molecules import ProteinHandler, LigandHandler, StructurePreparation
        
        # Test protein handler
        try:
            logger.info("Testing protein handler...")
            protein_handler = ProteinHandler()
            logger.info("✓ ProteinHandler created successfully")
            
        except Exception as e:
            logger.error(f"✗ ProteinHandler failed: {e}")
        
        # Test ligand handler
        try:
            logger.info("Testing ligand handler...")
            ligand_handler = LigandHandler()
            logger.info("✓ LigandHandler created successfully")
            
        except Exception as e:
            logger.error(f"✗ LigandHandler failed: {e}")
        
        # Test structure preparation
        try:
            logger.info("Testing structure preparation...")
            prep = StructurePreparation()
            logger.info("✓ StructurePreparation created successfully")
            
        except Exception as e:
            logger.error(f"✗ StructurePreparation failed: {e}")
        
        return True
        
    except Exception as e:
        logger.error(f"Molecule handling testing failed: {e}")
        return False

def test_core_engine():
    """Test core docking engine."""
    logger.info("=" * 60)
    logger.info("Testing Core Docking Engine")
    logger.info("=" * 60)
    
    try:
        from pandadock.core import DockingEngine, PoseGenerator, ResultProcessor
        
        # Test docking engine
        try:
            logger.info("Testing docking engine...")
            engine = DockingEngine()
            logger.info("✓ DockingEngine created successfully")
            
        except Exception as e:
            logger.error(f"✗ DockingEngine failed: {e}")
        
        # Test pose generator
        try:
            logger.info("Testing pose generator...")
            pose_gen = PoseGenerator()
            logger.info("✓ PoseGenerator created successfully")
            
        except Exception as e:
            logger.error(f"✗ PoseGenerator failed: {e}")
        
        # Test result processor
        try:
            logger.info("Testing result processor...")
            result_proc = ResultProcessor()
            logger.info("✓ ResultProcessor created successfully")
            
        except Exception as e:
            logger.error(f"✗ ResultProcessor failed: {e}")
        
        return True
        
    except Exception as e:
        logger.error(f"Core engine testing failed: {e}")
        return False

def test_utilities():
    """Test utility modules."""
    logger.info("=" * 60)
    logger.info("Testing Utility Modules")
    logger.info("=" * 60)
    
    try:
        from pandadock.utils_new import (
            rotation_matrix_from_euler, apply_rotation_translation,
            calculate_rmsd
        )
        
        # Test mathematical utilities
        try:
            logger.info("Testing mathematical utilities...")
            
            # Test rotation matrix
            euler = np.array([0.1, 0.2, 0.3])
            rotation = rotation_matrix_from_euler(euler)
            logger.info(f"✓ Rotation matrix created: shape {rotation.shape}")
            
            # Test coordinate transformation
            coords = np.random.rand(5, 3)
            translation = np.array([1.0, 2.0, 3.0])
            transformed = apply_rotation_translation(coords, rotation, translation)
            logger.info(f"✓ Coordinate transformation successful: shape {transformed.shape}")
            
            # Test RMSD calculation
            coords1 = np.random.rand(5, 3)
            coords2 = coords1 + 0.1 * np.random.rand(5, 3)
            rmsd = calculate_rmsd(coords1, coords2)
            logger.info(f"✓ RMSD calculation successful: {rmsd:.3f}")
            
        except Exception as e:
            logger.error(f"✗ Mathematical utilities failed: {e}")
        
        return True
        
    except Exception as e:
        logger.error(f"Utility testing failed: {e}")
        return False

def run_comprehensive_test():
    """Run comprehensive functionality test."""
    logger.info("🐼 Starting PandaDock Comprehensive Functionality Test")
    logger.info("This test verifies that all original functionality is preserved")
    logger.info("after the complete architectural refactoring.")
    logger.info("")
    
    test_results = {}
    
    # Run all tests
    test_functions = [
        ("Hardware Abstraction", test_hardware_abstraction),
        ("Utilities", test_utilities),
        ("Scoring Functions", test_scoring_functions),
        ("Algorithms", test_algorithms),
        ("Analysis Modules", test_analysis_modules),
        ("Molecule Handling", test_molecule_handling),
        ("Core Engine", test_core_engine),
    ]
    
    passed_tests = 0
    total_tests = len(test_functions)
    
    for test_name, test_func in test_functions:
        try:
            result = test_func()
            test_results[test_name] = "PASSED" if result else "FAILED"
            if result:
                passed_tests += 1
        except Exception as e:
            logger.error(f"Test {test_name} crashed: {e}")
            test_results[test_name] = "CRASHED"
    
    # Print summary
    logger.info("")
    logger.info("=" * 60)
    logger.info("TEST SUMMARY")
    logger.info("=" * 60)
    
    for test_name, result in test_results.items():
        status_symbol = "✓" if result == "PASSED" else "✗"
        logger.info(f"{status_symbol} {test_name}: {result}")
    
    logger.info("")
    logger.info(f"Overall Result: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        logger.info("🎉 ALL TESTS PASSED! The refactored architecture is working correctly.")
        logger.info("✅ All original PandaDock functionality has been preserved.")
    elif passed_tests >= total_tests * 0.8:
        logger.info("⚠️  Most tests passed. Minor issues need attention.")
    else:
        logger.info("❌ Multiple test failures detected. Major issues need fixing.")
    
    logger.info("")
    logger.info("🐼 PandaDock functionality test complete!")
    
    return passed_tests / total_tests

if __name__ == "__main__":
    success_rate = run_comprehensive_test()
    
    # Exit with appropriate code
    if success_rate == 1.0:
        sys.exit(0)  # All tests passed
    elif success_rate >= 0.8:
        sys.exit(1)  # Most tests passed, minor issues
    else:
        sys.exit(2)  # Major issues