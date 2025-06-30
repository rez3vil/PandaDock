"""
Main docking engine for PandaDock.

This module orchestrates the entire docking workflow, from input validation
to result generation, with clean separation of concerns.
"""

import logging
import time
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path

from ..hardware import DeviceManager, ComputeBackend, PerformanceMonitor
from .pose_generator import PoseGenerator
from .result_processor import ResultProcessor


class DockingEngineError(Exception):
    """Exception raised by the docking engine."""
    pass


class DockingEngine:
    """
    Main orchestrator for molecular docking operations.
    
    Coordinates all aspects of the docking process including hardware management,
    algorithm selection, scoring, and result processing.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize the docking engine.
        
        Args:
            config: Configuration dictionary with docking parameters (optional)
        """
        self.config = config or {}
        self.logger = logging.getLogger(__name__)
        
        # Initialize hardware and performance monitoring
        self.device_manager = self._setup_device_manager()
        self.compute_backend = ComputeBackend(self.device_manager)
        self.performance_monitor = PerformanceMonitor(self.device_manager)
        
        # Initialize core components
        self.pose_generator = PoseGenerator(self.compute_backend, self.performance_monitor)
        self.result_processor = ResultProcessor(self.performance_monitor)
        
        # State tracking
        self.is_initialized = True
        self.current_protein = None
        self.current_ligand = None
        self.results = []
        
        self.logger.info(f"DockingEngine initialized with {self.device_manager.selected_device.name}")
    
    def _setup_device_manager(self) -> DeviceManager:
        """Set up device manager based on configuration."""
        prefer_gpu = self.config.get('use_gpu', True)
        gpu_id = self.config.get('gpu_id', None)
        
        return DeviceManager(prefer_gpu=prefer_gpu, gpu_id=gpu_id)
    
    def run_docking(self, protein_path: str, ligand_path: str, 
                   output_dir: str, **kwargs) -> Dict[str, Any]:
        """
        Run complete docking workflow.
        
        Args:
            protein_path: Path to protein structure file
            ligand_path: Path to ligand structure file
            output_dir: Output directory for results
            **kwargs: Additional docking parameters
            
        Returns:
            Dictionary with docking results and metadata
        """
        if not self.is_initialized:
            raise DockingEngineError("DockingEngine not properly initialized")
        
        start_time = time.time()
        self.logger.info(f"Starting docking: {Path(protein_path).name} + {Path(ligand_path).name}")
        
        try:
            # Step 1: Load and validate structures
            with self.performance_monitor.start_operation("load_structures"):
                protein, ligand = self._load_structures(protein_path, ligand_path)
            
            # Step 2: Prepare molecules if requested
            if self.config.get('prepare_molecules', False):
                with self.performance_monitor.start_operation("prepare_molecules"):
                    protein, ligand = self._prepare_molecules(protein, ligand)
            
            # Step 3: Define binding site
            with self.performance_monitor.start_operation("define_binding_site"):
                self._define_binding_site(protein)
            
            # Step 4: Set up flexible residues if requested
            if self.config.get('flexible_residues') or self.config.get('auto_flex'):
                with self.performance_monitor.start_operation("setup_flexible_residues"):
                    protein = self._setup_flexible_residues(protein)
            
            # Step 5: Initialize scoring function
            with self.performance_monitor.start_operation("initialize_scoring"):
                scoring_function = self._create_scoring_function()
            
            # Step 6: Initialize search algorithm
            with self.performance_monitor.start_operation("initialize_algorithm"):
                search_algorithm = self._create_search_algorithm(scoring_function)
            
            # Step 7: Run docking search
            with self.performance_monitor.start_operation("docking_search"):
                raw_results = self._run_search(search_algorithm, protein, ligand)
            
            # Step 8: Post-process results
            with self.performance_monitor.start_operation("post_processing"):
                processed_results = self._post_process_results(raw_results, protein, ligand)
            
            # Step 9: Save results
            with self.performance_monitor.start_operation("save_results"):
                self._save_results(processed_results, output_dir)
            
            # Generate final report
            elapsed_time = time.time() - start_time
            report = self._generate_final_report(processed_results, elapsed_time)
            
            self.logger.info(f"Docking completed successfully in {elapsed_time:.1f} seconds")
            
            return {
                'success': True,
                'results': processed_results,
                'elapsed_time': elapsed_time,
                'report': report,
                'performance_report': self.performance_monitor.generate_performance_report()
            }
            
        except Exception as e:
            elapsed_time = time.time() - start_time
            error_msg = f"Docking failed after {elapsed_time:.1f} seconds: {str(e)}"
            self.logger.error(error_msg)
            
            return {
                'success': False,
                'error': str(e),
                'elapsed_time': elapsed_time,
                'performance_report': self.performance_monitor.generate_performance_report()
            }
    
    def _load_structures(self, protein_path: str, ligand_path: str) -> Tuple[Any, Any]:
        """Load protein and ligand structures."""
        from ..molecules import ProteinHandler, LigandHandler
        
        try:
            protein_handler = ProteinHandler()
            ligand_handler = LigandHandler()
            
            protein = protein_handler.load_protein(protein_path)
            ligand = ligand_handler.load_ligand(ligand_path)
            
            self.current_protein = protein
            self.current_ligand = ligand
            
            self.logger.info(f"Loaded protein: {protein.n_atoms} atoms")
            self.logger.info(f"Loaded ligand: {ligand.n_atoms} atoms")
            
            return protein, ligand
            
        except Exception as e:
            raise DockingEngineError(f"Failed to load structures: {e}")
    
    def _prepare_molecules(self, protein: Any, ligand: Any) -> Tuple[Any, Any]:
        """Prepare protein and ligand for docking."""
        from ..molecules import StructurePreparation
        
        try:
            preparation = StructurePreparation()
            
            # Prepare protein
            if self.config.get('prepare_protein', True):
                prepared_protein = preparation.prepare_protein(
                    protein, 
                    ph=self.config.get('ph', 7.4),
                    add_hydrogens=True
                )
                self.logger.info("Protein preparation completed")
            else:
                prepared_protein = protein
            
            # Prepare ligand
            if self.config.get('prepare_ligand', True):
                prepared_ligand = preparation.prepare_ligand(
                    ligand,
                    add_hydrogens=True,
                    generate_conformers=self.config.get('generate_conformers', True)
                )
                self.logger.info("Ligand preparation completed")
            else:
                prepared_ligand = ligand
            
            return prepared_protein, prepared_ligand
            
        except Exception as e:
            raise DockingEngineError(f"Molecule preparation failed: {e}")
    
    def _define_binding_site(self, protein: Any) -> None:
        """Define the binding site for docking."""
        try:
            if self.config.get('binding_site_center'):
                # Use provided center
                center = self.config['binding_site_center']
                radius = self.config.get('binding_site_radius', 10.0)
                protein.define_active_site(center, radius)
                self.logger.info(f"Using provided binding site: {center}, radius {radius}")
                
            elif self.config.get('detect_pockets', False):
                # Auto-detect binding pockets
                pockets = protein.detect_pockets()
                if pockets:
                    protein.define_active_site(pockets[0]['center'], pockets[0]['radius'])
                    self.logger.info(f"Auto-detected binding site: {len(pockets)} pockets found")
                else:
                    self.logger.warning("No pockets detected, using whole protein")
                    
            else:
                self.logger.info("No binding site specified, using whole protein")
                
        except Exception as e:
            raise DockingEngineError(f"Binding site definition failed: {e}")
    
    def _setup_flexible_residues(self, protein: Any) -> Any:
        """Set up flexible residues for protein flexibility."""
        try:
            flexible_residues = []
            
            if self.config.get('flexible_residues'):
                # Use specified residues
                flexible_residues = self.config['flexible_residues']
                self.logger.info(f"Using specified flexible residues: {len(flexible_residues)}")
                
            elif self.config.get('auto_flex', False):
                # Auto-detect flexible residues
                from ..molecules import FlexibleResidueDetector
                detector = FlexibleResidueDetector()
                flexible_residues = detector.detect_flexible_residues(
                    protein,
                    max_residues=self.config.get('max_flex_residues', 5)
                )
                self.logger.info(f"Auto-detected flexible residues: {len(flexible_residues)}")
            
            if flexible_residues:
                protein.define_flexible_residues(
                    flexible_residues,
                    max_rotatable_bonds=self.config.get('max_flex_bonds', 3)
                )
                self.logger.info(f"Set up {len(flexible_residues)} flexible residues")
            
            return protein
            
        except Exception as e:
            raise DockingEngineError(f"Flexible residue setup failed: {e}")
    
    def _create_scoring_function(self) -> Any:
        """Create and configure the scoring function."""
        try:
            from ..scoring import ScoringFunctionFactory
            
            factory = ScoringFunctionFactory(self.compute_backend)
            
            scoring_type = self._determine_scoring_type()
            scoring_function = factory.create_scoring_function(
                scoring_type=scoring_type,
                enhanced=self.config.get('enhanced_scoring', False),
                physics_based=self.config.get('physics_based', False)
            )
            
            self.logger.info(f"Created {scoring_type} scoring function")
            return scoring_function
            
        except Exception as e:
            raise DockingEngineError(f"Scoring function creation failed: {e}")
    
    def _determine_scoring_type(self) -> str:
        """Determine the appropriate scoring function type."""
        if self.config.get('physics_based', False):
            return 'physics'
        elif self.config.get('enhanced_scoring', False):
            return 'enhanced'
        else:
            return 'standard'
    
    def _create_search_algorithm(self, scoring_function: Any) -> Any:
        """Create and configure the search algorithm."""
        try:
            from ..algorithms import AlgorithmFactory
            
            factory = AlgorithmFactory(self.compute_backend)
            
            algorithm_type = self._determine_algorithm_type()
            algorithm_params = self._get_algorithm_parameters()
            
            search_algorithm = factory.create_algorithm(
                algorithm_type=algorithm_type,
                scoring_function=scoring_function,
                **algorithm_params
            )
            
            self.logger.info(f"Created {algorithm_type} search algorithm")
            return search_algorithm
            
        except Exception as e:
            raise DockingEngineError(f"Search algorithm creation failed: {e}")
    
    def _determine_algorithm_type(self) -> str:
        """Determine the appropriate search algorithm."""
        if self.config.get('monte_carlo', False):
            return 'monte-carlo'
        else:
            return self.config.get('algorithm', 'genetic')
    
    def _get_algorithm_parameters(self) -> Dict[str, Any]:
        """Get algorithm-specific parameters."""
        algorithm_type = self._determine_algorithm_type()
        params = {}
        
        # Common parameters
        params['max_iterations'] = self.config.get('iterations', 100)
        params['output_dir'] = self.config.get('output_dir')
        
        # Algorithm-specific parameters
        if algorithm_type == 'genetic':
            params.update({
                'population_size': self.config.get('population_size', 100),
                'mutation_rate': self.config.get('mutation_rate', 0.2),
                'local_optimization': self.config.get('local_opt', False)
            })
        elif algorithm_type == 'monte-carlo':
            params.update({
                'n_steps': self.config.get('mc_steps', 1000),
                'temperature': self.config.get('temperature', 300.0),
                'cooling_factor': self.config.get('cooling_factor', 0.95)
            })
        elif algorithm_type == 'pandadock':
            params.update({
                'high_temp': self.config.get('high_temp', 1000.0),
                'target_temp': self.config.get('target_temp', 300.0),
                'num_conformers': self.config.get('num_conformers', 10),
                'md_steps': self.config.get('md_steps', 1000)
            })
        
        return params
    
    def _run_search(self, algorithm: Any, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Run the docking search."""
        try:
            # Check for reference-based docking
            if self.config.get('reference_ligand'):
                return self._run_reference_docking(algorithm, protein, ligand)
            elif self.config.get('exhaustiveness', 1) > 1:
                return self._run_ensemble_docking(algorithm, protein, ligand)
            else:
                return self._run_standard_docking(algorithm, protein, ligand)
                
        except Exception as e:
            raise DockingEngineError(f"Docking search failed: {e}")
    
    def _run_standard_docking(self, algorithm: Any, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Run standard docking search."""
        self.logger.info("Running standard docking search")
        return algorithm.search(protein, ligand)
    
    def _run_reference_docking(self, algorithm: Any, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Run reference-guided docking."""
        reference_ligand = self.config['reference_ligand']
        
        if self.config.get('exact_alignment', False):
            self.logger.info("Running exact reference alignment")
            return algorithm.exact_reference_docking(protein, ligand, reference_ligand)
        else:
            self.logger.info("Running reference-guided docking")
            return algorithm.reference_guided_docking(protein, ligand, reference_ligand)
    
    def _run_ensemble_docking(self, algorithm: Any, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Run ensemble docking with multiple runs."""
        n_runs = self.config['exhaustiveness']
        self.logger.info(f"Running ensemble docking with {n_runs} runs")
        
        all_results = []
        for run_id in range(n_runs):
            self.logger.info(f"Running ensemble iteration {run_id + 1}/{n_runs}")
            run_results = algorithm.search(protein, ligand)
            all_results.extend(run_results)
        
        return all_results
    
    def _post_process_results(self, raw_results: List[Tuple[Any, float]], 
                            protein: Any, ligand: Any) -> Dict[str, Any]:
        """Post-process docking results."""
        try:
            return self.result_processor.process_results(
                raw_results, protein, ligand, self.config
            )
            
        except Exception as e:
            raise DockingEngineError(f"Result post-processing failed: {e}")
    
    def _save_results(self, results: Dict[str, Any], output_dir: str) -> None:
        """Save docking results to files."""
        try:
            from ..io import ResultWriters
            
            writers = ResultWriters(output_dir)
            writers.save_all_results(results, self.config)
            
            self.logger.info(f"Results saved to {output_dir}")
            
        except Exception as e:
            raise DockingEngineError(f"Result saving failed: {e}")
    
    def _generate_final_report(self, results: Dict[str, Any], elapsed_time: float) -> Dict[str, Any]:
        """Generate final docking report."""
        best_score = None
        n_poses = 0
        
        if results.get('poses'):
            n_poses = len(results['poses'])
            if n_poses > 0:
                best_score = min(pose['score'] for pose in results['poses'])
        
        return {
            'total_poses': n_poses,
            'best_score': best_score,
            'elapsed_time': elapsed_time,
            'device_used': self.device_manager.selected_device.name,
            'algorithm_used': self._determine_algorithm_type(),
            'scoring_used': self._determine_scoring_type()
        }
    
    def cleanup(self) -> None:
        """Clean up resources."""
        try:
            if hasattr(self, 'compute_backend'):
                self.compute_backend.cleanup()
            
            if hasattr(self, 'device_manager'):
                self.device_manager.cleanup()
            
            self.logger.info("DockingEngine cleanup completed")
            
        except Exception as e:
            self.logger.warning(f"Error during cleanup: {e}")
    
    def get_status(self) -> Dict[str, Any]:
        """Get current engine status."""
        return {
            'initialized': self.is_initialized,
            'device': self.device_manager.selected_device.name if self.device_manager else None,
            'protein_loaded': self.current_protein is not None,
            'ligand_loaded': self.current_ligand is not None,
            'memory_info': self.device_manager.get_memory_info() if self.device_manager else None
        }