"""
Batch screening module for PandaDock.

This module provides functionality for high-throughput batch processing
of ligand libraries with advanced parallel processing capabilities.
"""

import os
import time
import json
import csv
import glob
import multiprocessing as mp
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any, Union
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from .virtual_screening import VirtualScreeningManager
from ..core.docking_engine import DockingEngine
from ..algorithms.algorithm_factory import AlgorithmFactory
from ..scoring.scoring_factory import ScoringFunctionFactory
from ..hardware import DeviceManager, PerformanceMonitor

# Mock molecule classes for now
class Protein:
    def __init__(self, file_path):
        self.coords = np.random.randn(10, 3) * 5
        self.atoms = [{'element': 'C', 'coords': coord} for coord in self.coords]
        self.active_site = None

class Ligand:
    def __init__(self, file_path):
        self.coords = np.random.randn(5, 3) * 2
        self.atoms = [{'element': 'C', 'coords': coord} for coord in self.coords]


class BatchScreeningManager:
    """
    Manager for high-throughput batch screening of ligand libraries.
    
    Provides advanced features including:
    - Parallel processing with load balancing
    - Progress tracking and resumption
    - Error handling and recovery
    - Comprehensive reporting and analysis
    """
    
    def __init__(self, protein_file: Union[str, Path],
                 ligand_library: Union[str, Path],
                 output_dir: Union[str, Path],
                 screening_params: Optional[Dict[str, Any]] = None,
                 n_processes: Optional[int] = None):
        """
        Initialize batch screening manager.
        
        Args:
            protein_file: Path to protein PDB file
            ligand_library: Path to directory containing ligand files
            output_dir: Output directory for results
            screening_params: Dictionary of screening parameters
            n_processes: Number of parallel processes (None = auto-detect)
        """
        self.protein_file = Path(protein_file)
        self.ligand_library = Path(ligand_library)
        self.output_dir = Path(output_dir)
        self.screening_params = screening_params or {}
        self.n_processes = n_processes or mp.cpu_count()
        
        # Validate inputs
        self._validate_inputs()
        
        # Setup hardware and monitoring
        self.device_manager = DeviceManager()
        self.performance_monitor = PerformanceMonitor()
        
        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()
        
        # Initialize screening statistics
        self.screening_stats = {
            'start_time': None,
            'end_time': None,
            'total_ligands': 0,
            'processed_ligands': 0,
            'successful_ligands': 0,
            'failed_ligands': 0,
            'average_time_per_ligand': 0.0
        }
    
    def _validate_inputs(self) -> None:
        """Validate input files and directories."""
        if not self.protein_file.exists():
            raise FileNotFoundError(f"Protein file not found: {self.protein_file}")
        
        if not self.ligand_library.exists():
            raise FileNotFoundError(f"Ligand library not found: {self.ligand_library}")
        
        if not self.ligand_library.is_dir():
            raise ValueError(f"Ligand library must be a directory: {self.ligand_library}")
    
    def _setup_logging(self) -> None:
        """Setup logging for batch screening."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        log_file = self.output_dir / "batch_screening.log"
        
        # Create file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(formatter)
        
        # Add handler to logger
        self.logger.addHandler(file_handler)
    
    def run_batch_screening(self) -> Dict[str, Any]:
        """
        Run complete batch screening process.
        
        Returns:
            Dictionary containing screening results and statistics
        """
        self.screening_stats['start_time'] = datetime.now()
        self.logger.info("Starting batch screening process")
        
        try:
            # Prepare protein
            protein = self._prepare_protein()
            
            # Get ligand files
            ligand_files = self._get_ligand_files()
            self.screening_stats['total_ligands'] = len(ligand_files)
            
            self.logger.info(f"Found {len(ligand_files)} ligands for screening")
            self.logger.info(f"Using {self.n_processes} parallel processes")
            
            # Create output structure
            self._create_output_structure()
            
            # Initialize summary files
            self._initialize_summary_files()
            
            # Run screening
            if self.n_processes > 1:
                results = self._run_parallel_batch_screening(protein, ligand_files)
            else:
                results = self._run_sequential_batch_screening(protein, ligand_files)
            
            # Generate final reports
            self._generate_final_reports(results)
            
            # Update statistics
            self.screening_stats['end_time'] = datetime.now()
            self.screening_stats['processed_ligands'] = len(results)
            self.screening_stats['successful_ligands'] = len([r for r in results.values() 
                                                            if not r.get('error')])
            self.screening_stats['failed_ligands'] = (self.screening_stats['processed_ligands'] - 
                                                    self.screening_stats['successful_ligands'])
            
            total_time = (self.screening_stats['end_time'] - 
                         self.screening_stats['start_time']).total_seconds()
            
            if self.screening_stats['processed_ligands'] > 0:
                self.screening_stats['average_time_per_ligand'] = (
                    total_time / self.screening_stats['processed_ligands']
                )
            
            self.logger.info("Batch screening completed successfully")
            
            return {
                'results': results,
                'statistics': self.screening_stats,
                'output_directory': str(self.output_dir)
            }
            
        except Exception as e:
            self.logger.error(f"Batch screening failed: {e}")
            raise
    
    def _prepare_protein(self) -> Protein:
        """Prepare protein for screening."""
        self.logger.info(f"Loading protein: {self.protein_file}")
        
        # Load protein
        protein = Protein(str(self.protein_file))
        
        # Setup active site
        if 'active_site_center' in self.screening_params:
            center = np.array(self.screening_params['active_site_center'])
            radius = self.screening_params.get('active_site_radius', 10.0)
            protein.active_site = {'center': center, 'radius': radius}
            self.logger.info(f"Using specified active site at {center} with radius {radius}Å")
        elif self.screening_params.get('detect_binding_sites', False):
            # Auto-detect binding sites
            self.logger.info("Auto-detecting binding sites")
            binding_sites = self._detect_binding_sites(protein)
            if binding_sites:
                protein.active_site = binding_sites[0]
                self.logger.info(f"Using detected binding site at {binding_sites[0]['center']}")
            else:
                # Fallback to protein center
                center = np.mean(protein.coords, axis=0)
                protein.active_site = {'center': center, 'radius': 10.0}
                self.logger.warning(f"No binding sites detected. Using protein center: {center}")
        else:
            # Use protein center as default
            center = np.mean(protein.coords, axis=0)
            protein.active_site = {'center': center, 'radius': 10.0}
            self.logger.info(f"Using protein center as active site: {center}")
        
        return protein
    
    def _detect_binding_sites(self, protein: Protein) -> List[Dict[str, Any]]:
        """Detect potential binding sites in protein."""
        # Simple binding site detection based on cavities
        # In a real implementation, this would use more sophisticated algorithms
        
        coords = protein.coords
        center = np.mean(coords, axis=0)
        
        # For now, just return the protein center as a binding site
        return [{'center': center, 'radius': 10.0}]
    
    def _get_ligand_files(self) -> List[Path]:
        """Get list of ligand files from library directory."""
        supported_extensions = ['*.mol', '*.mol2', '*.sdf', '*.pdb']
        ligand_files = []
        
        for extension in supported_extensions:
            pattern = str(self.ligand_library / extension)
            ligand_files.extend([Path(f) for f in glob.glob(pattern)])
        
        if not ligand_files:
            raise ValueError(f"No ligand files found in {self.ligand_library}")
        
        return sorted(ligand_files)
    
    def _create_output_structure(self) -> None:
        """Create output directory structure."""
        directories = [
            'ligands',      # Individual ligand results
            'complexes',    # Protein-ligand complexes
            'reports',      # Summary reports
            'logs',         # Log files
            'checkpoints'   # Checkpoint files for resumption
        ]
        
        for dir_name in directories:
            (self.output_dir / dir_name).mkdir(parents=True, exist_ok=True)
    
    def _initialize_summary_files(self) -> None:
        """Initialize summary CSV files."""
        summary_file = self.output_dir / "screening_summary.csv"
        
        with open(summary_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Ligand_Name', 'Best_Score', 'Processing_Time_s', 
                'Status', 'Error_Message', 'Timestamp'
            ])
    
    def _run_sequential_batch_screening(self, protein: Protein, 
                                      ligand_files: List[Path]) -> Dict[str, Any]:
        """Run batch screening sequentially."""
        results = {}
        
        for i, ligand_file in enumerate(ligand_files):
            ligand_name = ligand_file.stem
            self.logger.info(f"Processing ligand {i+1}/{len(ligand_files)}: {ligand_name}")
            
            try:
                # Process single ligand
                result = self._process_single_ligand_batch(protein, ligand_file)
                results[ligand_name] = result
                
                # Update summary file
                self._update_summary_file(ligand_name, result)
                
                # Save checkpoint
                if (i + 1) % 10 == 0:  # Save checkpoint every 10 ligands
                    self._save_checkpoint(results, i + 1)
                
            except Exception as e:
                self.logger.error(f"Error processing ligand {ligand_name}: {e}")
                result = {
                    'error': str(e),
                    'best_score': None,
                    'processing_time': 0.0,
                    'poses': []
                }
                results[ligand_name] = result
                self._update_summary_file(ligand_name, result)
        
        return results
    
    def _run_parallel_batch_screening(self, protein: Protein, 
                                    ligand_files: List[Path]) -> Dict[str, Any]:
        """Run batch screening in parallel."""
        results = {}
        
        # Prepare tasks for parallel execution
        tasks = []
        for ligand_file in ligand_files:
            task = {
                'protein_file': str(self.protein_file),
                'ligand_file': str(ligand_file),
                'ligand_name': ligand_file.stem,
                'screening_params': self.screening_params,
                'output_dir': str(self.output_dir / 'ligands' / ligand_file.stem)
            }
            tasks.append(task)
        
        # Process tasks in parallel
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            # Submit all tasks
            future_to_name = {
                executor.submit(_process_ligand_batch_worker, task): task['ligand_name']
                for task in tasks
            }
            
            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_name):
                ligand_name = future_to_name[future]
                completed += 1
                
                try:
                    result = future.result()
                    results[ligand_name] = result
                    self.logger.info(f"Completed {completed}/{len(ligand_files)}: {ligand_name}")
                    
                    # Update summary file
                    self._update_summary_file(ligand_name, result)
                    
                    # Save checkpoint periodically
                    if completed % 20 == 0:
                        self._save_checkpoint(results, completed)
                    
                except Exception as e:
                    self.logger.error(f"Error processing ligand {ligand_name}: {e}")
                    result = {
                        'error': str(e),
                        'best_score': None,
                        'processing_time': 0.0,
                        'poses': []
                    }
                    results[ligand_name] = result
                    self._update_summary_file(ligand_name, result)
        
        return results
    
    def _process_single_ligand_batch(self, protein: Protein, ligand_file: Path) -> Dict[str, Any]:
        """Process a single ligand in batch mode."""
        start_time = time.time()
        ligand_name = ligand_file.stem
        
        try:
            # Load ligand
            ligand = Ligand(str(ligand_file))
            
            # Create virtual screening manager for this ligand
            vs_manager = VirtualScreeningManager(
                scoring_function_type=self.screening_params.get('scoring_function', 'physics_based'),
                algorithm_type=self.screening_params.get('algorithm', 'genetic'),
                output_dir=self.output_dir / 'ligands' / ligand_name,
                n_workers=1,  # Single worker for batch processing
                exhaustiveness=self.screening_params.get('exhaustiveness', 8),
                num_modes=self.screening_params.get('num_modes', 9)
            )
            
            # Run screening for single ligand
            screening_result = vs_manager.run_screening(protein, [ligand], [ligand_name])
            
            # Extract results
            ligand_results = screening_result['results'].get(ligand_name, {})
            
            processing_time = time.time() - start_time
            
            return {
                'best_score': ligand_results.get('best_score'),
                'poses': ligand_results.get('poses', []),
                'processing_time': processing_time,
                'statistics': ligand_results.get('statistics', {}),
                'analysis': ligand_results.get('analysis', {})
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            raise Exception(f"Failed to process ligand {ligand_name}: {e}")
    
    def _update_summary_file(self, ligand_name: str, result: Dict[str, Any]) -> None:
        """Update summary CSV file with ligand result."""
        summary_file = self.output_dir / "screening_summary.csv"
        
        with open(summary_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                ligand_name,
                result.get('best_score', 'N/A'),
                f"{result.get('processing_time', 0.0):.2f}",
                'Success' if not result.get('error') else 'Failed',
                result.get('error', ''),
                datetime.now().isoformat()
            ])
    
    def _save_checkpoint(self, results: Dict[str, Any], completed: int) -> None:
        """Save checkpoint for resuming interrupted screening."""
        checkpoint_file = self.output_dir / 'checkpoints' / f"checkpoint_{completed}.json"
        
        checkpoint_data = {
            'completed_ligands': completed,
            'total_ligands': self.screening_stats['total_ligands'],
            'results': results,
            'timestamp': datetime.now().isoformat(),
            'screening_params': self.screening_params
        }
        
        with open(checkpoint_file, 'w') as f:
            json.dump(checkpoint_data, f, indent=2, default=str)
        
        self.logger.info(f"Checkpoint saved: {completed} ligands completed")
    
    def _generate_final_reports(self, results: Dict[str, Any]) -> None:
        """Generate comprehensive final reports."""
        # Generate detailed statistics report
        self._generate_statistics_report(results)
        
        # Generate ranked results
        self._generate_ranked_results(results)
        
        # Generate performance analysis
        self._generate_performance_analysis(results)
        
        # Generate visualizations if possible
        try:
            self._generate_visualizations(results)
        except ImportError:
            self.logger.warning("Matplotlib not available. Skipping visualizations.")
    
    def _generate_statistics_report(self, results: Dict[str, Any]) -> None:
        """Generate detailed statistics report."""
        report_file = self.output_dir / 'reports' / 'statistics_report.txt'
        
        # Calculate statistics
        successful_results = {k: v for k, v in results.items() if not v.get('error')}
        failed_results = {k: v for k, v in results.items() if v.get('error')}
        
        if successful_results:
            scores = [r['best_score'] for r in successful_results.values() if r.get('best_score') is not None]
            processing_times = [r['processing_time'] for r in successful_results.values()]
        else:
            scores = []
            processing_times = []
        
        with open(report_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("              PandaDock Batch Screening Statistics\n")
            f.write("=" * 70 + "\n\n")
            
            # General statistics
            f.write("GENERAL STATISTICS\n")
            f.write("-" * 18 + "\n")
            f.write(f"Total Ligands Processed: {len(results)}\n")
            f.write(f"Successful: {len(successful_results)}\n")
            f.write(f"Failed: {len(failed_results)}\n")
            f.write(f"Success Rate: {len(successful_results)/len(results)*100:.1f}%\n\n")
            
            # Timing statistics
            total_time = (self.screening_stats['end_time'] - 
                         self.screening_stats['start_time']).total_seconds()
            f.write("TIMING STATISTICS\n")
            f.write("-" * 17 + "\n")
            f.write(f"Total Runtime: {total_time:.2f} seconds ({total_time/3600:.2f} hours)\n")
            f.write(f"Average Time per Ligand: {total_time/len(results):.2f} seconds\n")
            
            if processing_times:
                f.write(f"Fastest Ligand: {min(processing_times):.2f} seconds\n")
                f.write(f"Slowest Ligand: {max(processing_times):.2f} seconds\n")
                f.write(f"Median Processing Time: {np.median(processing_times):.2f} seconds\n\n")
            
            # Score statistics
            if scores:
                f.write("SCORING STATISTICS\n")
                f.write("-" * 18 + "\n")
                f.write(f"Best Score: {min(scores):.3f}\n")
                f.write(f"Worst Score: {max(scores):.3f}\n")
                f.write(f"Mean Score: {np.mean(scores):.3f}\n")
                f.write(f"Median Score: {np.median(scores):.3f}\n")
                f.write(f"Standard Deviation: {np.std(scores):.3f}\n\n")
            
            # Error analysis
            if failed_results:
                f.write("ERROR ANALYSIS\n")
                f.write("-" * 14 + "\n")
                error_types = {}
                for result in failed_results.values():
                    error = result.get('error', 'Unknown error')
                    error_type = error.split(':')[0] if ':' in error else error
                    error_types[error_type] = error_types.get(error_type, 0) + 1
                
                for error_type, count in sorted(error_types.items(), key=lambda x: x[1], reverse=True):
                    f.write(f"{error_type}: {count} occurrences\n")
    
    def _generate_ranked_results(self, results: Dict[str, Any]) -> None:
        """Generate ranked results table."""
        ranked_file = self.output_dir / 'reports' / 'ranked_results.csv'
        
        # Filter and sort successful results
        successful_results = [(name, data) for name, data in results.items() 
                            if not data.get('error') and data.get('best_score') is not None]
        
        sorted_results = sorted(successful_results, key=lambda x: x[1]['best_score'])
        
        with open(ranked_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Rank', 'Ligand_Name', 'Best_Score', 'Processing_Time_s', 
                'Number_of_Poses', 'Score_Range'
            ])
            
            for i, (name, data) in enumerate(sorted_results):
                poses = data.get('poses', [])
                pose_scores = [p.get('score', 0) for p in poses if 'score' in p]
                score_range = max(pose_scores) - min(pose_scores) if len(pose_scores) > 1 else 0
                
                writer.writerow([
                    i + 1,
                    name,
                    f"{data['best_score']:.3f}",
                    f"{data.get('processing_time', 0):.2f}",
                    len(poses),
                    f"{score_range:.3f}"
                ])
    
    def _generate_performance_analysis(self, results: Dict[str, Any]) -> None:
        """Generate performance analysis report."""
        perf_file = self.output_dir / 'reports' / 'performance_analysis.json'
        
        # Calculate performance metrics
        processing_times = [r.get('processing_time', 0) for r in results.values()]
        
        performance_data = {
            'hardware_info': {
                'cpu_count': mp.cpu_count(),
                'processes_used': self.n_processes,
                'device_info': self.device_manager.get_device_info()
            },
            'throughput_metrics': {
                'ligands_per_hour': len(results) / (
                    (self.screening_stats['end_time'] - 
                     self.screening_stats['start_time']).total_seconds() / 3600
                ),
                'average_processing_time': np.mean(processing_times) if processing_times else 0,
                'processing_time_std': np.std(processing_times) if processing_times else 0,
                'fastest_ligand_time': min(processing_times) if processing_times else 0,
                'slowest_ligand_time': max(processing_times) if processing_times else 0
            },
            'screening_parameters': self.screening_params,
            'resource_utilization': {
                'parallel_efficiency': self._calculate_parallel_efficiency(processing_times),
                'estimated_sequential_time': sum(processing_times) if processing_times else 0,
                'actual_wall_time': (
                    self.screening_stats['end_time'] - 
                    self.screening_stats['start_time']
                ).total_seconds()
            }
        }
        
        with open(perf_file, 'w') as f:
            json.dump(performance_data, f, indent=2, default=str)
    
    def _calculate_parallel_efficiency(self, processing_times: List[float]) -> float:
        """Calculate parallel processing efficiency."""
        if not processing_times:
            return 0.0
        
        estimated_sequential_time = sum(processing_times)
        actual_wall_time = (
            self.screening_stats['end_time'] - 
            self.screening_stats['start_time']
        ).total_seconds()
        
        if actual_wall_time > 0:
            speedup = estimated_sequential_time / actual_wall_time
            efficiency = speedup / self.n_processes
            return min(efficiency, 1.0)  # Cap at 100% efficiency
        
        return 0.0
    
    def _generate_visualizations(self, results: Dict[str, Any]) -> None:
        """Generate visualization plots."""
        import matplotlib.pyplot as plt
        
        plots_dir = self.output_dir / 'reports' / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        # Filter successful results
        successful_results = {k: v for k, v in results.items() if not v.get('error')}
        
        if not successful_results:
            return
        
        scores = [r['best_score'] for r in successful_results.values() if r.get('best_score') is not None]
        processing_times = [r['processing_time'] for r in successful_results.values()]
        
        # Score distribution
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        plt.hist(scores, bins=min(30, len(scores)//5 + 1), alpha=0.7, color='skyblue')
        plt.axvline(np.mean(scores), color='red', linestyle='--', label=f'Mean: {np.mean(scores):.2f}')
        plt.xlabel('Docking Score')
        plt.ylabel('Frequency')
        plt.title('Distribution of Best Docking Scores')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Processing time distribution
        plt.subplot(1, 2, 2)
        plt.hist(processing_times, bins=min(30, len(processing_times)//5 + 1), alpha=0.7, color='lightcoral')
        plt.axvline(np.mean(processing_times), color='blue', linestyle='--', 
                   label=f'Mean: {np.mean(processing_times):.2f}s')
        plt.xlabel('Processing Time (seconds)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Processing Times')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'distributions.png', dpi=300)
        plt.close()
        
        # Score vs processing time scatter plot
        plt.figure(figsize=(10, 6))
        plt.scatter(processing_times, scores, alpha=0.6, s=50)
        plt.xlabel('Processing Time (seconds)')
        plt.ylabel('Best Docking Score')
        plt.title('Docking Score vs Processing Time')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(plots_dir / 'score_vs_time.png', dpi=300)
        plt.close()


def _process_ligand_batch_worker(task: Dict[str, Any]) -> Dict[str, Any]:
    """Worker function for parallel batch ligand processing."""
    try:
        # Extract task parameters
        protein_file = task['protein_file']
        ligand_file = task['ligand_file']
        ligand_name = task['ligand_name']
        screening_params = task['screening_params']
        output_dir = task['output_dir']
        
        start_time = time.time()
        
        # Load protein and ligand
        protein = Protein(protein_file)
        ligand = Ligand(ligand_file)
        
        # Setup active site for protein
        if 'active_site_center' in screening_params:
            center = np.array(screening_params['active_site_center'])
            radius = screening_params.get('active_site_radius', 10.0)
            protein.active_site = {'center': center, 'radius': radius}
        else:
            center = np.mean(protein.coords, axis=0)
            protein.active_site = {'center': center, 'radius': 10.0}
        
        # Create virtual screening manager for this worker
        vs_manager = VirtualScreeningManager(
            scoring_function_type=screening_params.get('scoring_function', 'physics_based'),
            algorithm_type=screening_params.get('algorithm', 'genetic'),
            output_dir=output_dir,
            n_workers=1,
            exhaustiveness=screening_params.get('exhaustiveness', 8),
            num_modes=screening_params.get('num_modes', 9)
        )
        
        # Run screening
        screening_result = vs_manager.run_screening(protein, [ligand], [ligand_name])
        
        # Extract results
        ligand_results = screening_result['results'].get(ligand_name, {})
        
        processing_time = time.time() - start_time
        
        return {
            'best_score': ligand_results.get('best_score'),
            'poses': ligand_results.get('poses', []),
            'processing_time': processing_time,
            'statistics': ligand_results.get('statistics', {}),
            'analysis': ligand_results.get('analysis', {})
        }
        
    except Exception as e:
        processing_time = time.time() - start_time if 'start_time' in locals() else 0.0
        return {
            'error': str(e),
            'best_score': None,
            'processing_time': processing_time,
            'poses': []
        }