"""
Virtual screening module for PandaDock.

This module provides efficient screening of multiple ligands against a protein target,
with support for parallel processing, pose clustering, and results analysis.
"""

import os
import csv
import time
import json
import numpy as np
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, Union
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

from ..core.docking_engine import DockingEngine
from ..algorithms.algorithm_factory import AlgorithmFactory
from ..scoring.scoring_factory import ScoringFunctionFactory
from ..hardware import DeviceManager, PerformanceMonitor
from ..io.pdb_writer import PDBWriter


class VirtualScreeningManager:
    """
    Manager for virtual screening of multiple ligands against a protein target.
    Coordinates the screening process, distributes work, and manages results.
    """
    
    def __init__(self, scoring_function_type: str = 'physics_based',
                 algorithm_type: str = 'genetic',
                 output_dir: Optional[Union[str, Path]] = None,
                 n_workers: Optional[int] = None,
                 exhaustiveness: int = 8,
                 num_modes: int = 9,
                 max_evals: int = 10000,
                 rmsd_threshold: float = 2.0):
        """
        Initialize virtual screening manager.
        
        Args:
            scoring_function_type: Type of scoring function to use
            algorithm_type: Type of search algorithm to use
            output_dir: Directory for output files
            n_workers: Number of parallel workers (None = auto-detect)
            exhaustiveness: Search exhaustiveness (higher = more thorough)
            num_modes: Number of binding modes to generate per ligand
            max_evals: Maximum number of pose evaluations per ligand
            rmsd_threshold: RMSD threshold for clustering poses (Å)
        """
        self.scoring_function_type = scoring_function_type
        self.algorithm_type = algorithm_type
        self.output_dir = Path(output_dir) if output_dir else None
        
        # Hardware setup
        self.device_manager = DeviceManager()
        self.performance_monitor = PerformanceMonitor()
        self.n_workers = n_workers or mp.cpu_count()
        
        # Screening parameters
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.max_evals = max_evals
        self.rmsd_threshold = rmsd_threshold
        
        # Setup logging
        self.logger = logging.getLogger(__name__)
        if self.output_dir:
            self._setup_logging()
    
    def _setup_logging(self) -> None:
        """Setup logging for screening process."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        log_file = self.output_dir / "virtual_screening.log"
        
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
    
    def run_screening(self, protein: Any, ligands: List[Any], 
                     ligand_names: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Run virtual screening on multiple ligands.
        
        Args:
            protein: Protein target object
            ligands: List of ligand objects
            ligand_names: Optional list of ligand names/identifiers
            
        Returns:
            Dictionary of screening results
        """
        start_time = time.time()
        
        # Setup output directory structure
        if self.output_dir:
            self._create_output_structure()
        
        # Assign names to ligands if not provided
        if ligand_names is None:
            ligand_names = [f"ligand_{i+1}" for i in range(len(ligands))]
        
        # Setup active site if not defined
        self._setup_active_site(protein)
        
        # Initialize components
        self._initialize_screening_components()
        
        # Progress tracking
        total_ligands = len(ligands)
        self.logger.info(f"Starting virtual screening of {total_ligands} ligands")
        self.logger.info(f"Using {self.n_workers} workers")
        self.logger.info(f"Algorithm: {self.algorithm_type}, Scoring: {self.scoring_function_type}")
        
        # Process ligands
        if self.n_workers > 1:
            results = self._run_parallel_screening(protein, ligands, ligand_names)
        else:
            results = self._run_sequential_screening(protein, ligands, ligand_names)
        
        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        
        # Generate summary report
        if self.output_dir:
            self._generate_summary_report(results, elapsed_time)
        
        self.logger.info(f"Virtual screening completed in {elapsed_time:.2f} seconds")
        
        return {
            'results': results,
            'statistics': self._calculate_screening_statistics(results, elapsed_time),
            'output_directory': str(self.output_dir) if self.output_dir else None
        }
    
    def _create_output_structure(self) -> None:
        """Create output directory structure."""
        directories = ['poses', 'complexes', 'logs', 'reports']
        for dir_name in directories:
            (self.output_dir / dir_name).mkdir(parents=True, exist_ok=True)
    
    def _setup_active_site(self, protein: Any) -> None:
        """Setup active site for protein if not already defined."""
        if not hasattr(protein, 'active_site') or not protein.active_site:
            if hasattr(protein, 'coords'):
                center = np.mean(protein.coords, axis=0)
                radius = 10.0
                protein.active_site = {'center': center, 'radius': radius}
                self.logger.info(f"Auto-detected active site at {center} with radius {radius}Å")
            else:
                # Default active site
                protein.active_site = {'center': np.array([0., 0., 0.]), 'radius': 10.0}
                self.logger.warning("Using default active site at origin")
    
    def _initialize_screening_components(self) -> None:
        """Initialize scoring function and algorithm factories."""
        self.scoring_factory = ScoringFunctionFactory()
        self.algorithm_factory = AlgorithmFactory()
    
    def _run_sequential_screening(self, protein: Any, ligands: List[Any], 
                                ligand_names: List[str]) -> Dict[str, Any]:
        """Run screening sequentially."""
        results = {}
        
        for i, (ligand, ligand_name) in enumerate(zip(ligands, ligand_names)):
            self.logger.info(f"Processing ligand {i+1}/{len(ligands)}: {ligand_name}")
            
            try:
                # Process single ligand
                ligand_results = self._process_single_ligand(protein, ligand, ligand_name)
                results[ligand_name] = ligand_results
                
                # Update progress
                if self.output_dir:
                    self._update_progress(i+1, len(ligands), ligand_name)
                    
            except Exception as e:
                self.logger.error(f"Error processing ligand {ligand_name}: {e}")
                results[ligand_name] = {
                    'error': str(e),
                    'poses': [],
                    'best_score': None,
                    'processing_time': 0.0
                }
        
        return results
    
    def _run_parallel_screening(self, protein: Any, ligands: List[Any], 
                              ligand_names: List[str]) -> Dict[str, Any]:
        """Run screening in parallel using multiple processes."""
        results = {}
        
        # Prepare tasks for parallel execution
        tasks = []
        for ligand, ligand_name in zip(ligands, ligand_names):
            task = {
                'protein': protein,
                'ligand': ligand,
                'ligand_name': ligand_name,
                'scoring_function_type': self.scoring_function_type,
                'algorithm_type': self.algorithm_type,
                'output_dir': self.output_dir / 'poses' / ligand_name if self.output_dir else None
            }
            tasks.append(task)
        
        # Process tasks in parallel
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all tasks
            future_to_name = {
                executor.submit(_process_ligand_worker, task): task['ligand_name']
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
                    self.logger.info(f"Completed {completed}/{len(ligands)}: {ligand_name}")
                    
                except Exception as e:
                    self.logger.error(f"Error processing ligand {ligand_name}: {e}")
                    results[ligand_name] = {
                        'error': str(e),
                        'poses': [],
                        'best_score': None,
                        'processing_time': 0.0
                    }
                
                # Update progress
                if self.output_dir:
                    self._update_progress(completed, len(ligands), ligand_name)
        
        return results
    
    def _process_single_ligand(self, protein: Any, ligand: Any, 
                             ligand_name: str) -> Dict[str, Any]:
        """Process a single ligand through docking."""
        start_time = time.time()
        
        try:
            # Create docking engine
            docking_engine = DockingEngine(
                performance_monitor=self.performance_monitor
            )
            
            # Create scoring function
            scoring_function = self.scoring_factory.create_scoring_function(
                self.scoring_function_type
            )
            
            # Create search algorithm
            algorithm = self.algorithm_factory.create_algorithm(
                self.algorithm_type,
                scoring_function,
                max_iterations=self.max_evals // 10  # Adjust iterations based on max_evals
            )
            
            # Run docking
            with self.performance_monitor.start_operation(f"docking_{ligand_name}"):
                docking_results = docking_engine.run_docking(
                    protein=protein,
                    ligand=ligand,
                    algorithm=algorithm,
                    config={
                        'max_poses': self.num_modes,
                        'rmsd_threshold': self.rmsd_threshold,
                        'local_opt': False,
                        'analyze_results': True
                    }
                )
            
            # Extract poses and scores
            poses = docking_results.get('poses', [])
            best_score = poses[0]['score'] if poses else None
            
            # Save results if output directory specified
            if self.output_dir:
                self._save_ligand_results(protein, ligand_name, poses)
            
            processing_time = time.time() - start_time
            
            return {
                'poses': poses,
                'best_score': best_score,
                'processing_time': processing_time,
                'statistics': docking_results.get('statistics', {}),
                'analysis': docking_results.get('analysis', {})
            }
            
        except Exception as e:
            processing_time = time.time() - start_time
            self.logger.error(f"Failed to process ligand {ligand_name}: {e}")
            
            return {
                'error': str(e),
                'poses': [],
                'best_score': None,
                'processing_time': processing_time
            }
    
    def _save_ligand_results(self, protein: Any, ligand_name: str, 
                           poses: List[Dict[str, Any]]) -> None:
        """Save docking results for a single ligand."""
        # Create ligand directory
        ligand_dir = self.output_dir / "poses" / ligand_name
        ligand_dir.mkdir(parents=True, exist_ok=True)
        
        # Save pose files and create complexes
        for i, pose_data in enumerate(poses[:10]):  # Save top 10 poses
            score = pose_data['score']
            pose_coords = pose_data.get('coordinates', [])
            atom_types = pose_data.get('atom_types', None)
            atom_names = pose_data.get('atom_names', None)
            
            # Save ligand pose
            pose_file = ligand_dir / f"pose_{i+1}_score_{score:.2f}.pdb"
            self._write_pose_pdb(pose_file, pose_coords, ligand_name, i+1, score, atom_types)
            
            # Save protein-ligand complex for top 3 poses
            if i < 3:
                complex_file = self.output_dir / "complexes" / f"{ligand_name}_pose_{i+1}_score_{score:.2f}.pdb"
                self._write_complex_pdb(complex_file, protein, pose_coords, ligand_name, score)
        
        # Save scores CSV
        scores_file = ligand_dir / "scores.csv"
        with open(scores_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Pose', 'Score', 'Rank'])
            for i, pose_data in enumerate(poses):
                writer.writerow([i+1, pose_data['score'], pose_data.get('rank', i+1)])
    
    def _write_pose_pdb(self, filename: Path, coordinates: List[List[float]], 
                       ligand_name: str, pose_num: int, score: float, 
                       atom_types: Optional[List[str]] = None) -> None:
        """Write ligand pose to PDB file using proper PDB writer."""
        try:
            coords_array = np.array(coordinates)
            PDBWriter.write_ligand_pose(
                filename=filename,
                coordinates=coords_array,
                ligand_name=ligand_name,
                score=score,
                pose_rank=pose_num,
                atom_types=atom_types
            )
        except Exception as e:
            self.logger.warning(f"Failed to write pose PDB {filename}: {e}")
    
    def _write_complex_pdb(self, filename: Path, protein: Any, ligand_coords: List[List[float]], 
                          ligand_name: str, score: float) -> None:
        """Write protein-ligand complex to PDB file using proper PDB writer."""
        try:
            if hasattr(protein, 'coords'):
                protein_coords = protein.coords
                ligand_coords_array = np.array(ligand_coords)
                
                # Prepare protein atom information if available
                protein_atoms = None
                if hasattr(protein, 'atoms'):
                    protein_atoms = protein.atoms
                
                PDBWriter.write_protein_ligand_complex(
                    filename=filename,
                    protein_coords=protein_coords,
                    ligand_coords=ligand_coords_array,
                    protein_atoms=protein_atoms,
                    ligand_name=ligand_name,
                    score=score
                )
            else:
                self.logger.warning(f"Protein missing coordinates for complex {filename}")
        except Exception as e:
            self.logger.warning(f"Failed to write complex PDB {filename}: {e}")
    
    def _update_progress(self, completed: int, total: int, current_ligand: str) -> None:
        """Update progress tracking file."""
        progress_file = self.output_dir / "progress.json"
        progress_data = {
            'completed': completed,
            'total': total,
            'progress_percent': (completed / total) * 100,
            'current_ligand': current_ligand,
            'timestamp': datetime.now().isoformat()
        }
        
        with open(progress_file, 'w') as f:
            json.dump(progress_data, f, indent=2)
    
    def _generate_summary_report(self, results: Dict[str, Any], elapsed_time: float) -> None:
        """Generate comprehensive summary report."""
        # Sort ligands by best score
        valid_results = {name: data for name, data in results.items() 
                        if data.get('best_score') is not None}
        
        sorted_ligands = sorted(
            valid_results.items(),
            key=lambda x: x[1]['best_score']
        )
        
        # Create CSV summary
        csv_file = self.output_dir / "reports" / "screening_results.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Rank', 'Ligand', 'Best_Score', 'Processing_Time_s', 'Status'])
            
            for i, (name, data) in enumerate(sorted_ligands):
                writer.writerow([
                    i+1, name, data['best_score'], 
                    data.get('processing_time', 0.0), 'Success'
                ])
            
            # Add failed ligands
            for name, data in results.items():
                if data.get('error'):
                    writer.writerow([
                        '-', name, 'N/A', data.get('processing_time', 0.0), 
                        f"Error: {data['error']}"
                    ])
        
        # Create detailed text report
        self._create_text_report(results, sorted_ligands, elapsed_time)
        
        # Generate plots if possible
        try:
            self._generate_plots(sorted_ligands)
        except ImportError:
            self.logger.warning("Matplotlib not available. Skipping plots.")
    
    def _create_text_report(self, results: Dict[str, Any], sorted_ligands: List[Tuple[str, Dict]], 
                           elapsed_time: float) -> None:
        """Create detailed text report."""
        report_file = self.output_dir / "reports" / "screening_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("         PandaDock Virtual Screening Results\n")
            f.write("=" * 60 + "\n\n")
            
            # Screening information
            f.write("SCREENING INFORMATION\n")
            f.write("-" * 20 + "\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Ligands: {len(results)}\n")
            f.write(f"Successful: {len(sorted_ligands)}\n")
            f.write(f"Failed: {len(results) - len(sorted_ligands)}\n")
            f.write(f"Workers: {self.n_workers}\n")
            f.write(f"Algorithm: {self.algorithm_type}\n")
            f.write(f"Scoring Function: {self.scoring_function_type}\n")
            f.write(f"Total Runtime: {elapsed_time:.2f} seconds\n")
            f.write(f"Average Time per Ligand: {elapsed_time/len(results):.2f} seconds\n\n")
            
            # Top ligands
            f.write("TOP 20 LIGANDS\n")
            f.write("-" * 14 + "\n")
            f.write("Rank  Ligand                   Best Score    Time (s)\n")
            f.write("-" * 55 + "\n")
            
            for i, (name, data) in enumerate(sorted_ligands[:20]):
                f.write(f"{i+1:4d}  {name:22s}   {data['best_score']:8.3f}    {data.get('processing_time', 0):.1f}\n")
            
            f.write(f"\nFull results available in screening_results.csv\n")
            f.write("=" * 60 + "\n")
    
    def _generate_plots(self, sorted_ligands: List[Tuple[str, Dict]]) -> None:
        """Generate summary plots."""
        import matplotlib.pyplot as plt
        
        plots_dir = self.output_dir / "reports" / "plots"
        plots_dir.mkdir(exist_ok=True)
        
        if not sorted_ligands:
            return
        
        scores = [data['best_score'] for _, data in sorted_ligands]
        names = [name for name, _ in sorted_ligands]
        
        # Score distribution histogram
        plt.figure(figsize=(10, 6))
        plt.hist(scores, bins=min(20, len(scores)//3 + 1), alpha=0.7, color='skyblue')
        plt.axvline(np.mean(scores), color='red', linestyle='--', label=f'Mean: {np.mean(scores):.2f}')
        plt.axvline(np.median(scores), color='green', linestyle=':', label=f'Median: {np.median(scores):.2f}')
        plt.xlabel('Docking Score')
        plt.ylabel('Frequency')
        plt.title('Distribution of Docking Scores')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(plots_dir / "score_distribution.png", dpi=300)
        plt.close()
        
        # Top compounds bar chart
        top_n = min(20, len(sorted_ligands))
        plt.figure(figsize=(12, 8))
        top_scores = scores[:top_n]
        top_names = names[:top_n]
        
        bars = plt.bar(range(top_n), top_scores, color='lightcoral')
        plt.xticks(range(top_n), top_names, rotation=45, ha='right')
        plt.xlabel('Ligand')
        plt.ylabel('Docking Score')
        plt.title(f'Top {top_n} Ligands by Docking Score')
        plt.grid(True, axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(plots_dir / "top_ligands.png", dpi=300)
        plt.close()
    
    def _calculate_screening_statistics(self, results: Dict[str, Any], 
                                      elapsed_time: float) -> Dict[str, Any]:
        """Calculate screening statistics."""
        total_ligands = len(results)
        successful = len([r for r in results.values() if r.get('best_score') is not None])
        failed = total_ligands - successful
        
        if successful > 0:
            scores = [r['best_score'] for r in results.values() if r.get('best_score') is not None]
            processing_times = [r.get('processing_time', 0) for r in results.values()]
            
            statistics = {
                'total_ligands': total_ligands,
                'successful': successful,
                'failed': failed,
                'success_rate': successful / total_ligands,
                'total_runtime': elapsed_time,
                'average_runtime_per_ligand': elapsed_time / total_ligands,
                'best_score': min(scores),
                'worst_score': max(scores),
                'mean_score': np.mean(scores),
                'median_score': np.median(scores),
                'std_score': np.std(scores),
                'total_processing_time': sum(processing_times),
                'average_processing_time': np.mean(processing_times)
            }
        else:
            statistics = {
                'total_ligands': total_ligands,
                'successful': 0,
                'failed': total_ligands,
                'success_rate': 0.0,
                'total_runtime': elapsed_time,
                'message': 'No ligands were successfully processed'
            }
        
        return statistics


def _process_ligand_worker(task: Dict[str, Any]) -> Dict[str, Any]:
    """Worker function for parallel ligand processing."""
    try:
        # Extract task parameters
        protein = task['protein']
        ligand = task['ligand']
        ligand_name = task['ligand_name']
        scoring_function_type = task['scoring_function_type']
        algorithm_type = task['algorithm_type']
        output_dir = task.get('output_dir')
        
        # Create screening manager for this worker
        worker_manager = VirtualScreeningManager(
            scoring_function_type=scoring_function_type,
            algorithm_type=algorithm_type,
            output_dir=output_dir,
            n_workers=1  # Single worker mode
        )
        
        # Process the ligand
        result = worker_manager._process_single_ligand(protein, ligand, ligand_name)
        return result
        
    except Exception as e:
        return {
            'error': str(e),
            'poses': [],
            'best_score': None,
            'processing_time': 0.0
        }