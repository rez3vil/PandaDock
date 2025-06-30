"""
Result writing utilities for PandaDock.

This module handles writing docking results to various file formats.
"""

import json
import csv
import logging
import numpy as np
from typing import Dict, List, Any, Optional
from pathlib import Path
from datetime import datetime

from .pdb_writer import PDBWriter


class ResultWriters:
    """Handles writing docking results to various file formats."""
    
    def __init__(self, output_dir: str):
        """
        Initialize result writers.
        
        Args:
            output_dir: Output directory for results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
    
    def save_all_results(self, results: Dict[str, Any], config: Dict[str, Any], protein=None) -> Dict[str, str]:
        """
        Save results in all available formats.
        
        Args:
            results: Processed docking results
            config: Configuration parameters
            protein: Protein structure object (optional)
            
        Returns:
            Dictionary mapping format names to file paths
        """
        saved_files = {}
        
        # Save poses in PDB format
        poses_dir = self.output_dir / "poses"
        poses_dir.mkdir(exist_ok=True)
        
        if results.get('poses'):
            for i, pose_data in enumerate(results['poses']):
                pose_file = poses_dir / f"pose_{i+1}_score_{pose_data['score']:.2f}.pdb"
                self._save_pose_pdb(pose_data, pose_file)
                if i == 0:  # Track first pose file as example
                    saved_files['best_pose'] = str(pose_file)
        
        # Save summary in JSON format
        json_file = self.output_dir / "results_summary.json"
        saved_files['json_summary'] = str(json_file)
        self._save_json_summary(results, config, json_file)
        
        # Save detailed CSV results
        csv_file = self.output_dir / "results_detailed.csv"
        saved_files['csv_detailed'] = str(csv_file)
        self._save_csv_results(results, csv_file)
        
        # Save text report
        text_file = self.output_dir / "docking_report.txt"
        saved_files['text_report'] = str(text_file)
        self._save_text_report(results, config, text_file)
        
        # Generate plots
        plots_dir = self.output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        plot_files = self._generate_plots(results, plots_dir)
        saved_files.update(plot_files)
        
        # Save binding site sphere if protein provided
        if protein and hasattr(protein, 'active_site') and protein.active_site:
            sphere_file = self.output_dir / "sphere.pdb"
            saved_files['sphere'] = str(sphere_file)
            self._save_binding_site_sphere(protein, sphere_file)
        
        # Save configuration file
        config_file = self.output_dir / "docking_config.json"
        saved_files['config'] = str(config_file)
        self._save_config(config, config_file)
        
        self.logger.info(f"Saved results to {len(saved_files)} files in {self.output_dir}")
        
        return saved_files
    
    def _save_binding_site_sphere(self, protein, output_file: Path) -> None:
        """Save binding site sphere for visualization."""
        try:
            from ..molecules import ProteinHandler
            
            handler = ProteinHandler()
            handler.save_binding_site_sphere(protein, output_file)
            
        except Exception as e:
            self.logger.warning(f"Failed to save binding site sphere: {e}")
    
    def _generate_plots(self, results: Dict[str, Any], plots_dir: Path) -> Dict[str, str]:
        """Generate plots for docking results."""
        plot_files = {}
        
        try:
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import matplotlib.pyplot as plt
            
            poses = results.get('poses', [])
            if not poses:
                return plot_files
            
            scores = [pose['score'] for pose in poses]
            
            # 1. Score distribution histogram
            plt.figure(figsize=(10, 6))
            plt.hist(scores, bins=min(20, len(scores)//2), alpha=0.7, edgecolor='black')
            plt.title('Docking Score Distribution', fontsize=14, fontweight='bold')
            plt.xlabel('Docking Score')
            plt.ylabel('Frequency')
            plt.grid(True, alpha=0.3)
            
            score_dist_file = plots_dir / 'score_distribution.png'
            plt.savefig(score_dist_file, dpi=300, bbox_inches='tight')
            plt.close()
            plot_files['score_distribution'] = str(score_dist_file)
            
            # 2. Score vs Rank plot
            plt.figure(figsize=(10, 6))
            ranks = [pose['rank'] for pose in poses]
            plt.plot(ranks, scores, 'o-', linewidth=2, markersize=4)
            plt.title('Docking Scores by Rank', fontsize=14, fontweight='bold')
            plt.xlabel('Pose Rank')
            plt.ylabel('Docking Score')
            plt.grid(True, alpha=0.3)
            
            # Highlight best poses
            if len(scores) > 0:
                best_score_idx = np.argmin(scores)
                plt.plot(ranks[best_score_idx], scores[best_score_idx], 
                        'ro', markersize=10, label=f'Best Score: {scores[best_score_idx]:.3f}')
                plt.legend()
            
            score_rank_file = plots_dir / 'score_vs_rank.png'
            plt.savefig(score_rank_file, dpi=300, bbox_inches='tight')
            plt.close()
            plot_files['score_vs_rank'] = str(score_rank_file)
            
            self.logger.info(f"Generated {len(plot_files)} plots")
            
        except ImportError:
            self.logger.warning("Matplotlib not available, skipping plot generation")
        except Exception as e:
            self.logger.warning(f"Failed to generate plots: {e}")
        
        return plot_files
    
    def _save_pose_pdb(self, pose_data: Dict[str, Any], output_file: Path) -> None:
        """Save a single pose in PDB format using proper PDB writer."""
        
        try:
            if 'coordinates' in pose_data:
                coords = np.array(pose_data['coordinates'])
                
                # Extract atom types if available
                atom_types = pose_data.get('atom_types', None)
                atom_names = pose_data.get('atom_names', None)
                
                # Use the proper PDB writer
                PDBWriter.write_ligand_pose(
                    filename=output_file,
                    coordinates=coords,
                    ligand_name="LIG",
                    score=pose_data.get('score'),
                    pose_rank=pose_data.get('rank'),
                    atom_types=atom_types
                )
                
                self.logger.debug(f"Saved pose PDB: {output_file}")
            else:
                self.logger.warning(f"No coordinates found in pose data for {output_file}")
                
        except Exception as e:
            self.logger.warning(f"Failed to save pose PDB {output_file}: {e}")
    
    def _save_json_summary(self, results: Dict[str, Any], config: Dict[str, Any], 
                          output_file: Path) -> None:
        """Save results summary in JSON format."""
        
        try:
            summary = {
                'metadata': {
                    'generated_by': 'PandaDock',
                    'timestamp': datetime.now().isoformat(),
                    'version': '2.0.0-refactored'
                },
                'configuration': self._serialize_config(config),
                'statistics': results.get('statistics', {}),
                'analysis': results.get('analysis', {}),
                'poses': []
            }
            
            # Add pose summaries (without full coordinate data)
            if results.get('poses'):
                for pose in results['poses']:
                    pose_summary = {
                        'rank': pose['rank'],
                        'score': pose['score'],
                        'n_atoms': pose.get('n_atoms', 0),
                        'energy_components': pose.get('energy_components', {})
                    }
                    summary['poses'].append(pose_summary)
            
            with open(output_file, 'w') as f:
                json.dump(summary, f, indent=2, default=self._json_serializer)
                
        except Exception as e:
            self.logger.warning(f"Failed to save JSON summary {output_file}: {e}")
    
    def _save_csv_results(self, results: Dict[str, Any], output_file: Path) -> None:
        """Save detailed results in CSV format."""
        
        try:
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                
                # Header
                writer.writerow([
                    'Rank', 'Score', 'N_Atoms', 'VdW_Energy', 'Electrostatic_Energy',
                    'HBond_Energy', 'Solvation_Energy', 'Center_X', 'Center_Y', 'Center_Z'
                ])
                
                # Data rows
                if results.get('poses'):
                    for pose in results['poses']:
                        energy_components = pose.get('energy_components', {})
                        
                        # Calculate center of mass if coordinates available
                        center = [0.0, 0.0, 0.0]
                        if 'coordinates' in pose:
                            coords = np.array(pose['coordinates'])
                            center = np.mean(coords, axis=0).tolist()
                        
                        writer.writerow([
                            pose['rank'],
                            pose['score'],
                            pose.get('n_atoms', 0),
                            energy_components.get('van_der_waals', 0.0),
                            energy_components.get('electrostatic', 0.0),
                            energy_components.get('hydrogen_bonds', 0.0),
                            energy_components.get('solvation', 0.0),
                            center[0], center[1], center[2]
                        ])
                        
        except Exception as e:
            self.logger.warning(f"Failed to save CSV results {output_file}: {e}")
    
    def _save_text_report(self, results: Dict[str, Any], config: Dict[str, Any], 
                         output_file: Path) -> None:
        """Save comprehensive text report."""
        
        try:
            with open(output_file, 'w') as f:
                # Header
                f.write("=" * 80 + "\n")
                f.write("  PandaDock Molecular Docking Results\n")
                f.write("  Refactored Architecture - Easy to Debug & Extend\n")
                f.write("=" * 80 + "\n\n")
                
                # Metadata
                f.write("RUN INFORMATION\n")
                f.write("-" * 20 + "\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Algorithm: {config.get('algorithm', 'Unknown')}\n")
                f.write(f"Scoring: {self._get_scoring_description(config)}\n")
                f.write(f"Hardware: {config.get('device_used', 'Unknown')}\n\n")
                
                # Statistics
                stats = results.get('statistics', {})
                if stats:
                    f.write("RESULTS STATISTICS\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"Total poses: {stats.get('total_poses', 0)}\n")
                    if stats.get('best_score') is not None:
                        f.write(f"Best score: {stats['best_score']:.4f}\n")
                        f.write(f"Worst score: {stats.get('worst_score', 0):.4f}\n")
                        f.write(f"Mean score: {stats.get('mean_score', 0):.4f}\n")
                        f.write(f"Score std: {stats.get('std_score', 0):.4f}\n")
                    f.write("\n")
                
                # Top poses
                poses = results.get('poses', [])
                if poses:
                    f.write("TOP POSES\n")
                    f.write("-" * 20 + "\n")
                    f.write(f"{'Rank':<6}{'Score':<12}{'N_Atoms':<10}{'File'}\n")
                    f.write("-" * 50 + "\n")
                    
                    for pose in poses[:10]:  # Top 10
                        filename = f"pose_{pose['rank']}_score_{pose['score']:.2f}.pdb"
                        f.write(
                            f"{pose['rank']:<6}{pose['score']:<12.4f}"
                            f"{pose.get('n_atoms', 0):<10}{filename}\n"
                        )
                    f.write("\n")
                
                # Analysis results
                analysis = results.get('analysis', {})
                if analysis:
                    f.write("ANALYSIS RESULTS\n")
                    f.write("-" * 20 + "\n")
                    
                    # Clustering
                    clustering = analysis.get('clustering', {})
                    if clustering:
                        f.write(f"Pose clusters: {clustering.get('n_clusters', 0)}\n")
                        f.write(f"Largest cluster: {clustering.get('largest_cluster_size', 0)} poses\n")
                    
                    # Interactions
                    interactions = analysis.get('interactions', {})
                    if interactions:
                        f.write("Key interactions detected for top poses\n")
                    
                    f.write("\n")
                
                # Footer
                f.write("=" * 80 + "\n")
                f.write("For detailed pose structures, see the poses/ directory\n")
                f.write("For machine-readable data, see results_summary.json\n")
                f.write("=" * 80 + "\n")
                
        except Exception as e:
            self.logger.warning(f"Failed to save text report {output_file}: {e}")
    
    def _save_config(self, config: Dict[str, Any], output_file: Path) -> None:
        """Save configuration parameters."""
        
        try:
            serializable_config = self._serialize_config(config)
            
            with open(output_file, 'w') as f:
                json.dump(serializable_config, f, indent=2, default=self._json_serializer)
                
        except Exception as e:
            self.logger.warning(f"Failed to save config {output_file}: {e}")
    
    def _serialize_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Convert config to JSON-serializable format."""
        
        serializable = {}
        
        for key, value in config.items():
            if isinstance(value, (str, int, float, bool, list, dict)):
                serializable[key] = value
            elif isinstance(value, np.ndarray):
                serializable[key] = value.tolist()
            elif value is None:
                serializable[key] = None
            else:
                serializable[key] = str(value)
        
        return serializable
    
    def _json_serializer(self, obj):
        """Custom JSON serializer for numpy types."""
        
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        else:
            return str(obj)
    
    def _get_scoring_description(self, config: Dict[str, Any]) -> str:
        """Get human-readable scoring description."""
        
        if config.get('physics_based', False):
            return "Physics-based (MM-GBSA)"
        elif config.get('enhanced_scoring', False):
            return "Enhanced (with electrostatics)"
        else:
            return "Standard composite"
    
    def save_complex_structures(self, results: Dict[str, Any], protein, 
                              output_dir: Optional[Path] = None) -> List[str]:
        """
        Save protein-ligand complex structures.
        
        Args:
            results: Docking results
            protein: Protein structure object
            output_dir: Optional custom output directory
            
        Returns:
            List of saved complex file paths
        """
        if output_dir is None:
            output_dir = self.output_dir
        
        complex_dir = output_dir / "complexes"
        complex_dir.mkdir(exist_ok=True)
        
        saved_files = []
        
        poses = results.get('poses', [])
        for pose_data in poses[:5]:  # Save top 5 complexes
            complex_file = complex_dir / f"complex_rank_{pose_data['rank']}.pdb"
            
            try:
                # Extract coordinates and atom information
                if hasattr(protein, 'coords') and 'coordinates' in pose_data:
                    protein_coords = protein.coords
                    ligand_coords = np.array(pose_data['coordinates'])
                    
                    # Prepare protein atom information if available
                    protein_atoms = None
                    if hasattr(protein, 'atoms'):
                        protein_atoms = protein.atoms
                    
                    # Use proper PDB writer for complex
                    PDBWriter.write_protein_ligand_complex(
                        filename=complex_file,
                        protein_coords=protein_coords,
                        ligand_coords=ligand_coords,
                        protein_atoms=protein_atoms,
                        ligand_name="LIG",
                        score=pose_data.get('score')
                    )
                    
                    saved_files.append(str(complex_file))
                else:
                    self.logger.warning(f"Missing coordinates for complex {complex_file}")
                    
                
            except Exception as e:
                self.logger.warning(f"Failed to save complex {complex_file}: {e}")
        
        if saved_files:
            self.logger.info(f"Saved {len(saved_files)} protein-ligand complexes")
        
        return saved_files