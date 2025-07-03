"""
Result processing for molecular docking.

This module handles post-processing of docking results including clustering,
analysis, optimization, and validation.
"""

import logging
import numpy as np
from typing import List, Tuple, Dict, Any, Optional
import time

from ..hardware import PerformanceMonitor


class ResultProcessor:
    """
    Processes and analyzes docking results.
    
    Handles clustering, local optimization, energy analysis,
    and result validation.
    """
    
    def __init__(self, performance_monitor: PerformanceMonitor = None):
        """
        Initialize result processor.
        
        Args:
            performance_monitor: Performance monitoring system (optional)
        """
        self.performance_monitor = performance_monitor
        self.logger = logging.getLogger(__name__)
    
    def process_results(self, raw_results: List[Tuple[Any, float]], 
                       protein: Any, ligand: Any, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process raw docking results.
        
        Args:
            raw_results: List of (pose, score) tuples from search
            protein: Protein molecule object
            ligand: Ligand molecule object
            config: Configuration dictionary
            
        Returns:
            Dictionary with processed results
        """
        self.logger.info(f"Processing {len(raw_results)} raw results")
        
        if not raw_results:
            return {
                'poses': [],
                'statistics': {'total_poses': 0},
                'analysis': {'message': 'No poses generated'}
            }
        
        processed_results = {
            'poses': [],
            'statistics': {},
            'analysis': {}
        }
        
        # Sort results by score (best first)
        sorted_results = sorted(raw_results, key=lambda x: x[1])
        
        # Apply local optimization if requested
        if config.get('local_opt', False):
            sorted_results = self._apply_local_optimization(
                sorted_results, protein, config
            )
        
        # Remove duplicates and limit results
        unique_results = self._remove_duplicates(sorted_results)
        final_results = unique_results[:config.get('max_poses', 20)]
        
        # Convert to structured format
        processed_results['poses'] = self._format_poses(final_results, ligand)
        
        # Generate statistics
        processed_results['statistics'] = self._calculate_statistics(final_results)
        
        # Perform analysis if requested
        if config.get('analyze_results', True):
            processed_results['analysis'] = self._analyze_results(
                final_results, protein, ligand, config
            )
        
        self.logger.info(f"Processed results: {len(processed_results['poses'])} final poses")
        
        return processed_results
    
    def _apply_local_optimization(self, results: List[Tuple[Any, float]], 
                                 protein: Any, config: Dict[str, Any]) -> List[Tuple[Any, float]]:
        """Apply local optimization to top poses."""
        
        self.logger.info("Applying local optimization to top poses")
        
        optimized_results = []
        max_optimize = min(config.get('max_optimize', 10), len(results))
        
        with self.performance_monitor.start_operation("local_optimization"):
            for i, (pose, score) in enumerate(results):
                if i < max_optimize:
                    # Apply optimization
                    try:
                        optimized_pose = self._optimize_single_pose(pose, protein, config)
                        # Re-score optimized pose
                        optimized_score = self._score_pose(optimized_pose, protein, config)
                        optimized_results.append((optimized_pose, optimized_score))
                        
                        self.logger.debug(f"Optimized pose {i+1}: {score:.3f} -> {optimized_score:.3f}")
                        
                    except Exception as e:
                        self.logger.warning(f"Optimization failed for pose {i+1}: {e}")
                        optimized_results.append((pose, score))
                else:
                    # Keep remaining poses as-is
                    optimized_results.append((pose, score))
        
        return optimized_results
    
    def _optimize_single_pose(self, pose: Any, protein: Any, config: Dict[str, Any]) -> Any:
        """Optimize a single pose using local optimization."""
        
        # Simple gradient-based optimization
        try:
            # Extract coordinates
            if hasattr(pose, 'coords'):
                coords = pose.coords.copy()
            else:
                coords = np.array(pose)
            
            # Apply small random perturbations and keep best
            best_coords = coords.copy()
            best_score = self._score_pose_coords(coords, protein, config)
            
            for _ in range(10):  # Limited optimization steps
                # Small random perturbation
                perturbed_coords = coords + np.random.normal(0, 0.1, coords.shape)
                score = self._score_pose_coords(perturbed_coords, protein, config)
                
                if score < best_score:
                    best_coords = perturbed_coords
                    best_score = score
            
            # Create optimized pose object
            if hasattr(pose, 'coords'):
                optimized_pose = pose.__class__(pose)
                optimized_pose.coords = best_coords
                return optimized_pose
            else:
                return best_coords
                
        except Exception as e:
            self.logger.warning(f"Single pose optimization failed: {e}")
            return pose
    
    def _score_pose(self, pose: Any, protein: Any, config: Dict[str, Any]) -> float:
        """Score a pose using the configured scoring function."""
        try:
            # Simple distance-based scoring as fallback
            return self._score_pose_coords(
                pose.coords if hasattr(pose, 'coords') else pose,
                protein, config
            )
        except Exception:
            return 0.0
    
    def _score_pose_coords(self, coords: np.ndarray, protein: Any, config: Dict[str, Any]) -> float:
        """Score pose coordinates."""
        try:
            # Simple distance-based scoring
            protein_coords = protein.coords if hasattr(protein, 'coords') else np.array([[0, 0, 0]])
            
            # Calculate minimum distance
            distances = np.linalg.norm(
                protein_coords[:, np.newaxis, :] - coords[np.newaxis, :, :],
                axis=2
            )
            min_distance = np.min(distances)
            
            # Simple scoring function
            if min_distance < 1.0:
                return 100.0  # Clash penalty
            elif min_distance > 20.0:
                return 25.0   # Too far penalty (reduced from 50.0)
            elif min_distance > 15.0:
                return 10.0   # Moderately far penalty
            else:
                return -min_distance  # Better when closer
                
        except Exception:
            return 0.0
    
    def _remove_duplicates(self, results: List[Tuple[Any, float]], 
                          rmsd_threshold: float = 1.0) -> List[Tuple[Any, float]]:
        """Remove duplicate poses based on RMSD."""
        
        if len(results) <= 1:
            return results
        
        unique_results = [results[0]]  # Keep first pose
        
        for pose, score in results[1:]:
            is_duplicate = False
            
            for unique_pose, _ in unique_results:
                try:
                    rmsd = self._calculate_rmsd(pose, unique_pose)
                    if rmsd < rmsd_threshold:
                        is_duplicate = True
                        break
                except Exception:
                    # If RMSD calculation fails, assume not duplicate
                    continue
            
            if not is_duplicate:
                unique_results.append((pose, score))
        
        self.logger.info(f"Removed {len(results) - len(unique_results)} duplicate poses")
        return unique_results
    
    def _calculate_rmsd(self, pose1: Any, pose2: Any) -> float:
        """Calculate RMSD between two poses."""
        try:
            coords1 = pose1.coords if hasattr(pose1, 'coords') else np.array(pose1)
            coords2 = pose2.coords if hasattr(pose2, 'coords') else np.array(pose2)
            
            if coords1.shape != coords2.shape:
                return float('inf')
            
            diff = coords1 - coords2
            msd = np.mean(np.sum(diff**2, axis=1))
            return np.sqrt(msd)
            
        except Exception:
            return float('inf')
    
    def _format_poses(self, results: List[Tuple[Any, float]], ligand: Any = None) -> List[Dict[str, Any]]:
        """Format poses into structured dictionary format."""
        
        formatted_poses = []
        
        for i, (pose, score) in enumerate(results):
            # Extract energy components if available
            energy_components = {}
            if hasattr(pose, 'energy_components') and pose.energy_components:
                energy_components = pose.energy_components
            
            pose_dict = {
                'rank': i + 1,
                'score': float(score),
                'pose_data': pose,
                'energy_components': energy_components
            }
            
            # Extract coordinates if available
            if hasattr(pose, 'coords'):
                pose_dict['coordinates'] = pose.coords.tolist()
                pose_dict['n_atoms'] = len(pose.coords)
            elif isinstance(pose, np.ndarray):
                pose_dict['coordinates'] = pose.tolist()
                pose_dict['n_atoms'] = len(pose)
            
            # Add ligand atom information for proper PDB writing
            if hasattr(ligand, 'atom_types'):
                pose_dict['atom_types'] = ligand.atom_types
            if hasattr(ligand, 'atom_names'):
                pose_dict['atom_names'] = ligand.atom_names
            
            formatted_poses.append(pose_dict)
        
        return formatted_poses
    
    def _calculate_statistics(self, results: List[Tuple[Any, float]]) -> Dict[str, Any]:
        """Calculate statistics for the results."""
        
        if not results:
            return {'total_poses': 0}
        
        scores = [score for _, score in results]
        
        statistics = {
            'total_poses': len(results),
            'best_score': min(scores),
            'worst_score': max(scores),
            'mean_score': np.mean(scores),
            'std_score': np.std(scores),
            'score_range': max(scores) - min(scores)
        }
        
        return statistics
    
    def _analyze_results(self, results: List[Tuple[Any, float]], 
                        protein: Any, ligand: Any, config: Dict[str, Any]) -> Dict[str, Any]:
        """Perform analysis on the results."""
        
        analysis = {
            'clustering': {},
            'interactions': {},
            'binding_modes': {},
            'energy_analysis': {}
        }
        
        # Basic clustering analysis
        if config.get('cluster_poses', False):
            analysis['clustering'] = self._cluster_poses(results, config)
        
        # Interaction analysis for top poses
        if config.get('analyze_interactions', False):
            analysis['interactions'] = self._analyze_interactions(
                results[:5], protein, config
            )
        
        # Energy component analysis
        if config.get('detailed_energy', False):
            analysis['energy_analysis'] = self._analyze_energy_components(
                results[:3], protein, config
            )
        
        return analysis
    
    def _cluster_poses(self, results: List[Tuple[Any, float]], 
                      config: Dict[str, Any]) -> Dict[str, Any]:
        """Perform pose clustering."""
        
        if len(results) < 2:
            return {'clusters': [], 'n_clusters': 0}
        
        # Simple clustering based on RMSD
        rmsd_cutoff = config.get('rmsd_cutoff', 2.0)
        clusters = []
        unclustered = list(results)
        
        while unclustered:
            # Start new cluster with first unclustered pose
            cluster_center = unclustered.pop(0)
            cluster = [cluster_center]
            
            # Find poses similar to cluster center
            remaining = []
            for pose, score in unclustered:
                try:
                    rmsd = self._calculate_rmsd(cluster_center[0], pose)
                    if rmsd < rmsd_cutoff:
                        cluster.append((pose, score))
                    else:
                        remaining.append((pose, score))
                except Exception:
                    remaining.append((pose, score))
            
            unclustered = remaining
            clusters.append({
                'center': cluster_center,
                'members': cluster,
                'size': len(cluster),
                'best_score': min(score for _, score in cluster)
            })
        
        return {
            'clusters': clusters,
            'n_clusters': len(clusters),
            'largest_cluster_size': max(len(c['members']) for c in clusters) if clusters else 0
        }
    
    def _analyze_interactions(self, top_results: List[Tuple[Any, float]], 
                            protein: Any, config: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze protein-ligand interactions for top poses."""
        
        interactions = {
            'hydrogen_bonds': [],
            'hydrophobic_contacts': [],
            'electrostatic_interactions': []
        }
        
        for i, (pose, score) in enumerate(top_results):
            pose_interactions = {
                'pose_rank': i + 1,
                'score': score,
                'interactions': self._find_interactions(pose, protein)
            }
            
            # Categorize interactions
            for interaction_type in interactions:
                interactions[interaction_type].append(pose_interactions)
        
        return interactions
    
    def _find_interactions(self, pose: Any, protein: Any) -> List[str]:
        """Find interactions between pose and protein."""
        
        # Simplified interaction detection
        interactions = []
        
        try:
            pose_coords = pose.coords if hasattr(pose, 'coords') else np.array(pose)
            protein_coords = protein.coords if hasattr(protein, 'coords') else np.array([[0, 0, 0]])
            
            # Calculate distances
            distances = np.linalg.norm(
                protein_coords[:, np.newaxis, :] - pose_coords[np.newaxis, :, :],
                axis=2
            )
            
            # Find close contacts
            close_contacts = np.where(distances < 4.0)
            
            for prot_idx, lig_idx in zip(close_contacts[0], close_contacts[1]):
                distance = distances[prot_idx, lig_idx]
                interactions.append(f"Contact: protein_{prot_idx} - ligand_{lig_idx} ({distance:.2f}Å)")
        
        except Exception as e:
            interactions.append(f"Interaction analysis failed: {e}")
        
        return interactions
    
    def _analyze_energy_components(self, top_results: List[Tuple[Any, float]], 
                                  protein: Any, config: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze energy components for top poses."""
        
        energy_analysis = {
            'per_pose_breakdown': [],
            'average_components': {},
            'component_trends': {}
        }
        
        for i, (pose, total_score) in enumerate(top_results):
            # Simplified energy component estimation
            components = {
                'van_der_waals': total_score * 0.4,
                'electrostatic': total_score * 0.3,
                'hydrogen_bonds': total_score * 0.2,
                'solvation': total_score * 0.1
            }
            
            energy_analysis['per_pose_breakdown'].append({
                'pose_rank': i + 1,
                'total_score': total_score,
                'components': components
            })
        
        return energy_analysis