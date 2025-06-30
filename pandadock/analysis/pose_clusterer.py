"""
Pose clustering analysis for molecular docking results.

This module provides hierarchical and density-based clustering methods
for analyzing docking poses based on RMSD and other structural features.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Tuple
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from ..utils_new import calculate_rmsd


class PoseClusterer:
    """Clustering of docking poses using various methods."""
    
    def __init__(self, method: str = 'hierarchical', rmsd_cutoff: float = 2.0, 
                 min_cluster_size: int = 3):
        """
        Initialize pose clusterer.
        
        Args:
            method: Clustering method ('hierarchical' or 'dbscan')
            rmsd_cutoff: RMSD cutoff for clustering (Angstroms)
            min_cluster_size: Minimum number of poses for a valid cluster
        """
        self.method = method
        self.rmsd_cutoff = rmsd_cutoff
        self.min_cluster_size = min_cluster_size
        self.logger = logging.getLogger(__name__)
    
    def cluster_poses(self, poses: List[Any], scores: Optional[List[float]] = None) -> Dict[str, Any]:
        """
        Cluster ligand poses based on RMSD.
        
        Args:
            poses: List of pose objects with coordinate data
            scores: Corresponding scores for each pose
            
        Returns:
            Clustering results with cluster assignments and representatives
        """
        if not poses:
            return {
                'clusters': [],
                'n_clusters': 0,
                'rmsd_cutoff': self.rmsd_cutoff,
                'rmsd_matrix': np.array([])
            }
        
        if scores is None:
            scores = [0.0] * len(poses)
        
        self.logger.info(f"Clustering {len(poses)} poses using {self.method} method")
        
        if self.method == 'hierarchical':
            return self._hierarchical_clustering(poses, scores)
        elif self.method == 'dbscan':
            return self._dbscan_clustering(poses, scores)
        else:
            raise ValueError(f"Clustering method '{self.method}' not implemented")
    
    def _hierarchical_clustering(self, poses: List[Any], scores: List[float]) -> Dict[str, Any]:
        """
        Perform hierarchical clustering based on RMSD.
        
        Args:
            poses: List of pose objects
            scores: Corresponding scores for each pose
            
        Returns:
            Clustering results
        """
        n_poses = len(poses)
        self.logger.info(f"Computing RMSD matrix for {n_poses} poses...")
        
        # Compute pairwise RMSD matrix
        rmsd_matrix = np.zeros((n_poses, n_poses))
        
        for i in range(n_poses):
            for j in range(i+1, n_poses):
                try:
                    # Extract coordinates from poses
                    coords1 = self._extract_coordinates(poses[i])
                    coords2 = self._extract_coordinates(poses[j])
                    
                    if coords1.shape != coords2.shape:
                        rmsd = float('inf')  # Different sized molecules
                    else:
                        rmsd = calculate_rmsd(coords1, coords2)
                    
                    rmsd_matrix[i, j] = rmsd
                    rmsd_matrix[j, i] = rmsd
                    
                except Exception as e:
                    self.logger.warning(f"RMSD calculation failed for poses {i}, {j}: {e}")
                    rmsd_matrix[i, j] = float('inf')
                    rmsd_matrix[j, i] = float('inf')
        
        # Convert distance matrix to condensed form for scipy
        condensed_matrix = squareform(rmsd_matrix)
        
        # Perform hierarchical clustering
        self.logger.info("Performing hierarchical clustering...")
        try:
            Z = linkage(condensed_matrix, method='average')
            
            # Form flat clusters based on RMSD cutoff
            cluster_indices = fcluster(Z, self.rmsd_cutoff, criterion='distance')
            
        except Exception as e:
            self.logger.error(f"Hierarchical clustering failed: {e}")
            # Fallback: each pose is its own cluster
            cluster_indices = list(range(1, n_poses + 1))
        
        # Process clusters
        clusters = {}
        for i, cluster_idx in enumerate(cluster_indices):
            if cluster_idx not in clusters:
                clusters[cluster_idx] = []
            clusters[cluster_idx].append({
                'pose_idx': i,
                'pose': poses[i],
                'score': scores[i]
            })
        
        # Format results
        result_clusters = []
        for cluster_idx, members in clusters.items():
            # Sort cluster members by score (best first)
            members.sort(key=lambda x: x['score'])
            
            # Skip small clusters
            if len(members) < self.min_cluster_size:
                continue
                
            # Add cluster to results
            result_clusters.append({
                'cluster_id': cluster_idx,
                'members': members,
                'size': len(members),
                'representative': members[0]['pose_idx'],
                'best_score': members[0]['score'],
                'worst_score': members[-1]['score'],
                'avg_score': np.mean([m['score'] for m in members])
            })
        
        # Sort clusters by size (largest first)
        result_clusters.sort(key=lambda x: x['size'], reverse=True)
        
        self.logger.info(f"Found {len(result_clusters)} clusters (min size: {self.min_cluster_size})")
        
        return {
            'clusters': result_clusters,
            'n_clusters': len(result_clusters),
            'rmsd_cutoff': self.rmsd_cutoff,
            'rmsd_matrix': rmsd_matrix,
            'method': 'hierarchical'
        }
    
    def _dbscan_clustering(self, poses: List[Any], scores: List[float]) -> Dict[str, Any]:
        """
        Perform DBSCAN clustering based on RMSD.
        
        Args:
            poses: List of pose objects
            scores: Corresponding scores
            
        Returns:
            Clustering results
        """
        try:
            from sklearn.cluster import DBSCAN
        except ImportError:
            self.logger.error("scikit-learn not available for DBSCAN clustering")
            # Fallback to hierarchical clustering
            return self._hierarchical_clustering(poses, scores)
        
        n_poses = len(poses)
        self.logger.info(f"Computing coordinate matrix for {n_poses} poses...")
        
        # Extract coordinates and flatten for DBSCAN
        coordinate_matrix = []
        valid_indices = []
        
        for i, pose in enumerate(poses):
            try:
                coords = self._extract_coordinates(pose)
                if coords.size > 0:
                    coordinate_matrix.append(coords.flatten())
                    valid_indices.append(i)
            except Exception as e:
                self.logger.warning(f"Failed to extract coordinates from pose {i}: {e}")
        
        if not coordinate_matrix:
            return {
                'clusters': [],
                'n_clusters': 0,
                'rmsd_cutoff': self.rmsd_cutoff,
                'rmsd_matrix': np.array([])
            }
        
        # Convert to numpy array
        coordinate_matrix = np.array(coordinate_matrix)
        
        # Run DBSCAN
        self.logger.info("Performing DBSCAN clustering...")
        
        # Convert RMSD cutoff to appropriate eps for coordinate space
        eps = self.rmsd_cutoff * np.sqrt(coordinate_matrix.shape[1])
        
        dbscan = DBSCAN(eps=eps, min_samples=self.min_cluster_size, metric='euclidean')
        cluster_labels = dbscan.fit_predict(coordinate_matrix)
        
        # Process clusters
        clusters = {}
        for i, cluster_label in enumerate(cluster_labels):
            if cluster_label == -1:  # Noise point
                continue
                
            original_idx = valid_indices[i]
            
            if cluster_label not in clusters:
                clusters[cluster_label] = []
            clusters[cluster_label].append({
                'pose_idx': original_idx,
                'pose': poses[original_idx],
                'score': scores[original_idx]
            })
        
        # Format results
        result_clusters = []
        for cluster_idx, members in clusters.items():
            # Sort cluster members by score
            members.sort(key=lambda x: x['score'])
            
            result_clusters.append({
                'cluster_id': cluster_idx,
                'members': members,
                'size': len(members),
                'representative': members[0]['pose_idx'],
                'best_score': members[0]['score'],
                'worst_score': members[-1]['score'],
                'avg_score': np.mean([m['score'] for m in members])
            })
        
        # Sort clusters by size
        result_clusters.sort(key=lambda x: x['size'], reverse=True)
        
        self.logger.info(f"DBSCAN found {len(result_clusters)} clusters")
        
        return {
            'clusters': result_clusters,
            'n_clusters': len(result_clusters),
            'rmsd_cutoff': self.rmsd_cutoff,
            'rmsd_matrix': np.array([]),  # Not computed for DBSCAN
            'method': 'dbscan'
        }
    
    def _extract_coordinates(self, pose: Any) -> np.ndarray:
        """Extract coordinates from a pose object."""
        
        # Try different attribute names for coordinates
        coord_attrs = ['coords', 'coordinates', 'xyz', 'positions']
        
        for attr in coord_attrs:
            if hasattr(pose, attr):
                coords = getattr(pose, attr)
                if isinstance(coords, np.ndarray):
                    return coords
                elif isinstance(coords, (list, tuple)):
                    return np.array(coords)
        
        # If pose is already a numpy array
        if isinstance(pose, np.ndarray):
            return pose
        
        # If pose is a list/tuple of coordinates
        if isinstance(pose, (list, tuple)):
            return np.array(pose)
        
        raise ValueError(f"Could not extract coordinates from pose of type {type(pose)}")
    
    def visualize_clusters(self, clustering_results: Dict[str, Any], 
                          output_file: Optional[str] = None) -> Optional[Any]:
        """
        Generate visualization of clustering results.
        
        Args:
            clustering_results: Results from cluster_poses method
            output_file: Path to save the visualization
            
        Returns:
            Matplotlib figure object (if matplotlib available)
        """
        if not MATPLOTLIB_AVAILABLE:
            self.logger.warning("matplotlib not available for cluster visualization")
            return None
        
        clusters = clustering_results['clusters']
        rmsd_matrix = clustering_results.get('rmsd_matrix', np.array([]))
        
        # Create figure
        if rmsd_matrix.size > 0:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Plot 1: RMSD heatmap
            im = ax1.imshow(rmsd_matrix, cmap='viridis', interpolation='nearest')
            ax1.set_title('Pairwise RMSD Matrix')
            ax1.set_xlabel('Pose Index')
            ax1.set_ylabel('Pose Index')
            fig.colorbar(im, ax=ax1, label='RMSD (Å)')
        else:
            fig, ax2 = plt.subplots(1, 1, figsize=(10, 6))
        
        # Plot cluster sizes and best scores
        if clusters:
            cluster_sizes = [c['size'] for c in clusters]
            cluster_scores = [c['best_score'] for c in clusters]
            
            # Create bar plot for cluster sizes
            bars = ax2.bar(range(len(clusters)), cluster_sizes, alpha=0.7)
            ax2.set_title('Cluster Sizes and Best Scores')
            ax2.set_xlabel('Cluster Index')
            ax2.set_ylabel('Number of Poses')
            ax2.set_xticks(range(len(clusters)))
            
            # Add score labels
            for i, (bar, score) in enumerate(zip(bars, cluster_scores)):
                ax2.text(
                    bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.5,
                    f'{score:.2f}',
                    ha='center',
                    va='bottom',
                    rotation=90,
                    fontsize=9
                )
            
            # Add second y-axis for scores
            ax2_score = ax2.twinx()
            ax2_score.plot(range(len(clusters)), cluster_scores, 'ro-', alpha=0.7)
            ax2_score.set_ylabel('Best Score', color='r')
            ax2_score.tick_params(axis='y', labelcolor='r')
        else:
            ax2.text(0.5, 0.5, 'No clusters found', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('No Clustering Results')
        
        plt.tight_layout()
        
        # Save figure if output file specified
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            self.logger.info(f"Cluster visualization saved to {output_file}")
            
        return fig
    
    def get_cluster_representatives(self, clustering_results: Dict[str, Any]) -> List[Tuple[int, Any, float]]:
        """
        Get representative poses from each cluster.
        
        Args:
            clustering_results: Results from cluster_poses method
            
        Returns:
            List of (pose_idx, pose, score) tuples for cluster representatives
        """
        representatives = []
        
        for cluster in clustering_results['clusters']:
            best_member = cluster['members'][0]  # Already sorted by score
            representatives.append((
                best_member['pose_idx'],
                best_member['pose'],
                best_member['score']
            ))
        
        return representatives
    
    def cluster_statistics(self, clustering_results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compute clustering statistics.
        
        Args:
            clustering_results: Results from cluster_poses method
            
        Returns:
            Dictionary with clustering statistics
        """
        clusters = clustering_results['clusters']
        
        if not clusters:
            return {
                'n_clusters': 0,
                'total_clustered_poses': 0,
                'largest_cluster_size': 0,
                'average_cluster_size': 0.0,
                'cluster_score_range': (0.0, 0.0)
            }
        
        cluster_sizes = [c['size'] for c in clusters]
        cluster_scores = [c['best_score'] for c in clusters]
        
        return {
            'n_clusters': len(clusters),
            'total_clustered_poses': sum(cluster_sizes),
            'largest_cluster_size': max(cluster_sizes),
            'smallest_cluster_size': min(cluster_sizes),
            'average_cluster_size': np.mean(cluster_sizes),
            'median_cluster_size': np.median(cluster_sizes),
            'best_cluster_score': min(cluster_scores),
            'worst_cluster_score': max(cluster_scores),
            'cluster_score_range': (min(cluster_scores), max(cluster_scores)),
            'rmsd_cutoff_used': clustering_results['rmsd_cutoff'],
            'method_used': clustering_results.get('method', 'unknown')
        }