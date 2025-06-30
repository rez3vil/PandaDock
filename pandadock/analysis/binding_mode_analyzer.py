"""
Binding mode classification for molecular docking results.

This module provides methods for classifying binding modes of docking poses
based on interaction fingerprints and similarity analysis.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Union

from .interaction_analyzer import InteractionFingerprinter

logger = logging.getLogger(__name__)


class BindingModeClassifier:
    """Classify binding modes of docking poses."""
    
    def __init__(self, reference_modes: Optional[Dict] = None, 
                 similarity_threshold: float = 0.7):
        """
        Initialize binding mode classifier.
        
        Args:
            reference_modes: Dictionary of reference binding modes
            similarity_threshold: Threshold for similarity to classify a pose into a mode
        """
        self.reference_modes = reference_modes or {}
        self.similarity_threshold = similarity_threshold
        self.fingerprinter = InteractionFingerprinter()
        self.logger = logging.getLogger(__name__)
    
    def classify_pose(self, protein: Any, pose: Any) -> str:
        """
        Classify a docking pose into a binding mode category.
        
        Args:
            protein: Protein object
            pose: Ligand pose
            
        Returns:
            Binding mode classification
        """
        try:
            # Generate fingerprint for pose
            pose_fp = self.fingerprinter.generate_fingerprint(protein, pose)
            
            # Compare with reference modes
            best_mode = "Unknown"
            best_similarity = 0.0
            
            for mode_name, mode_data in self.reference_modes.items():
                if 'fingerprint' in mode_data:
                    similarity = self.fingerprinter.compare_fingerprints(
                        pose_fp, mode_data['fingerprint']
                    )
                    
                    if similarity > best_similarity:
                        best_similarity = similarity
                        best_mode = mode_name
            
            # Check if similarity exceeds threshold
            if best_similarity >= self.similarity_threshold:
                self.logger.debug(f"Classified pose as '{best_mode}' with similarity {best_similarity:.3f}")
                return best_mode
            else:
                self.logger.debug(f"Pose classified as 'Novel' (best similarity: {best_similarity:.3f})")
                return "Novel"
                
        except Exception as e:
            self.logger.error(f"Error classifying pose: {e}")
            return "Unknown"
    
    def define_reference_mode(self, name: str, protein: Any, reference_pose: Any) -> None:
        """
        Define a new reference binding mode.
        
        Args:
            name: Name for the binding mode
            protein: Protein object
            reference_pose: Reference ligand pose for this binding mode
        """
        try:
            # Generate fingerprint for reference pose
            fp = self.fingerprinter.generate_fingerprint(protein, reference_pose)
            
            # Store reference mode
            self.reference_modes[name] = {
                'fingerprint': fp,
                'pose': reference_pose
            }
            
            total_interactions = sum(fp[k] for k in fp if k != 'interactions')
            self.logger.info(f"Defined reference binding mode '{name}' with {total_interactions} interactions")
            
        except Exception as e:
            self.logger.error(f"Error defining reference mode '{name}': {e}")
    
    def discover_modes(self, protein: Any, poses: List[Any], 
                      n_modes: int = 5) -> List[Dict[str, Any]]:
        """
        Discover binding modes from a set of poses.
        
        Args:
            protein: Protein object
            poses: List of ligand poses
            n_modes: Number of binding modes to discover
            
        Returns:
            List of discovered binding modes
        """
        if not poses:
            return []
        
        try:
            self.logger.info(f"Discovering binding modes from {len(poses)} poses")
            
            # Generate fingerprints for all poses
            fingerprints = []
            valid_poses = []
            
            for i, pose in enumerate(poses):
                try:
                    fp = self.fingerprinter.generate_fingerprint(protein, pose)
                    fingerprints.append(fp)
                    valid_poses.append(pose)
                except Exception as e:
                    self.logger.warning(f"Failed to generate fingerprint for pose {i}: {e}")
            
            if not fingerprints:
                self.logger.warning("No valid fingerprints generated")
                return []
            
            # Calculate pairwise similarity matrix
            n_poses = len(valid_poses)
            similarity_matrix = np.zeros((n_poses, n_poses))
            
            for i in range(n_poses):
                for j in range(i+1, n_poses):
                    similarity = self.fingerprinter.compare_fingerprints(
                        fingerprints[i], fingerprints[j]
                    )
                    similarity_matrix[i, j] = similarity
                    similarity_matrix[j, i] = similarity
            
            # Try to use scikit-learn for clustering
            modes = self._cluster_with_sklearn(
                similarity_matrix, fingerprints, valid_poses, n_modes
            )
            
            if not modes:
                # Fallback to simple mode discovery
                modes = self._simple_mode_discovery(
                    fingerprints, valid_poses, n_modes
                )
            
            self.logger.info(f"Discovered {len(modes)} binding modes")
            return modes
            
        except Exception as e:
            self.logger.error(f"Error discovering modes: {e}")
            return []
    
    def _cluster_with_sklearn(self, similarity_matrix: np.ndarray, 
                            fingerprints: List[Dict], poses: List[Any],
                            n_modes: int) -> List[Dict[str, Any]]:
        """Use scikit-learn for clustering-based mode discovery."""
        try:
            from sklearn.cluster import AgglomerativeClustering
            
            n_poses = len(poses)
            
            # Convert similarity to distance
            distance_matrix = 1.0 - similarity_matrix
            
            # Apply hierarchical clustering
            clustering = AgglomerativeClustering(
                n_clusters=min(n_modes, n_poses),
                affinity='precomputed',
                linkage='average'
            )
            
            cluster_labels = clustering.fit_predict(distance_matrix)
            
            # Create binding mode data
            modes = []
            clusters = {}
            
            for i, label in enumerate(cluster_labels):
                if label not in clusters:
                    clusters[label] = []
                clusters[label].append(i)
            
            # Process each cluster
            for label, indices in clusters.items():
                if len(indices) == 0:
                    continue
                
                # Find most central pose (highest average similarity to other poses in cluster)
                if len(indices) > 1:
                    avg_similarities = []
                    for i in indices:
                        avg_sim = np.mean([similarity_matrix[i, j] for j in indices if j != i])
                        avg_similarities.append(avg_sim)
                    central_idx = indices[np.argmax(avg_similarities)]
                else:
                    central_idx = indices[0]
                
                # Create mode data
                mode_data = {
                    'name': f"Mode {label+1}",
                    'representative': central_idx,
                    'members': indices,
                    'count': len(indices),
                    'fingerprint': fingerprints[central_idx],
                    'pose': poses[central_idx],
                    'best_score': min([self._get_pose_score(poses[i]) for i in indices])
                }
                
                modes.append(mode_data)
            
            # Sort modes by count (largest first)
            modes.sort(key=lambda x: x['count'], reverse=True)
            
            return modes
            
        except ImportError:
            self.logger.warning("scikit-learn not available for clustering")
            return []
        except Exception as e:
            self.logger.warning(f"Clustering failed: {e}")
            return []
    
    def _simple_mode_discovery(self, fingerprints: List[Dict], 
                             poses: List[Any], n_modes: int) -> List[Dict[str, Any]]:
        """Simple approach: take top n_modes poses as representatives."""
        modes = []
        
        for i in range(min(n_modes, len(poses))):
            mode_data = {
                'name': f"Mode {i+1}",
                'representative': i,
                'members': [i],
                'count': 1,
                'fingerprint': fingerprints[i],
                'pose': poses[i],
                'best_score': self._get_pose_score(poses[i])
            }
            modes.append(mode_data)
        
        return modes
    
    def _get_pose_score(self, pose: Any) -> float:
        """Helper to extract score from pose if available."""
        if hasattr(pose, 'score'):
            return pose.score
        elif isinstance(pose, dict) and 'score' in pose:
            return pose['score']
        else:
            return 0.0
    
    def classify_poses(self, protein: Any, poses: List[Any]) -> List[str]:
        """
        Classify multiple poses into binding modes.
        
        Args:
            protein: Protein object
            poses: List of ligand poses
            
        Returns:
            List of binding mode classifications for each pose
        """
        classifications = []
        
        for i, pose in enumerate(poses):
            try:
                classification = self.classify_pose(protein, pose)
                classifications.append(classification)
            except Exception as e:
                self.logger.warning(f"Failed to classify pose {i}: {e}")
                classifications.append("Unknown")
        
        return classifications
    
    def get_mode_statistics(self, classifications: List[str]) -> Dict[str, Any]:
        """
        Get statistics about binding mode classifications.
        
        Args:
            classifications: List of classifications from classify_poses
            
        Returns:
            Dictionary with mode statistics
        """
        from collections import Counter
        
        mode_counts = Counter(classifications)
        total_poses = len(classifications)
        
        stats = {
            'total_poses': total_poses,
            'unique_modes': len(mode_counts),
            'mode_counts': dict(mode_counts),
            'mode_percentages': {
                mode: count/total_poses * 100 
                for mode, count in mode_counts.items()
            }
        }
        
        return stats
    
    def export_reference_modes(self, filename: str) -> None:
        """
        Export reference modes to a file.
        
        Args:
            filename: Output filename
        """
        try:
            import json
            
            # Convert modes to serializable format
            export_data = {}
            for name, mode_data in self.reference_modes.items():
                export_data[name] = {
                    'fingerprint': mode_data['fingerprint'],
                    'similarity_threshold': self.similarity_threshold
                }
            
            with open(filename, 'w') as f:
                json.dump(export_data, f, indent=2)
            
            self.logger.info(f"Exported {len(self.reference_modes)} reference modes to {filename}")
            
        except Exception as e:
            self.logger.error(f"Error exporting reference modes: {e}")
    
    def import_reference_modes(self, filename: str) -> None:
        """
        Import reference modes from a file.
        
        Args:
            filename: Input filename
        """
        try:
            import json
            
            with open(filename, 'r') as f:
                import_data = json.load(f)
            
            for name, mode_data in import_data.items():
                self.reference_modes[name] = {
                    'fingerprint': mode_data['fingerprint']
                }
                
                if 'similarity_threshold' in mode_data:
                    self.similarity_threshold = mode_data['similarity_threshold']
            
            self.logger.info(f"Imported {len(import_data)} reference modes from {filename}")
            
        except Exception as e:
            self.logger.error(f"Error importing reference modes: {e}")