# -*- coding: utf-8 -*-
"""
Pose filtering and validation module
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import Pose
from utils.math_utils import distance_matrix, calculate_rmsd


class PoseFiltering:
    """
    Pose filtering and validation system
    
    Features:
    - Clash detection using van der Waals radii
    - Energy-based filtering
    - RMSD-based clustering
    - Binding site occupancy validation
    - Physicochemical property filtering
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Filtering thresholds
        self.clash_threshold = 2.0  # Angstroms
        self.energy_threshold = config.docking.energy_range  # kcal/mol
        self.rmsd_threshold = 2.0  # Angstroms for clustering
        self.binding_site_threshold = 0.5  # Minimum fraction of binding site occupied
        
        # Van der Waals radii (Angstroms)
        self.vdw_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80,
            'S': 1.80, 'F': 1.47, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
        }
        
        # Physicochemical limits
        self.max_molecular_weight = 800.0  # Da
        self.max_logp = 5.0
        self.max_hbd = 5  # Hydrogen bond donors
        self.max_hba = 10  # Hydrogen bond acceptors
        
        self.logger.info("Initialized PoseFiltering")
    
    def filter_poses(self, poses: List[Pose]) -> List[Pose]:
        """
        Apply all filtering criteria to poses
        
        Args:
            poses: List of poses to filter
            
        Returns:
            Filtered and clustered poses
        """
        self.logger.info(f"Filtering {len(poses)} poses")
        
        if not poses:
            return poses
        
        # Step 1: Basic validation
        valid_poses = self.filter_by_validity(poses)
        self.logger.info(f"After validity filter: {len(valid_poses)} poses")
        
        # Step 2: Energy filtering
        energy_filtered = self.filter_by_energy(valid_poses, self.energy_threshold)
        self.logger.info(f"After energy filter: {len(energy_filtered)} poses")
        
        # Step 3: Clash detection
        clash_filtered = self.filter_by_clash(energy_filtered)
        self.logger.info(f"After clash filter: {len(clash_filtered)} poses")
        
        # Step 4: Binding site occupancy
        occupancy_filtered = self.filter_by_binding_site_occupancy(clash_filtered)
        self.logger.info(f"After occupancy filter: {len(occupancy_filtered)} poses")
        
        # Step 5: Physicochemical properties
        property_filtered = self.filter_by_physicochemical_properties(occupancy_filtered)
        self.logger.info(f"After property filter: {len(property_filtered)} poses")
        
        # Step 6: RMSD clustering
        clustered_poses = self.cluster_poses_by_rmsd(property_filtered)
        self.logger.info(f"After clustering: {len(clustered_poses)} poses")
        
        return clustered_poses
    
    def filter_by_validity(self, poses: List[Pose]) -> List[Pose]:
        """Filter poses by basic validity criteria"""
        valid_poses = []
        
        for pose in poses:
            if self.is_valid_pose(pose):
                valid_poses.append(pose)
        
        return valid_poses
    
    def is_valid_pose(self, pose: Pose) -> bool:
        """Check if pose passes basic validity tests"""
        # Check for NaN or infinite coordinates
        if np.any(~np.isfinite(pose.coordinates)):
            return False
        
        # Check if pose has reasonable number of atoms
        if len(pose.coordinates) < 3 or len(pose.coordinates) > 200:
            return False
        
        # Check if pose is within reasonable distance from origin
        center = np.mean(pose.coordinates, axis=0)
        if np.linalg.norm(center) > 100.0:  # 100 Angstroms from origin
            return False
        
        # Check for reasonable bond lengths
        if not self.check_reasonable_geometry(pose.coordinates):
            return False
        
        return True
    
    def check_reasonable_geometry(self, coordinates: np.ndarray) -> bool:
        """Check if molecular geometry is reasonable"""
        # Calculate all pairwise distances
        distances = distance_matrix(coordinates, coordinates)
        
        # Remove diagonal (self-distances)
        np.fill_diagonal(distances, np.inf)
        
        # Check for atoms that are too close
        min_distance = np.min(distances)
        if min_distance < 0.5:  # Less than 0.5 Angstroms
            return False
        
        # Check for atoms that are too far from all others
        for i in range(len(coordinates)):
            min_dist_to_others = np.min(distances[i])
            if min_dist_to_others > 10.0:  # More than 10 Angstroms
                return False
        
        return True
    
    def filter_by_energy(self, poses: List[Pose], threshold: float) -> List[Pose]:
        """Filter poses by energy threshold"""
        if not poses:
            return poses
        
        # Filter out poses with None energy and collect valid energies
        valid_poses = []
        energies = []
        
        for pose in poses:
            if pose.energy is not None:
                valid_poses.append(pose)
                energies.append(pose.energy)
            else:
                self.logger.warning(f"Pose {pose.pose_id} has None energy, skipping in energy filter")
        
        if not energies:
            self.logger.warning("No poses with valid energy values found")
            return poses  # Return original poses if none have valid energies
        
        # Find reference energy (best pose)
        min_energy = min(energies)
        max_energy = max(energies)
        
        self.logger.info(f"Energy range: {min_energy:.2f} to {max_energy:.2f} kcal/mol")
        self.logger.info(f"Energy threshold: {threshold:.2f} kcal/mol")
        
        # Filter poses within threshold of best energy
        filtered_poses = []
        for pose in valid_poses:
            if pose.energy <= min_energy + threshold:
                filtered_poses.append(pose)
        
        self.logger.info(f"Energy filter: {len(filtered_poses)}/{len(poses)} poses passed")
        return filtered_poses
    
    def filter_by_clash(self, poses: List[Pose]) -> List[Pose]:
        """Filter poses by clash detection"""
        if not poses:
            return poses
            
        filtered_poses = []
        clash_scores = []
        
        for pose in poses:
            clash_score = self.calculate_clash_score(pose)
            pose.clash_score = clash_score
            clash_scores.append(clash_score)
            
            if clash_score < 5.0:  # Threshold for acceptable clash
                filtered_poses.append(pose)
        
        if clash_scores:
            min_clash = min(clash_scores)
            max_clash = max(clash_scores)
            self.logger.info(f"Clash score range: {min_clash:.2f} to {max_clash:.2f}")
        
        self.logger.info(f"Clash filter: {len(filtered_poses)}/{len(poses)} poses passed")
        return filtered_poses
    
    def calculate_clash_score(self, pose: Pose) -> float:
        """
        Calculate clash score for a pose
        
        Only counts severe clashes, not normal bonded interactions
        """
        coordinates = pose.coordinates
        distances = distance_matrix(coordinates, coordinates)
        
        # Remove diagonal
        np.fill_diagonal(distances, np.inf)
        
        clash_score = 0.0
        
        # Check for severe clashes between atoms
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                distance = distances[i, j]
                
                # Only consider severe clashes (atoms much closer than they should be)
                # Normal bonds are 1.2-1.8 Å, so only flag distances < 1.0 Å as clashes
                if distance < 1.0:  # Very close atoms
                    clash_score += (1.0 - distance) ** 2  # Quadratic penalty
                
                # Check for moderately close non-bonded atoms
                # Skip 1-2 (bonded) and 1-3 (angle) interactions
                if j - i > 3 and distance < 2.0:  # Non-bonded atoms too close
                    clash_score += 0.1 * (2.0 - distance)  # Small penalty
        
        return clash_score
    
    def filter_by_binding_site_occupancy(self, poses: List[Pose]) -> List[Pose]:
        """Filter poses by binding site occupancy"""
        filtered_poses = []
        
        for pose in poses:
            occupancy = self.calculate_binding_site_occupancy(pose)
            pose.binding_site_coverage = occupancy
            
            if occupancy >= self.binding_site_threshold:
                filtered_poses.append(pose)
        
        return filtered_poses
    
    def calculate_binding_site_occupancy(self, pose: Pose) -> float:
        """
        Calculate how well the pose occupies the binding site
        
        Returns fraction of binding site volume occupied
        """
        # This is a simplified implementation
        # In reality, this would use detailed binding site definition
        
        # Check if pose center is within binding site
        pose_center = np.mean(pose.coordinates, axis=0)
        binding_site_center = np.array([0.0, 0.0, 0.0])  # Placeholder
        
        distance_to_center = np.linalg.norm(pose_center - binding_site_center)
        
        # Simple occupancy score based on distance to center
        if distance_to_center > 10.0:
            return 0.0
        
        occupancy = max(0.0, 1.0 - distance_to_center / 10.0)
        return occupancy
    
    def filter_by_physicochemical_properties(self, poses: List[Pose]) -> List[Pose]:
        """Filter poses by physicochemical properties (Lipinski's Rule of Five)"""
        filtered_poses = []
        
        for pose in poses:
            if self.passes_physicochemical_filter(pose):
                filtered_poses.append(pose)
        
        return filtered_poses
    
    def passes_physicochemical_filter(self, pose: Pose) -> bool:
        """Check if pose passes physicochemical property filters"""
        # This is a simplified implementation
        # In reality, this would calculate actual molecular properties
        
        # Estimate molecular weight from number of atoms
        estimated_mw = len(pose.coordinates) * 12.0  # Rough estimate
        
        if estimated_mw > self.max_molecular_weight:
            return False
        
        # Other properties would be calculated from actual molecular structure
        # For now, assume all poses pass
        return True
    
    def cluster_poses_by_rmsd(self, poses: List[Pose]) -> List[Pose]:
        """
        Cluster poses by RMSD and return representative poses
        
        Uses a simple greedy clustering algorithm
        """
        if not poses:
            return poses
        
        # Sort poses by score
        sorted_poses = sorted(poses, key=lambda x: x.score)
        
        clustered_poses = []
        
        for pose in sorted_poses:
            # Check if this pose is similar to any existing cluster representative
            is_new_cluster = True
            
            for representative in clustered_poses:
                rmsd = self.calculate_rmsd(pose, representative)
                pose.rmsd = rmsd
                
                if rmsd < self.rmsd_threshold:
                    is_new_cluster = False
                    break
            
            if is_new_cluster:
                clustered_poses.append(pose)
        
        return clustered_poses
    
    def calculate_rmsd(self, pose1: Pose, pose2: Pose) -> float:
        """Calculate RMSD between two poses"""
        if pose1.coordinates.shape != pose2.coordinates.shape:
            return float('inf')
        
        return calculate_rmsd(pose1.coordinates, pose2.coordinates)
    
    def filter_by_diversity(self, poses: List[Pose], max_poses: int = 100) -> List[Pose]:
        """
        Filter poses to maintain diversity
        
        Uses maximum diversity selection
        """
        if len(poses) <= max_poses:
            return poses
        
        # Start with best pose
        selected_poses = [poses[0]]
        remaining_poses = poses[1:]
        
        while len(selected_poses) < max_poses and remaining_poses:
            # Find pose with maximum minimum distance to selected poses
            best_pose = None
            best_min_distance = 0.0
            
            for pose in remaining_poses:
                min_distance = float('inf')
                
                for selected_pose in selected_poses:
                    distance = self.calculate_rmsd(pose, selected_pose)
                    min_distance = min(min_distance, distance)
                
                if min_distance > best_min_distance:
                    best_min_distance = min_distance
                    best_pose = pose
            
            if best_pose:
                selected_poses.append(best_pose)
                remaining_poses.remove(best_pose)
            else:
                break
        
        return selected_poses
    
    def analyze_pose_quality(self, pose: Pose) -> Dict[str, Any]:
        """
        Analyze the quality of a pose
        
        Returns detailed quality metrics
        """
        quality_metrics = {
            'clash_score': self.calculate_clash_score(pose),
            'binding_site_occupancy': self.calculate_binding_site_occupancy(pose),
            'geometry_score': self.calculate_geometry_score(pose),
            'compactness_score': self.calculate_compactness_score(pose),
            'strain_energy': self.calculate_strain_energy(pose)
        }
        
        # Overall quality score
        quality_score = (
            -quality_metrics['clash_score'] * 0.3 +
            quality_metrics['binding_site_occupancy'] * 0.3 +
            quality_metrics['geometry_score'] * 0.2 +
            quality_metrics['compactness_score'] * 0.1 +
            -quality_metrics['strain_energy'] * 0.1
        )
        
        quality_metrics['overall_quality'] = quality_score
        
        return quality_metrics
    
    def calculate_geometry_score(self, pose: Pose) -> float:
        """Calculate geometric quality score"""
        # Check for reasonable bond lengths and angles
        # This is a simplified implementation
        
        coordinates = pose.coordinates
        if len(coordinates) < 3:
            return 0.0
        
        # Calculate distances
        distances = distance_matrix(coordinates, coordinates)
        np.fill_diagonal(distances, np.inf)
        
        # Check for reasonable distance distribution
        min_distances = np.min(distances, axis=1)
        
        # Good geometry has consistent minimum distances
        mean_min_dist = np.mean(min_distances)
        std_min_dist = np.std(min_distances)
        
        # Higher score for more consistent geometry
        geometry_score = max(0.0, 1.0 - std_min_dist / mean_min_dist)
        
        return geometry_score
    
    def calculate_compactness_score(self, pose: Pose) -> float:
        """Calculate compactness score"""
        coordinates = pose.coordinates
        
        # Calculate radius of gyration
        center = np.mean(coordinates, axis=0)
        distances_to_center = np.linalg.norm(coordinates - center, axis=1)
        radius_of_gyration = np.sqrt(np.mean(distances_to_center**2))
        
        # Normalize by number of atoms
        normalized_rg = radius_of_gyration / np.sqrt(len(coordinates))
        
        # Higher score for more compact structures
        compactness_score = max(0.0, 1.0 - normalized_rg / 5.0)
        
        return compactness_score
    
    def calculate_strain_energy(self, pose: Pose) -> float:
        """Calculate internal strain energy"""
        # This would calculate conformational strain
        # Simplified implementation
        
        coordinates = pose.coordinates
        if len(coordinates) < 3:
            return 0.0
        
        # Calculate distances
        distances = distance_matrix(coordinates, coordinates)
        np.fill_diagonal(distances, np.inf)
        
        # Find very short distances (indicating strain)
        strain_energy = 0.0
        threshold = 2.0  # Angstroms
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                distance = distances[i, j]
                if distance < threshold:
                    strain_energy += (threshold - distance)**2
        
        return strain_energy
    
    def get_filtering_statistics(self, original_poses: List[Pose], filtered_poses: List[Pose]) -> Dict[str, Any]:
        """Get statistics about the filtering process"""
        stats = {
            'original_count': len(original_poses),
            'filtered_count': len(filtered_poses),
            'rejection_rate': 1.0 - len(filtered_poses) / max(1, len(original_poses)),
            'energy_range': None,
            'clash_score_range': None,
            'rmsd_range': None
        }
        
        if filtered_poses:
            energies = [pose.energy for pose in filtered_poses]
            clash_scores = [pose.clash_score for pose in filtered_poses]
            
            stats['energy_range'] = (min(energies), max(energies))
            stats['clash_score_range'] = (min(clash_scores), max(clash_scores))
            
            if len(filtered_poses) > 1:
                rmsds = []
                for i in range(len(filtered_poses)):
                    for j in range(i + 1, len(filtered_poses)):
                        rmsd = self.calculate_rmsd(filtered_poses[i], filtered_poses[j])
                        rmsds.append(rmsd)
                
                if rmsds:
                    stats['rmsd_range'] = (min(rmsds), max(rmsds))
        
        return stats