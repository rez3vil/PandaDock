"""
Physics-based docking engine (Glide-style)
Implements detailed molecular mechanics scoring and flexible docking
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from scipy.optimize import minimize
from tqdm import tqdm

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import DockingEngine, Pose
from scoring.scoring_functions import ScoringFunctions
from utils.math_utils import rotation_matrix, quaternion_to_matrix
from docking.flexible_docking import FlexibleDocking
from docking.pose_filtering import PoseFiltering


class PhysicsEngine(DockingEngine):
    """
    Physics-based docking engine similar to Glide
    
    Features:
    - Detailed molecular mechanics scoring
    - Flexible side chain sampling
    - Energy minimization
    - Clash detection and resolution
    - Systematic conformer generation
    """
    
    def __init__(self, config):
        super().__init__(config)
        self.scoring = ScoringFunctions(config)
        self.flexible_docking = FlexibleDocking(config) if config.docking.flexible_residues else None
        self.pose_filter = PoseFiltering(config)
        
        # Physics engine specific parameters
        self.energy_threshold = config.docking.energy_range * 1.4  # kcal/mol
        self.minimization_steps = config.docking.minimization_steps
        self.vdw_scale = config.docking.vdw_scale
        self.electrostatic_scale = config.docking.electrostatic_scale
        
        # Conformer generation parameters
        self.num_conformers = 1000
        self.torsion_increment = 30.0  # degrees
        
        self.logger.info("Initialized PhysicsEngine (Glide-style)")
    
    def dock(self, protein_file: str, ligand_file: str) -> List[Pose]:
        """
        Main docking method using physics-based approach
        
        Steps:
        1. Prepare receptor and ligand
        2. Generate ligand conformers
        3. Sample binding poses
        4. Optimize poses with flexible sidechains
        5. Score and filter poses
        6. Return top poses
        """
        self.logger.info(f"Starting physics-based docking: {protein_file} + {ligand_file}")
        
        # Prepare structures
        self.prepare_receptor(protein_file)
        self.prepare_ligand(ligand_file)
        
        # Generate conformers systematically
        conformers = self.generate_systematic_conformers()
        self.logger.info(f"Generated {len(conformers)} conformers")
        
        # Sample binding poses for each conformer with progress bar
        all_poses = []
        conformer_pbar = tqdm(conformers, desc="🧬 Processing conformers", unit="conformer", 
                             bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
        for conformer in conformer_pbar:
            poses = self.sample_binding_poses(conformer)
            all_poses.extend(poses)
            conformer_pbar.set_postfix({"poses": len(all_poses)})
        conformer_pbar.close()
        
        self.logger.info(f"Generated {len(all_poses)} initial poses")
        
        # Filter poses by basic criteria
        print("🔍 Filtering poses by energy and clash criteria...")
        
        # First calculate energy for all poses if not already done
        print("⚡ Computing initial energies for pose filtering...")
        energy_pbar = tqdm(all_poses, desc="🔬 Initial energy calc", unit="pose",
                          bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
        for pose in energy_pbar:
            if pose.energy == 0.0:  # Only calculate if not already done
                pose.energy = self.scoring.calculate_total_energy(pose.coordinates)
        energy_pbar.close()
        
        filtered_poses = self.pose_filter.filter_by_energy(all_poses, self.energy_threshold)
        filtered_poses = self.pose_filter.filter_by_clash(filtered_poses)
        
        self.logger.info(f"After filtering: {len(filtered_poses)} poses")
        
        # Optimize poses with flexible sidechains
        if self.flexible_docking:
            print("🔧 Optimizing side-chain flexibility...")
            optimized_poses = []
            sidechain_pbar = tqdm(filtered_poses, desc="⚡ Side-chain optimization", unit="pose",
                                 bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
            for pose in sidechain_pbar:
                optimized_pose = self.optimize_with_flexibility(pose)
                optimized_poses.append(optimized_pose)
            sidechain_pbar.close()
            filtered_poses = optimized_poses
        
        # Final scoring and ranking
        print("📊 Computing final physics-based scores...")
        scoring_pbar = tqdm(filtered_poses, desc="🎯 Final scoring", unit="pose",
                           bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
        for pose in scoring_pbar:
            pose.score = self.score(pose)
        scoring_pbar.close()
        
        # Sort by score and cluster
        filtered_poses.sort(key=lambda x: x.score)
        final_poses = self.cluster_poses(filtered_poses)
        
        # Take top N poses
        final_poses = final_poses[:self.config.docking.num_poses]
        
        self.logger.info(f"Final result: {len(final_poses)} poses")
        return final_poses
    
    def generate_systematic_conformers(self) -> List[np.ndarray]:
        """
        Generate conformers using systematic torsion sampling
        Similar to Glide's conformer generation
        """
        if not self.ligand:
            raise ValueError("Ligand must be prepared before generating conformers")
        
        conformers = []
        
        # Start with the original ligand structure
        base_coords = self.ligand['coordinates'].copy()
        conformers.append(base_coords)
        
        # Generate additional conformers with systematic variations
        num_torsions = min(4, len(base_coords) // 3)  # Reasonable number of torsions
        angles = np.arange(0, 360, self.torsion_increment)
        
        target_conformers = min(self.num_conformers, 50)  # Limit for efficiency
        
        for i in range(1, target_conformers):
            # Generate torsion combination
            torsion_angles = np.random.choice(angles, num_torsions)
            
            # Build conformer coordinates
            conformer = self.build_conformer_from_torsions(torsion_angles)
            
            # Quick energy filter
            if self.quick_energy_check(conformer):
                conformers.append(conformer)
        
        self.logger.info(f"Generated {len(conformers)} valid conformers from {target_conformers} attempts")
        return conformers
    
    def build_conformer_from_torsions(self, torsion_angles: np.ndarray) -> np.ndarray:
        """Build 3D coordinates from torsion angles"""
        if not self.ligand:
            raise ValueError("Ligand must be prepared before building conformers")
        
        # Use actual ligand coordinates as base
        base_coords = self.ligand['coordinates'].copy()
        
        # Apply torsional changes (simplified implementation)
        # In a real implementation, this would properly rotate around torsion bonds
        conformer = base_coords.copy()
        
        # Apply small perturbations based on torsion angles
        for i, angle in enumerate(torsion_angles):
            if i < len(conformer):
                # Simple rotation around z-axis for demonstration
                rad = np.radians(angle)
                cos_a, sin_a = np.cos(rad), np.sin(rad)
                rotation = np.array([[cos_a, -sin_a, 0], [sin_a, cos_a, 0], [0, 0, 1]])
                
                # Center molecule
                center = np.mean(conformer, axis=0)
                centered = conformer - center
                
                # Apply rotation
                rotated = np.dot(centered, rotation.T)
                conformer = rotated + center
                
                # Apply small displacement
                conformer += np.random.normal(0, 0.1, conformer.shape)
        
        return conformer
    
    def quick_energy_check(self, conformer: np.ndarray) -> bool:
        """Quick energy check to filter bad conformers"""
        # Placeholder implementation
        # Check for severe clashes
        distances = np.linalg.norm(conformer[:, np.newaxis] - conformer[np.newaxis, :], axis=2)
        np.fill_diagonal(distances, np.inf)
        
        # Check if any atoms are too close
        min_distance = np.min(distances)
        return min_distance > 1.0  # Angstroms
    
    def sample_binding_poses(self, conformer: np.ndarray) -> List[Pose]:
        """
        Sample binding poses for a given conformer
        Uses systematic rotation and translation
        """
        poses = []
        
        # Sample translations within grid box
        num_translations = 50
        translations = self.sample_translations(num_translations)
        
        # Sample rotations
        num_rotations = 20
        rotations = self.sample_rotations(num_rotations)
        
        for translation in translations:
            for rotation in rotations:
                # Apply transformation
                transformed_coords = self.apply_transformation(conformer, rotation, translation)
                
                # Create pose
                pose = Pose(
                    coordinates=transformed_coords,
                    score=0.0,
                    energy=0.0,
                    ligand_name=self.ligand['name'] if self.ligand else 'unknown',
                    pose_id=f"pose_{len(poses)}"
                )
                
                # Calculate initial energy for the pose
                pose.energy = self.scoring.calculate_total_energy(transformed_coords)
                
                # Quick check if pose is reasonable
                if self.validate_pose(pose):
                    poses.append(pose)
        
        return poses
    
    def sample_translations(self, num_samples: int) -> List[np.ndarray]:
        """Sample translation vectors within grid box"""
        translations = []
        min_bounds, max_bounds = self.grid_box.get_bounds()
        
        for _ in range(num_samples):
            translation = np.random.uniform(min_bounds, max_bounds)
            translations.append(translation)
        
        return translations
    
    def sample_rotations(self, num_samples: int) -> List[np.ndarray]:
        """Sample rotation matrices"""
        rotations = []
        
        for _ in range(num_samples):
            # Random rotation using Euler angles
            angles = np.random.uniform(0, 2*np.pi, 3)
            rotation = rotation_matrix(angles)
            rotations.append(rotation)
        
        return rotations
    
    def apply_transformation(self, coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
        """Apply rotation and translation to coordinates while preserving molecular geometry"""
        # Center coordinates around their geometric center
        center = np.mean(coords, axis=0)
        centered = coords - center
        
        # Apply rotation around center
        rotated = np.dot(centered, rotation.T)
        
        # Translate to new position within grid box
        grid_center = self.grid_box.center
        final_coords = rotated + grid_center + translation
        
        # Ensure coordinates stay within grid box bounds
        min_bounds, max_bounds = self.grid_box.get_bounds()
        
        # Check if any atoms are outside bounds and adjust if needed
        coord_center = np.mean(final_coords, axis=0)
        for i in range(3):
            if coord_center[i] < min_bounds[i] + 2.0:
                final_coords[:, i] += (min_bounds[i] + 2.0 - coord_center[i])
            elif coord_center[i] > max_bounds[i] - 2.0:
                final_coords[:, i] -= (coord_center[i] - (max_bounds[i] - 2.0))
        
        return final_coords
    
    def optimize_with_flexibility(self, pose: Pose) -> Pose:
        """
        Optimize pose with flexible sidechains
        Uses the FlexibleDocking module
        """
        if not self.flexible_docking:
            return self.local_minimize(pose)
        
        # Optimize with flexible sidechains
        optimized_pose = self.flexible_docking.optimize_with_sidechains(pose)
        
        # Additional local minimization
        optimized_pose = self.local_minimize(optimized_pose)
        
        return optimized_pose
    
    def local_minimize(self, pose: Pose) -> Pose:
        """
        Local energy minimization of a pose
        Uses gradient-based optimization
        """
        def energy_function(coords_flat):
            coords = coords_flat.reshape(-1, 3)
            temp_pose = Pose(coordinates=coords, score=0.0, energy=0.0)
            return self.calculate_energy(temp_pose)
        
        # Flatten coordinates
        initial_coords = pose.coordinates.flatten()
        
        # Minimize
        result = minimize(
            energy_function,
            initial_coords,
            method='L-BFGS-B',
            options={'maxiter': self.minimization_steps}
        )
        
        # Create optimized pose
        optimized_coords = result.x.reshape(-1, 3)
        optimized_pose = Pose(
            coordinates=optimized_coords,
            score=pose.score,
            energy=result.fun,
            ligand_name=pose.ligand_name,
            pose_id=pose.pose_id,
            flexible_residues=pose.flexible_residues
        )
        
        return optimized_pose
    
    def score(self, pose: Pose) -> float:
        """
        Score a pose using physics-based scoring function
        """
        # Calculate individual energy terms
        pose.vdw_energy = self.scoring.calculate_vdw_energy(pose.coordinates)
        pose.electrostatic_energy = self.scoring.calculate_electrostatic_energy(pose.coordinates)
        pose.hbond_energy = self.scoring.calculate_hbond_energy(pose.coordinates)
        pose.hydrophobic_energy = self.scoring.calculate_hydrophobic_energy(pose.coordinates)
        pose.solvation_energy = self.scoring.calculate_solvation_energy(pose.coordinates)
        pose.entropy_energy = self.scoring.calculate_entropy_penalty(pose.coordinates)
        
        # Calculate total energy using the scoring function
        pose.energy = self.scoring.calculate_total_energy(pose.coordinates)
        
        # Calculate interactions
        pose.hbond_interactions = self.scoring.find_hbond_interactions(pose.coordinates)
        pose.hydrophobic_interactions = self.scoring.find_hydrophobic_interactions(pose.coordinates)
        pose.salt_bridge_interactions = self.scoring.find_salt_bridge_interactions(pose.coordinates)
        
        # Calculate clash score
        pose.clash_score = self.scoring.calculate_clash_score(pose.coordinates)
        
        # Final score equals energy (lower is better for both)
        # Add small clash penalty
        score = pose.energy + pose.clash_score * 0.1
        
        return score
    
    def calculate_energy(self, pose: Pose) -> float:
        """Calculate total energy for a pose"""
        return self.scoring.calculate_total_energy(pose.coordinates)
    
    def get_engine_info(self) -> Dict[str, Any]:
        """Get information about the physics engine"""
        info = super().get_engine_info()
        info.update({
            'engine_type': 'PhysicsEngine',
            'description': 'Glide-style physics-based docking',
            'features': [
                'Systematic conformer generation',
                'Flexible side chain sampling',
                'Energy minimization',
                'Detailed molecular mechanics scoring',
                'Clash detection and resolution'
            ],
            'parameters': {
                'energy_threshold': self.energy_threshold,
                'minimization_steps': self.minimization_steps,
                'vdw_scale': self.vdw_scale,
                'electrostatic_scale': self.electrostatic_scale,
                'num_conformers': self.num_conformers,
                'torsion_increment': self.torsion_increment
            }
        })
        return info