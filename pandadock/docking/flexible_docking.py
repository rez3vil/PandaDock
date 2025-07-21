# -*- coding: utf-8 -*-
"""
Flexible docking module for side-chain sampling and optimization
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from scipy.optimize import minimize

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import Pose
from utils.rotamer_lib import RotamerLibrary
from utils.math_utils import rotation_matrix, distance_matrix


class FlexibleResidue:
    """Represents a flexible residue with rotamer states"""
    
    def __init__(self, residue_id: str, residue_type: str, backbone_coords: np.ndarray):
        self.residue_id = residue_id
        self.residue_type = residue_type
        self.backbone_coords = backbone_coords
        self.rotamers = []
        self.current_rotamer_idx = 0
        self.original_coords = None
        
    def add_rotamer(self, sidechain_coords: np.ndarray, probability: float = 1.0):
        """Add a rotamer conformation"""
        self.rotamers.append({
            'coords': sidechain_coords,
            'probability': probability
        })
        
    def get_current_coords(self) -> np.ndarray:
        """Get current rotamer coordinates"""
        if not self.rotamers:
            return self.backbone_coords
        
        rotamer = self.rotamers[self.current_rotamer_idx]
        return np.vstack([self.backbone_coords, rotamer['coords']])
    
    def set_rotamer(self, idx: int):
        """Set current rotamer"""
        if 0 <= idx < len(self.rotamers):
            self.current_rotamer_idx = idx
    
    def get_rotamer_probability(self) -> float:
        """Get probability of current rotamer"""
        if not self.rotamers:
            return 1.0
        return self.rotamers[self.current_rotamer_idx]['probability']


class FlexibleDocking:
    """
    Flexible docking implementation with side-chain sampling
    
    Features:
    - Rotamer library-based side-chain sampling
    - Dead-end elimination
    - Monte Carlo optimization
    - Clash detection and resolution
    - Energy minimization
    """
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Initialize rotamer library
        self.rotamer_lib = RotamerLibrary()
        
        # Flexible residues
        self.flexible_residues = {}
        self.flexible_residue_names = config.docking.flexible_residues
        
        # Optimization parameters
        self.max_rotamer_combinations = 10000
        self.clash_threshold = 2.0  # Angstroms
        self.monte_carlo_steps = 1000
        self.temperature = 1.0
        
        # Energy weights
        self.sidechain_weight = 0.5
        self.ligand_sidechain_weight = 1.0
        self.sidechain_sidechain_weight = 0.3
        
        self.logger.info(f"Initialized FlexibleDocking with {len(self.flexible_residue_names)} flexible residues")
    
    def setup_flexible_residues(self, protein_coords: np.ndarray, residue_info: List[Dict]):
        """Setup flexible residues from protein structure"""
        self.logger.info("Setting up flexible residues")
        
        for residue in residue_info:
            residue_id = residue['id']
            residue_type = residue['type']
            
            if residue_id in self.flexible_residue_names:
                # Extract backbone coordinates
                backbone_coords = residue['backbone_coords']
                
                # Create flexible residue
                flex_residue = FlexibleResidue(residue_id, residue_type, backbone_coords)
                
                # Generate rotamers
                rotamers = self.rotamer_lib.get_rotamers(residue_type)
                
                for rotamer_data in rotamers:
                    sidechain_coords = self._build_sidechain_coords(
                        backbone_coords, rotamer_data
                    )
                    flex_residue.add_rotamer(
                        sidechain_coords, 
                        rotamer_data['probability']
                    )
                
                self.flexible_residues[residue_id] = flex_residue
                self.logger.debug(f"Added flexible residue {residue_id} with {len(rotamers)} rotamers")
    
    def _build_sidechain_coords(self, backbone_coords: np.ndarray, rotamer_data: Dict) -> np.ndarray:
        """Build side-chain coordinates from rotamer data"""
        # This is a simplified implementation
        # In reality, this would use proper molecular mechanics
        
        # Get chi angles from rotamer data
        chi_angles = rotamer_data.get('chi_angles', [])
        
        # Build side-chain coordinates (placeholder)
        num_sidechain_atoms = len(chi_angles) + 2  # Simplified
        sidechain_coords = np.random.randn(num_sidechain_atoms, 3) * 1.5
        
        # Position relative to backbone
        ca_pos = backbone_coords[1]  # Assuming CA is second atom
        sidechain_coords += ca_pos
        
        return sidechain_coords
    
    def optimize_with_sidechains(self, pose: Pose) -> Pose:
        """
        Optimize pose with flexible side-chains
        
        Uses dead-end elimination followed by Monte Carlo optimization
        """
        self.logger.debug(f"Optimizing pose {pose.pose_id} with flexible sidechains")
        
        if not self.flexible_residues:
            return pose
        
        # Initialize side-chain conformations
        self._initialize_sidechain_conformations(pose)
        
        # Dead-end elimination
        self._dead_end_elimination(pose)
        
        # Monte Carlo optimization
        optimized_pose = self._monte_carlo_optimization(pose)
        
        # Final energy minimization
        final_pose = self._minimize_with_sidechains(optimized_pose)
        
        return final_pose
    
    def _initialize_sidechain_conformations(self, pose: Pose):
        """Initialize side-chain conformations"""
        for residue_id, flex_residue in self.flexible_residues.items():
            # Find best rotamer based on ligand interaction
            best_rotamer_idx = self._find_best_initial_rotamer(flex_residue, pose)
            flex_residue.set_rotamer(best_rotamer_idx)
    
    def _find_best_initial_rotamer(self, flex_residue: FlexibleResidue, pose: Pose) -> int:
        """Find best initial rotamer for a residue"""
        best_idx = 0
        best_energy = float('inf')
        
        for i in range(len(flex_residue.rotamers)):
            flex_residue.set_rotamer(i)
            
            # Calculate interaction energy with ligand
            sidechain_coords = flex_residue.get_current_coords()
            energy = self._calculate_ligand_sidechain_energy(pose.coordinates, sidechain_coords)
            
            if energy < best_energy:
                best_energy = energy
                best_idx = i
        
        return best_idx
    
    def _dead_end_elimination(self, pose: Pose):
        """
        Dead-end elimination to remove impossible rotamer combinations
        """
        self.logger.debug("Performing dead-end elimination")
        
        # Calculate pairwise energies
        pairwise_energies = self._calculate_pairwise_energies(pose)
        
        # Eliminate dead-end rotamers
        eliminated = True
        while eliminated:
            eliminated = False
            
            for residue_id, flex_residue in self.flexible_residues.items():
                for i in range(len(flex_residue.rotamers)):
                    if self._is_dead_end_rotamer(residue_id, i, pairwise_energies):
                        # Mark rotamer as eliminated
                        flex_residue.rotamers[i]['eliminated'] = True
                        eliminated = True
    
    def _calculate_pairwise_energies(self, pose: Pose) -> Dict:
        """Calculate pairwise energies between all rotamer combinations"""
        pairwise_energies = {}
        
        residue_ids = list(self.flexible_residues.keys())
        
        for i, res_id1 in enumerate(residue_ids):
            for j, res_id2 in enumerate(residue_ids):
                if i <= j:
                    continue
                
                res1 = self.flexible_residues[res_id1]
                res2 = self.flexible_residues[res_id2]
                
                # Calculate energy for all rotamer combinations
                energies = np.zeros((len(res1.rotamers), len(res2.rotamers)))
                
                for r1_idx in range(len(res1.rotamers)):
                    for r2_idx in range(len(res2.rotamers)):
                        res1.set_rotamer(r1_idx)
                        res2.set_rotamer(r2_idx)
                        
                        coords1 = res1.get_current_coords()
                        coords2 = res2.get_current_coords()
                        
                        energy = self._calculate_sidechain_sidechain_energy(coords1, coords2)
                        energies[r1_idx, r2_idx] = energy
                
                pairwise_energies[(res_id1, res_id2)] = energies
        
        return pairwise_energies
    
    def _is_dead_end_rotamer(self, residue_id: str, rotamer_idx: int, pairwise_energies: Dict) -> bool:
        """Check if a rotamer is a dead-end (always worse than alternatives)"""
        # Simplified dead-end elimination
        # In reality, this would be more sophisticated
        
        flex_residue = self.flexible_residues[residue_id]
        
        # Check if rotamer has very high energy
        if len(flex_residue.rotamers) > rotamer_idx:
            rotamer = flex_residue.rotamers[rotamer_idx]
            if rotamer.get('eliminated', False):
                return True
        
        return False
    
    def _monte_carlo_optimization(self, pose: Pose) -> Pose:
        """
        Monte Carlo optimization of side-chain conformations
        """
        self.logger.debug("Performing Monte Carlo optimization")
        
        current_energy = self._calculate_total_flexible_energy(pose)
        best_energy = current_energy
        best_conformations = self._get_current_conformations()
        
        for step in range(self.monte_carlo_steps):
            # Propose move
            self._propose_rotamer_move()
            
            # Calculate new energy
            new_energy = self._calculate_total_flexible_energy(pose)
            
            # Accept or reject move
            if self._accept_move(current_energy, new_energy, step):
                current_energy = new_energy
                
                # Update best
                if new_energy < best_energy:
                    best_energy = new_energy
                    best_conformations = self._get_current_conformations()
            else:
                # Reject move
                self._revert_last_move()
        
        # Set best conformations
        self._set_conformations(best_conformations)
        
        return pose
    
    def _propose_rotamer_move(self):
        """Propose a random rotamer move"""
        # Select random flexible residue
        residue_ids = list(self.flexible_residues.keys())
        residue_id = np.random.choice(residue_ids)
        
        flex_residue = self.flexible_residues[residue_id]
        
        # Select random rotamer
        available_rotamers = [i for i, rot in enumerate(flex_residue.rotamers) 
                            if not rot.get('eliminated', False)]
        
        if available_rotamers:
            new_rotamer_idx = np.random.choice(available_rotamers)
            self.last_move = (residue_id, flex_residue.current_rotamer_idx, new_rotamer_idx)
            flex_residue.set_rotamer(new_rotamer_idx)
    
    def _accept_move(self, current_energy: float, new_energy: float, step: int) -> bool:
        """Accept or reject Monte Carlo move"""
        if new_energy < current_energy:
            return True
        
        # Calculate temperature for simulated annealing
        temperature = self.temperature * (1.0 - step / self.monte_carlo_steps)
        
        # Metropolis criterion
        if temperature > 0:
            probability = np.exp(-(new_energy - current_energy) / temperature)
            return np.random.random() < probability
        
        return False
    
    def _revert_last_move(self):
        """Revert last Monte Carlo move"""
        if hasattr(self, 'last_move'):
            residue_id, old_rotamer_idx, _ = self.last_move
            self.flexible_residues[residue_id].set_rotamer(old_rotamer_idx)
    
    def _get_current_conformations(self) -> Dict:
        """Get current side-chain conformations"""
        conformations = {}
        for residue_id, flex_residue in self.flexible_residues.items():
            conformations[residue_id] = flex_residue.current_rotamer_idx
        return conformations
    
    def _set_conformations(self, conformations: Dict):
        """Set side-chain conformations"""
        for residue_id, rotamer_idx in conformations.items():
            if residue_id in self.flexible_residues:
                self.flexible_residues[residue_id].set_rotamer(rotamer_idx)
    
    def _calculate_total_flexible_energy(self, pose: Pose) -> float:
        """Calculate total energy including flexible side-chains"""
        total_energy = 0.0
        
        # Ligand-sidechain interactions
        for flex_residue in self.flexible_residues.values():
            sidechain_coords = flex_residue.get_current_coords()
            energy = self._calculate_ligand_sidechain_energy(pose.coordinates, sidechain_coords)
            total_energy += energy * self.ligand_sidechain_weight
        
        # Sidechain-sidechain interactions
        residue_ids = list(self.flexible_residues.keys())
        for i, res_id1 in enumerate(residue_ids):
            for j, res_id2 in enumerate(residue_ids):
                if i < j:
                    coords1 = self.flexible_residues[res_id1].get_current_coords()
                    coords2 = self.flexible_residues[res_id2].get_current_coords()
                    energy = self._calculate_sidechain_sidechain_energy(coords1, coords2)
                    total_energy += energy * self.sidechain_sidechain_weight
        
        # Rotamer probabilities
        for flex_residue in self.flexible_residues.values():
            prob = flex_residue.get_rotamer_probability()
            if prob > 0:
                total_energy -= np.log(prob) * 0.1  # Entropy term
        
        return total_energy
    
    def _calculate_ligand_sidechain_energy(self, ligand_coords: np.ndarray, sidechain_coords: np.ndarray) -> float:
        """Calculate interaction energy between ligand and side-chain"""
        # Calculate distances
        distances = distance_matrix(ligand_coords, sidechain_coords)
        
        # Simple Lennard-Jones potential
        r_min = 3.0  # Angstroms
        energy = 0.0
        
        for i in range(distances.shape[0]):
            for j in range(distances.shape[1]):
                r = distances[i, j]
                if r < 10.0:  # Cutoff
                    # Lennard-Jones potential
                    lj_energy = 4.0 * ((r_min/r)**12 - (r_min/r)**6)
                    energy += lj_energy
        
        return energy
    
    def _calculate_sidechain_sidechain_energy(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Calculate interaction energy between two side-chains"""
        # Calculate distances
        distances = distance_matrix(coords1, coords2)
        
        # Simple repulsive potential
        energy = 0.0
        clash_threshold = 2.5  # Angstroms
        
        for i in range(distances.shape[0]):
            for j in range(distances.shape[1]):
                r = distances[i, j]
                if r < clash_threshold:
                    # Repulsive potential
                    energy += 100.0 * (clash_threshold - r)**2
        
        return energy
    
    def _minimize_with_sidechains(self, pose: Pose) -> Pose:
        """Final energy minimization with side-chains"""
        self.logger.debug("Performing final minimization with side-chains")
        
        # Extract current coordinates
        all_coords = [pose.coordinates]
        
        for flex_residue in self.flexible_residues.values():
            sidechain_coords = flex_residue.get_current_coords()
            all_coords.append(sidechain_coords)
        
        # Flatten coordinates
        initial_coords = np.vstack(all_coords).flatten()
        
        # Define energy function
        def energy_function(coords_flat):
            # Reshape coordinates
            coords_3d = coords_flat.reshape(-1, 3)
            
            # Split into ligand and side-chain coordinates
            ligand_coords = coords_3d[:len(pose.coordinates)]
            sidechain_start = len(pose.coordinates)
            
            energy = 0.0
            
            # Ligand-sidechain energies
            for flex_residue in self.flexible_residues.values():
                sidechain_len = len(flex_residue.get_current_coords())
                sidechain_coords = coords_3d[sidechain_start:sidechain_start + sidechain_len]
                energy += self._calculate_ligand_sidechain_energy(ligand_coords, sidechain_coords)
                sidechain_start += sidechain_len
            
            return energy
        
        # Minimize
        result = minimize(
            energy_function,
            initial_coords,
            method='L-BFGS-B',
            options={'maxiter': 100}
        )
        
        # Create optimized pose
        optimized_coords = result.x.reshape(-1, 3)
        ligand_coords = optimized_coords[:len(pose.coordinates)]
        
        optimized_pose = Pose(
            coordinates=ligand_coords,
            score=pose.score,
            energy=result.fun,
            ligand_name=pose.ligand_name,
            pose_id=pose.pose_id,
            flexible_residues=list(self.flexible_residues.keys())
        )
        
        return optimized_pose
    
    def get_flexible_residue_info(self) -> Dict[str, Any]:
        """Get information about flexible residues"""
        info = {}
        
        for residue_id, flex_residue in self.flexible_residues.items():
            info[residue_id] = {
                'type': flex_residue.residue_type,
                'num_rotamers': len(flex_residue.rotamers),
                'current_rotamer': flex_residue.current_rotamer_idx,
                'current_probability': flex_residue.get_rotamer_probability()
            }
        
        return info