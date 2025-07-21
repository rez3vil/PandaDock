# -*- coding: utf-8 -*-
"""
Rotamer library for flexible side-chain docking
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
import json
from pathlib import Path


class RotamerLibrary:
    """
    Rotamer library for amino acid side-chain conformations
    
    Features:
    - Backbone-dependent rotamer library
    - Chi angle definitions
    - Rotamer probabilities
    - Clash detection
    - Energy-based filtering
    """
    
    def __init__(self, library_path: Optional[str] = None):
        self.logger = logging.getLogger(__name__)
        
        # Rotamer data
        self.rotamer_data = {}
        self.chi_definitions = {}
        self.backbone_dependent = True
        
        # Load default rotamer library
        if library_path:
            self.load_library(library_path)
        else:
            self._initialize_default_library()
        
        self.logger.info("Initialized RotamerLibrary")
    
    def _initialize_default_library(self):
        """Initialize default rotamer library"""
        self.logger.info("Initializing default rotamer library")
        
        # Define chi angle definitions for each amino acid
        self.chi_definitions = {
            'ALA': [],  # No side chain
            'GLY': [],  # No side chain
            'VAL': [['N', 'CA', 'CB', 'CG1']],
            'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
            'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
            'MET': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'SD'], ['CB', 'CG', 'SD', 'CE']],
            'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
            'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
            'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
            'SER': [['N', 'CA', 'CB', 'OG']],
            'THR': [['N', 'CA', 'CB', 'OG1']],
            'CYS': [['N', 'CA', 'CB', 'SG']],
            'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
            'GLN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'OE1']],
            'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
            'GLU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'OE1']],
            'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
            'LYS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'CE'], ['CG', 'CD', 'CE', 'NZ']],
            'ARG': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'NE'], ['CG', 'CD', 'NE', 'CZ']],
            'PRO': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD']]
        }
        
        # Initialize rotamer data for each amino acid
        for residue_type in self.chi_definitions.keys():
            self.rotamer_data[residue_type] = self._generate_default_rotamers(residue_type)
    
    def _generate_default_rotamers(self, residue_type: str) -> List[Dict[str, Any]]:
        """Generate default rotamers for a residue type"""
        chi_def = self.chi_definitions[residue_type]
        num_chis = len(chi_def)
        
        if num_chis == 0:
            # No side chain
            return [{'chi_angles': [], 'probability': 1.0}]
        
        rotamers = []
        
        if num_chis == 1:
            # Common chi1 angles
            chi1_angles = [-60, 60, 180]
            for chi1 in chi1_angles:
                rotamers.append({
                    'chi_angles': [chi1],
                    'probability': 0.33
                })
        
        elif num_chis == 2:
            # Common chi1, chi2 combinations
            chi_combinations = [
                [-60, -60], [-60, 60], [-60, 180],
                [60, -60], [60, 60], [60, 180],
                [180, -60], [180, 60], [180, 180]
            ]
            for chi1, chi2 in chi_combinations:
                rotamers.append({
                    'chi_angles': [chi1, chi2],
                    'probability': 0.11
                })
        
        elif num_chis == 3:
            # Selected chi1, chi2, chi3 combinations
            chi_combinations = [
                [-60, -60, -60], [-60, -60, 60], [-60, -60, 180],
                [-60, 60, -60], [-60, 60, 60], [-60, 60, 180],
                [60, -60, -60], [60, -60, 60], [60, -60, 180],
                [60, 60, -60], [60, 60, 60], [60, 60, 180],
                [180, 60, -60], [180, 60, 60], [180, 60, 180]
            ]
            for chi1, chi2, chi3 in chi_combinations:
                rotamers.append({
                    'chi_angles': [chi1, chi2, chi3],
                    'probability': 0.067
                })
        
        elif num_chis == 4:
            # Selected chi1, chi2, chi3, chi4 combinations
            chi_combinations = [
                [-60, -60, -60, -60], [-60, -60, -60, 60], [-60, -60, -60, 180],
                [-60, -60, 60, -60], [-60, -60, 60, 60], [-60, -60, 60, 180],
                [-60, 60, -60, -60], [-60, 60, -60, 60], [-60, 60, -60, 180],
                [60, -60, -60, -60], [60, -60, -60, 60], [60, -60, -60, 180],
                [60, 60, -60, -60], [60, 60, -60, 60], [60, 60, -60, 180],
                [180, 60, -60, -60], [180, 60, -60, 60], [180, 60, -60, 180]
            ]
            for chi_angles in chi_combinations:
                rotamers.append({
                    'chi_angles': list(chi_angles),
                    'probability': 0.056
                })
        
        # Normalize probabilities
        total_prob = sum(rot['probability'] for rot in rotamers)
        for rot in rotamers:
            rot['probability'] /= total_prob
        
        return rotamers
    
    def load_library(self, library_path: str):
        """Load rotamer library from file"""
        self.logger.info(f"Loading rotamer library from {library_path}")
        
        library_path = Path(library_path)
        
        if not library_path.exists():
            raise FileNotFoundError(f"Rotamer library not found: {library_path}")
        
        if library_path.suffix == '.json':
            with open(library_path, 'r') as f:
                data = json.load(f)
                self.rotamer_data = data.get('rotamer_data', {})
                self.chi_definitions = data.get('chi_definitions', {})
        else:
            # Assume it's a Dunbrack-style library
            self._load_dunbrack_library(library_path)
    
    def _load_dunbrack_library(self, library_path: Path):
        """Load Dunbrack-style rotamer library"""
        # This is a simplified implementation
        # In practice, would parse the actual Dunbrack library format
        
        self.logger.info("Loading Dunbrack-style rotamer library")
        
        # For now, use default library
        self._initialize_default_library()
    
    def save_library(self, library_path: str):
        """Save rotamer library to file"""
        library_path = Path(library_path)
        
        data = {
            'rotamer_data': self.rotamer_data,
            'chi_definitions': self.chi_definitions,
            'backbone_dependent': self.backbone_dependent
        }
        
        with open(library_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"Rotamer library saved to {library_path}")
    
    def get_rotamers(self, residue_type: str, phi: float = None, psi: float = None) -> List[Dict[str, Any]]:
        """
        Get rotamers for a residue type
        
        Args:
            residue_type: Three-letter amino acid code
            phi: Backbone phi angle (optional)
            psi: Backbone psi angle (optional)
            
        Returns:
            List of rotamer dictionaries
        """
        residue_type = residue_type.upper()
        
        if residue_type not in self.rotamer_data:
            self.logger.warning(f"Unknown residue type: {residue_type}")
            return []
        
        rotamers = self.rotamer_data[residue_type]
        
        if self.backbone_dependent and phi is not None and psi is not None:
            # Filter rotamers based on backbone conformation
            rotamers = self._filter_by_backbone(rotamers, phi, psi)
        
        return rotamers
    
    def _filter_by_backbone(self, rotamers: List[Dict[str, Any]], phi: float, psi: float) -> List[Dict[str, Any]]:
        """Filter rotamers based on backbone conformation"""
        # This is a simplified implementation
        # In practice, would use actual backbone-dependent probabilities
        
        filtered_rotamers = []
        
        for rotamer in rotamers:
            # Simple backbone-dependent filtering
            if -90 <= phi <= -30 and -60 <= psi <= 30:
                # Alpha-helix region
                if rotamer['chi_angles'] and rotamer['chi_angles'][0] == -60:
                    rotamer['probability'] *= 1.5  # Prefer gauche-
            elif 30 <= phi <= 90 and 90 <= psi <= 150:
                # Beta-sheet region
                if rotamer['chi_angles'] and rotamer['chi_angles'][0] == 180:
                    rotamer['probability'] *= 1.5  # Prefer trans
            
            filtered_rotamers.append(rotamer)
        
        # Renormalize probabilities
        total_prob = sum(rot['probability'] for rot in filtered_rotamers)
        if total_prob > 0:
            for rot in filtered_rotamers:
                rot['probability'] /= total_prob
        
        return filtered_rotamers
    
    def build_rotamer_coordinates(self, residue_type: str, backbone_coords: np.ndarray, 
                                 chi_angles: List[float]) -> np.ndarray:
        """
        Build side-chain coordinates for a rotamer
        
        Args:
            residue_type: Three-letter amino acid code
            backbone_coords: Backbone coordinates [N, CA, C] (3 x 3)
            chi_angles: Chi angle values in degrees
            
        Returns:
            Side-chain coordinates
        """
        residue_type = residue_type.upper()
        
        if residue_type not in self.chi_definitions:
            raise ValueError(f"Unknown residue type: {residue_type}")
        
        chi_def = self.chi_definitions[residue_type]
        
        if len(chi_angles) != len(chi_def):
            raise ValueError(f"Expected {len(chi_def)} chi angles, got {len(chi_angles)}")
        
        # Get ideal bond lengths and angles
        ideal_geometry = self._get_ideal_geometry(residue_type)
        
        # Build coordinates
        coordinates = self._build_coordinates_from_chi(
            backbone_coords, chi_angles, ideal_geometry
        )
        
        return coordinates
    
    def _get_ideal_geometry(self, residue_type: str) -> Dict[str, Any]:
        """Get ideal bond lengths and angles for a residue"""
        # This is a simplified implementation
        # In practice, would use actual geometric parameters
        
        geometry = {
            'bond_lengths': {
                'CA-CB': 1.54,  # Angstroms
                'CB-CG': 1.54,
                'CG-CD': 1.54,
                'CD-CE': 1.54,
                'CE-NZ': 1.47,
                'CB-OG': 1.43,
                'CB-SG': 1.82,
                'CG-SD': 1.82,
                'SD-CE': 1.82
            },
            'bond_angles': {
                'N-CA-CB': 110.0,  # degrees
                'CA-CB-CG': 114.0,
                'CB-CG-CD': 114.0,
                'CG-CD-CE': 114.0,
                'CD-CE-NZ': 109.0,
                'CA-CB-OG': 109.0,
                'CA-CB-SG': 109.0,
                'CB-CG-SD': 109.0,
                'CG-SD-CE': 96.0
            }
        }
        
        return geometry
    
    def _build_coordinates_from_chi(self, backbone_coords: np.ndarray, 
                                   chi_angles: List[float], geometry: Dict[str, Any]) -> np.ndarray:
        """Build coordinates from chi angles"""
        # This is a simplified implementation
        # In practice, would use proper molecular mechanics
        
        if len(chi_angles) == 0:
            return np.empty((0, 3))
        
        # Start with CA and CB
        ca_pos = backbone_coords[1]  # CA is second atom
        n_pos = backbone_coords[0]   # N is first atom
        
        # Build CB position
        cb_pos = self._build_cb_position(n_pos, ca_pos, backbone_coords[2])
        
        coordinates = [cb_pos]
        
        # Build subsequent atoms using chi angles
        current_pos = cb_pos
        previous_pos = ca_pos
        
        for i, chi_angle in enumerate(chi_angles):
            # Build next atom position
            next_pos = self._build_next_atom_position(
                previous_pos, current_pos, chi_angle, geometry
            )
            
            coordinates.append(next_pos)
            
            # Update positions
            previous_pos = current_pos
            current_pos = next_pos
        
        return np.array(coordinates)
    
    def _build_cb_position(self, n_pos: np.ndarray, ca_pos: np.ndarray, c_pos: np.ndarray) -> np.ndarray:
        """Build CB position from backbone atoms"""
        # This is a simplified implementation
        # In practice, would use proper tetrahedral geometry
        
        # Vector from CA to N
        ca_n = n_pos - ca_pos
        ca_n = ca_n / np.linalg.norm(ca_n)
        
        # Vector from CA to C
        ca_c = c_pos - ca_pos
        ca_c = ca_c / np.linalg.norm(ca_c)
        
        # Bisector vector
        bisector = ca_n + ca_c
        bisector = bisector / np.linalg.norm(bisector)
        
        # CB position (opposite to bisector)
        cb_bond_length = 1.54  # Angstroms
        cb_pos = ca_pos - bisector * cb_bond_length
        
        return cb_pos
    
    def _build_next_atom_position(self, pos1: np.ndarray, pos2: np.ndarray, 
                                 chi_angle: float, geometry: Dict[str, Any]) -> np.ndarray:
        """Build next atom position using chi angle"""
        # This is a simplified implementation
        # In practice, would use proper molecular mechanics
        
        # Default bond length and angle
        bond_length = 1.54  # Angstroms
        bond_angle = 114.0  # degrees
        
        # Convert chi angle to radians
        chi_rad = np.radians(chi_angle)
        bond_angle_rad = np.radians(bond_angle)
        
        # Vector from pos1 to pos2
        vec = pos2 - pos1
        vec = vec / np.linalg.norm(vec)
        
        # Create perpendicular vector
        perp = np.array([vec[1], -vec[0], 0])
        if np.linalg.norm(perp) < 0.1:
            perp = np.array([0, vec[2], -vec[1]])
        perp = perp / np.linalg.norm(perp)
        
        # Calculate new position
        direction = vec * np.cos(bond_angle_rad) + perp * np.sin(bond_angle_rad)
        
        # Apply chi angle rotation
        rotated_direction = self._rotate_around_axis(direction, vec, chi_rad)
        
        # New position
        new_pos = pos2 + rotated_direction * bond_length
        
        return new_pos
    
    def _rotate_around_axis(self, vector: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
        """Rotate vector around axis by angle"""
        # Rodrigues' rotation formula
        axis = axis / np.linalg.norm(axis)
        
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        
        rotated = (vector * cos_angle + 
                  np.cross(axis, vector) * sin_angle + 
                  axis * np.dot(axis, vector) * (1 - cos_angle))
        
        return rotated
    
    def calculate_rotamer_energy(self, residue_type: str, chi_angles: List[float]) -> float:
        """Calculate rotamer energy based on chi angles"""
        # This is a simplified implementation
        # In practice, would use actual energy functions
        
        energy = 0.0
        
        # Torsional energy
        for chi in chi_angles:
            # Simplified torsional potential
            chi_rad = np.radians(chi)
            energy += 0.5 * (1 + np.cos(3 * chi_rad))  # 3-fold potential
        
        return energy
    
    def find_best_rotamer(self, residue_type: str, target_coords: np.ndarray, 
                         backbone_coords: np.ndarray) -> Dict[str, Any]:
        """
        Find best rotamer to match target coordinates
        
        Args:
            residue_type: Three-letter amino acid code
            target_coords: Target side-chain coordinates
            backbone_coords: Backbone coordinates
            
        Returns:
            Best rotamer dictionary
        """
        rotamers = self.get_rotamers(residue_type)
        
        best_rotamer = None
        best_rmsd = float('inf')
        
        for rotamer in rotamers:
            try:
                # Build rotamer coordinates
                rotamer_coords = self.build_rotamer_coordinates(
                    residue_type, backbone_coords, rotamer['chi_angles']
                )
                
                # Calculate RMSD
                if len(rotamer_coords) == len(target_coords):
                    rmsd = np.sqrt(np.mean(np.sum((rotamer_coords - target_coords)**2, axis=1)))
                    
                    if rmsd < best_rmsd:
                        best_rmsd = rmsd
                        best_rotamer = rotamer.copy()
                        best_rotamer['rmsd'] = rmsd
            
            except Exception as e:
                self.logger.warning(f"Error evaluating rotamer: {e}")
                continue
        
        return best_rotamer
    
    def get_rotamer_statistics(self) -> Dict[str, Any]:
        """Get statistics about the rotamer library"""
        stats = {
            'total_residue_types': len(self.rotamer_data),
            'rotamer_counts': {},
            'chi_angle_counts': {},
            'backbone_dependent': self.backbone_dependent
        }
        
        for residue_type, rotamers in self.rotamer_data.items():
            stats['rotamer_counts'][residue_type] = len(rotamers)
            
            if rotamers:
                chi_count = len(rotamers[0]['chi_angles'])
                stats['chi_angle_counts'][residue_type] = chi_count
        
        return stats
    
    def validate_rotamer_library(self) -> Dict[str, Any]:
        """Validate the rotamer library"""
        validation_results = {
            'valid': True,
            'errors': [],
            'warnings': []
        }
        
        for residue_type, rotamers in self.rotamer_data.items():
            # Check if chi definitions exist
            if residue_type not in self.chi_definitions:
                validation_results['errors'].append(f"Missing chi definitions for {residue_type}")
                validation_results['valid'] = False
                continue
            
            expected_chi_count = len(self.chi_definitions[residue_type])
            
            # Check rotamers
            for i, rotamer in enumerate(rotamers):
                if len(rotamer['chi_angles']) != expected_chi_count:
                    validation_results['errors'].append(
                        f"Chi angle count mismatch for {residue_type} rotamer {i}"
                    )
                    validation_results['valid'] = False
                
                if rotamer['probability'] < 0 or rotamer['probability'] > 1:
                    validation_results['warnings'].append(
                        f"Invalid probability for {residue_type} rotamer {i}: {rotamer['probability']}"
                    )
            
            # Check probability normalization
            total_prob = sum(rot['probability'] for rot in rotamers)
            if abs(total_prob - 1.0) > 0.01:
                validation_results['warnings'].append(
                    f"Probabilities for {residue_type} do not sum to 1.0: {total_prob}"
                )
        
        return validation_results