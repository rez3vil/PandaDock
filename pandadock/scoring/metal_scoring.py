# -*- coding: utf-8 -*-
"""
Metal-Specific Scoring Functions

This module provides specialized scoring functions for metal-ligand interactions
in metalloproteins, including coordination energy, geometric constraints, and
electrostatic effects specific to metal centers.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass
import logging

from scoring.scoring_functions import ScoringFunctions
from docking.metal_coordination import MetalCenter, CoordinationGeometry, MetalType, MetalConstraintValidator

logger = logging.getLogger(__name__)


@dataclass
class MetalScoringParameters:
    """Parameters for metal-specific scoring"""
    # Coordination energy parameters
    coordination_strength: Dict[str, float]  # Element-specific coordination strengths
    geometric_penalty_weight: float = 2.0   # Weight for geometric violations
    distance_tolerance: float = 0.2          # Distance tolerance in Angstroms
    angle_tolerance: float = 10.0            # Angle tolerance in degrees
    
    # Electrostatic parameters
    metal_charge_scaling: float = 0.8        # Scale metal charge for softer interactions
    polarization_factor: float = 1.2         # Account for metal polarization effects
    charge_transfer_factor: float = 0.3      # Charge transfer contribution
    
    # Solvation parameters
    metal_desolvation_penalty: float = 2.0   # Penalty for metal desolvation
    coordination_solvation_bonus: float = 1.0 # Bonus for coordinated solvation
    
    # Entropy parameters
    coordination_entropy_penalty: float = 1.5  # Entropy loss upon coordination
    geometric_entropy_factor: float = 0.5     # Geometric constraint entropy
    
    def __post_init__(self):
        """Set default coordination strengths if not provided"""
        if not hasattr(self, 'coordination_strength') or not self.coordination_strength:
            self.coordination_strength = {
                'N': -4.0,  # Strong coordination (histidine, etc.)
                'O': -3.5,  # Moderate to strong (carboxylates, water)
                'S': -3.0,  # Moderate (cysteine, methionine)
                'P': -2.5,  # Weaker (phosphates)
                'C': -1.5,  # Weak (carbonyls, aromatics)
                'Cl': -2.0, # Halides
                'F': -2.5,
                'Br': -1.8,
                'I': -1.5
            }


class MetalScoringFunction(ScoringFunctions):
    """Enhanced scoring function for metal-containing systems"""
    
    def __init__(self, metal_centers: List[MetalCenter], 
                 parameters: Optional[MetalScoringParameters] = None):
        super().__init__()
        self.metal_centers = metal_centers
        self.parameters = parameters or MetalScoringParameters()
        self.validator = MetalConstraintValidator()
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Metal-specific interaction parameters
        self.metal_vdw_radii = self._initialize_metal_vdw_radii()
        self.metal_electronegativities = self._initialize_electronegativities()
        
    def _initialize_metal_vdw_radii(self) -> Dict[MetalType, float]:
        """Initialize van der Waals radii for metal ions"""
        return {
            MetalType.ZN: 1.39,
            MetalType.FE: 1.52,
            MetalType.MG: 1.73,
            MetalType.CA: 2.31,
            MetalType.MN: 1.61,
            MetalType.CU: 1.40,
            MetalType.NI: 1.24,
            MetalType.CO: 1.26,
            MetalType.MO: 1.54,
            MetalType.W: 1.62,
            MetalType.V: 1.53,
            MetalType.CR: 1.39
        }
    
    def _initialize_electronegativities(self) -> Dict[MetalType, float]:
        """Initialize Pauling electronegativities for metals"""
        return {
            MetalType.ZN: 1.65,
            MetalType.FE: 1.83,
            MetalType.MG: 1.31,
            MetalType.CA: 1.00,
            MetalType.MN: 1.55,
            MetalType.CU: 1.90,
            MetalType.NI: 1.91,
            MetalType.CO: 1.88,
            MetalType.MO: 2.16,
            MetalType.W: 2.36,
            MetalType.V: 1.63,
            MetalType.CR: 1.66
        }
    
    def calculate_total_energy(self, ligand_coordinates: np.ndarray, 
                             protein_coordinates: Optional[np.ndarray] = None,
                             ligand_atom_types: Optional[List[str]] = None) -> float:
        """Calculate total energy including metal-specific terms"""
        # Base scoring function energy
        base_energy = super().calculate_total_energy(ligand_coordinates, protein_coordinates, ligand_atom_types)
        
        # Metal-specific energy contributions
        metal_energy = 0.0
        
        for metal_center in self.metal_centers:
            # Coordination energy
            coord_energy = self.calculate_coordination_energy(
                metal_center, ligand_coordinates, ligand_atom_types
            )
            
            # Geometric penalty
            geom_penalty = self.calculate_geometric_penalty(
                metal_center, ligand_coordinates, ligand_atom_types
            )
            
            # Metal-specific electrostatics
            electrostatic_energy = self.calculate_metal_electrostatics(
                metal_center, ligand_coordinates, ligand_atom_types
            )
            
            # Metal solvation effects
            solvation_energy = self.calculate_metal_solvation(
                metal_center, ligand_coordinates, ligand_atom_types
            )
            
            # Coordination entropy
            entropy_penalty = self.calculate_coordination_entropy(
                metal_center, ligand_coordinates, ligand_atom_types
            )
            
            metal_energy += (coord_energy + geom_penalty + 
                           electrostatic_energy + solvation_energy + entropy_penalty)
        
        total_energy = base_energy + metal_energy
        
        self.logger.debug(f"Energy breakdown - Base: {base_energy:.2f}, "
                         f"Metal: {metal_energy:.2f}, Total: {total_energy:.2f}")
        
        return total_energy
    
    def calculate_coordination_energy(self, metal_center: MetalCenter, 
                                    ligand_coordinates: np.ndarray,
                                    ligand_atom_types: Optional[List[str]] = None) -> float:
        """Calculate energy from metal-ligand coordination bonds"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        coordination_energy = 0.0
        
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            
            # Check if atom can coordinate
            if atom_type in self.parameters.coordination_strength and distance < 4.0:
                # Calculate coordination strength based on distance and geometry
                strength = self._calculate_coordination_strength(
                    metal_center, distance, atom_type
                )
                
                # Apply distance-dependent energy function
                energy = self._coordination_energy_function(distance, strength, metal_center.metal_type)
                coordination_energy += energy
        
        return coordination_energy
    
    def _calculate_coordination_strength(self, metal_center: MetalCenter, 
                                       distance: float, atom_type: str) -> float:
        """Calculate coordination strength based on metal type and geometry"""
        base_strength = self.parameters.coordination_strength.get(atom_type, -1.0)
        
        # Modify strength based on metal type
        metal_factor = self._get_metal_coordination_factor(metal_center.metal_type)
        
        # Modify based on coordination number (more coordination = weaker individual bonds)
        coord_number_factor = 1.0 / np.sqrt(metal_center.coordination_number)
        
        # Modify based on oxidation state (higher charge = stronger binding)
        if metal_center.oxidation_state:
            charge_factor = 1.0 + 0.1 * metal_center.oxidation_state
        else:
            charge_factor = 1.0
        
        return base_strength * metal_factor * coord_number_factor * charge_factor
    
    def _get_metal_coordination_factor(self, metal_type: MetalType) -> float:
        """Get metal-specific coordination strength factors"""
        factors = {
            MetalType.ZN: 1.0,    # Reference
            MetalType.FE: 1.2,    # Strong field metal
            MetalType.MG: 0.8,    # Weaker than Zn
            MetalType.CA: 0.6,    # Much weaker
            MetalType.MN: 1.0,    # Similar to Zn
            MetalType.CU: 1.3,    # Strong coordination
            MetalType.NI: 1.2,    # Strong field
            MetalType.CO: 1.1,    # Moderate
            MetalType.MO: 1.4,    # Very strong
            MetalType.W: 1.4,     # Very strong
            MetalType.V: 1.1,     # Moderate
            MetalType.CR: 1.1     # Moderate
        }
        return factors.get(metal_type, 1.0)
    
    def _coordination_energy_function(self, distance: float, strength: float, 
                                    metal_type: MetalType) -> float:
        """Calculate distance-dependent coordination energy"""
        # Get optimal coordination distance
        optimal_distance = self._get_optimal_coordination_distance(metal_type)
        
        # Use a Morse-like potential for coordination
        # E = strength * (1 - exp(-alpha*(r - r0)))^2 - strength
        alpha = 2.0  # Controls the width of the potential well
        r0 = optimal_distance
        
        if distance > 4.0:  # Beyond coordination range
            return 0.0
        
        # Morse potential
        exp_term = np.exp(-alpha * (distance - r0))
        energy = strength * ((1 - exp_term) ** 2 - 1)
        
        # Add repulsive term for very short distances
        if distance < r0 * 0.8:
            repulsion = 10.0 * ((r0 * 0.8 / distance) ** 12)
            energy += repulsion
        
        return energy
    
    def _get_optimal_coordination_distance(self, metal_type: MetalType) -> float:
        """Get optimal coordination distance for metal type"""
        distances = {
            MetalType.ZN: 2.1,
            MetalType.FE: 2.0,
            MetalType.MG: 2.1,
            MetalType.CA: 2.4,
            MetalType.MN: 2.2,
            MetalType.CU: 2.0,
            MetalType.NI: 2.0,
            MetalType.CO: 2.0,
            MetalType.MO: 2.3,
            MetalType.W: 2.3,
            MetalType.V: 2.1,
            MetalType.CR: 2.1
        }
        return distances.get(metal_type, 2.1)
    
    def calculate_geometric_penalty(self, metal_center: MetalCenter,
                                  ligand_coordinates: np.ndarray,
                                  ligand_atom_types: Optional[List[str]] = None) -> float:
        """Calculate penalty for geometric constraint violations"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        # Find coordinating atoms in ligand
        ligand_atoms = []
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            if atom_type in ['N', 'O', 'S', 'P'] and distance < 3.5:
                ligand_atoms.append({
                    'coordinates': coord,
                    'element': atom_type,
                    'atom_id': i
                })
        
        if not ligand_atoms:
            return 0.0
        
        # Validate coordination geometry
        validation = self.validator.validate_coordination(metal_center, ligand_atoms)
        
        # Calculate penalty based on violations
        penalty = 0.0
        
        # Distance violations
        for violation in validation['distance_violations']:
            distance_error = min(
                abs(violation['actual_distance'] - violation['expected_range'][0]),
                abs(violation['actual_distance'] - violation['expected_range'][1])
            )
            penalty += self.parameters.geometric_penalty_weight * distance_error
        
        # Angle violations
        for violation in validation['angle_violations']:
            expected_angles = violation['expected_angles']
            actual_angle = violation['actual_angle']
            angle_error = min(abs(actual_angle - exp_angle) for exp_angle in expected_angles)
            if angle_error > self.parameters.angle_tolerance:
                penalty += self.parameters.geometric_penalty_weight * (angle_error / 180.0)
        
        # Geometric score penalty (lower score = higher penalty)
        if validation['geometric_score'] < 0.8:
            penalty += self.parameters.geometric_penalty_weight * (0.8 - validation['geometric_score'])
        
        return penalty
    
    def calculate_metal_electrostatics(self, metal_center: MetalCenter,
                                     ligand_coordinates: np.ndarray,
                                     ligand_atom_types: Optional[List[str]] = None,
                                     ligand_charges: Optional[np.ndarray] = None) -> float:
        """Calculate metal-specific electrostatic interactions"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        if ligand_charges is None:
            # Estimate partial charges based on atom types
            ligand_charges = self._estimate_partial_charges(ligand_atom_types)
        
        electrostatic_energy = 0.0
        
        # Metal charge (scaled for softer interactions)
        metal_charge = metal_center.charge * self.parameters.metal_charge_scaling
        
        for i, (coord, charge) in enumerate(zip(ligand_coordinates, ligand_charges)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            
            if distance > 0.1:  # Avoid division by zero
                # Coulomb interaction with distance-dependent dielectric
                dielectric = self._calculate_dielectric(distance, metal_center)
                coulomb_energy = 332.0 * metal_charge * charge / (dielectric * distance)
                
                # Add polarization effects
                polarization_energy = self._calculate_polarization_energy(
                    metal_center, coord, ligand_atom_types[i], distance
                )
                
                # Add charge transfer effects
                charge_transfer_energy = self._calculate_charge_transfer_energy(
                    metal_center, ligand_atom_types[i], distance
                )
                
                electrostatic_energy += (coulomb_energy + polarization_energy + 
                                       charge_transfer_energy)
        
        return electrostatic_energy
    
    def _estimate_partial_charges(self, atom_types: List[str]) -> np.ndarray:
        """Estimate partial charges for ligand atoms"""
        charge_map = {
            'C': 0.0,
            'N': -0.3,
            'O': -0.4,
            'S': -0.2,
            'P': 0.3,
            'F': -0.4,
            'Cl': -0.3,
            'Br': -0.2,
            'I': -0.1,
            'H': 0.1
        }
        
        return np.array([charge_map.get(atom_type, 0.0) for atom_type in atom_types])
    
    def _calculate_dielectric(self, distance: float, metal_center: MetalCenter) -> float:
        """Calculate distance-dependent dielectric constant"""
        # Near metal: low dielectric (tight binding)
        # Far from metal: higher dielectric (bulk solvent)
        min_dielectric = 4.0
        max_dielectric = 80.0
        transition_distance = 5.0
        
        if distance < transition_distance:
            # Sigmoidal transition
            ratio = distance / transition_distance
            dielectric = min_dielectric + (max_dielectric - min_dielectric) * (ratio ** 2)
        else:
            dielectric = max_dielectric
        
        return dielectric
    
    def _calculate_polarization_energy(self, metal_center: MetalCenter, 
                                     ligand_coord: np.ndarray, atom_type: str,
                                     distance: float) -> float:
        """Calculate metal-induced polarization energy"""
        if distance > 5.0:
            return 0.0
        
        # Metal polarizability (approximate values in Ų)
        metal_polarizability = {
            MetalType.ZN: 2.0,
            MetalType.FE: 2.5,
            MetalType.MG: 1.0,
            MetalType.CA: 3.0,
            MetalType.CU: 2.2,
            MetalType.MN: 2.3,
            MetalType.NI: 2.1,
            MetalType.CO: 2.0
        }.get(metal_center.metal_type, 2.0)
        
        # Ligand atom polarizability
        atom_polarizability = {
            'C': 1.5,
            'N': 1.0,
            'O': 0.8,
            'S': 2.5,
            'P': 3.0
        }.get(atom_type, 1.0)
        
        # Metal field strength
        field_strength = abs(metal_center.charge) / (distance ** 2)
        
        # Induced dipole interaction
        polarization_energy = -0.5 * self.parameters.polarization_factor * \
                             metal_polarizability * atom_polarizability * \
                             field_strength / (distance ** 4)
        
        return polarization_energy
    
    def _calculate_charge_transfer_energy(self, metal_center: MetalCenter,
                                        atom_type: str, distance: float) -> float:
        """Calculate charge transfer contribution"""
        if distance > 3.0:
            return 0.0
        
        # Only significant for direct coordination
        if atom_type not in ['N', 'O', 'S', 'P']:
            return 0.0
        
        # Electronegativity difference drives charge transfer
        metal_en = self.metal_electronegativities.get(metal_center.metal_type, 1.5)
        ligand_en = {'N': 3.04, 'O': 3.44, 'S': 2.58, 'P': 2.19}.get(atom_type, 2.5)
        
        en_diff = abs(ligand_en - metal_en)
        
        # Exponential decay with distance
        distance_factor = np.exp(-distance / 1.0)
        
        charge_transfer_energy = -self.parameters.charge_transfer_factor * \
                               en_diff * distance_factor
        
        return charge_transfer_energy
    
    def calculate_metal_solvation(self, metal_center: MetalCenter,
                                ligand_coordinates: np.ndarray,
                                ligand_atom_types: Optional[List[str]] = None) -> float:
        """Calculate metal solvation effects"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        solvation_energy = 0.0
        
        # Metal desolvation penalty (metal loses water coordination)
        coordinating_atoms = 0
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            if atom_type in ['N', 'O', 'S', 'P'] and distance < 3.0:
                coordinating_atoms += 1
        
        # Penalty for replacing water coordination
        if coordinating_atoms > 0:
            desolvation_penalty = self.parameters.metal_desolvation_penalty * coordinating_atoms
            solvation_energy += desolvation_penalty
            
            # Bonus for forming stable coordination (coordination sphere)
            coordination_bonus = -self.parameters.coordination_solvation_bonus * \
                               np.sqrt(coordinating_atoms)
            solvation_energy += coordination_bonus
        
        return solvation_energy
    
    def calculate_coordination_entropy(self, metal_center: MetalCenter,
                                     ligand_coordinates: np.ndarray,
                                     ligand_atom_types: Optional[List[str]] = None) -> float:
        """Calculate entropy effects of coordination"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        entropy_penalty = 0.0
        
        # Count coordinating atoms
        coordinating_atoms = 0
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            if atom_type in ['N', 'O', 'S', 'P'] and distance < 3.0:
                coordinating_atoms += 1
        
        if coordinating_atoms > 0:
            # Entropy loss from coordination (freezing degrees of freedom)
            coordination_entropy = self.parameters.coordination_entropy_penalty * coordinating_atoms
            
            # Additional penalty for geometric constraints
            if metal_center.geometry != CoordinationGeometry.TETRAHEDRAL:
                # More restrictive geometries have higher entropy penalties
                geometry_penalties = {
                    CoordinationGeometry.LINEAR: 2.0,
                    CoordinationGeometry.SQUARE_PLANAR: 1.5,
                    CoordinationGeometry.OCTAHEDRAL: 1.3,
                    CoordinationGeometry.TRIGONAL_PLANAR: 1.2,
                    CoordinationGeometry.TRIGONAL_BIPYRAMIDAL: 1.4
                }
                geometry_penalty = geometry_penalties.get(metal_center.geometry, 1.0)
                coordination_entropy *= geometry_penalty
            
            entropy_penalty += coordination_entropy
        
        return entropy_penalty
    
    def find_metal_interactions(self, metal_center: MetalCenter,
                              ligand_coordinates: np.ndarray,
                              ligand_atom_types: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """Find and classify metal-ligand interactions"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        interactions = []
        
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            distance = np.linalg.norm(coord - metal_center.coordinates)
            
            # Classify interaction type
            if distance < 3.5:  # Within coordination range
                interaction_type = self._classify_metal_interaction(
                    metal_center, atom_type, distance
                )
                
                if interaction_type != 'none':
                    interaction = {
                        'type': 'metal_coordination',
                        'subtype': interaction_type,
                        'metal_type': metal_center.metal_type.value,
                        'ligand_atom': i,
                        'ligand_element': atom_type,
                        'distance': distance,
                        'energy': self._coordination_energy_function(
                            distance, 
                            self._calculate_coordination_strength(metal_center, distance, atom_type),
                            metal_center.metal_type
                        ),
                        'geometry': metal_center.geometry.value,
                        'coordination_number': metal_center.coordination_number
                    }
                    interactions.append(interaction)
        
        return interactions
    
    def _classify_metal_interaction(self, metal_center: MetalCenter, 
                                  atom_type: str, distance: float) -> str:
        """Classify the type of metal-ligand interaction"""
        if atom_type in ['N', 'O', 'S', 'P']:
            if distance < 2.8:
                return 'direct_coordination'
            elif distance < 3.5:
                return 'weak_coordination'
        elif atom_type == 'C':
            if distance < 2.5:
                return 'organometallic'
            elif distance < 3.5:
                return 'weak_interaction'
        
        return 'none'
    
    def get_metal_scoring_report(self, metal_center: MetalCenter,
                               ligand_coordinates: np.ndarray,
                               ligand_atom_types: Optional[List[str]] = None) -> Dict[str, Any]:
        """Generate comprehensive metal scoring report"""
        if ligand_atom_types is None:
            ligand_atom_types = ['C'] * len(ligand_coordinates)
        
        # Calculate all energy components
        coord_energy = self.calculate_coordination_energy(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        geom_penalty = self.calculate_geometric_penalty(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        electrostatic_energy = self.calculate_metal_electrostatics(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        solvation_energy = self.calculate_metal_solvation(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        entropy_penalty = self.calculate_coordination_entropy(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        
        # Find interactions
        interactions = self.find_metal_interactions(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        
        # Validation
        ligand_atoms = []
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            ligand_atoms.append({
                'coordinates': coord,
                'element': atom_type,
                'atom_id': i
            })
        
        validation = self.validator.validate_coordination(metal_center, ligand_atoms)
        
        total_metal_energy = (coord_energy + geom_penalty + electrostatic_energy + 
                            solvation_energy + entropy_penalty)
        
        report = {
            'metal_center': {
                'type': metal_center.metal_type.value,
                'coordinates': metal_center.coordinates.tolist(),
                'geometry': metal_center.geometry.value,
                'coordination_number': metal_center.coordination_number,
                'oxidation_state': metal_center.oxidation_state,
                'charge': metal_center.charge
            },
            'energy_breakdown': {
                'coordination_energy': coord_energy,
                'geometric_penalty': geom_penalty,
                'electrostatic_energy': electrostatic_energy,
                'solvation_energy': solvation_energy,
                'entropy_penalty': entropy_penalty,
                'total_metal_energy': total_metal_energy
            },
            'interactions': interactions,
            'validation': validation,
            'coordination_quality': {
                'coordination_score': validation.get('coordination_score', 0.0),
                'geometric_score': validation.get('geometric_score', 0.0),
                'overall_validity': validation.get('valid', False),
                'violation_count': len(validation.get('violations', []))
            }
        }
        
        return report