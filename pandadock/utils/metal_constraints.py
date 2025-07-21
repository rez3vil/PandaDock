# -*- coding: utf-8 -*-
"""
Metal Constraint Handling Utilities

This module provides utilities for handling geometric and chemical constraints
specific to metal coordination in metalloproteins during docking.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from enum import Enum
import logging

from ..docking.metal_coordination import MetalCenter, CoordinationGeometry, MetalType

logger = logging.getLogger(__name__)


class ConstraintType(Enum):
    """Types of metal coordination constraints"""
    DISTANCE = "distance"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"
    COORDINATION_NUMBER = "coordination_number"
    GEOMETRY = "geometry"
    OXIDATION_STATE = "oxidation_state"
    CHARGE_NEUTRALITY = "charge_neutrality"


@dataclass
class MetalConstraint:
    """Individual metal coordination constraint"""
    constraint_type: ConstraintType
    target_value: float
    tolerance: float
    weight: float = 1.0
    atoms_involved: List[int] = field(default_factory=list)
    metal_center_id: Optional[int] = None
    description: str = ""
    
    def __post_init__(self):
        if not self.description:
            self.description = f"{self.constraint_type.value} constraint: {self.target_value} ± {self.tolerance}"


@dataclass
class ConstraintViolation:
    """Represents a constraint violation"""
    constraint: MetalConstraint
    actual_value: float
    violation_magnitude: float
    penalty: float
    atoms_involved: List[int]
    description: str = ""
    
    def __post_init__(self):
        if not self.description:
            self.description = (
                f"{self.constraint.constraint_type.value} violation: "
                f"expected {self.constraint.target_value:.2f} ± {self.constraint.tolerance:.2f}, "
                f"got {self.actual_value:.2f}"
            )


class MetalConstraintManager:
    """Manages and enforces metal coordination constraints"""
    
    def __init__(self, metal_centers: List[MetalCenter]):
        self.metal_centers = metal_centers
        self.constraints: List[MetalConstraint] = []
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Automatically generate standard constraints
        self._generate_standard_constraints()
    
    def _generate_standard_constraints(self):
        """Generate standard constraints for all metal centers"""
        for i, metal_center in enumerate(self.metal_centers):
            self._add_distance_constraints(metal_center, i)
            self._add_angle_constraints(metal_center, i)
            self._add_coordination_number_constraints(metal_center, i)
            self._add_geometry_constraints(metal_center, i)
    
    def _add_distance_constraints(self, metal_center: MetalCenter, center_id: int):
        """Add distance constraints for metal center"""
        # Get typical distances for this metal type
        distance_ranges = metal_center.distance_constraints
        
        for element, (min_dist, max_dist) in distance_ranges.items():
            ideal_distance = (min_dist + max_dist) / 2
            tolerance = (max_dist - min_dist) / 2
            
            constraint = MetalConstraint(
                constraint_type=ConstraintType.DISTANCE,
                target_value=ideal_distance,
                tolerance=tolerance,
                weight=2.0,  # Distance constraints are important
                metal_center_id=center_id,
                description=f"Metal-{element} distance constraint"
            )
            self.constraints.append(constraint)
    
    def _add_angle_constraints(self, metal_center: MetalCenter, center_id: int):
        """Add angle constraints based on coordination geometry"""
        geometry = metal_center.geometry
        expected_angles = self._get_ideal_angles(geometry)
        
        for angle in expected_angles:
            constraint = MetalConstraint(
                constraint_type=ConstraintType.ANGLE,
                target_value=angle,
                tolerance=15.0,  # 15 degree tolerance
                weight=1.5,
                metal_center_id=center_id,
                description=f"{geometry.value} L-M-L angle constraint"
            )
            self.constraints.append(constraint)
    
    def _add_coordination_number_constraints(self, metal_center: MetalCenter, center_id: int):
        """Add coordination number constraints"""
        target_coord_number = metal_center.coordination_number
        
        constraint = MetalConstraint(
            constraint_type=ConstraintType.COORDINATION_NUMBER,
            target_value=float(target_coord_number),
            tolerance=1.0,  # Allow ±1 coordination
            weight=3.0,  # Very important constraint
            metal_center_id=center_id,
            description=f"Coordination number constraint: {target_coord_number}"
        )
        self.constraints.append(constraint)
    
    def _add_geometry_constraints(self, metal_center: MetalCenter, center_id: int):
        """Add geometric constraints for coordination geometry"""
        # This is more of a composite constraint that checks overall geometry
        constraint = MetalConstraint(
            constraint_type=ConstraintType.GEOMETRY,
            target_value=1.0,  # Perfect geometry match
            tolerance=0.3,     # Allow some deviation
            weight=2.0,
            metal_center_id=center_id,
            description=f"{metal_center.geometry.value} geometry constraint"
        )
        self.constraints.append(constraint)
    
    def _get_ideal_angles(self, geometry: CoordinationGeometry) -> List[float]:
        """Get ideal angles for coordination geometry"""
        angle_map = {
            CoordinationGeometry.LINEAR: [180.0],
            CoordinationGeometry.TRIGONAL_PLANAR: [120.0],
            CoordinationGeometry.TETRAHEDRAL: [109.47],
            CoordinationGeometry.SQUARE_PLANAR: [90.0, 180.0],
            CoordinationGeometry.OCTAHEDRAL: [90.0, 180.0],
            CoordinationGeometry.TRIGONAL_BIPYRAMIDAL: [90.0, 120.0, 180.0],
            CoordinationGeometry.SQUARE_PYRAMIDAL: [90.0, 180.0],
            CoordinationGeometry.PENTAGONAL_BIPYRAMIDAL: [72.0, 90.0, 180.0]
        }
        return angle_map.get(geometry, [109.47])
    
    def add_custom_constraint(self, constraint: MetalConstraint):
        """Add a custom constraint"""
        self.constraints.append(constraint)
        self.logger.info(f"Added custom constraint: {constraint.description}")
    
    def evaluate_constraints(self, ligand_coordinates: np.ndarray,
                           ligand_atom_types: List[str]) -> Tuple[List[ConstraintViolation], float]:
        """Evaluate all constraints and return violations and total penalty"""
        violations = []
        total_penalty = 0.0
        
        for constraint in self.constraints:
            try:
                violation = self._evaluate_single_constraint(
                    constraint, ligand_coordinates, ligand_atom_types
                )
                if violation:
                    violations.append(violation)
                    total_penalty += violation.penalty
            except Exception as e:
                self.logger.warning(f"Error evaluating constraint {constraint.description}: {e}")
        
        self.logger.debug(f"Evaluated {len(self.constraints)} constraints, "
                         f"found {len(violations)} violations, "
                         f"total penalty: {total_penalty:.3f}")
        
        return violations, total_penalty
    
    def _evaluate_single_constraint(self, constraint: MetalConstraint,
                                   ligand_coordinates: np.ndarray,
                                   ligand_atom_types: List[str]) -> Optional[ConstraintViolation]:
        """Evaluate a single constraint"""
        
        if constraint.constraint_type == ConstraintType.DISTANCE:
            return self._evaluate_distance_constraint(constraint, ligand_coordinates, ligand_atom_types)
        elif constraint.constraint_type == ConstraintType.ANGLE:
            return self._evaluate_angle_constraint(constraint, ligand_coordinates, ligand_atom_types)
        elif constraint.constraint_type == ConstraintType.COORDINATION_NUMBER:
            return self._evaluate_coordination_number_constraint(constraint, ligand_coordinates, ligand_atom_types)
        elif constraint.constraint_type == ConstraintType.GEOMETRY:
            return self._evaluate_geometry_constraint(constraint, ligand_coordinates, ligand_atom_types)
        else:
            self.logger.warning(f"Unknown constraint type: {constraint.constraint_type}")
            return None
    
    def _evaluate_distance_constraint(self, constraint: MetalConstraint,
                                    ligand_coordinates: np.ndarray,
                                    ligand_atom_types: List[str]) -> Optional[ConstraintViolation]:
        """Evaluate distance constraints"""
        metal_center = self.metal_centers[constraint.metal_center_id]
        metal_pos = metal_center.coordinates
        
        # Find coordinating atoms
        coordinating_atoms = []
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            if atom_type in ['N', 'O', 'S', 'P']:  # Potential coordinators
                distance = np.linalg.norm(coord - metal_pos)
                if distance < 4.0:  # Within coordination range
                    coordinating_atoms.append((i, distance, atom_type))
        
        if not coordinating_atoms:
            return None  # No coordinating atoms to check
        
        # Check each coordinating atom
        violations = []
        for atom_idx, distance, atom_type in coordinating_atoms:
            # Check if this atom type matches the constraint
            expected_distance = constraint.target_value
            tolerance = constraint.tolerance
            
            if abs(distance - expected_distance) > tolerance:
                violation_magnitude = abs(distance - expected_distance) - tolerance
                penalty = constraint.weight * violation_magnitude
                
                violation = ConstraintViolation(
                    constraint=constraint,
                    actual_value=distance,
                    violation_magnitude=violation_magnitude,
                    penalty=penalty,
                    atoms_involved=[atom_idx]
                )
                violations.append(violation)
        
        # Return the worst violation
        if violations:
            return max(violations, key=lambda v: v.penalty)
        return None
    
    def _evaluate_angle_constraint(self, constraint: MetalConstraint,
                                 ligand_coordinates: np.ndarray,
                                 ligand_atom_types: List[str]) -> Optional[ConstraintViolation]:
        """Evaluate angle constraints"""
        metal_center = self.metal_centers[constraint.metal_center_id]
        metal_pos = metal_center.coordinates
        
        # Find coordinating atoms
        coordinating_atoms = []
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            if atom_type in ['N', 'O', 'S', 'P']:
                distance = np.linalg.norm(coord - metal_pos)
                if distance < 3.5:
                    coordinating_atoms.append((i, coord))
        
        if len(coordinating_atoms) < 2:
            return None  # Need at least 2 atoms for angle
        
        # Calculate all L-M-L angles
        worst_violation = None
        expected_angle = constraint.target_value
        tolerance = constraint.tolerance
        
        for i in range(len(coordinating_atoms)):
            for j in range(i + 1, len(coordinating_atoms)):
                idx1, coord1 = coordinating_atoms[i]
                idx2, coord2 = coordinating_atoms[j]
                
                # Calculate angle
                v1 = coord1 - metal_pos
                v2 = coord2 - metal_pos
                angle = self._calculate_angle(v1, v2)
                
                # Check against expected angle
                angle_diff = min(
                    abs(angle - expected_angle),
                    abs(angle - (360 - expected_angle))  # Handle wraparound
                )
                
                if angle_diff > tolerance:
                    violation_magnitude = angle_diff - tolerance
                    penalty = constraint.weight * violation_magnitude / 180.0  # Normalize by 180°
                    
                    violation = ConstraintViolation(
                        constraint=constraint,
                        actual_value=angle,
                        violation_magnitude=violation_magnitude,
                        penalty=penalty,
                        atoms_involved=[idx1, idx2]
                    )
                    
                    if worst_violation is None or violation.penalty > worst_violation.penalty:
                        worst_violation = violation
        
        return worst_violation
    
    def _evaluate_coordination_number_constraint(self, constraint: MetalConstraint,
                                               ligand_coordinates: np.ndarray,
                                               ligand_atom_types: List[str]) -> Optional[ConstraintViolation]:
        """Evaluate coordination number constraints"""
        metal_center = self.metal_centers[constraint.metal_center_id]
        metal_pos = metal_center.coordinates
        
        # Count coordinating atoms from ligand
        ligand_coordinators = 0
        coordinating_atom_indices = []
        
        for i, (coord, atom_type) in enumerate(zip(ligand_coordinates, ligand_atom_types)):
            if atom_type in ['N', 'O', 'S', 'P']:
                distance = np.linalg.norm(coord - metal_pos)
                if distance < 3.0:  # Direct coordination
                    ligand_coordinators += 1
                    coordinating_atom_indices.append(i)
        
        # Total coordination (protein + ligand)
        protein_coordinators = len(metal_center.coordinating_atoms)
        total_coordination = protein_coordinators + ligand_coordinators
        
        expected_coordination = constraint.target_value
        tolerance = constraint.tolerance
        
        if abs(total_coordination - expected_coordination) > tolerance:
            violation_magnitude = abs(total_coordination - expected_coordination) - tolerance
            penalty = constraint.weight * violation_magnitude
            
            return ConstraintViolation(
                constraint=constraint,
                actual_value=float(total_coordination),
                violation_magnitude=violation_magnitude,
                penalty=penalty,
                atoms_involved=coordinating_atom_indices,
                description=f"Coordination number violation: expected {expected_coordination}, got {total_coordination}"
            )
        
        return None
    
    def _evaluate_geometry_constraint(self, constraint: MetalConstraint,
                                    ligand_coordinates: np.ndarray,
                                    ligand_atom_types: List[str]) -> Optional[ConstraintViolation]:
        """Evaluate overall geometry constraints"""
        metal_center = self.metal_centers[constraint.metal_center_id]
        
        # This is a composite constraint that evaluates how well the
        # overall coordination geometry matches the expected pattern
        
        geometry_score = self._calculate_geometry_match_score(
            metal_center, ligand_coordinates, ligand_atom_types
        )
        
        expected_score = constraint.target_value
        tolerance = constraint.tolerance
        
        if geometry_score < (expected_score - tolerance):
            violation_magnitude = (expected_score - tolerance) - geometry_score
            penalty = constraint.weight * violation_magnitude
            
            return ConstraintViolation(
                constraint=constraint,
                actual_value=geometry_score,
                violation_magnitude=violation_magnitude,
                penalty=penalty,
                atoms_involved=[],  # Involves all coordinating atoms
                description=f"Geometry constraint violation: score {geometry_score:.3f} < {expected_score - tolerance:.3f}"
            )
        
        return None
    
    def _calculate_geometry_match_score(self, metal_center: MetalCenter,
                                      ligand_coordinates: np.ndarray,
                                      ligand_atom_types: List[str]) -> float:
        """Calculate how well the geometry matches the expected pattern"""
        metal_pos = metal_center.coordinates
        expected_geometry = metal_center.geometry
        
        # Find all coordinating atoms (protein + ligand)
        all_coordinators = []
        
        # Add protein coordinators
        for coord_atom in metal_center.coordinating_atoms:
            all_coordinators.append(coord_atom['coordinates'])
        
        # Add ligand coordinators
        for coord, atom_type in zip(ligand_coordinates, ligand_atom_types):
            if atom_type in ['N', 'O', 'S', 'P']:
                distance = np.linalg.norm(coord - metal_pos)
                if distance < 3.5:
                    all_coordinators.append(coord)
        
        if len(all_coordinators) < 2:
            return 0.0  # Cannot evaluate geometry with < 2 coordinators
        
        # Calculate actual angles
        actual_angles = []
        for i in range(len(all_coordinators)):
            for j in range(i + 1, len(all_coordinators)):
                v1 = all_coordinators[i] - metal_pos
                v2 = all_coordinators[j] - metal_pos
                angle = self._calculate_angle(v1, v2)
                actual_angles.append(angle)
        
        # Compare with ideal angles for this geometry
        ideal_angles = self._get_ideal_angles(expected_geometry)
        
        # Calculate RMS deviation from ideal angles
        if not ideal_angles:
            return 0.5  # Default score for unknown geometry
        
        min_rmsd = float('inf')
        
        # Try to match actual angles with ideal angles
        for ideal_angle in ideal_angles:
            deviations = []
            for actual_angle in actual_angles:
                min_dev = min(
                    abs(actual_angle - ideal_angle),
                    abs(actual_angle - (360 - ideal_angle))
                )
                deviations.append(min_dev)
            
            if deviations:
                rmsd = np.sqrt(np.mean(np.array(deviations) ** 2))
                min_rmsd = min(min_rmsd, rmsd)
        
        # Convert RMSD to score (lower RMSD = higher score)
        if min_rmsd == float('inf'):
            return 0.0
        
        # Score between 0 and 1, where 1 is perfect match
        max_acceptable_rmsd = 30.0  # degrees
        score = max(0.0, 1.0 - (min_rmsd / max_acceptable_rmsd))
        
        return score
    
    def _calculate_angle(self, v1: np.ndarray, v2: np.ndarray) -> float:
        """Calculate angle between two vectors in degrees"""
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
    
    def get_constraint_summary(self) -> Dict[str, Any]:
        """Get summary of all constraints"""
        constraint_counts = {}
        for constraint in self.constraints:
            constraint_type = constraint.constraint_type.value
            constraint_counts[constraint_type] = constraint_counts.get(constraint_type, 0) + 1
        
        return {
            'total_constraints': len(self.constraints),
            'constraint_types': constraint_counts,
            'metal_centers': len(self.metal_centers),
            'constraints_per_metal': len(self.constraints) / len(self.metal_centers) if self.metal_centers else 0
        }
    
    def optimize_pose_with_constraints(self, ligand_coordinates: np.ndarray,
                                     ligand_atom_types: List[str],
                                     max_iterations: int = 100,
                                     step_size: float = 0.1) -> Tuple[np.ndarray, float]:
        """Optimize ligand pose to satisfy constraints"""
        
        best_coords = ligand_coordinates.copy()
        best_penalty = self.evaluate_constraints(ligand_coordinates, ligand_atom_types)[1]
        
        current_coords = ligand_coordinates.copy()
        
        for iteration in range(max_iterations):
            # Calculate gradients for constraint violations
            gradients = self._calculate_constraint_gradients(current_coords, ligand_atom_types)
            
            # Apply gradient descent step
            current_coords -= step_size * gradients
            
            # Evaluate new position
            violations, penalty = self.evaluate_constraints(current_coords, ligand_atom_types)
            
            if penalty < best_penalty:
                best_coords = current_coords.copy()
                best_penalty = penalty
                
                # Early termination if constraints are satisfied
                if penalty < 0.01:
                    break
            
            # Adaptive step size
            if iteration % 10 == 0 and penalty > best_penalty * 1.1:
                step_size *= 0.9  # Reduce step size if not improving
        
        self.logger.debug(f"Constraint optimization: {max_iterations} iterations, "
                         f"penalty {best_penalty:.3f}")
        
        return best_coords, best_penalty
    
    def _calculate_constraint_gradients(self, ligand_coordinates: np.ndarray,
                                      ligand_atom_types: List[str]) -> np.ndarray:
        """Calculate gradients for constraint optimization"""
        gradients = np.zeros_like(ligand_coordinates)
        epsilon = 0.01  # Small step for numerical gradients
        
        # Calculate numerical gradients
        for i in range(len(ligand_coordinates)):
            for j in range(3):  # x, y, z components
                # Positive step
                coords_plus = ligand_coordinates.copy()
                coords_plus[i, j] += epsilon
                _, penalty_plus = self.evaluate_constraints(coords_plus, ligand_atom_types)
                
                # Negative step
                coords_minus = ligand_coordinates.copy()
                coords_minus[i, j] -= epsilon
                _, penalty_minus = self.evaluate_constraints(coords_minus, ligand_atom_types)
                
                # Numerical gradient
                gradients[i, j] = (penalty_plus - penalty_minus) / (2 * epsilon)
        
        return gradients


class ConstraintSetPresets:
    """Predefined constraint sets for common scenarios"""
    
    @staticmethod
    def create_strict_coordination_constraints(metal_centers: List[MetalCenter]) -> MetalConstraintManager:
        """Create strict constraints for precise coordination geometry"""
        manager = MetalConstraintManager(metal_centers)
        
        # Increase weights for strict enforcement
        for constraint in manager.constraints:
            if constraint.constraint_type == ConstraintType.DISTANCE:
                constraint.weight = 5.0
                constraint.tolerance *= 0.5  # Tighter tolerance
            elif constraint.constraint_type == ConstraintType.ANGLE:
                constraint.weight = 3.0
                constraint.tolerance *= 0.7
            elif constraint.constraint_type == ConstraintType.COORDINATION_NUMBER:
                constraint.weight = 10.0
                constraint.tolerance = 0.5  # Very strict
        
        return manager
    
    @staticmethod
    def create_flexible_coordination_constraints(metal_centers: List[MetalCenter]) -> MetalConstraintManager:
        """Create flexible constraints allowing more deviation"""
        manager = MetalConstraintManager(metal_centers)
        
        # Reduce weights and increase tolerances
        for constraint in manager.constraints:
            constraint.weight *= 0.5
            constraint.tolerance *= 1.5
        
        return manager
    
    @staticmethod
    def create_distance_only_constraints(metal_centers: List[MetalCenter]) -> MetalConstraintManager:
        """Create constraints focusing only on distances"""
        manager = MetalConstraintManager(metal_centers)
        
        # Remove non-distance constraints
        distance_constraints = [
            c for c in manager.constraints 
            if c.constraint_type == ConstraintType.DISTANCE
        ]
        manager.constraints = distance_constraints
        
        return manager


def apply_metal_constraints_to_pose(pose_coordinates: np.ndarray,
                                  ligand_atom_types: List[str],
                                  metal_centers: List[MetalCenter],
                                  constraint_preset: str = "standard") -> Tuple[np.ndarray, Dict[str, Any]]:
    """Apply metal constraints to optimize a pose"""
    
    # Select constraint preset
    if constraint_preset == "strict":
        manager = ConstraintSetPresets.create_strict_coordination_constraints(metal_centers)
    elif constraint_preset == "flexible":
        manager = ConstraintSetPresets.create_flexible_coordination_constraints(metal_centers)
    elif constraint_preset == "distance_only":
        manager = ConstraintSetPresets.create_distance_only_constraints(metal_centers)
    else:
        manager = MetalConstraintManager(metal_centers)
    
    # Optimize pose
    optimized_coords, final_penalty = manager.optimize_pose_with_constraints(
        pose_coordinates, ligand_atom_types
    )
    
    # Evaluate final constraints
    violations, penalty = manager.evaluate_constraints(optimized_coords, ligand_atom_types)
    
    # Create results summary
    results = {
        'optimized_coordinates': optimized_coords,
        'final_penalty': final_penalty,
        'violations': [
            {
                'type': v.constraint.constraint_type.value,
                'actual_value': v.actual_value,
                'expected_value': v.constraint.target_value,
                'tolerance': v.constraint.tolerance,
                'penalty': v.penalty,
                'description': v.description
            }
            for v in violations
        ],
        'constraint_summary': manager.get_constraint_summary(),
        'total_violations': len(violations),
        'constraint_satisfaction': 1.0 - min(1.0, penalty / 10.0)  # Rough satisfaction score
    }
    
    return optimized_coords, results