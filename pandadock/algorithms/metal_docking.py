"""
Metal-based docking algorithm for PandaDock.

This module implements specialized docking algorithms and scoring functions
for metal-containing systems with coordination constraints.
"""

import numpy as np
import logging
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass, field
import copy

from .base_algorithm import BaseAlgorithm
from ..scoring.base_scoring import BaseScoringFunction


# Metal coordination geometry definitions
METAL_COORDINATION_GEOMETRIES = {
    'octahedral': {
        'coordination_number': 6,
        'ideal_angles': [90.0, 180.0],
        'tolerance': 15.0,
        'metals': ['Fe', 'Co', 'Ni', 'Cr', 'Mn', 'Mo', 'W', 'Ru', 'Os', 'Rh', 'Ir', 'Pd', 'Pt']
    },
    'tetrahedral': {
        'coordination_number': 4,
        'ideal_angles': [109.47],
        'tolerance': 15.0,
        'metals': ['Zn', 'Cd', 'Hg', 'Cu', 'Ag', 'Au', 'Be', 'Mg', 'Al', 'Ga', 'In']
    },
    'square_planar': {
        'coordination_number': 4,
        'ideal_angles': [90.0, 180.0],
        'tolerance': 15.0,
        'metals': ['Pt', 'Pd', 'Ni', 'Cu', 'Au']
    },
    'trigonal_bipyramidal': {
        'coordination_number': 5,
        'ideal_angles': [90.0, 120.0, 180.0],
        'tolerance': 15.0,
        'metals': ['Fe', 'Co', 'Ni', 'Cu', 'Zn']
    },
    'linear': {
        'coordination_number': 2,
        'ideal_angles': [180.0],
        'tolerance': 10.0,
        'metals': ['Ag', 'Au', 'Cu', 'Hg']
    }
}

# Metal-ligand bond lengths (Angstroms)
METAL_BOND_LENGTHS = {
    'Fe': {'N': 2.0, 'O': 1.9, 'S': 2.3, 'C': 2.0, 'P': 2.4},
    'Zn': {'N': 2.1, 'O': 2.0, 'S': 2.3, 'C': 2.0, 'P': 2.4},
    'Cu': {'N': 2.0, 'O': 1.9, 'S': 2.3, 'C': 1.9, 'P': 2.3},
    'Ni': {'N': 2.0, 'O': 2.0, 'S': 2.3, 'C': 1.9, 'P': 2.3},
    'Co': {'N': 2.0, 'O': 2.0, 'S': 2.3, 'C': 1.9, 'P': 2.3},
    'Mn': {'N': 2.2, 'O': 2.1, 'S': 2.4, 'C': 2.1, 'P': 2.4},
    'Mg': {'N': 2.1, 'O': 2.0, 'S': 2.5, 'C': 2.2, 'P': 2.5},
    'Ca': {'N': 2.4, 'O': 2.3, 'S': 2.8, 'C': 2.5, 'P': 2.8},
    'Pt': {'N': 2.0, 'O': 2.0, 'S': 2.3, 'C': 2.0, 'P': 2.3},
    'Pd': {'N': 2.0, 'O': 2.0, 'S': 2.3, 'C': 2.0, 'P': 2.3}
}


@dataclass
class MetalCenter:
    """Represents a metal center in a protein."""
    
    coords: np.ndarray
    element: str
    coordination_geometry: str = 'auto'
    max_coordination: int = 6
    coordinating_atoms: List[Dict[str, Any]] = field(default_factory=list)
    coordination_number: int = 0
    
    def __post_init__(self):
        """Initialize metal center after dataclass creation."""
        self.logger = logging.getLogger(f"MetalCenter_{self.element}")
        
        if self.coordination_geometry == 'auto':
            self.preferred_geometry = self._determine_preferred_geometry()
        else:
            self.preferred_geometry = self.coordination_geometry
    
    def _determine_preferred_geometry(self) -> str:
        """Determine preferred coordination geometry based on metal type."""
        for geometry, info in METAL_COORDINATION_GEOMETRIES.items():
            if self.element in info['metals']:
                return geometry
        return 'octahedral'
    
    def add_coordinating_atom(self, atom: Dict[str, Any], bond_type: str = 'coordinate') -> None:
        """Add a coordinating atom to this metal center."""
        atom_coords = np.array(atom['coords'])
        distance = np.linalg.norm(atom_coords - self.coords)
        
        coordination_info = {
            'atom': atom,
            'coords': atom_coords,
            'distance': distance,
            'bond_type': bond_type,
            'element': atom.get('element', atom.get('name', ''))[:1]
        }
        
        self.coordinating_atoms.append(coordination_info)
        self.coordination_number = len(self.coordinating_atoms)
    
    def get_coordination_score(self) -> float:
        """Calculate coordination geometry score (lower is better)."""
        if self.coordination_number == 0:
            return 1000.0
        
        geometry_info = METAL_COORDINATION_GEOMETRIES.get(self.preferred_geometry, {})
        ideal_coordination = geometry_info.get('coordination_number', 6)
        ideal_angles = geometry_info.get('ideal_angles', [90.0])
        tolerance = geometry_info.get('tolerance', 15.0)
        
        score = 0.0
        
        # Coordination number penalty
        coord_diff = abs(self.coordination_number - ideal_coordination)
        score += coord_diff * 10.0
        
        # Bond length penalties
        for coord_info in self.coordinating_atoms:
            element = coord_info['element']
            distance = coord_info['distance']
            
            if self.element in METAL_BOND_LENGTHS and element in METAL_BOND_LENGTHS[self.element]:
                ideal_distance = METAL_BOND_LENGTHS[self.element][element]
                distance_penalty = abs(distance - ideal_distance) * 5.0
                score += distance_penalty
            else:
                if distance < 1.5 or distance > 3.0:
                    score += 20.0
        
        # Angular penalties
        if self.coordination_number >= 2:
            angles = self._calculate_coordination_angles()
            angle_score = self._score_coordination_angles(angles, ideal_angles, tolerance)
            score += angle_score
        
        return score
    
    def _calculate_coordination_angles(self) -> List[float]:
        """Calculate angles between coordinating atoms."""
        if self.coordination_number < 2:
            return []
        
        angles = []
        coords_list = [info['coords'] for info in self.coordinating_atoms]
        
        for i in range(len(coords_list)):
            for j in range(i + 1, len(coords_list)):
                vec1 = coords_list[i] - self.coords
                vec2 = coords_list[j] - self.coords
                
                cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_angle))
                angles.append(angle)
        
        return angles
    
    def _score_coordination_angles(self, actual_angles: List[float], 
                                 ideal_angles: List[float], tolerance: float) -> float:
        """Score coordination angles against ideal geometry."""
        if not actual_angles:
            return 0.0
        
        score = 0.0
        
        for actual_angle in actual_angles:
            min_deviation = float('inf')
            for ideal_angle in ideal_angles:
                deviation = abs(actual_angle - ideal_angle)
                min_deviation = min(min_deviation, deviation)
            
            if min_deviation > tolerance:
                score += (min_deviation - tolerance) * 2.0
        
        return score


@dataclass
class MetalConstraint:
    """Represents a constraint for metal coordination during docking."""
    
    metal_center: MetalCenter
    constraint_type: str = 'coordination'
    strength: float = 1.0
    target_atoms: List[Dict[str, Any]] = field(default_factory=list)
    required_coordination: Optional[int] = None
    
    def evaluate(self, ligand_pose: Any) -> float:
        """Evaluate constraint violation penalty."""
        if self.constraint_type == 'coordination':
            return self._evaluate_coordination_constraint(ligand_pose)
        elif self.constraint_type == 'distance':
            return self._evaluate_distance_constraint(ligand_pose)
        elif self.constraint_type == 'angle':
            return self._evaluate_angle_constraint(ligand_pose)
        else:
            return 0.0
    
    def _evaluate_coordination_constraint(self, ligand_pose: Any) -> float:
        """Evaluate coordination constraint."""
        coordinating_atoms = self._find_coordinating_atoms(ligand_pose)
        
        temp_metal_center = copy.deepcopy(self.metal_center)
        for atom in coordinating_atoms:
            temp_metal_center.add_coordinating_atom(atom)
        
        coord_score = temp_metal_center.get_coordination_score()
        return coord_score * self.strength
    
    def _find_coordinating_atoms(self, ligand_pose: Any, max_distance: float = 3.5) -> List[Dict[str, Any]]:
        """Find atoms in ligand that could coordinate to metal."""
        coordinating_elements = ['N', 'O', 'S', 'P', 'C']
        coordinating_atoms = []
        
        ligand_atoms = getattr(ligand_pose, 'atoms', [])
        if hasattr(ligand_pose, 'coords'):
            # Handle coordinate arrays
            coords = ligand_pose.coords
            for i, coord in enumerate(coords):
                atom = {
                    'coords': coord,
                    'element': 'C',  # Default element
                    'index': i
                }
                distance = np.linalg.norm(coord - self.metal_center.coords)
                if distance <= max_distance:
                    coordinating_atoms.append(atom)
        else:
            # Handle atom lists
            for atom in ligand_atoms:
                element = atom.get('element', atom.get('name', ''))[:1]
                if element in coordinating_elements:
                    atom_coords = np.array(atom['coords'])
                    distance = np.linalg.norm(atom_coords - self.metal_center.coords)
                    
                    if distance <= max_distance:
                        coordinating_atoms.append(atom)
        
        return coordinating_atoms
    
    def _evaluate_distance_constraint(self, ligand_pose: Any) -> float:
        """Evaluate distance constraint."""
        penalty = 0.0
        
        for target_atom_info in self.target_atoms:
            ligand_atom = self._find_ligand_atom(ligand_pose, target_atom_info)
            if ligand_atom:
                distance = np.linalg.norm(
                    np.array(ligand_atom['coords']) - self.metal_center.coords
                )
                ideal_distance = target_atom_info.get('ideal_distance', 2.0)
                tolerance = target_atom_info.get('tolerance', 0.3)
                
                if abs(distance - ideal_distance) > tolerance:
                    penalty += abs(distance - ideal_distance) * self.strength
        
        return penalty
    
    def _evaluate_angle_constraint(self, ligand_pose: Any) -> float:
        """Evaluate angle constraint."""
        return 0.0
    
    def _find_ligand_atom(self, ligand_pose: Any, target_info: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Find specific atom in ligand based on target information."""
        element = target_info.get('element', '')
        atom_name = target_info.get('name', '')
        
        ligand_atoms = getattr(ligand_pose, 'atoms', [])
        for atom in ligand_atoms:
            if element and atom.get('element', '').startswith(element):
                return atom
            if atom_name and atom.get('name', '') == atom_name:
                return atom
        
        return None


class MetalDockingScorer(BaseScoringFunction):
    """Specialized scoring function for metal-containing systems."""
    
    def __init__(self, base_scoring_function: BaseScoringFunction, 
                 metal_centers: Optional[List[MetalCenter]] = None,
                 metal_constraints: Optional[List[MetalConstraint]] = None,
                 metal_weight: float = 10.0):
        """Initialize metal docking scorer."""
        super().__init__()
        self.base_scoring_function = base_scoring_function
        self.metal_centers = metal_centers or []
        self.metal_constraints = metal_constraints or []
        self.metal_weight = metal_weight
        self.logger = logging.getLogger("MetalDockingScorer")
    
    def score(self, protein: Any, ligand: Any) -> float:
        """Score protein-ligand complex including metal coordination terms."""
        # Get base score
        base_score = self.base_scoring_function.score(protein, ligand)
        
        # Calculate metal coordination score
        if hasattr(ligand, 'coords'):
            metal_score = self._calculate_metal_score(ligand.coords)
        else:
            metal_score = 0.0
        
        # Combine scores
        total_score = base_score + (metal_score * self.metal_weight)
        
        self.logger.debug(f"Base: {base_score:.2f}, Metal: {metal_score:.2f}, Total: {total_score:.2f}")
        
        return total_score
    
    def score_pose(self, ligand_coords: np.ndarray, protein_coords: np.ndarray) -> float:
        """Score ligand pose including metal coordination terms."""
        # Get base score
        base_score = self.base_scoring_function.score_pose(ligand_coords, protein_coords)
        
        # Calculate metal coordination score
        metal_score = self._calculate_metal_score(ligand_coords)
        
        # Combine scores
        total_score = base_score + (metal_score * self.metal_weight)
        
        self.logger.debug(f"Base: {base_score:.2f}, Metal: {metal_score:.2f}, Total: {total_score:.2f}")
        
        return total_score
    
    def _calculate_metal_score(self, ligand_coords: np.ndarray) -> float:
        """Calculate metal coordination contribution to score."""
        total_metal_score = 0.0
        
        # Create a temporary ligand object for constraint evaluation
        temp_ligand = type('TempLigand', (), {
            'coords': ligand_coords,
            'atoms': [{'coords': coord, 'element': 'C'} for coord in ligand_coords]
        })()
        
        # Evaluate each metal constraint
        for constraint in self.metal_constraints:
            constraint_penalty = constraint.evaluate(temp_ligand)
            total_metal_score += constraint_penalty
        
        # Add coordination bonuses for well-coordinated metals
        for metal_center in self.metal_centers:
            coordinating_atoms = self._find_coordinating_atoms(ligand_coords, metal_center)
            
            temp_center = copy.deepcopy(metal_center)
            for atom in coordinating_atoms:
                temp_center.add_coordinating_atom(atom)
            
            coord_penalty = temp_center.get_coordination_score()
            total_metal_score += coord_penalty
        
        return total_metal_score
    
    def _find_coordinating_atoms(self, ligand_coords: np.ndarray, 
                               metal_center: MetalCenter, max_distance: float = 3.5) -> List[Dict[str, Any]]:
        """Find potential coordinating atoms."""
        coordinating_atoms = []
        
        for i, coord in enumerate(ligand_coords):
            distance = np.linalg.norm(coord - metal_center.coords)
            if distance <= max_distance:
                atom = {
                    'coords': coord,
                    'element': 'C',  # Default element
                    'index': i,
                    'distance': distance
                }
                coordinating_atoms.append(atom)
        
        return coordinating_atoms


class MetalDockingAlgorithm(BaseAlgorithm):
    """Specialized search algorithm for metal-containing ligands."""
    
    def __init__(self, base_algorithm: BaseAlgorithm, 
                 metal_centers: List[MetalCenter],
                 metal_constraints: List[MetalConstraint],
                 refinement_iterations: int = 20):
        """Initialize metal docking search algorithm."""
        super().__init__(base_algorithm.scoring_function)
        self.base_algorithm = base_algorithm
        self.metal_centers = metal_centers
        self.metal_constraints = metal_constraints
        self.refinement_iterations = refinement_iterations
        self.logger = logging.getLogger("MetalDockingAlgorithm")
    
    def search(self, protein: Any, ligand: Any) -> List[Tuple[Any, float]]:
        """Perform metal-aware docking search."""
        self.logger.info(f"Starting metal-aware docking with {len(self.metal_centers)} metal centers")
        
        # Perform base search
        base_results = self.base_algorithm.search(protein, ligand)
        
        # Apply metal-specific refinement
        refined_results = []
        
        for pose, score in base_results:
            # Refine pose for better metal coordination
            refined_pose = self._refine_metal_coordination(protein, pose)
            
            # Re-score with metal terms
            if isinstance(self.scoring_function, MetalDockingScorer):
                if hasattr(refined_pose, 'coords'):
                    refined_score = self.scoring_function.score_pose(refined_pose.coords, protein.coords)
                else:
                    refined_score = self.scoring_function.score_pose(refined_pose, protein.coords)
            else:
                refined_score = score
            
            refined_results.append((refined_pose, refined_score))
        
        # Sort by score
        refined_results.sort(key=lambda x: x[1])
        
        self.logger.info(f"Metal-aware docking completed. Best score: {refined_results[0][1]:.2f}")
        
        return refined_results
    
    def _refine_metal_coordination(self, protein: Any, pose: Any) -> Any:
        """Refine ligand pose to improve metal coordination."""
        refined_pose = copy.deepcopy(pose)
        
        for iteration in range(self.refinement_iterations):
            improved = False
            
            for metal_center in self.metal_centers:
                # Find coordinating atoms
                if hasattr(refined_pose, 'coords'):
                    coordinating_atoms = self._find_ligand_coordinating_atoms(refined_pose.coords, metal_center)
                else:
                    coordinating_atoms = self._find_ligand_coordinating_atoms(refined_pose, metal_center)
                
                if coordinating_atoms:
                    # Try to optimize coordination geometry
                    optimized_pose = self._optimize_coordination_geometry(
                        refined_pose, metal_center, coordinating_atoms
                    )
                    
                    # Check if improvement was made
                    if self._is_better_coordination(optimized_pose, refined_pose, metal_center):
                        refined_pose = optimized_pose
                        improved = True
            
            if not improved:
                break
        
        return refined_pose
    
    def _find_ligand_coordinating_atoms(self, ligand_coords: np.ndarray, 
                                      metal_center: MetalCenter, max_distance: float = 3.5) -> List[Dict[str, Any]]:
        """Find atoms in ligand that could coordinate to metal."""
        coordinating_atoms = []
        
        if isinstance(ligand_coords, np.ndarray):
            coords_array = ligand_coords
        else:
            coords_array = getattr(ligand_coords, 'coords', ligand_coords)
        
        for i, coord in enumerate(coords_array):
            distance = np.linalg.norm(coord - metal_center.coords)
            if distance <= max_distance:
                coordinating_atoms.append({
                    'index': i,
                    'coords': coord,
                    'distance': distance,
                    'element': 'C'  # Default element
                })
        
        return coordinating_atoms
    
    def _optimize_coordination_geometry(self, pose: Any, metal_center: MetalCenter, 
                                      coordinating_atoms: List[Dict[str, Any]]) -> Any:
        """Optimize ligand pose for better metal coordination geometry."""
        if not coordinating_atoms:
            return pose
        
        optimized_pose = copy.deepcopy(pose)
        
        # Simple optimization: adjust distances to ideal values
        for coord_info in coordinating_atoms:
            element = coord_info['element']
            current_distance = coord_info['distance']
            
            if (metal_center.element in METAL_BOND_LENGTHS and 
                element in METAL_BOND_LENGTHS[metal_center.element]):
                
                ideal_distance = METAL_BOND_LENGTHS[metal_center.element][element]
                
                if abs(current_distance - ideal_distance) > 0.2:
                    # Calculate adjustment
                    atom_coords = coord_info['coords']
                    metal_coords = metal_center.coords
                    
                    vector = atom_coords - metal_coords
                    vector_norm = np.linalg.norm(vector)
                    
                    if vector_norm > 0:
                        scaling_factor = ideal_distance / vector_norm
                        adjustment = vector * (scaling_factor - 1.0)
                        
                        # Apply small adjustment
                        if hasattr(optimized_pose, 'coords'):
                            optimized_pose.coords += adjustment * 0.1
                        elif hasattr(optimized_pose, 'translate'):
                            optimized_pose.translate(adjustment * 0.1)
        
        return optimized_pose
    
    def _is_better_coordination(self, new_pose: Any, old_pose: Any, metal_center: MetalCenter) -> bool:
        """Check if new pose has better metal coordination than old pose."""
        old_score = self._calculate_coordination_score(old_pose, metal_center)
        new_score = self._calculate_coordination_score(new_pose, metal_center)
        
        return new_score < old_score
    
    def _calculate_coordination_score(self, pose: Any, metal_center: MetalCenter) -> float:
        """Calculate coordination score for a pose."""
        if hasattr(pose, 'coords'):
            coordinating_atoms = self._find_ligand_coordinating_atoms(pose.coords, metal_center)
        else:
            coordinating_atoms = self._find_ligand_coordinating_atoms(pose, metal_center)
        
        temp_center = copy.deepcopy(metal_center)
        for coord_info in coordinating_atoms:
            atom = {
                'coords': coord_info['coords'],
                'element': coord_info['element']
            }
            temp_center.add_coordinating_atom(atom)
        
        return temp_center.get_coordination_score()


class MetalDockingPreparation:
    """Utilities for preparing metal-containing systems for docking."""
    
    @staticmethod
    def detect_metal_centers(protein: Any, metal_elements: Optional[List[str]] = None) -> List[MetalCenter]:
        """Detect metal centers in protein structure."""
        if metal_elements is None:
            metal_elements = ['FE', 'ZN', 'CU', 'MG', 'CA', 'MN', 'CO', 
                            'NI', 'PT', 'PD', 'MO', 'W', 'CD', 'HG']
        
        metal_centers = []
        
        protein_atoms = getattr(protein, 'atoms', [])
        for atom in protein_atoms:
            element = atom.get('element', atom.get('name', ''))[:2].strip().upper()
            
            if element in metal_elements:
                metal_center = MetalCenter(
                    coords=np.array(atom['coords']),
                    element=element
                )
                
                # Find existing coordinating atoms in protein
                coordinating_atoms = MetalDockingPreparation._find_protein_coordinating_atoms(
                    protein, metal_center
                )
                
                for coord_atom in coordinating_atoms:
                    metal_center.add_coordinating_atom(coord_atom, 'protein_coordinate')
                
                metal_centers.append(metal_center)
        
        return metal_centers
    
    @staticmethod
    def _find_protein_coordinating_atoms(protein: Any, metal_center: MetalCenter, 
                                       max_distance: float = 3.0) -> List[Dict[str, Any]]:
        """Find atoms in protein already coordinating to metal."""
        coordinating_elements = ['N', 'O', 'S', 'P']
        coordinating_atoms = []
        
        protein_atoms = getattr(protein, 'atoms', [])
        for atom in protein_atoms:
            if atom.get('coords') == metal_center.coords.tolist():
                continue
            
            element = atom.get('element', atom.get('name', ''))[:1]
            if element in coordinating_elements:
                distance = np.linalg.norm(
                    np.array(atom['coords']) - metal_center.coords
                )
                if distance <= max_distance:
                    coordinating_atoms.append(atom)
        
        return coordinating_atoms
    
    @staticmethod
    def create_metal_constraints(metal_centers: List[MetalCenter], 
                               constraint_types: Optional[List[str]] = None,
                               strengths: Optional[List[float]] = None) -> List[MetalConstraint]:
        """Create metal constraints for docking."""
        if constraint_types is None:
            constraint_types = ['coordination']
        if strengths is None:
            strengths = [1.0] * len(constraint_types)
        
        constraints = []
        
        for metal_center in metal_centers:
            for constraint_type, strength in zip(constraint_types, strengths):
                constraint = MetalConstraint(
                    metal_center=metal_center,
                    constraint_type=constraint_type,
                    strength=strength
                )
                constraints.append(constraint)
        
        return constraints
    
    @staticmethod
    def analyze_ligand_metal_binding_sites(ligand: Any) -> Dict[str, Any]:
        """Analyze potential metal binding sites in ligand."""
        coordinating_elements = ['N', 'O', 'S', 'P']
        binding_sites = []
        
        ligand_atoms = getattr(ligand, 'atoms', [])
        for i, atom in enumerate(ligand_atoms):
            element = atom.get('element', atom.get('name', ''))[:1]
            if element in coordinating_elements:
                binding_site = {
                    'atom_index': i,
                    'atom': atom,
                    'element': element,
                    'coords': np.array(atom['coords']),
                    'binding_strength': MetalDockingPreparation._estimate_binding_strength(atom)
                }
                binding_sites.append(binding_site)
        
        # Identify potential chelating groups
        chelating_groups = MetalDockingPreparation._identify_chelating_groups(binding_sites)
        
        return {
            'binding_sites': binding_sites,
            'chelating_groups': chelating_groups,
            'total_coordinating_atoms': len(binding_sites)
        }
    
    @staticmethod
    def _estimate_binding_strength(atom: Dict[str, Any]) -> float:
        """Estimate relative metal binding strength of an atom."""
        element = atom.get('element', atom.get('name', ''))[:1]
        
        strength_map = {
            'N': 0.8,  # Nitrogen
            'O': 0.6,  # Oxygen  
            'S': 0.9,  # Sulfur
            'P': 0.7   # Phosphorus
        }
        
        return strength_map.get(element, 0.3)
    
    @staticmethod
    def _identify_chelating_groups(binding_sites: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Identify potential chelating groups."""
        chelating_groups = []
        
        for i, site1 in enumerate(binding_sites):
            for j, site2 in enumerate(binding_sites[i+1:], i+1):
                distance = np.linalg.norm(site1['coords'] - site2['coords'])
                
                if 2.5 <= distance <= 4.5:
                    chelating_groups.append({
                        'type': 'bidentate',
                        'atoms': [site1['atom_index'], site2['atom_index']],
                        'distance': distance
                    })
        
        return chelating_groups