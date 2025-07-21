# -*- coding: utf-8 -*-
"""
Metal Coordination Geometry Module

This module handles the geometric constraints and coordination patterns
for metal ions in metalloproteins during docking.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Union, Any
from dataclasses import dataclass, field
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class CoordinationGeometry(Enum):
    """Supported metal coordination geometries"""
    LINEAR = "linear"
    TRIGONAL_PLANAR = "trigonal_planar"
    TETRAHEDRAL = "tetrahedral"
    SQUARE_PLANAR = "square_planar"
    TRIGONAL_BIPYRAMIDAL = "trigonal_bipyramidal"
    OCTAHEDRAL = "octahedral"
    SQUARE_PYRAMIDAL = "square_pyramidal"
    PENTAGONAL_BIPYRAMIDAL = "pentagonal_bipyramidal"


class MetalType(Enum):
    """Common metal types in metalloproteins"""
    ZN = "Zn"
    FE = "Fe"
    MG = "Mg"
    CA = "Ca"
    MN = "Mn"
    CU = "Cu"
    NI = "Ni"
    CO = "Co"
    MO = "Mo"
    W = "W"
    V = "V"
    CR = "Cr"


@dataclass
class MetalCenter:
    """Represents a metal center with its coordination environment"""
    metal_type: MetalType
    coordinates: np.ndarray  # 3D coordinates of metal ion
    coordination_number: int
    geometry: CoordinationGeometry
    coordinating_atoms: List[Dict] = field(default_factory=list)
    distance_constraints: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    angle_constraints: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    oxidation_state: Optional[int] = None
    charge: float = 0.0
    
    def __post_init__(self):
        """Initialize default distance and angle constraints"""
        if not self.distance_constraints:
            self.distance_constraints = self._get_default_distances()
        if not self.angle_constraints:
            self.angle_constraints = self._get_default_angles()
    
    def _get_default_distances(self) -> Dict[str, Tuple[float, float]]:
        """Get default metal-ligand distance constraints"""
        # Distance ranges (min, max) in Angstroms
        distances = {
            MetalType.ZN: {
                'N': (1.9, 2.4),
                'O': (1.8, 2.3),
                'S': (2.2, 2.8),
                'C': (1.9, 2.5),
                'P': (2.3, 2.9)
            },
            MetalType.FE: {
                'N': (1.8, 2.5),
                'O': (1.7, 2.4),
                'S': (2.1, 2.9),
                'C': (1.8, 2.6),
                'P': (2.2, 3.0)
            },
            MetalType.MG: {
                'N': (2.0, 2.6),
                'O': (1.9, 2.5),
                'S': (2.4, 3.0),
                'C': (2.1, 2.7),
                'P': (2.5, 3.1)
            },
            MetalType.CA: {
                'N': (2.3, 2.9),
                'O': (2.2, 2.8),
                'S': (2.7, 3.3),
                'C': (2.4, 3.0),
                'P': (2.8, 3.4)
            },
            MetalType.CU: {
                'N': (1.8, 2.4),
                'O': (1.7, 2.3),
                'S': (2.1, 2.7),
                'C': (1.8, 2.5),
                'P': (2.2, 2.8)
            },
            MetalType.MN: {
                'N': (1.9, 2.5),
                'O': (1.8, 2.4),
                'S': (2.2, 2.8),
                'C': (1.9, 2.6),
                'P': (2.3, 2.9)
            }
        }
        
        return distances.get(self.metal_type, distances[MetalType.ZN])
    
    def _get_default_angles(self) -> Dict[str, Tuple[float, float]]:
        """Get default coordination angle constraints"""
        # Angle ranges (min, max) in degrees
        angles = {
            CoordinationGeometry.LINEAR: {
                'L-M-L': (170, 190)
            },
            CoordinationGeometry.TRIGONAL_PLANAR: {
                'L-M-L': (110, 130)
            },
            CoordinationGeometry.TETRAHEDRAL: {
                'L-M-L': (100, 120)
            },
            CoordinationGeometry.SQUARE_PLANAR: {
                'L-M-L_cis': (85, 95),
                'L-M-L_trans': (175, 185)
            },
            CoordinationGeometry.OCTAHEDRAL: {
                'L-M-L_cis': (85, 95),
                'L-M-L_trans': (175, 185)
            },
            CoordinationGeometry.TRIGONAL_BIPYRAMIDAL: {
                'L-M-L_equatorial': (115, 125),
                'L-M-L_axial': (175, 185),
                'L_eq-M-L_ax': (85, 95)
            }
        }
        
        return angles.get(self.geometry, angles[CoordinationGeometry.TETRAHEDRAL])


class MetalCoordinationAnalyzer:
    """Analyzes and validates metal coordination environments"""
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Standard coordination geometries and their ideal angles
        self.ideal_geometries = {
            CoordinationGeometry.LINEAR: {
                'coordination_number': 2,
                'angles': [180.0]
            },
            CoordinationGeometry.TRIGONAL_PLANAR: {
                'coordination_number': 3,
                'angles': [120.0, 120.0, 120.0]
            },
            CoordinationGeometry.TETRAHEDRAL: {
                'coordination_number': 4,
                'angles': [109.47] * 6  # All L-M-L angles
            },
            CoordinationGeometry.SQUARE_PLANAR: {
                'coordination_number': 4,
                'angles': [90.0] * 4 + [180.0] * 2  # cis + trans angles
            },
            CoordinationGeometry.OCTAHEDRAL: {
                'coordination_number': 6,
                'angles': [90.0] * 12 + [180.0] * 3  # cis + trans angles
            }
        }
    
    def detect_metal_centers(self, protein_structure: Dict) -> List[MetalCenter]:
        """Detect metal centers in protein structure"""
        metal_centers = []
        
        # Look for metal atoms in the structure
        for atom in protein_structure.get('atoms', []):
            if self._is_metal_atom(atom):
                metal_center = self._analyze_metal_environment(atom, protein_structure)
                if metal_center:
                    metal_centers.append(metal_center)
        
        self.logger.info(f"Detected {len(metal_centers)} metal centers")
        return metal_centers
    
    def _is_metal_atom(self, atom: Dict) -> bool:
        """Check if an atom is a metal"""
        element = atom.get('element', '').upper()
        metal_elements = {mt.value for mt in MetalType}
        return element in metal_elements
    
    def _analyze_metal_environment(self, metal_atom: Dict, protein_structure: Dict) -> Optional[MetalCenter]:
        """Analyze the coordination environment of a metal atom"""
        metal_coords = np.array([metal_atom['x'], metal_atom['y'], metal_atom['z']])
        metal_type = MetalType(metal_atom['element'])
        
        # Find coordinating atoms within reasonable distance
        coordinating_atoms = self._find_coordinating_atoms(metal_coords, protein_structure)
        
        if not coordinating_atoms:
            return None
        
        # Determine coordination geometry
        coordination_number = len(coordinating_atoms)
        geometry = self._determine_geometry(metal_coords, coordinating_atoms, coordination_number)
        
        return MetalCenter(
            metal_type=metal_type,
            coordinates=metal_coords,
            coordination_number=coordination_number,
            geometry=geometry,
            coordinating_atoms=coordinating_atoms,
            oxidation_state=self._estimate_oxidation_state(metal_type, coordinating_atoms),
            charge=self._calculate_formal_charge(metal_type, coordinating_atoms)
        )
    
    def _find_coordinating_atoms(self, metal_coords: np.ndarray, protein_structure: Dict, 
                                max_distance: float = 3.5) -> List[Dict]:
        """Find atoms coordinating to the metal center"""
        coordinating_atoms = []
        
        for atom in protein_structure.get('atoms', []):
            if atom.get('element') in ['H']:  # Skip hydrogens
                continue
            
            atom_coords = np.array([atom['x'], atom['y'], atom['z']])
            distance = np.linalg.norm(metal_coords - atom_coords)
            
            if 0.1 < distance <= max_distance:  # Avoid self and very close atoms
                coordinating_atoms.append({
                    'atom_id': atom.get('atom_id'),
                    'element': atom.get('element'),
                    'residue': atom.get('residue'),
                    'coordinates': atom_coords,
                    'distance': distance,
                    'atom_type': self._classify_coordination_atom(atom)
                })
        
        # Sort by distance and take closest atoms
        coordinating_atoms.sort(key=lambda x: x['distance'])
        return coordinating_atoms[:8]  # Max 8 coordination
    
    def _classify_coordination_atom(self, atom: Dict) -> str:
        """Classify the type of coordinating atom"""
        element = atom.get('element', '')
        residue = atom.get('residue', '')
        atom_name = atom.get('atom_name', '')
        
        # Classify based on element and chemical environment
        if element == 'N':
            if residue in ['HIS', 'CYS', 'MET']:
                return 'protein_sidechain'
            elif 'N' in atom_name and residue in ['ARG', 'LYS']:
                return 'protein_basic'
            else:
                return 'protein_backbone'
        elif element == 'O':
            if residue in ['ASP', 'GLU']:
                return 'protein_acidic'
            elif residue in ['SER', 'THR', 'TYR']:
                return 'protein_polar'
            elif 'OW' in atom_name:
                return 'water'
            else:
                return 'protein_backbone'
        elif element == 'S':
            if residue in ['CYS', 'MET']:
                return 'protein_sulfur'
            else:
                return 'ligand_sulfur'
        else:
            return 'ligand_other'
    
    def _determine_geometry(self, metal_coords: np.ndarray, coordinating_atoms: List[Dict], 
                          coordination_number: int) -> CoordinationGeometry:
        """Determine the coordination geometry based on atom positions"""
        if coordination_number == 2:
            return CoordinationGeometry.LINEAR
        elif coordination_number == 3:
            return self._analyze_three_coordinate(metal_coords, coordinating_atoms)
        elif coordination_number == 4:
            return self._analyze_four_coordinate(metal_coords, coordinating_atoms)
        elif coordination_number == 5:
            return self._analyze_five_coordinate(metal_coords, coordinating_atoms)
        elif coordination_number == 6:
            return CoordinationGeometry.OCTAHEDRAL
        else:
            return CoordinationGeometry.OCTAHEDRAL  # Default for high coordination
    
    def _analyze_three_coordinate(self, metal_coords: np.ndarray, 
                                coordinating_atoms: List[Dict]) -> CoordinationGeometry:
        """Analyze 3-coordinate geometry"""
        # Calculate angles between ligands
        angles = []
        coords = [atom['coordinates'] for atom in coordinating_atoms]
        
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                v1 = coords[i] - metal_coords
                v2 = coords[j] - metal_coords
                angle = self._calculate_angle(v1, v2)
                angles.append(angle)
        
        # Check if angles are close to 120° (trigonal planar)
        avg_angle = np.mean(angles)
        if 110 < avg_angle < 130:
            return CoordinationGeometry.TRIGONAL_PLANAR
        else:
            return CoordinationGeometry.TRIGONAL_PLANAR  # Default for 3-coordinate
    
    def _analyze_four_coordinate(self, metal_coords: np.ndarray, 
                               coordinating_atoms: List[Dict]) -> CoordinationGeometry:
        """Analyze 4-coordinate geometry (tetrahedral vs square planar)"""
        coords = [atom['coordinates'] for atom in coordinating_atoms]
        
        # Calculate all L-M-L angles
        angles = []
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                v1 = coords[i] - metal_coords
                v2 = coords[j] - metal_coords
                angle = self._calculate_angle(v1, v2)
                angles.append(angle)
        
        # Tetrahedral: all angles ~109.5°
        # Square planar: 4 angles ~90°, 2 angles ~180°
        angles_90 = sum(1 for angle in angles if 85 < angle < 95)
        angles_180 = sum(1 for angle in angles if 175 < angle < 185)
        angles_109 = sum(1 for angle in angles if 105 < angle < 115)
        
        if angles_90 >= 4 and angles_180 >= 2:
            return CoordinationGeometry.SQUARE_PLANAR
        elif angles_109 >= 4:
            return CoordinationGeometry.TETRAHEDRAL
        else:
            # Use dihedral angle test for ambiguous cases
            return self._dihedral_geometry_test(metal_coords, coords)
    
    def _analyze_five_coordinate(self, metal_coords: np.ndarray, 
                               coordinating_atoms: List[Dict]) -> CoordinationGeometry:
        """Analyze 5-coordinate geometry"""
        # For simplicity, default to trigonal bipyramidal
        # Could be extended to distinguish from square pyramidal
        return CoordinationGeometry.TRIGONAL_BIPYRAMIDAL
    
    def _dihedral_geometry_test(self, metal_coords: np.ndarray, 
                              ligand_coords: List[np.ndarray]) -> CoordinationGeometry:
        """Use dihedral angles to distinguish tetrahedral from square planar"""
        if len(ligand_coords) != 4:
            return CoordinationGeometry.TETRAHEDRAL
        
        # Calculate dihedral angles for all combinations
        dihedrals = []
        for i in range(4):
            for j in range(i + 1, 4):
                for k in range(j + 1, 4):
                    for l in range(k + 1, 4):
                        dihedral = self._calculate_dihedral(
                            ligand_coords[i], metal_coords, 
                            ligand_coords[j], ligand_coords[k]
                        )
                        dihedrals.append(abs(dihedral))
        
        # Square planar should have dihedrals close to 0° or 180°
        planar_dihedrals = sum(1 for d in dihedrals if d < 20 or d > 160)
        
        if planar_dihedrals > len(dihedrals) * 0.7:
            return CoordinationGeometry.SQUARE_PLANAR
        else:
            return CoordinationGeometry.TETRAHEDRAL
    
    def _calculate_angle(self, v1: np.ndarray, v2: np.ndarray) -> float:
        """Calculate angle between two vectors in degrees"""
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
    
    def _calculate_dihedral(self, p1: np.ndarray, p2: np.ndarray, 
                          p3: np.ndarray, p4: np.ndarray) -> float:
        """Calculate dihedral angle between four points"""
        v1 = p2 - p1
        v2 = p3 - p2
        v3 = p4 - p3
        
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        
        cos_dihedral = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
        cos_dihedral = np.clip(cos_dihedral, -1.0, 1.0)
        
        return np.degrees(np.arccos(cos_dihedral))
    
    def _estimate_oxidation_state(self, metal_type: MetalType, 
                                 coordinating_atoms: List[Dict]) -> int:
        """Estimate metal oxidation state from coordination environment"""
        # Simplified estimation based on common oxidation states
        oxidation_states = {
            MetalType.ZN: 2,
            MetalType.FE: 2,  # Could be 2 or 3
            MetalType.MG: 2,
            MetalType.CA: 2,
            MetalType.MN: 2,
            MetalType.CU: 2,
            MetalType.NI: 2,
            MetalType.CO: 2
        }
        
        return oxidation_states.get(metal_type, 2)
    
    def _calculate_formal_charge(self, metal_type: MetalType, 
                               coordinating_atoms: List[Dict]) -> float:
        """Calculate formal charge on metal center"""
        # Simplified charge calculation
        oxidation_state = self._estimate_oxidation_state(metal_type, coordinating_atoms)
        
        # Account for charge contribution from ligands
        ligand_charge = 0.0
        for atom in coordinating_atoms:
            if atom['atom_type'] == 'protein_acidic':
                ligand_charge -= 0.5  # Partial negative charge
            elif atom['atom_type'] == 'protein_basic':
                ligand_charge += 0.2  # Partial positive charge
        
        return float(oxidation_state) + ligand_charge


class MetalConstraintValidator:
    """Validates metal coordination constraints during docking"""
    
    def __init__(self, tolerance_factor: float = 1.2):
        self.tolerance_factor = tolerance_factor
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def validate_coordination(self, metal_center: MetalCenter, 
                            ligand_atoms: List[Dict]) -> Dict[str, Any]:
        """Validate ligand coordination to metal center"""
        validation_results = {
            'valid': True,
            'violations': [],
            'coordination_score': 0.0,
            'geometric_score': 0.0,
            'distance_violations': [],
            'angle_violations': []
        }
        
        # Find potential coordinating atoms in ligand
        coordinating_ligand_atoms = self._find_ligand_coordinators(
            metal_center, ligand_atoms
        )
        
        if not coordinating_ligand_atoms:
            validation_results['valid'] = False
            validation_results['violations'].append("No coordinating atoms found in ligand")
            return validation_results
        
        # Validate distances
        distance_score = self._validate_distances(
            metal_center, coordinating_ligand_atoms, validation_results
        )
        
        # Validate angles
        angle_score = self._validate_angles(
            metal_center, coordinating_ligand_atoms, validation_results
        )
        
        # Calculate overall scores
        validation_results['coordination_score'] = distance_score
        validation_results['geometric_score'] = angle_score
        
        # Overall validity
        validation_results['valid'] = (
            len(validation_results['violations']) == 0 and
            distance_score > 0.5 and
            angle_score > 0.5
        )
        
        return validation_results
    
    def _find_ligand_coordinators(self, metal_center: MetalCenter, 
                                ligand_atoms: List[Dict]) -> List[Dict]:
        """Find atoms in ligand that can coordinate to metal"""
        coordinating_atoms = []
        max_distance = 4.0  # Angstroms
        
        for atom in ligand_atoms:
            if atom.get('element') in ['N', 'O', 'S', 'P']:  # Common coordinating elements
                distance = np.linalg.norm(
                    metal_center.coordinates - atom['coordinates']
                )
                
                if distance <= max_distance:
                    coordinating_atoms.append({
                        **atom,
                        'distance': distance
                    })
        
        return coordinating_atoms
    
    def _validate_distances(self, metal_center: MetalCenter, 
                          coordinating_atoms: List[Dict], 
                          validation_results: Dict) -> float:
        """Validate metal-ligand distances"""
        distance_scores = []
        
        for atom in coordinating_atoms:
            element = atom.get('element', 'C')
            distance = atom['distance']
            
            # Get expected distance range
            distance_range = metal_center.distance_constraints.get(element, (1.5, 3.5))
            min_dist, max_dist = distance_range
            
            # Apply tolerance
            min_dist /= self.tolerance_factor
            max_dist *= self.tolerance_factor
            
            if min_dist <= distance <= max_dist:
                # Calculate score based on how close to ideal distance
                ideal_dist = (min_dist + max_dist) / 2
                score = 1.0 - abs(distance - ideal_dist) / (max_dist - min_dist)
                distance_scores.append(max(0.0, score))
            else:
                # Distance violation
                violation = {
                    'type': 'distance',
                    'atom': atom.get('atom_id'),
                    'element': element,
                    'actual_distance': distance,
                    'expected_range': (min_dist, max_dist)
                }
                validation_results['distance_violations'].append(violation)
                validation_results['violations'].append(
                    f"Distance violation for {element}: {distance:.2f} Å "
                    f"(expected {min_dist:.2f}-{max_dist:.2f} Å)"
                )
                distance_scores.append(0.0)
        
        return np.mean(distance_scores) if distance_scores else 0.0
    
    def _validate_angles(self, metal_center: MetalCenter, 
                       coordinating_atoms: List[Dict], 
                       validation_results: Dict) -> float:
        """Validate coordination angles"""
        if len(coordinating_atoms) < 2:
            return 1.0  # No angle constraints for single coordination
        
        angle_scores = []
        
        # Calculate all L-M-L angles
        for i in range(len(coordinating_atoms)):
            for j in range(i + 1, len(coordinating_atoms)):
                atom1 = coordinating_atoms[i]
                atom2 = coordinating_atoms[j]
                
                # Calculate angle
                v1 = atom1['coordinates'] - metal_center.coordinates
                v2 = atom2['coordinates'] - metal_center.coordinates
                angle = self._calculate_angle(v1, v2)
                
                # Validate against expected angles for geometry
                angle_score = self._score_coordination_angle(
                    angle, metal_center.geometry, validation_results
                )
                angle_scores.append(angle_score)
        
        return np.mean(angle_scores) if angle_scores else 0.0
    
    def _score_coordination_angle(self, angle: float, geometry: CoordinationGeometry, 
                                validation_results: Dict) -> float:
        """Score an individual coordination angle"""
        expected_angles = self._get_expected_angles(geometry)
        
        # Find best matching expected angle
        best_score = 0.0
        best_expected = None
        
        for expected_angle in expected_angles:
            tolerance = 15.0  # degrees
            if abs(angle - expected_angle) <= tolerance:
                score = 1.0 - abs(angle - expected_angle) / tolerance
                if score > best_score:
                    best_score = score
                    best_expected = expected_angle
        
        if best_score == 0.0:
            # Angle violation
            violation = {
                'type': 'angle',
                'actual_angle': angle,
                'expected_angles': expected_angles,
                'geometry': geometry.value
            }
            validation_results['angle_violations'].append(violation)
            validation_results['violations'].append(
                f"Angle violation: {angle:.1f}° does not match {geometry.value} geometry"
            )
        
        return best_score
    
    def _get_expected_angles(self, geometry: CoordinationGeometry) -> List[float]:
        """Get expected angles for a coordination geometry"""
        angle_map = {
            CoordinationGeometry.LINEAR: [180.0],
            CoordinationGeometry.TRIGONAL_PLANAR: [120.0],
            CoordinationGeometry.TETRAHEDRAL: [109.47],
            CoordinationGeometry.SQUARE_PLANAR: [90.0, 180.0],
            CoordinationGeometry.OCTAHEDRAL: [90.0, 180.0],
            CoordinationGeometry.TRIGONAL_BIPYRAMIDAL: [90.0, 120.0, 180.0]
        }
        
        return angle_map.get(geometry, [109.47])  # Default to tetrahedral
    
    def _calculate_angle(self, v1: np.ndarray, v2: np.ndarray) -> float:
        """Calculate angle between two vectors in degrees"""
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))


def create_metal_center_from_pdb(pdb_line: str, protein_structure: Dict) -> Optional[MetalCenter]:
    """Create a MetalCenter from a PDB HETATM line"""
    try:
        # Parse PDB line
        element = pdb_line[76:78].strip()
        if not element:
            element = pdb_line[12:16].strip()[0]
        
        # Check if it's a metal
        try:
            metal_type = MetalType(element.upper())
        except ValueError:
            return None
        
        # Extract coordinates
        x = float(pdb_line[30:38])
        y = float(pdb_line[38:46])
        z = float(pdb_line[46:54])
        coordinates = np.array([x, y, z])
        
        # Analyze coordination environment
        analyzer = MetalCoordinationAnalyzer()
        metal_atom = {
            'element': element.upper(),
            'x': x, 'y': y, 'z': z
        }
        
        return analyzer._analyze_metal_environment(metal_atom, protein_structure)
        
    except Exception as e:
        logger.error(f"Error parsing metal center from PDB line: {e}")
        return None