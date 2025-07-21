# -*- coding: utf-8 -*-
"""
Metal Docking Engine

Advanced docking engine specialized for metalloproteins with metal-containing
active sites. Handles coordination constraints, geometric requirements, and
metal-specific scoring functions.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple, Any, Union
import logging
import copy
from dataclasses import dataclass, field

from docking.base_engine import DockingEngine, Pose, GridBox
from docking.metal_coordination import (
    MetalCenter, MetalCoordinationAnalyzer, MetalConstraintValidator,
    CoordinationGeometry, MetalType, create_metal_center_from_pdb
)
from scoring.metal_scoring import MetalScoringFunction, MetalScoringParameters
from utils.math_utils import rotation_matrix, quaternion_to_matrix

logger = logging.getLogger(__name__)


@dataclass
class MetalPose(Pose):
    """Extended pose class for metal docking results"""
    # Metal-specific attributes
    metal_interactions: List[Dict[str, Any]] = field(default_factory=list)
    coordination_quality: Dict[str, float] = field(default_factory=dict)
    geometric_violations: List[Dict[str, Any]] = field(default_factory=list)
    metal_energy_breakdown: Dict[str, float] = field(default_factory=dict)
    coordinating_atoms: List[int] = field(default_factory=list)
    coordination_distances: List[float] = field(default_factory=list)
    coordination_angles: List[float] = field(default_factory=list)
    
    def get_metal_binding_affinity(self) -> float:
        """Get metal-corrected binding affinity"""
        base_affinity = super().get_binding_affinity()
        
        # Apply metal-specific corrections
        metal_correction = self.metal_energy_breakdown.get('total_metal_energy', 0.0)
        coordination_bonus = -2.0 * len(self.coordinating_atoms)  # Coordination stabilization
        
        return base_affinity + metal_correction + coordination_bonus


@dataclass  
class MetalDockingConfig:
    """Configuration for metal docking"""
    # Metal detection and analysis
    detect_metals_automatically: bool = True
    metal_detection_distance: float = 3.5  # Angstroms
    require_metal_coordination: bool = True
    min_coordinating_atoms: int = 1
    max_coordinating_atoms: int = 6
    
    # Coordination constraints
    enforce_geometric_constraints: bool = True
    geometric_constraint_weight: float = 2.0
    distance_tolerance: float = 0.3  # Angstroms
    angle_tolerance: float = 15.0    # degrees
    
    # Sampling parameters
    coordination_focused_sampling: bool = True
    metal_focused_exhaustiveness: int = 20
    coordination_cone_angle: float = 30.0  # degrees
    metal_proximity_bias: float = 2.0
    
    # Scoring parameters
    use_metal_scoring: bool = True
    metal_scoring_weight: float = 1.0
    coordination_energy_scaling: float = 1.0
    geometric_penalty_scaling: float = 1.0
    
    # Pose filtering
    filter_non_coordinating_poses: bool = True
    require_geometric_validity: bool = False
    min_coordination_score: float = 0.3
    max_geometric_violations: int = 2


class MetalDockingEngine(DockingEngine):
    """Specialized docking engine for metalloproteins"""
    
    def __init__(self, config, metal_config: Optional[MetalDockingConfig] = None):
        super().__init__(config)
        self.metal_config = metal_config or MetalDockingConfig()
        self.metal_centers: List[MetalCenter] = []
        self.metal_analyzer = MetalCoordinationAnalyzer()
        self.metal_validator = MetalConstraintValidator()
        self.metal_scorer: Optional[MetalScoringFunction] = None
        
        self.logger = logging.getLogger(self.__class__.__name__)
        
    def prepare_receptor(self, protein_file: str):
        """Prepare receptor and detect metal centers"""
        super().prepare_receptor(protein_file)
        
        # Parse protein structure to find metals
        self.logger.info("Analyzing protein structure for metal centers...")
        protein_structure = self._parse_protein_structure(protein_file)
        
        if self.metal_config.detect_metals_automatically:
            self.metal_centers = self.metal_analyzer.detect_metal_centers(protein_structure)
            self.logger.info(f"Detected {len(self.metal_centers)} metal centers")
            
            for i, metal_center in enumerate(self.metal_centers):
                self.logger.info(
                    f"Metal {i+1}: {metal_center.metal_type.value} at "
                    f"{metal_center.coordinates}, "
                    f"coordination: {metal_center.coordination_number}, "
                    f"geometry: {metal_center.geometry.value}"
                )
        
        # Initialize metal-specific scoring function
        if self.metal_centers and self.metal_config.use_metal_scoring:
            metal_params = MetalScoringParameters()
            metal_params.geometric_penalty_weight = self.metal_config.geometric_penalty_scaling
            self.metal_scorer = MetalScoringFunction(self.metal_centers, metal_params)
            self.logger.info("Initialized metal-specific scoring function")
    
    def _parse_protein_structure(self, protein_file: str) -> Dict[str, Any]:
        """Parse protein structure to extract atomic information"""
        atoms = []
        
        with open(protein_file, 'r') as f:
            for line_num, line in enumerate(f):
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        # Parse PDB line
                        atom = {
                            'line_number': line_num,
                            'record_type': line[:6].strip(),
                            'atom_id': int(line[6:11].strip()),
                            'atom_name': line[12:16].strip(),
                            'residue': line[17:20].strip(),
                            'chain': line[21:22].strip(),
                            'residue_number': int(line[22:26].strip()),
                            'x': float(line[30:38].strip()),
                            'y': float(line[38:46].strip()),
                            'z': float(line[46:54].strip()),
                            'occupancy': float(line[54:60].strip()) if line[54:60].strip() else 1.0,
                            'b_factor': float(line[60:66].strip()) if line[60:66].strip() else 20.0,
                            'element': line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0]
                        }
                        atoms.append(atom)
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"Error parsing line {line_num}: {e}")
                        continue
        
        return {'atoms': atoms}
    
    def dock(self, protein_file: str, ligand_file: str) -> List[MetalPose]:
        """Main metal docking method"""
        self.logger.info(f"Starting metal docking: {protein_file} + {ligand_file}")
        
        # Prepare receptor and ligand
        self.prepare_receptor(protein_file)
        self.prepare_ligand(ligand_file)
        
        if not self.metal_centers:
            self.logger.warning("No metal centers found, using standard docking")
            return self._fallback_to_standard_docking()
        
        # Generate metal-focused conformers
        if self.metal_config.coordination_focused_sampling:
            conformers = self._generate_metal_focused_conformers()
        else:
            conformers = self.generate_conformers(self.config.docking.num_poses * 5)
        
        self.logger.info(f"Generated {len(conformers)} conformers for metal docking")
        
        # Score and evaluate each conformer
        poses = []
        for i, conformer in enumerate(conformers):
            try:
                pose = self._evaluate_metal_conformer(conformer, i)
                if pose and self._is_valid_metal_pose(pose):
                    poses.append(pose)
            except Exception as e:
                self.logger.warning(f"Error evaluating conformer {i}: {e}")
                continue
        
        self.logger.info(f"Generated {len(poses)} valid metal poses")
        
        # Filter and rank poses
        if self.metal_config.filter_non_coordinating_poses:
            poses = self._filter_metal_poses(poses)
        
        # Sort by score and return top poses
        poses.sort(key=lambda x: x.score)
        final_poses = poses[:self.config.docking.num_poses]
        
        self.logger.info(f"Returning {len(final_poses)} final metal poses")
        return final_poses
    
    def _generate_metal_focused_conformers(self) -> List[np.ndarray]:
        """Generate conformers focused on metal coordination"""
        conformers = []
        base_coords = self.ligand['coordinates'].copy()
        
        # Identify potential coordinating atoms in ligand
        coordinating_indices = self._find_potential_coordinators()
        
        if not coordinating_indices:
            self.logger.warning("No potential coordinating atoms found in ligand")
            return self.generate_conformers(100)
        
        self.logger.info(f"Found {len(coordinating_indices)} potential coordinating atoms")
        
        # For each metal center, generate focused conformers
        for metal_center in self.metal_centers:
            metal_conformers = self._generate_conformers_for_metal(
                base_coords, metal_center, coordinating_indices
            )
            conformers.extend(metal_conformers)
        
        # Add some random conformers for diversity
        random_conformers = self.generate_conformers(50)
        conformers.extend(random_conformers)
        
        return conformers
    
    def _find_potential_coordinators(self) -> List[int]:
        """Find atoms in ligand that could coordinate to metals"""
        coordinating_indices = []
        atom_types = self.ligand.get('atom_types', [])
        
        for i, atom_type in enumerate(atom_types):
            if atom_type in ['N', 'O', 'S', 'P']:  # Common coordinating elements
                coordinating_indices.append(i)
        
        return coordinating_indices
    
    def _generate_conformers_for_metal(self, base_coords: np.ndarray, 
                                     metal_center: MetalCenter,
                                     coordinating_indices: List[int]) -> List[np.ndarray]:
        """Generate conformers targeting a specific metal center"""
        conformers = []
        metal_pos = metal_center.coordinates
        
        num_conformers = self.metal_config.metal_focused_exhaustiveness
        
        for _ in range(num_conformers):
            conformer = base_coords.copy()
            
            # Choose a random coordinating atom to target the metal
            target_atom_idx = np.random.choice(coordinating_indices)
            target_pos = conformer[target_atom_idx]
            
            # Calculate desired position for coordination
            desired_distance = self._get_target_coordination_distance(
                metal_center.metal_type, 
                self.ligand['atom_types'][target_atom_idx]
            )
            
            # Generate random direction within coordination cone
            direction = self._sample_coordination_direction(
                metal_center, target_pos
            )
            
            desired_pos = metal_pos + direction * desired_distance
            
            # Translate ligand to place coordinating atom at desired position
            translation = desired_pos - target_pos
            conformer += translation
            
            # Apply random rotation around the coordination bond
            rotation_angle = np.random.uniform(0, 2 * np.pi)
            rotation_axis = direction / np.linalg.norm(direction)
            rot_matrix = self._rotation_matrix_axis_angle(rotation_axis, rotation_angle)
            
            # Rotate around the metal-ligand bond
            centered_coords = conformer - metal_pos
            rotated_coords = np.dot(centered_coords, rot_matrix.T)
            conformer = rotated_coords + metal_pos
            
            # Add small random perturbations
            noise = np.random.normal(0, 0.1, conformer.shape)
            conformer += noise
            
            conformers.append(conformer)
        
        return conformers
    
    def _get_target_coordination_distance(self, metal_type: MetalType, 
                                        ligand_atom_type: str) -> float:
        """Get target coordination distance for metal-ligand pair"""
        distance_map = {
            (MetalType.ZN, 'N'): 2.1,
            (MetalType.ZN, 'O'): 2.0,
            (MetalType.ZN, 'S'): 2.3,
            (MetalType.FE, 'N'): 2.0,
            (MetalType.FE, 'O'): 1.9,
            (MetalType.FE, 'S'): 2.3,
            (MetalType.MG, 'N'): 2.2,
            (MetalType.MG, 'O'): 2.1,
            (MetalType.CA, 'O'): 2.4,
            (MetalType.CU, 'N'): 2.0,
            (MetalType.CU, 'O'): 1.9,
            (MetalType.CU, 'S'): 2.2
        }
        
        return distance_map.get((metal_type, ligand_atom_type), 2.1)
    
    def _sample_coordination_direction(self, metal_center: MetalCenter, 
                                     ligand_pos: np.ndarray) -> np.ndarray:
        """Sample direction for coordination within geometric constraints"""
        
        # Get existing coordination directions to avoid clashes
        existing_directions = []
        for coord_atom in metal_center.coordinating_atoms:
            direction = coord_atom['coordinates'] - metal_center.coordinates
            existing_directions.append(direction / np.linalg.norm(direction))
        
        # Sample random direction
        max_attempts = 100
        for _ in range(max_attempts):
            # Random unit vector
            direction = np.random.normal(0, 1, 3)
            direction = direction / np.linalg.norm(direction)
            
            # Check if direction is compatible with coordination geometry
            if self._is_valid_coordination_direction(
                direction, existing_directions, metal_center.geometry
            ):
                return direction
        
        # Fallback: return random direction if no valid one found
        direction = np.random.normal(0, 1, 3)
        return direction / np.linalg.norm(direction)
    
    def _is_valid_coordination_direction(self, new_direction: np.ndarray,
                                       existing_directions: List[np.ndarray],
                                       geometry: CoordinationGeometry) -> bool:
        """Check if coordination direction is geometrically valid"""
        
        cone_angle = self.metal_config.coordination_cone_angle
        
        for existing_dir in existing_directions:
            angle = np.degrees(np.arccos(
                np.clip(np.dot(new_direction, existing_dir), -1, 1)
            ))
            
            # Check minimum angle based on geometry
            min_angle = self._get_minimum_coordination_angle(geometry)
            
            if angle < min_angle - cone_angle or angle > min_angle + cone_angle:
                continue  # This direction might work
            else:
                return False  # Too close to existing direction
        
        return True
    
    def _get_minimum_coordination_angle(self, geometry: CoordinationGeometry) -> float:
        """Get minimum angle between coordination bonds for geometry"""
        angle_map = {
            CoordinationGeometry.LINEAR: 180.0,
            CoordinationGeometry.TRIGONAL_PLANAR: 120.0,
            CoordinationGeometry.TETRAHEDRAL: 109.5,
            CoordinationGeometry.SQUARE_PLANAR: 90.0,
            CoordinationGeometry.OCTAHEDRAL: 90.0,
            CoordinationGeometry.TRIGONAL_BIPYRAMIDAL: 90.0
        }
        return angle_map.get(geometry, 109.5)
    
    def _rotation_matrix_axis_angle(self, axis: np.ndarray, angle: float) -> np.ndarray:
        """Create rotation matrix for rotation around axis by angle"""
        axis = axis / np.linalg.norm(axis)
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        
        # Rodrigues' rotation formula
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        
        R = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
        return R
    
    def _evaluate_metal_conformer(self, conformer: np.ndarray, 
                                conformer_id: int) -> Optional[MetalPose]:
        """Evaluate a conformer for metal docking"""
        
        # Create base pose
        pose = MetalPose(
            coordinates=conformer,
            score=0.0,
            energy=0.0,
            pose_id=f"metal_pose_{conformer_id}",
            ligand_name=self.ligand.get('name', 'unknown')
        )
        
        # Calculate metal-specific score
        if self.metal_scorer:
            total_energy = self.metal_scorer.calculate_total_energy(
                conformer, 
                ligand_atom_types=self.ligand.get('atom_types')
            )
            pose.energy = total_energy
            pose.score = total_energy  # Lower is better
            
            # Get detailed metal analysis for all metal centers
            for metal_center in self.metal_centers:
                metal_report = self.metal_scorer.get_metal_scoring_report(
                    metal_center, conformer, self.ligand.get('atom_types')
                )
                
                # Extract metal-specific information
                pose.metal_interactions.extend(metal_report['interactions'])
                pose.coordination_quality.update(metal_report['coordination_quality'])
                pose.geometric_violations.extend(
                    metal_report['validation'].get('violations', [])
                )
                pose.metal_energy_breakdown.update(metal_report['energy_breakdown'])
                
                # Extract coordination details
                for interaction in metal_report['interactions']:
                    if interaction['type'] == 'metal_coordination':
                        pose.coordinating_atoms.append(interaction['ligand_atom'])
                        pose.coordination_distances.append(interaction['distance'])
        
        else:
            # Fallback to standard scoring
            pose.score = self._calculate_basic_score(conformer)
            pose.energy = pose.score
        
        # Calculate additional pose metrics
        pose.confidence = self._calculate_pose_confidence(pose)
        
        return pose
    
    def _calculate_basic_score(self, conformer: np.ndarray) -> float:
        """Basic scoring when metal scorer is not available"""
        score = 0.0
        
        # Simple distance-based scoring to metal centers
        for metal_center in self.metal_centers:
            min_distance = float('inf')
            for coord in conformer:
                distance = np.linalg.norm(coord - metal_center.coordinates)
                min_distance = min(min_distance, distance)
            
            # Favorable score for close approach to metal
            if min_distance < 3.0:
                score -= 2.0 * (3.0 - min_distance)
            elif min_distance < 5.0:
                score -= 1.0 * (5.0 - min_distance)
        
        return score
    
    def _calculate_pose_confidence(self, pose: MetalPose) -> float:
        """Calculate confidence score for metal pose"""
        confidence = 0.0
        
        # Coordination quality component
        coord_score = pose.coordination_quality.get('coordination_score', 0.0)
        geom_score = pose.coordination_quality.get('geometric_score', 0.0)
        confidence += 0.4 * coord_score + 0.3 * geom_score
        
        # Interaction count component
        if pose.coordinating_atoms:
            interaction_score = min(1.0, len(pose.coordinating_atoms) / 3.0)
            confidence += 0.2 * interaction_score
        
        # Violation penalty
        violation_penalty = min(0.1, len(pose.geometric_violations) * 0.05)
        confidence -= violation_penalty
        
        return max(0.0, min(1.0, confidence))
    
    def _is_valid_metal_pose(self, pose: MetalPose) -> bool:
        """Check if pose meets metal docking validity criteria"""
        
        # Must have some coordination if required
        if (self.metal_config.require_metal_coordination and 
            len(pose.coordinating_atoms) < self.metal_config.min_coordinating_atoms):
            return False
        
        # Check coordination count limits
        if len(pose.coordinating_atoms) > self.metal_config.max_coordinating_atoms:
            return False
        
        # Check coordination quality
        coord_score = pose.coordination_quality.get('coordination_score', 0.0)
        if coord_score < self.metal_config.min_coordination_score:
            return False
        
        # Check geometric violations
        if (self.metal_config.require_geometric_validity and
            len(pose.geometric_violations) > self.metal_config.max_geometric_violations):
            return False
        
        # Check if pose is within grid box
        if not self.grid_box.contains_point(np.mean(pose.coordinates, axis=0)):
            return False
        
        return True
    
    def _filter_metal_poses(self, poses: List[MetalPose]) -> List[MetalPose]:
        """Filter poses based on metal-specific criteria"""
        filtered_poses = []
        
        for pose in poses:
            # Apply all validity checks
            if self._is_valid_metal_pose(pose):
                # Additional quality checks
                if self._passes_quality_filters(pose):
                    filtered_poses.append(pose)
        
        self.logger.info(f"Filtered {len(poses)} poses to {len(filtered_poses)} quality poses")
        return filtered_poses
    
    def _passes_quality_filters(self, pose: MetalPose) -> bool:
        """Apply additional quality filters for metal poses"""
        
        # Minimum energy threshold
        if pose.energy > 10.0:  # Very unfavorable
            return False
        
        # Coordination distance checks
        for distance in pose.coordination_distances:
            if distance < 1.5 or distance > 4.0:  # Unrealistic distances
                return False
        
        # Must have reasonable confidence
        if pose.confidence < 0.2:
            return False
        
        return True
    
    def _fallback_to_standard_docking(self) -> List[MetalPose]:
        """Fallback to standard docking when no metals found"""
        self.logger.info("Performing standard docking fallback")
        
        # Generate conformers normally
        conformers = self.generate_conformers(100)
        
        # Create standard poses
        poses = []
        for i, conformer in enumerate(conformers):
            pose = MetalPose(
                coordinates=conformer,
                score=self._calculate_basic_score(conformer),
                energy=self._calculate_basic_score(conformer),
                pose_id=f"standard_pose_{i}",
                ligand_name=self.ligand.get('name', 'unknown'),
                confidence=0.5  # Default confidence
            )
            poses.append(pose)
        
        # Sort and return
        poses.sort(key=lambda x: x.score)
        return poses[:self.config.docking.num_poses]
    
    def score(self, pose: Pose) -> float:
        """Score a pose using metal-aware scoring"""
        if isinstance(pose, MetalPose) and self.metal_scorer:
            # Use metal-specific scoring
            return self.metal_scorer.calculate_total_energy(
                pose.coordinates,
                ligand_atom_types=self.ligand.get('atom_types')
            )
        else:
            # Basic scoring
            return self._calculate_basic_score(pose.coordinates)
    
    def get_metal_docking_report(self, poses: List[MetalPose]) -> Dict[str, Any]:
        """Generate comprehensive metal docking report"""
        
        if not poses:
            return {'error': 'No poses to analyze'}
        
        # Overall statistics
        coordination_counts = [len(pose.coordinating_atoms) for pose in poses]
        coordination_scores = [
            pose.coordination_quality.get('coordination_score', 0.0) 
            for pose in poses
        ]
        
        # Metal center analysis
        metal_analysis = []
        for i, metal_center in enumerate(self.metal_centers):
            analysis = {
                'metal_id': i,
                'metal_type': metal_center.metal_type.value,
                'coordinates': metal_center.coordinates.tolist(),
                'geometry': metal_center.geometry.value,
                'coordination_number': metal_center.coordination_number,
                'poses_coordinating': sum(
                    1 for pose in poses 
                    if any(
                        interaction.get('metal_type') == metal_center.metal_type.value
                        for interaction in pose.metal_interactions
                    )
                )
            }
            metal_analysis.append(analysis)
        
        # Best pose detailed analysis
        best_pose = poses[0] if poses else None
        best_pose_analysis = None
        
        if best_pose and self.metal_scorer:
            best_pose_analysis = {
                'pose_id': best_pose.pose_id,
                'score': best_pose.score,
                'energy_breakdown': best_pose.metal_energy_breakdown,
                'coordinating_atoms': best_pose.coordinating_atoms,
                'coordination_distances': best_pose.coordination_distances,
                'interactions': best_pose.metal_interactions,
                'geometric_violations': best_pose.geometric_violations,
                'coordination_quality': best_pose.coordination_quality
            }
        
        report = {
            'docking_summary': {
                'total_poses': len(poses),
                'metal_centers_detected': len(self.metal_centers),
                'poses_with_coordination': sum(
                    1 for pose in poses if pose.coordinating_atoms
                ),
                'average_coordination_count': np.mean(coordination_counts) if coordination_counts else 0,
                'average_coordination_score': np.mean(coordination_scores) if coordination_scores else 0,
                'best_score': poses[0].score if poses else None,
                'best_coordination_score': max(coordination_scores) if coordination_scores else 0
            },
            'metal_centers': metal_analysis,
            'best_pose': best_pose_analysis,
            'configuration': {
                'detect_metals_automatically': self.metal_config.detect_metals_automatically,
                'enforce_geometric_constraints': self.metal_config.enforce_geometric_constraints,
                'coordination_focused_sampling': self.metal_config.coordination_focused_sampling,
                'use_metal_scoring': self.metal_config.use_metal_scoring,
                'filter_non_coordinating_poses': self.metal_config.filter_non_coordinating_poses
            }
        }
        
        return report
    
    def save_metal_poses(self, poses: List[MetalPose], output_dir: str):
        """Save metal poses with additional metal-specific information"""
        import os
        
        # Use base class method for standard formats
        super().save_poses(poses, output_dir)
        
        # Save metal-specific analysis
        metal_report = self.get_metal_docking_report(poses)
        
        # Save metal docking report
        import json
        report_file = os.path.join(output_dir, "metal_docking_report.json")
        with open(report_file, 'w') as f:
            json.dump(metal_report, f, indent=2, default=str)
        
        # Save metal interactions summary
        interactions_file = os.path.join(output_dir, "metal_interactions.csv")
        with open(interactions_file, 'w') as f:
            f.write("Pose_ID,Metal_Type,Ligand_Atom,Element,Distance,Energy,Interaction_Type\n")
            
            for pose in poses:
                for interaction in pose.metal_interactions:
                    f.write(f"{pose.pose_id},"
                           f"{interaction.get('metal_type', '')},"
                           f"{interaction.get('ligand_atom', '')},"
                           f"{interaction.get('ligand_element', '')},"
                           f"{interaction.get('distance', 0.0):.3f},"
                           f"{interaction.get('energy', 0.0):.3f},"
                           f"{interaction.get('subtype', '')}\n")
        
        self.logger.info(f"Saved metal docking results to {output_dir}")
        self.logger.info(f"Metal report: {report_file}")
        self.logger.info(f"Interactions: {interactions_file}")


def create_metal_docking_engine(config, metal_config: Optional[MetalDockingConfig] = None) -> MetalDockingEngine:
    """Factory function to create metal docking engine"""
    return MetalDockingEngine(config, metal_config)