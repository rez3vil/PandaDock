"""
Base classes for docking engines
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field
import numpy as np
import logging


@dataclass
class Pose:
    """Represents a docked pose"""
    coordinates: np.ndarray  # 3D coordinates of ligand atoms
    score: float
    energy: float
    rmsd: Optional[float] = None
    
    # Detailed scoring breakdown
    vdw_energy: float = 0.0
    electrostatic_energy: float = 0.0
    hbond_energy: float = 0.0
    hydrophobic_energy: float = 0.0
    solvation_energy: float = 0.0
    entropy_energy: float = 0.0
    
    # Pose quality metrics
    clash_score: float = 0.0
    binding_site_coverage: float = 0.0
    ligand_efficiency: float = 0.0
    
    # Interactions
    hbond_interactions: List[Dict[str, Any]] = field(default_factory=list)
    hydrophobic_interactions: List[Dict[str, Any]] = field(default_factory=list)
    salt_bridge_interactions: List[Dict[str, Any]] = field(default_factory=list)
    
    # Metadata
    pose_id: str = ""
    ligand_name: str = ""
    confidence: float = 0.0
    flexible_residues: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        """Post-initialization calculations"""
        if self.pose_id == "":
            self.pose_id = f"pose_{id(self)}"
    
    def get_binding_affinity(self) -> float:
        """Calculate binding affinity (ΔG) from score"""
        # Convert score to binding affinity (kcal/mol)
        # This is a simplified conversion; real implementations would be more complex
        return self.score * 1.36  # Rough conversion factor
    
    def get_ic50(self, temperature: float = 298.15) -> float:
        """Calculate IC50 from binding affinity"""
        # ΔG = -RTln(Ki), assuming IC50 ≈ Ki
        R = 1.987e-3  # kcal/mol/K
        delta_g = self.get_binding_affinity()
        
        if delta_g >= 0:
            return float('inf')
        
        ki = np.exp(-delta_g / (R * temperature))
        return ki * 1e9  # Convert to nM
    
    def calculate_ligand_efficiency(self, num_heavy_atoms: int) -> float:
        """Calculate ligand efficiency (LE)"""
        if num_heavy_atoms == 0:
            return 0.0
        
        delta_g = self.get_binding_affinity()
        self.ligand_efficiency = delta_g / num_heavy_atoms
        return self.ligand_efficiency
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert pose to dictionary for serialization"""
        return {
            'pose_id': self.pose_id,
            'ligand_name': self.ligand_name,
            'score': self.score,
            'energy': self.energy,
            'rmsd': self.rmsd,
            'binding_affinity': self.get_binding_affinity(),
            'ic50': self.get_ic50(),
            'ligand_efficiency': self.ligand_efficiency,
            'coordinates': self.coordinates.tolist(),
            'energy_breakdown': {
                'vdw': self.vdw_energy,
                'electrostatic': self.electrostatic_energy,
                'hbond': self.hbond_energy,
                'hydrophobic': self.hydrophobic_energy,
                'solvation': self.solvation_energy,
                'entropy': self.entropy_energy
            },
            'quality_metrics': {
                'clash_score': self.clash_score,
                'binding_site_coverage': self.binding_site_coverage,
                'confidence': self.confidence
            },
            'interactions': {
                'hbonds': self.hbond_interactions,
                'hydrophobic': self.hydrophobic_interactions,
                'salt_bridges': self.salt_bridge_interactions
            },
            'flexible_residues': self.flexible_residues
        }


@dataclass
class GridBox:
    """Represents a docking grid box"""
    center: np.ndarray  # [x, y, z] coordinates
    size: np.ndarray    # [x, y, z] dimensions
    resolution: float = 0.375
    
    def contains_point(self, point: np.ndarray) -> bool:
        """Check if point is within grid box"""
        half_size = self.size / 2
        return np.all(np.abs(point - self.center) <= half_size)
    
    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get min and max bounds of grid box"""
        half_size = self.size / 2
        min_bounds = self.center - half_size
        max_bounds = self.center + half_size
        return min_bounds, max_bounds


class DockingEngine(ABC):
    """Abstract base class for docking engines"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        self.receptor = None
        self.ligand = None
        self.grid_box = None
        self.poses = []
        
        # Initialize grid box from config
        self.grid_box = GridBox(
            center=np.array([config.io.center_x, config.io.center_y, config.io.center_z]),
            size=np.array([config.io.size_x, config.io.size_y, config.io.size_z])
        )
    
    @abstractmethod
    def dock(self, protein_file: str, ligand_file: str) -> List[Pose]:
        """
        Main docking method - must be implemented by subclasses
        
        Args:
            protein_file: Path to protein PDB file
            ligand_file: Path to ligand file (SDF, MOL2, etc.)
            
        Returns:
            List of Pose objects sorted by score
        """
        pass
    
    @abstractmethod
    def score(self, pose: Pose) -> float:
        """
        Score a given pose - must be implemented by subclasses
        
        Args:
            pose: Pose object to score
            
        Returns:
            Score value (lower is better)
        """
        pass
    
    def prepare_receptor(self, protein_file: str):
        """Prepare receptor for docking"""
        self.logger.info(f"Preparing receptor: {protein_file}")
        # Implementation would load and prepare protein structure
        # This is a placeholder for the actual implementation
        pass
    
    def prepare_ligand(self, ligand_file: str):
        """Prepare ligand for docking"""
        self.logger.info(f"Preparing ligand: {ligand_file}")
        # Implementation would load and prepare ligand structure
        # This is a placeholder for the actual implementation
        pass
    
    def generate_conformers(self, num_conformers: int = 100) -> List[np.ndarray]:
        """Generate ligand conformers"""
        self.logger.info(f"Generating {num_conformers} conformers")
        # Placeholder implementation
        conformers = []
        for i in range(num_conformers):
            # Generate random conformer (placeholder)
            conformer = np.random.randn(10, 3)  # 10 atoms, 3D coordinates
            conformers.append(conformer)
        return conformers
    
    def optimize_pose(self, pose: Pose) -> Pose:
        """Optimize a pose using local minimization"""
        self.logger.debug(f"Optimizing pose {pose.pose_id}")
        # Placeholder implementation
        # In real implementation, this would perform energy minimization
        return pose
    
    def filter_poses(self, poses: List[Pose]) -> List[Pose]:
        """Filter poses based on quality criteria"""
        self.logger.info(f"Filtering {len(poses)} poses")
        
        # Remove poses with high clash scores
        filtered_poses = [pose for pose in poses if pose.clash_score < 5.0]
        
        # Sort by score and take top N
        filtered_poses.sort(key=lambda x: x.score)
        return filtered_poses[:self.config.docking.num_poses]
    
    def calculate_rmsd(self, pose1: Pose, pose2: Pose) -> float:
        """Calculate RMSD between two poses"""
        if pose1.coordinates.shape != pose2.coordinates.shape:
            return float('inf')
        
        diff = pose1.coordinates - pose2.coordinates
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def cluster_poses(self, poses: List[Pose], rmsd_threshold: float = 2.0) -> List[Pose]:
        """Cluster poses by RMSD and return representative poses"""
        if not poses:
            return []
        
        clustered_poses = []
        for pose in poses:
            is_new_cluster = True
            for representative in clustered_poses:
                if self.calculate_rmsd(pose, representative) < rmsd_threshold:
                    is_new_cluster = False
                    break
            
            if is_new_cluster:
                clustered_poses.append(pose)
        
        return clustered_poses
    
    def validate_pose(self, pose: Pose) -> bool:
        """Validate that a pose is reasonable"""
        # Check if pose is within grid box
        if not self.grid_box.contains_point(np.mean(pose.coordinates, axis=0)):
            return False
        
        # Check for reasonable energy
        if pose.energy > 1000 or pose.energy < -1000:
            return False
        
        # Check for reasonable clash score
        if pose.clash_score > 10.0:
            return False
        
        return True
    
    def save_poses(self, poses: List[Pose], output_dir: str):
        """Save poses to output directory"""
        self.logger.info(f"Saving {len(poses)} poses to {output_dir}")
        # Implementation would save poses in various formats
        # This is a placeholder for the actual implementation
        pass
    
    def get_engine_info(self) -> Dict[str, Any]:
        """Get information about the docking engine"""
        return {
            'engine_type': self.__class__.__name__,
            'version': '1.0.0',
            'config': self.config.to_dict()
        }