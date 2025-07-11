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
        # Convert score to realistic binding affinity
        # Lower scores indicate better binding (more negative ΔG)
        
        # Handle different score ranges for different algorithms
        if self.score < 0:
            # Physics/precise mode: negative scores (e.g., -5 to -4)
            # More negative scores = better binding
            # Map score range -6 to -3 to binding affinity range -12 to -6 kcal/mol
            score_clamped = max(-7.0, min(-2.0, self.score))
            min_score, max_score = -6.0, -3.0
            min_affinity, max_affinity = -12.0, -6.0
            
            if score_clamped <= min_score:
                binding_affinity = min_affinity
            elif score_clamped >= max_score:
                binding_affinity = max_affinity
            else:
                # Linear interpolation: more negative score → more negative ΔG
                slope = (max_affinity - min_affinity) / (max_score - min_score)
                binding_affinity = min_affinity + slope * (score_clamped - min_score)
        else:
            # ML/balanced mode: positive scores (e.g., 0.1 to 0.4)
            # Lower positive scores = better binding
            # Map score range 0.1-0.4 to binding affinity range -12 to -6 kcal/mol
            score_clamped = max(0.05, min(0.5, self.score))
            min_score, max_score = 0.1, 0.4
            min_affinity, max_affinity = -12.0, -6.0
            
            if score_clamped <= min_score:
                binding_affinity = min_affinity
            elif score_clamped >= max_score:
                binding_affinity = max_affinity
            else:
                # Linear interpolation
                slope = (max_affinity - min_affinity) / (max_score - min_score)
                binding_affinity = min_affinity + slope * (score_clamped - min_score)
        
        return binding_affinity
    
    def get_ic50(self, temperature: float = 298.15) -> float:
        """Calculate IC50 from binding affinity"""
        # ΔG = RT ln(Kd), so Kd = exp(ΔG/RT)
        # For favorable binding, ΔG < 0, so Kd will be < 1 M (nanomolar range)
        # For competitive inhibition, IC50 ≈ Ki ≈ Kd
        
        R = 1.987e-3  # kcal/(mol·K)
        binding_affinity = self.get_binding_affinity()  # Negative for favorable binding
        
        if binding_affinity >= 0:
            return float('inf')  # No binding
        
        # Calculate Kd (dissociation constant) in M
        # ΔG = RT ln(Kd) → Kd = exp(ΔG/RT)
        kd = np.exp(binding_affinity / (R * temperature))
        
        # Convert to nM (IC50 ≈ Kd for competitive inhibition)
        ic50_nM = kd * 1e9
        
        return ic50_nM
    
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
        
        # Parse ligand structure from file
        if ligand_file.endswith('.pdb'):
            self.ligand = self._parse_pdb_ligand(ligand_file)
        elif ligand_file.endswith('.sdf') or ligand_file.endswith('.mol'):
            self.ligand = self._parse_sdf_ligand(ligand_file)
        else:
            raise ValueError(f"Unsupported ligand file format: {ligand_file}")
        
        self.logger.info(f"Loaded ligand with {len(self.ligand['coordinates'])} atoms")
    
    def _parse_pdb_ligand(self, pdb_file: str) -> Dict[str, Any]:
        """Parse ligand coordinates from PDB file"""
        coordinates = []
        atom_types = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    # Extract coordinates and atom type
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_type = line[76:78].strip() or line[12:16].strip()[0]
                    
                    coordinates.append([x, y, z])
                    atom_types.append(atom_type)
        
        return {
            'coordinates': np.array(coordinates),
            'atom_types': atom_types,
            'name': pdb_file.split('/')[-1].replace('.pdb', '')
        }
    
    def _parse_sdf_ligand(self, sdf_file: str) -> Dict[str, Any]:
        """Parse ligand coordinates from SDF file"""
        coordinates = []
        atom_types = []
        
        with open(sdf_file, 'r') as f:
            lines = f.readlines()
        
        # Find the molecule block
        for i, line in enumerate(lines):
            if 'V2000' in line or 'V3000' in line:
                # Parse number of atoms
                num_atoms = int(line[:3].strip())
                
                # Parse atom coordinates and types
                for j in range(i + 1, i + 1 + num_atoms):
                    if j < len(lines):
                        parts = lines[j].split()
                        if len(parts) >= 4:
                            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                            atom_type = parts[3]
                            
                            coordinates.append([x, y, z])
                            atom_types.append(atom_type)
                break
        
        if not coordinates:
            raise ValueError(f"No valid coordinates found in SDF file: {sdf_file}")
        
        return {
            'coordinates': np.array(coordinates),
            'atom_types': atom_types,
            'name': sdf_file.split('/')[-1].replace('.sdf', '').replace('.mol', '')
        }
    
    def generate_conformers(self, num_conformers: int = 100) -> List[np.ndarray]:
        """Generate ligand conformers"""
        if not self.ligand:
            raise ValueError("Ligand must be prepared before generating conformers")
        
        self.logger.info(f"Generating {num_conformers} conformers")
        conformers = []
        base_coords = self.ligand['coordinates'].copy()
        
        # First conformer is the original structure
        conformers.append(base_coords)
        
        # Generate additional conformers with small perturbations
        for i in range(1, num_conformers):
            # Apply small random rotations and bond twists
            perturbed_coords = base_coords.copy()
            
            # Add small random noise to simulate conformational changes
            noise = np.random.normal(0, 0.3, perturbed_coords.shape)
            perturbed_coords += noise
            
            # Apply a random overall rotation
            center = np.mean(perturbed_coords, axis=0)
            centered_coords = perturbed_coords - center
            
            # Random rotation matrix
            angles = np.random.uniform(0, 2*np.pi, 3)
            from utils.math_utils import rotation_matrix
            rot_matrix = rotation_matrix(angles)
            rotated_coords = np.dot(centered_coords, rot_matrix.T) + center
            
            conformers.append(rotated_coords)
        
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
    
    def _generate_basic_bonds(self, coordinates: np.ndarray) -> List[Tuple[int, int, str]]:
        """Generate basic bonds based on atomic distances"""
        bonds = []
        
        # Typical bond distances (in Angstroms)
        # C-C: 1.4-1.6, C-O: 1.2-1.5, C-N: 1.3-1.5
        max_bond_distance = 1.8  # Maximum distance to consider a bond
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                distance = np.linalg.norm(coordinates[i] - coordinates[j])
                
                # If atoms are close enough, consider them bonded
                if distance <= max_bond_distance:
                    bonds.append((i, j, 'single'))
        
        return bonds
    
    def save_poses(self, poses: List[Pose], output_dir: str):
        """Save poses to output directory in multiple formats"""
        import os
        
        self.logger.info(f"Saving {len(poses)} poses to {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
        
        # Save individual poses as PDB files
        for i, pose in enumerate(poses):
            # Save as PDB
            pdb_filename = os.path.join(output_dir, f"pose_{i+1}_{pose.pose_id}.pdb")
            self._save_pose_as_pdb(pose, pdb_filename)
            
            # Save as SDF
            sdf_filename = os.path.join(output_dir, f"pose_{i+1}_{pose.pose_id}.sdf")
            self._save_pose_as_sdf(pose, sdf_filename)
        
        # Save summary file with all poses and scores
        summary_filename = os.path.join(output_dir, "poses_summary.csv")
        self._save_poses_summary(poses, summary_filename)
        
        # Save multi-SDF file with all poses
        multi_sdf_filename = os.path.join(output_dir, "all.ligands.sdf")
        self._save_poses_as_multi_sdf(poses, multi_sdf_filename)
    
    def _save_pose_as_pdb(self, pose: Pose, filename: str):
        """Save a single pose as PDB file"""
        with open(filename, 'w') as f:
            f.write("REMARK PandaDock generated pose\n")
            f.write(f"REMARK Pose ID: {pose.pose_id}\n")
            f.write(f"REMARK Score: {pose.score:.3f}\n")
            f.write(f"REMARK Energy: {pose.energy:.2f} kcal/mol\n")
            f.write(f"REMARK Confidence: {pose.confidence:.3f}\n")
            f.write(f"REMARK Binding Affinity: {pose.get_binding_affinity():.2f} kcal/mol\n")
            f.write(f"REMARK IC50: {pose.get_ic50():.1f} nM\n")
            
            # Get atom types from ligand if available
            atom_types = []
            if hasattr(self, 'ligand') and self.ligand and 'atom_types' in self.ligand:
                atom_types = self.ligand['atom_types']
            
            # Write coordinates with proper atom types
            for i, coord in enumerate(pose.coordinates):
                # Use actual atom type if available, otherwise default to carbon
                atom_type = atom_types[i] if i < len(atom_types) else 'C'
                atom_symbol = atom_type[:1]  # First character for element symbol
                
                f.write(f"HETATM{i+1:5d}  {atom_symbol:<3s} LIG A   1    "
                       f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                       f"  1.00 20.00           {atom_symbol:<2s}\n")
            f.write("END\n")
    
    def _save_pose_as_sdf(self, pose: Pose, filename: str):
        """Save a single pose as SDF file"""
        with open(filename, 'w') as f:
            f.write(f"{pose.pose_id}\n")
            f.write("  PandaDock generated structure\n")
            f.write("\n")
            
            # Get atom types and bonds from ligand if available
            atom_types = []
            bonds = []
            if hasattr(self, 'ligand') and self.ligand:
                if 'atom_types' in self.ligand:
                    atom_types = self.ligand['atom_types']
                if 'bonds' in self.ligand:
                    bonds = self.ligand['bonds']
            
            # If no bonds available, generate basic connectivity based on distance
            if not bonds and len(pose.coordinates) > 1:
                bonds = self._generate_basic_bonds(pose.coordinates)
            
            # Molecule block
            num_atoms = len(pose.coordinates)
            num_bonds = len(bonds)
            f.write(f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n")
            
            # Atom block with proper atom types
            for i, coord in enumerate(pose.coordinates):
                # Use actual atom type if available, otherwise default to carbon
                atom_type = atom_types[i] if i < len(atom_types) else 'C'
                f.write(f"{coord[0]:10.4f}{coord[1]:10.4f}{coord[2]:10.4f} {atom_type:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
            
            # Bond block
            for atom1, atom2, bond_type in bonds:
                # Convert bond type to SDF format
                bond_order = '1'
                if bond_type in ['single', '1']:
                    bond_order = '1'
                elif bond_type in ['double', '2']:
                    bond_order = '2'
                elif bond_type in ['triple', '3']:
                    bond_order = '3'
                elif bond_type in ['aromatic', 'ar']:
                    bond_order = '4'
                else:
                    bond_order = '1'  # Default to single bond
                
                f.write(f"{atom1+1:3d}{atom2+1:3d}  {bond_order}  0  0  0  0\n")
            
            # Properties
            f.write("M  END\n")
            f.write(f">  <Score>\n{pose.score:.6f}\n\n")
            f.write(f">  <Energy>\n{pose.energy:.6f}\n\n")
            f.write(f">  <Confidence>\n{pose.confidence:.6f}\n\n")
            f.write(f">  <Binding_Affinity>\n{pose.get_binding_affinity():.6f}\n\n")
            f.write(f">  <IC50_nM>\n{pose.get_ic50():.1f}\n\n")
            f.write("$$$$\n")
    
    def _save_poses_summary(self, poses: List[Pose], filename: str):
        """Save poses summary as CSV file"""
        with open(filename, 'w') as f:
            f.write("Rank,Pose_ID,Score,Energy,Confidence,Binding_Affinity,IC50_nM,Ligand_Efficiency,Clash_Score\n")
            for i, pose in enumerate(poses):
                num_heavy_atoms = len(pose.coordinates) if hasattr(pose, 'coordinates') and pose.coordinates is not None else 13
                ligand_efficiency = pose.calculate_ligand_efficiency(num_heavy_atoms)
                f.write(f"{i+1},{pose.pose_id},{pose.score:.6f},{pose.energy:.6f},"
                       f"{pose.confidence:.6f},{pose.get_binding_affinity():.6f},"
                       f"{pose.get_ic50():.1f},{ligand_efficiency:.6f},{pose.clash_score:.6f}\n")
    
    def _save_poses_as_multi_sdf(self, poses: List[Pose], filename: str):
        """Save all poses in a single multi-structure SDF file"""
        with open(filename, 'w') as f:
            for pose in poses:
                f.write(f"{pose.pose_id}\n")
                f.write("  PandaDock generated structure\n")
                f.write("\n")
                
                # Get atom types and bonds from ligand if available
                atom_types = []
                bonds = []
                if hasattr(self, 'ligand') and self.ligand:
                    if 'atom_types' in self.ligand:
                        atom_types = self.ligand['atom_types']
                    if 'bonds' in self.ligand:
                        bonds = self.ligand['bonds']
                
                # If no bonds available, generate basic bonds based on distances
                if not bonds:
                    bonds = self._generate_basic_bonds(pose.coordinates)
                
                # Molecule block
                num_atoms = len(pose.coordinates)
                num_bonds = len(bonds)
                f.write(f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n")
                
                # Atom block with proper atom types
                for i, coord in enumerate(pose.coordinates):
                    # Use actual atom type if available, otherwise default to carbon
                    atom_type = atom_types[i] if i < len(atom_types) else 'C'
                    f.write(f"{coord[0]:10.4f}{coord[1]:10.4f}{coord[2]:10.4f} {atom_type:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
                
                # Bond block
                for atom1, atom2, bond_type in bonds:
                    # Convert bond type to SDF format
                    bond_order = '1'
                    if bond_type in ['single', '1']:
                        bond_order = '1'
                    elif bond_type in ['double', '2']:
                        bond_order = '2'
                    elif bond_type in ['triple', '3']:
                        bond_order = '3'
                    elif bond_type in ['aromatic', 'ar']:
                        bond_order = '4'
                    else:
                        bond_order = '1'  # Default to single bond
                    
                    f.write(f"{atom1+1:3d}{atom2+1:3d}  {bond_order}  0  0  0  0\n")
                
                # Properties
                f.write("M  END\n")
                f.write(f">  <Score>\n{pose.score:.6f}\n\n")
                f.write(f">  <Energy>\n{pose.energy:.6f}\n\n")
                f.write(f">  <Confidence>\n{pose.confidence:.6f}\n\n")
                f.write(f">  <Binding_Affinity>\n{pose.get_binding_affinity():.6f}\n\n")
                f.write(f">  <IC50_nM>\n{pose.get_ic50():.1f}\n\n")
                f.write("$$$$\n")
    
    def get_engine_info(self) -> Dict[str, Any]:
        """Get information about the docking engine"""
        return {
            'engine_type': self.__class__.__name__,
            'version': '1.0.0',
            'config': self.config.to_dict()
        }