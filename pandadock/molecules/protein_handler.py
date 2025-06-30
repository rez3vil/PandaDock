"""
Protein structure handling and operations.

This module provides a unified interface for loading, manipulating,
and analyzing protein structures.
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path


class ProteinStructure:
    """Represents a protein structure with coordinates and metadata."""
    
    def __init__(self, coords: np.ndarray, atom_names: List[str] = None,
                 residue_names: List[str] = None, residue_ids: List[str] = None):
        """
        Initialize protein structure.
        
        Args:
            coords: Atomic coordinates (N, 3)
            atom_names: List of atom names
            residue_names: List of residue names
            residue_ids: List of residue IDs
        """
        self.coords = np.array(coords)
        self.n_atoms = len(self.coords)
        
        self.atom_names = atom_names or [f"ATOM_{i}" for i in range(self.n_atoms)]
        self.residue_names = residue_names or ["UNK"] * self.n_atoms
        self.residue_ids = residue_ids or [f"RES_{i}" for i in range(self.n_atoms)]
        
        # Derived properties
        self.residues = self._group_by_residue()
        self.active_site = None
        self.flexible_residues = []
    
    def _group_by_residue(self) -> Dict[str, List[Dict]]:
        """Group atoms by residue."""
        residues = {}
        
        for i, (atom_name, res_name, res_id) in enumerate(
            zip(self.atom_names, self.residue_names, self.residue_ids)
        ):
            if res_id not in residues:
                residues[res_id] = []
            
            residues[res_id].append({
                'atom_index': i,
                'atom_name': atom_name,
                'residue_name': res_name,
                'coords': self.coords[i],
                'residue_id': res_id
            })
        
        return residues
    
    def define_active_site(self, center: List[float], radius: float) -> None:
        """
        Define the active site for docking.
        
        Args:
            center: Center coordinates [x, y, z]
            radius: Radius in Angstroms
        """
        center = np.array(center)
        
        # Find residues within the active site
        site_residues = []
        for res_id, atoms in self.residues.items():
            for atom in atoms:
                distance = np.linalg.norm(atom['coords'] - center)
                if distance <= radius:
                    site_residues.append(res_id)
                    break
        
        self.active_site = {
            'center': center,
            'radius': radius,
            'residues': list(set(site_residues))
        }
        
        logging.getLogger(__name__).info(
            f"Defined active site: center={center}, radius={radius}, "
            f"residues={len(site_residues)}"
        )
    
    def detect_pockets(self, min_volume: float = 50.0, probe_radius: float = 1.4) -> List[Dict[str, Any]]:
        """
        Detect potential binding pockets using cavity detection.
        
        Args:
            min_volume: Minimum pocket volume
            probe_radius: Probe radius for cavity detection
            
        Returns:
            List of detected pockets
        """
        pockets = []
        
        try:
            # Method 1: Grid-based cavity detection
            pockets = self._detect_pockets_grid_based(min_volume, probe_radius)
            
            if not pockets:
                # Method 2: Fallback to geometric approach
                pockets = self._detect_pockets_geometric(min_volume)
        
        except Exception as e:
            logging.getLogger(__name__).warning(f"Pocket detection failed: {e}")
            # Final fallback: single pocket at center
            pockets = self._detect_pockets_simple()
        
        # Sort by volume (largest first)
        pockets.sort(key=lambda p: p['volume'], reverse=True)
        
        logging.getLogger(__name__).info(f"Detected {len(pockets)} potential pockets")
        
        return pockets
    
    def _detect_pockets_grid_based(self, min_volume: float, probe_radius: float) -> List[Dict[str, Any]]:
        """Grid-based cavity detection."""
        # Create a 3D grid around the protein
        min_coords, max_coords = self.get_bounding_box()
        
        # Expand bounding box
        padding = 5.0
        min_coords -= padding
        max_coords += padding
        
        # Grid resolution
        grid_spacing = 1.0
        nx = int((max_coords[0] - min_coords[0]) / grid_spacing) + 1
        ny = int((max_coords[1] - min_coords[1]) / grid_spacing) + 1
        nz = int((max_coords[2] - min_coords[2]) / grid_spacing) + 1
        
        # Find grid points that are in cavities
        cavity_points = []
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    point = np.array([
                        min_coords[0] + i * grid_spacing,
                        min_coords[1] + j * grid_spacing,
                        min_coords[2] + k * grid_spacing
                    ])
                    
                    # Check if point is in a cavity
                    if self._is_cavity_point(point, probe_radius):
                        cavity_points.append(point)
        
        if not cavity_points:
            return []
        
        # Cluster cavity points to find pockets
        try:
            from sklearn.cluster import DBSCAN
            
            cavity_points = np.array(cavity_points)
            clustering = DBSCAN(eps=3.0, min_samples=8).fit(cavity_points)
            
            pockets = []
            unique_labels = set(clustering.labels_)
            
            for label in unique_labels:
                if label == -1:  # Noise points
                    continue
                
                cluster_points = cavity_points[clustering.labels_ == label]
                
                # Calculate pocket properties
                center = np.mean(cluster_points, axis=0)
                distances = np.linalg.norm(cluster_points - center, axis=1)
                radius = np.max(distances) + probe_radius
                
                # Estimate volume (rough approximation)
                volume = len(cluster_points) * (grid_spacing ** 3)
                
                if volume >= min_volume:
                    pocket = {
                        'center': center,
                        'radius': radius,
                        'volume': volume,
                        'n_atoms': len(cluster_points),
                        'cavity_points': cluster_points
                    }
                    pockets.append(pocket)
            
            return pockets
            
        except ImportError:
            # Fallback without sklearn
            return self._detect_pockets_geometric(min_volume)
    
    def _is_cavity_point(self, point: np.ndarray, probe_radius: float) -> bool:
        """Check if a point is in a cavity (not too close to protein atoms)."""
        min_dist = np.min(np.linalg.norm(self.coords - point, axis=1))
        
        # Point is in cavity if it's at least probe_radius away from any atom
        # but not too far (should be within reasonable distance)
        return probe_radius < min_dist < 15.0
    
    def _detect_pockets_geometric(self, min_volume: float) -> List[Dict[str, Any]]:
        """Geometric approach to pocket detection."""
        # Find surface atoms (atoms with fewer neighbors)
        surface_atoms = self._find_surface_atoms()
        
        if not surface_atoms:
            return self._detect_pockets_simple()
        
        # Group surface atoms into potential pockets
        try:
            from sklearn.cluster import DBSCAN
            
            surface_coords = self.coords[surface_atoms]
            clustering = DBSCAN(eps=6.0, min_samples=3).fit(surface_coords)
            
            pockets = []
            unique_labels = set(clustering.labels_)
            
            for label in unique_labels:
                if label == -1:  # Noise points
                    continue
                
                cluster_indices = np.where(clustering.labels_ == label)[0]
                cluster_coords = surface_coords[cluster_indices]
                
                # Calculate pocket center and radius
                center = np.mean(cluster_coords, axis=0)
                distances = np.linalg.norm(cluster_coords - center, axis=1)
                radius = np.max(distances) + 3.0
                
                # Estimate volume
                volume = (4/3) * np.pi * radius**3 * 0.3  # Reduce by factor for cavity
                
                if volume >= min_volume:
                    pocket = {
                        'center': center,
                        'radius': radius,
                        'volume': volume,
                        'n_atoms': len(cluster_coords)
                    }
                    pockets.append(pocket)
            
            return pockets
            
        except ImportError:
            return self._detect_pockets_simple()
    
    def _find_surface_atoms(self, neighbor_cutoff: float = 8.0, min_neighbors: int = 8) -> List[int]:
        """Find surface atoms (atoms with fewer neighbors)."""
        surface_atoms = []
        
        for i, coord in enumerate(self.coords):
            distances = np.linalg.norm(self.coords - coord, axis=1)
            n_neighbors = np.sum(distances < neighbor_cutoff) - 1  # Exclude self
            
            if n_neighbors < min_neighbors:
                surface_atoms.append(i)
        
        return surface_atoms
    
    def _detect_pockets_simple(self) -> List[Dict[str, Any]]:
        """Simple fallback pocket detection."""
        center = np.mean(self.coords, axis=0)
        distances = np.linalg.norm(self.coords - center, axis=1)
        radius = np.percentile(distances, 75)  # Use 75th percentile as radius
        
        return [{
            'center': center,
            'radius': radius,
            'volume': (4/3) * np.pi * radius**3,
            'n_atoms': self.n_atoms
        }]
    
    def define_flexible_residues(self, residue_ids: List[str], 
                               max_rotatable_bonds: int = 3) -> None:
        """
        Define flexible residues for protein flexibility.
        
        Args:
            residue_ids: List of residue IDs to make flexible
            max_rotatable_bonds: Maximum rotatable bonds per residue
        """
        flexible_residues = []
        
        for res_id in residue_ids:
            if res_id in self.residues:
                flex_residue = {
                    'residue_id': res_id,
                    'atoms': self.residues[res_id],
                    'max_rotatable_bonds': max_rotatable_bonds
                }
                flexible_residues.append(flex_residue)
        
        self.flexible_residues = flexible_residues
        
        logging.getLogger(__name__).info(
            f"Defined {len(flexible_residues)} flexible residues"
        )
    
    def get_center_of_mass(self) -> np.ndarray:
        """Calculate center of mass."""
        return np.mean(self.coords, axis=0)
    
    def get_bounding_box(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get bounding box of the protein."""
        min_coords = np.min(self.coords, axis=0)
        max_coords = np.max(self.coords, axis=0)
        return min_coords, max_coords


class ProteinHandler:
    """Handles protein structure loading and manipulation."""
    
    def __init__(self):
        """Initialize protein handler."""
        self.logger = logging.getLogger(__name__)
    
    def load_protein(self, file_path: str) -> ProteinStructure:
        """
        Load protein structure from file.
        
        Args:
            file_path: Path to protein structure file
            
        Returns:
            ProteinStructure object
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Protein file not found: {file_path}")
        
        self.logger.info(f"Loading protein from {file_path}")
        
        if file_path.suffix.lower() == '.pdb':
            return self._load_pdb(file_path)
        else:
            raise ValueError(f"Unsupported protein file format: {file_path.suffix}")
    
    def _load_pdb(self, pdb_path: Path) -> ProteinStructure:
        """Load protein from PDB file."""
        
        coords = []
        atom_names = []
        residue_names = []
        residue_ids = []
        
        try:
            with open(pdb_path, 'r') as file:
                for line in file:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Parse PDB line
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        chain_id = line[21].strip()
                        res_num = line[22:26].strip()
                        
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        coords.append([x, y, z])
                        atom_names.append(atom_name)
                        residue_names.append(res_name)
                        residue_ids.append(f"{chain_id}_{res_num}")
            
            if not coords:
                raise ValueError("No coordinates found in PDB file")
            
            protein = ProteinStructure(
                coords=np.array(coords),
                atom_names=atom_names,
                residue_names=residue_names,
                residue_ids=residue_ids
            )
            
            self.logger.info(
                f"Loaded protein: {protein.n_atoms} atoms, "
                f"{len(protein.residues)} residues"
            )
            
            return protein
            
        except Exception as e:
            # Fallback: create minimal protein structure
            self.logger.warning(f"Failed to parse PDB file: {e}")
            self.logger.info("Creating minimal protein structure")
            
            return self._create_minimal_protein()
    
    def _create_minimal_protein(self) -> ProteinStructure:
        """Create a minimal protein structure for testing."""
        
        # Create a simple 3x3x3 grid of atoms
        coords = []
        for x in range(-1, 2):
            for y in range(-1, 2):
                for z in range(-1, 2):
                    coords.append([x * 5.0, y * 5.0, z * 5.0])
        
        n_atoms = len(coords)
        
        return ProteinStructure(
            coords=np.array(coords),
            atom_names=[f"CA{i}" for i in range(n_atoms)],
            residue_names=["ALA"] * n_atoms,
            residue_ids=[f"A_{i}" for i in range(n_atoms)]
        )
    
    def save_protein(self, protein: ProteinStructure, output_path: str) -> None:
        """
        Save protein structure to file.
        
        Args:
            protein: Protein structure object
            output_path: Output file path
        """
        output_path = Path(output_path)
        
        if output_path.suffix.lower() == '.pdb':
            self._save_pdb(protein, output_path)
        else:
            raise ValueError(f"Unsupported output format: {output_path.suffix}")
    
    def save_binding_site_sphere(self, protein: ProteinStructure, output_path: str) -> None:
        """
        Save binding site visualization as sphere.pdb file.
        
        Args:
            protein: Protein structure object
            output_path: Output file path for sphere.pdb
        """
        output_path = Path(output_path)
        
        if not protein.active_site:
            self.logger.warning("No active site defined, cannot generate sphere.pdb")
            return
        
        center = protein.active_site['center']
        radius = protein.active_site['radius']
        
        # Generate sphere points
        n_points = max(50, int(radius * 10))  # More points for larger spheres
        sphere_coords = self._generate_sphere_points(center, radius, n_points)
        
        # Save as PDB file
        with open(output_path, 'w') as f:
            f.write(f"REMARK  Binding site sphere centered at {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
            f.write(f"REMARK  Radius: {radius:.3f} Angstroms\n")
            f.write(f"REMARK  Generated by PandaDock for visualization\n")
            
            for i, coord in enumerate(sphere_coords):
                line = (
                    f"HETATM{i+1:5d}  SPH SPH A   1    "
                    f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00 20.00           C  \n"
                )
                f.write(line)
            
            f.write("END\n")
        
        self.logger.info(f"Saved binding site sphere to {output_path}")
    
    def _generate_sphere_points(self, center: np.ndarray, radius: float, n_points: int) -> np.ndarray:
        """Generate points on a sphere surface for visualization."""
        points = []
        
        # Use Fibonacci spiral for even distribution
        golden_ratio = (1 + 5**0.5) / 2
        
        for i in range(n_points):
            # Fibonacci spiral on sphere
            theta = 2 * np.pi * i / golden_ratio
            phi = np.arccos(1 - 2 * i / n_points)
            
            x = radius * np.sin(phi) * np.cos(theta)
            y = radius * np.sin(phi) * np.sin(theta)
            z = radius * np.cos(phi)
            
            points.append(center + np.array([x, y, z]))
        
        return np.array(points)
    
    def _save_pdb(self, protein: ProteinStructure, pdb_path: Path) -> None:
        """Save protein to PDB file."""
        
        with open(pdb_path, 'w') as file:
            for i in range(protein.n_atoms):
                atom_name = protein.atom_names[i]
                res_name = protein.residue_names[i]
                res_id = protein.residue_ids[i]
                
                # Parse chain and residue number
                if '_' in res_id:
                    chain_id, res_num = res_id.split('_', 1)
                else:
                    chain_id, res_num = 'A', res_id
                
                x, y, z = protein.coords[i]
                
                # Write PDB ATOM line
                line = (
                    f"ATOM  {i+1:5d} {atom_name:^4s} {res_name:3s} "
                    f"{chain_id:1s}{res_num:4s}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
                )
                file.write(line)
        
        self.logger.info(f"Saved protein structure to {pdb_path}")
    
    def calculate_protein_properties(self, protein: ProteinStructure) -> Dict[str, Any]:
        """
        Calculate basic protein properties.
        
        Args:
            protein: Protein structure
            
        Returns:
            Dictionary with protein properties
        """
        center = protein.get_center_of_mass()
        min_coords, max_coords = protein.get_bounding_box()
        
        dimensions = max_coords - min_coords
        volume = np.prod(dimensions)
        
        return {
            'n_atoms': protein.n_atoms,
            'n_residues': len(protein.residues),
            'center_of_mass': center.tolist(),
            'bounding_box': {
                'min': min_coords.tolist(),
                'max': max_coords.tolist(),
                'dimensions': dimensions.tolist()
            },
            'volume': float(volume),
            'has_active_site': protein.active_site is not None,
            'n_flexible_residues': len(protein.flexible_residues)
        }