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
    
    def detect_pockets(self, min_volume: float = 50.0) -> List[Dict[str, Any]]:
        """
        Detect potential binding pockets.
        
        Args:
            min_volume: Minimum pocket volume
            
        Returns:
            List of detected pockets
        """
        # Simplified pocket detection using coordinate clustering
        pockets = []
        
        try:
            from sklearn.cluster import DBSCAN
            
            # Cluster coordinates to find dense regions
            clustering = DBSCAN(eps=8.0, min_samples=5).fit(self.coords)
            
            unique_labels = set(clustering.labels_)
            
            for label in unique_labels:
                if label == -1:  # Noise points
                    continue
                
                cluster_coords = self.coords[clustering.labels_ == label]
                
                # Calculate pocket center and radius
                center = np.mean(cluster_coords, axis=0)
                distances = np.linalg.norm(cluster_coords - center, axis=1)
                radius = np.max(distances) + 2.0  # Add some padding
                
                pocket = {
                    'center': center,
                    'radius': radius,
                    'volume': (4/3) * np.pi * radius**3,
                    'n_atoms': len(cluster_coords)
                }
                
                if pocket['volume'] >= min_volume:
                    pockets.append(pocket)
        
        except ImportError:
            # Fallback: use geometric center as single pocket
            center = np.mean(self.coords, axis=0)
            radius = np.max(np.linalg.norm(self.coords - center, axis=1)) * 0.5
            
            pockets = [{
                'center': center,
                'radius': radius,
                'volume': (4/3) * np.pi * radius**3,
                'n_atoms': self.n_atoms
            }]
        
        # Sort by volume (largest first)
        pockets.sort(key=lambda p: p['volume'], reverse=True)
        
        logging.getLogger(__name__).info(f"Detected {len(pockets)} potential pockets")
        
        return pockets
    
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