"""
Ligand structure handling and operations.

This module provides a unified interface for loading, manipulating,
and analyzing ligand structures.
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path


class LigandStructure:
    """Represents a ligand structure with coordinates and metadata."""
    
    def __init__(self, coords: np.ndarray, atom_names: List[str] = None,
                 atom_types: List[str] = None, bonds: List[Tuple[int, int]] = None):
        """
        Initialize ligand structure.
        
        Args:
            coords: Atomic coordinates (N, 3)
            atom_names: List of atom names
            atom_types: List of atom types
            bonds: List of bonds as (atom1_idx, atom2_idx) tuples
        """
        self.coords = np.array(coords)
        self.n_atoms = len(self.coords)
        
        self.atom_names = atom_names or [f"ATOM_{i}" for i in range(self.n_atoms)]
        self.atom_types = atom_types or ["C"] * self.n_atoms
        self.bonds = bonds or []
        
        # Conformers storage
        self.conformers = [self.coords.copy()]
        self.current_conformer = 0
    
    def add_conformer(self, coords: np.ndarray) -> None:
        """Add a new conformer."""
        if coords.shape == self.coords.shape:
            self.conformers.append(coords.copy())
        else:
            raise ValueError("Conformer coordinates shape mismatch")
    
    def set_conformer(self, conformer_index: int) -> None:
        """Set the active conformer."""
        if 0 <= conformer_index < len(self.conformers):
            self.current_conformer = conformer_index
            self.coords = self.conformers[conformer_index].copy()
        else:
            raise IndexError("Conformer index out of range")
    
    def get_center_of_mass(self) -> np.ndarray:
        """Calculate center of mass."""
        return np.mean(self.coords, axis=0)
    
    def get_bounding_box(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get bounding box of the ligand."""
        min_coords = np.min(self.coords, axis=0)
        max_coords = np.max(self.coords, axis=0)
        return min_coords, max_coords
    
    def translate(self, translation: np.ndarray) -> None:
        """Translate the ligand."""
        self.coords += translation
    
    def rotate(self, rotation_matrix: np.ndarray, center: Optional[np.ndarray] = None) -> None:
        """Rotate the ligand around a center point."""
        if center is None:
            center = self.get_center_of_mass()
        
        # Center coordinates
        centered_coords = self.coords - center
        
        # Apply rotation
        rotated_coords = np.dot(centered_coords, rotation_matrix.T)
        
        # Restore position
        self.coords = rotated_coords + center


class LigandHandler:
    """Handles ligand structure loading and manipulation."""
    
    def __init__(self):
        """Initialize ligand handler."""
        self.logger = logging.getLogger(__name__)
    
    def load_ligand(self, file_path: str) -> LigandStructure:
        """
        Load ligand structure from file.
        
        Args:
            file_path: Path to ligand structure file
            
        Returns:
            LigandStructure object
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Ligand file not found: {file_path}")
        
        self.logger.info(f"Loading ligand from {file_path}")
        
        if file_path.suffix.lower() in ['.sdf', '.mol']:
            return self._load_sdf_mol(file_path)
        elif file_path.suffix.lower() == '.pdb':
            return self._load_pdb(file_path)
        else:
            raise ValueError(f"Unsupported ligand file format: {file_path.suffix}")
    
    def _load_sdf_mol(self, file_path: Path) -> LigandStructure:
        """Load ligand from SDF/MOL file."""
        
        try:
            # Try using RDKit if available
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            if file_path.suffix.lower() == '.sdf':
                supplier = Chem.SDMolSupplier(str(file_path))
                mol = next(supplier)
            else:
                mol = Chem.MolFromMolFile(str(file_path))
            
            if mol is None:
                raise ValueError("Could not parse molecule with RDKit")
            
            # Add hydrogens if needed
            mol = Chem.AddHs(mol)
            
            # Get coordinates
            conf = mol.GetConformer()
            coords = []
            atom_names = []
            atom_types = []
            
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                coords.append([pos.x, pos.y, pos.z])
                atom_names.append(atom.GetSymbol() + str(atom.GetIdx()))
                atom_types.append(atom.GetSymbol())
            
            # Get bonds
            bonds = []
            for bond in mol.GetBonds():
                bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            
            ligand = LigandStructure(
                coords=np.array(coords),
                atom_names=atom_names,
                atom_types=atom_types,
                bonds=bonds
            )
            
            self.logger.info(
                f"Loaded ligand with RDKit: {ligand.n_atoms} atoms, "
                f"{len(bonds)} bonds"
            )
            
            return ligand
            
        except ImportError:
            self.logger.warning("RDKit not available, using simple parser")
            return self._load_simple_mol(file_path)
    
    def _load_simple_mol(self, file_path: Path) -> LigandStructure:
        """Simple MOL file parser without RDKit."""
        
        coords = []
        atom_names = []
        atom_types = []
        bonds = []
        
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            # Skip header lines
            if len(lines) < 4:
                raise ValueError("Invalid MOL file format")
            
            # Parse counts line
            counts_line = lines[3].strip()
            n_atoms = int(counts_line[:3])
            n_bonds = int(counts_line[3:6]) if len(counts_line) >= 6 else 0
            
            # Parse atoms
            for i in range(4, 4 + n_atoms):
                if i >= len(lines):
                    break
                
                line = lines[i]
                x = float(line[0:10])
                y = float(line[10:20])
                z = float(line[20:30])
                atom_type = line[31:34].strip()
                
                coords.append([x, y, z])
                atom_names.append(f"{atom_type}{i-3}")
                atom_types.append(atom_type)
            
            # Parse bonds
            for i in range(4 + n_atoms, 4 + n_atoms + n_bonds):
                if i >= len(lines):
                    break
                
                line = lines[i]
                atom1 = int(line[0:3]) - 1  # Convert to 0-based
                atom2 = int(line[3:6]) - 1
                bonds.append((atom1, atom2))
            
            if not coords:
                raise ValueError("No coordinates found in MOL file")
            
            ligand = LigandStructure(
                coords=np.array(coords),
                atom_names=atom_names,
                atom_types=atom_types,
                bonds=bonds
            )
            
            self.logger.info(
                f"Loaded ligand with simple parser: {ligand.n_atoms} atoms, "
                f"{len(bonds)} bonds"
            )
            
            return ligand
            
        except Exception as e:
            self.logger.warning(f"Failed to parse MOL file: {e}")
            return self._create_minimal_ligand()
    
    def _load_pdb(self, pdb_path: Path) -> LigandStructure:
        """Load ligand from PDB file."""
        
        coords = []
        atom_names = []
        atom_types = []
        
        try:
            with open(pdb_path, 'r') as file:
                for line in file:
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_name = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        element = line[76:78].strip() or atom_name[0]
                        
                        coords.append([x, y, z])
                        atom_names.append(atom_name)
                        atom_types.append(element)
            
            if not coords:
                raise ValueError("No coordinates found in PDB file")
            
            ligand = LigandStructure(
                coords=np.array(coords),
                atom_names=atom_names,
                atom_types=atom_types
            )
            
            self.logger.info(f"Loaded ligand from PDB: {ligand.n_atoms} atoms")
            
            return ligand
            
        except Exception as e:
            self.logger.warning(f"Failed to parse PDB file: {e}")
            return self._create_minimal_ligand()
    
    def _create_minimal_ligand(self) -> LigandStructure:
        """Create a minimal ligand structure for testing."""
        
        # Create a simple linear molecule
        coords = []
        for i in range(5):
            coords.append([i * 1.5, 0.0, 0.0])
        
        return LigandStructure(
            coords=np.array(coords),
            atom_names=[f"C{i}" for i in range(5)],
            atom_types=["C"] * 5,
            bonds=[(i, i+1) for i in range(4)]
        )
    
    def save_ligand(self, ligand: LigandStructure, output_path: str) -> None:
        """
        Save ligand structure to file.
        
        Args:
            ligand: Ligand structure object
            output_path: Output file path
        """
        output_path = Path(output_path)
        
        if output_path.suffix.lower() in ['.sdf', '.mol']:
            self._save_sdf(ligand, output_path)
        elif output_path.suffix.lower() == '.pdb':
            self._save_pdb(ligand, output_path)
        else:
            raise ValueError(f"Unsupported output format: {output_path.suffix}")
    
    def _save_sdf(self, ligand: LigandStructure, sdf_path: Path) -> None:
        """Save ligand to SDF file."""
        
        try:
            # Try using RDKit if available
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Create RDKit molecule
            mol = Chem.RWMol()
            
            # Add atoms
            for i, atom_type in enumerate(ligand.atom_types):
                atom = Chem.Atom(atom_type)
                mol.AddAtom(atom)
            
            # Add bonds
            for atom1, atom2 in ligand.bonds:
                mol.AddBond(atom1, atom2, Chem.BondType.SINGLE)
            
            # Set coordinates
            conf = Chem.Conformer(ligand.n_atoms)
            for i, coord in enumerate(ligand.coords):
                conf.SetAtomPosition(i, coord.tolist())
            
            mol.AddConformer(conf)
            
            # Write to file
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(mol)
            writer.close()
            
        except ImportError:
            # Fallback to simple SDF format
            self._save_simple_sdf(ligand, sdf_path)
        
        self.logger.info(f"Saved ligand structure to {sdf_path}")
    
    def _save_simple_sdf(self, ligand: LigandStructure, sdf_path: Path) -> None:
        """Save ligand to simple SDF format."""
        
        with open(sdf_path, 'w') as file:
            # Header
            file.write("Ligand\n")
            file.write("  Generated by PandaDock\n")
            file.write("\n")
            
            # Counts line
            file.write(f"{ligand.n_atoms:3d}{len(ligand.bonds):3d}  0  0  0  0  0  0  0  0999 V2000\n")
            
            # Atoms
            for i, (coord, atom_type) in enumerate(zip(ligand.coords, ligand.atom_types)):
                x, y, z = coord
                file.write(f"{x:10.4f}{y:10.4f}{z:10.4f} {atom_type:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
            
            # Bonds
            for atom1, atom2 in ligand.bonds:
                file.write(f"{atom1+1:3d}{atom2+1:3d}  1  0  0  0  0\n")
            
            # End
            file.write("M  END\n$$$$\n")
    
    def _save_pdb(self, ligand: LigandStructure, pdb_path: Path) -> None:
        """Save ligand to PDB file."""
        
        with open(pdb_path, 'w') as file:
            for i in range(ligand.n_atoms):
                atom_name = ligand.atom_names[i]
                x, y, z = ligand.coords[i]
                
                line = (
                    f"HETATM{i+1:5d} {atom_name:^4s} LIG A   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           "
                    f"{ligand.atom_types[i]:>2s}  \n"
                )
                file.write(line)
        
        self.logger.info(f"Saved ligand structure to {pdb_path}")
    
    def generate_conformers(self, ligand: LigandStructure, n_conformers: int = 10) -> LigandStructure:
        """
        Generate conformers for the ligand.
        
        Args:
            ligand: Input ligand structure
            n_conformers: Number of conformers to generate
            
        Returns:
            Ligand with multiple conformers
        """
        try:
            # Try using RDKit for conformer generation
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Create RDKit molecule
            mol = Chem.RWMol()
            
            # Add atoms
            for atom_type in ligand.atom_types:
                atom = Chem.Atom(atom_type)
                mol.AddAtom(atom)
            
            # Add bonds
            for atom1, atom2 in ligand.bonds:
                mol.AddBond(atom1, atom2, Chem.BondType.SINGLE)
            
            # Generate conformers
            AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers)
            
            # Extract conformer coordinates
            for conf_id in range(mol.GetNumConformers()):
                conf = mol.GetConformer(conf_id)
                coords = []
                for i in range(mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    coords.append([pos.x, pos.y, pos.z])
                
                if conf_id == 0:
                    ligand.coords = np.array(coords)
                    ligand.conformers[0] = np.array(coords)
                else:
                    ligand.add_conformer(np.array(coords))
            
            self.logger.info(f"Generated {mol.GetNumConformers()} conformers with RDKit")
            
        except ImportError:
            # Fallback to simple conformer generation
            self._generate_simple_conformers(ligand, n_conformers)
        
        return ligand
    
    def _generate_simple_conformers(self, ligand: LigandStructure, n_conformers: int) -> None:
        """Generate simple conformers without RDKit."""
        
        # Simple conformer generation using random perturbations
        for i in range(n_conformers - 1):  # -1 because we already have the original
            # Apply random rotation and small perturbations
            perturbed_coords = ligand.coords.copy()
            
            # Random rotation
            angles = np.random.uniform(0, 2*np.pi, 3)
            rotation_matrix = self._euler_to_rotation_matrix(angles)
            
            center = np.mean(perturbed_coords, axis=0)
            centered_coords = perturbed_coords - center
            rotated_coords = np.dot(centered_coords, rotation_matrix.T)
            perturbed_coords = rotated_coords + center
            
            # Small random perturbations
            perturbations = np.random.normal(0, 0.2, perturbed_coords.shape)
            perturbed_coords += perturbations
            
            ligand.add_conformer(perturbed_coords)
        
        self.logger.info(f"Generated {len(ligand.conformers)} conformers with simple method")
    
    def _euler_to_rotation_matrix(self, angles: np.ndarray) -> np.ndarray:
        """Convert Euler angles to rotation matrix."""
        
        cos_a, sin_a = np.cos(angles), np.sin(angles)
        
        # Rotation matrices for each axis
        Rx = np.array([[1, 0, 0],
                       [0, cos_a[0], -sin_a[0]],
                       [0, sin_a[0], cos_a[0]]])
        
        Ry = np.array([[cos_a[1], 0, sin_a[1]],
                       [0, 1, 0],
                       [-sin_a[1], 0, cos_a[1]]])
        
        Rz = np.array([[cos_a[2], -sin_a[2], 0],
                       [sin_a[2], cos_a[2], 0],
                       [0, 0, 1]])
        
        # Combined rotation
        return np.dot(Rz, np.dot(Ry, Rx))
    
    def calculate_ligand_properties(self, ligand: LigandStructure) -> Dict[str, Any]:
        """
        Calculate basic ligand properties.
        
        Args:
            ligand: Ligand structure
            
        Returns:
            Dictionary with ligand properties
        """
        center = ligand.get_center_of_mass()
        min_coords, max_coords = ligand.get_bounding_box()
        
        dimensions = max_coords - min_coords
        
        return {
            'n_atoms': ligand.n_atoms,
            'n_bonds': len(ligand.bonds),
            'n_conformers': len(ligand.conformers),
            'center_of_mass': center.tolist(),
            'bounding_box': {
                'min': min_coords.tolist(),
                'max': max_coords.tolist(),
                'dimensions': dimensions.tolist()
            },
            'molecular_weight': self._estimate_molecular_weight(ligand),
            'atom_types': list(set(ligand.atom_types))
        }
    
    def _estimate_molecular_weight(self, ligand: LigandStructure) -> float:
        """Estimate molecular weight based on atom types."""
        
        # Approximate atomic weights
        atomic_weights = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
            'Br': 79.904, 'I': 126.904
        }
        
        total_weight = 0.0
        for atom_type in ligand.atom_types:
            weight = atomic_weights.get(atom_type, 12.011)  # Default to carbon
            total_weight += weight
        
        return total_weight