# -*- coding: utf-8 -*-
"""
Ligand preparation module for PandaDock
Handles ligand preprocessing, conformer generation, and 3D structure preparation
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging
from pathlib import Path
import json
import subprocess
import tempfile
import os


class LigandPreparer:
    """
    Ligand preparation system for docking
    
    Features:
    - Multiple input format support (SMILES, SDF, MOL2, PDB)
    - 3D structure generation from SMILES
    - Torsion angle identification
    - Conformer generation
    - Protonation state assignment
    - Charge calculation
    - Molecular property calculation
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Ligand properties
        self.ligand_data = {}
        self.coordinates = None
        self.atom_types = []
        self.bonds = []
        self.partial_charges = []
        self.rotatable_bonds = []
        
        # Conformer generation parameters
        self.max_conformers = 1000
        self.energy_threshold = 10.0  # kcal/mol
        self.rmsd_threshold = 0.5  # Angstroms
        
        # Supported file formats
        self.supported_formats = {'.smi', '.smiles', '.sdf', '.mol', '.mol2', '.pdb'}
        
        self.logger.info("Initialized LigandPreparer")
    
    def prepare_ligand(self, ligand_input: str) -> Dict[str, Any]:
        """
        Main ligand preparation method
        
        Args:
            ligand_input: Path to ligand file or SMILES string
            
        Returns:
            Dictionary containing prepared ligand data
        """
        self.logger.info(f"Preparing ligand: {ligand_input}")
        
        # Determine input type
        if self.is_file_path(ligand_input):
            ligand_data = self.prepare_from_file(ligand_input)
        else:
            # Assume SMILES string
            ligand_data = self.prepare_from_smiles(ligand_input)
        
        # Calculate molecular properties
        ligand_data['properties'] = self.calculate_molecular_properties(ligand_data)
        
        # Validate ligand
        if not self.validate_ligand(ligand_data):
            raise ValueError("Ligand validation failed")
        
        self.ligand_data = ligand_data
        return ligand_data
    
    def is_file_path(self, input_str: str) -> bool:
        """Check if input is a file path"""
        return Path(input_str).exists() or '.' in input_str
    
    def prepare_from_file(self, file_path: str) -> Dict[str, Any]:
        """Prepare ligand from file"""
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Ligand file not found: {file_path}")
        
        suffix = file_path.suffix.lower()
        
        if suffix not in self.supported_formats:
            raise ValueError(f"Unsupported file format: {suffix}")
        
        # Read file based on format
        if suffix in {'.smi', '.smiles'}:
            return self.read_smiles_file(file_path)
        elif suffix == '.sdf':
            return self.read_sdf_file(file_path)
        elif suffix == '.mol2':
            return self.read_mol2_file(file_path)
        elif suffix == '.pdb':
            return self.read_pdb_file(file_path)
        else:
            raise ValueError(f"Format {suffix} not yet implemented")
    
    def prepare_from_smiles(self, smiles: str) -> Dict[str, Any]:
        """Prepare ligand from SMILES string"""
        self.logger.info(f"Preparing from SMILES: {smiles}")
        
        # Generate 3D structure from SMILES
        coordinates = self.generate_3d_from_smiles(smiles)
        
        # Extract molecular information
        atom_types = self.extract_atom_types_from_smiles(smiles)
        bonds = self.extract_bonds_from_smiles(smiles)
        
        # Calculate properties
        ligand_data = {
            'smiles': smiles,
            'coordinates': coordinates,
            'atom_types': atom_types,
            'bonds': bonds,
            'partial_charges': self.calculate_partial_charges(smiles),
            'rotatable_bonds': self.identify_rotatable_bonds(bonds),
            'num_atoms': len(coordinates),
            'num_heavy_atoms': len([atom for atom in atom_types if atom != 'H']),
            'molecular_weight': self.calculate_molecular_weight(atom_types),
            'format': 'smiles'
        }
        
        return ligand_data
    
    def generate_3d_from_smiles(self, smiles: str) -> np.ndarray:
        """
        Generate 3D coordinates from SMILES
        
        This is a simplified implementation. In practice, you would use:
        - RDKit for conformer generation
        - OpenEye OMEGA for systematic conformer generation
        - Molecular mechanics optimization
        """
        self.logger.debug(f"Generating 3D coordinates from SMILES: {smiles}")
        
        # Placeholder implementation
        # In real implementation, this would use RDKit or similar
        
        # Estimate number of atoms from SMILES
        num_atoms = self.estimate_atom_count_from_smiles(smiles)
        
        # Generate random 3D coordinates (placeholder)
        coordinates = np.random.randn(num_atoms, 3) * 3.0
        
        # Apply some basic geometric constraints
        coordinates = self.apply_geometric_constraints(coordinates)
        
        return coordinates
    
    def estimate_atom_count_from_smiles(self, smiles: str) -> int:
        """Estimate number of atoms from SMILES string"""
        # Very simplified estimation
        # Count characters that typically represent atoms
        atom_chars = set('CNOPSFClBrI')
        atom_count = sum(1 for char in smiles if char in atom_chars)
        
        # Add implicit hydrogens (rough estimate)
        hydrogen_count = max(1, atom_count // 2)
        
        return atom_count + hydrogen_count
    
    def extract_atom_types_from_smiles(self, smiles: str) -> List[str]:
        """Extract atom types from SMILES string"""
        # Simplified extraction
        # In practice, would use proper SMILES parser
        
        atom_types = []
        i = 0
        while i < len(smiles):
            char = smiles[i]
            
            # Check for two-letter atoms
            if i + 1 < len(smiles):
                two_char = smiles[i:i+2]
                if two_char in {'Cl', 'Br'}:
                    atom_types.append(two_char)
                    i += 2
                    continue
            
            # Single letter atoms
            if char in 'CNOPSFHI':
                atom_types.append(char)
            
            i += 1
        
        # Add implicit hydrogens
        for _ in range(len(atom_types) // 2):
            atom_types.append('H')
        
        return atom_types
    
    def extract_bonds_from_smiles(self, smiles: str) -> List[Tuple[int, int, str]]:
        """Extract bond information from SMILES string"""
        # Simplified bond extraction
        # Returns list of (atom1, atom2, bond_type) tuples
        
        bonds = []
        
        # This is a very simplified implementation
        # In practice, would use proper SMILES parser
        
        # Generate some bonds based on atom count
        num_atoms = len(self.extract_atom_types_from_smiles(smiles))
        
        for i in range(num_atoms - 1):
            bonds.append((i, i + 1, 'single'))
        
        # Add some branching
        if num_atoms > 5:
            bonds.append((0, num_atoms - 1, 'single'))
        
        return bonds
    
    def read_smiles_file(self, file_path: Path) -> Dict[str, Any]:
        """Read SMILES file"""
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Take first non-empty line
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#'):
                # SMILES can have optional name
                parts = line.split()
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else file_path.stem
                
                ligand_data = self.prepare_from_smiles(smiles)
                ligand_data['name'] = name
                return ligand_data
        
        raise ValueError(f"No valid SMILES found in {file_path}")
    
    def read_sdf_file(self, file_path: Path) -> Dict[str, Any]:
        """Read SDF file"""
        self.logger.debug(f"Reading SDF file: {file_path}")
        
        # Simplified SDF reader
        # In practice, would use proper SDF parser
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        if len(lines) < 4:
            raise ValueError(f"Invalid SDF file: {file_path}")
        
        name = lines[0].strip()
        
        # Parse counts line
        counts_line = lines[3].strip()
        try:
            num_atoms = int(counts_line[:3])
            num_bonds = int(counts_line[3:6])
        except (ValueError, IndexError):
            raise ValueError(f"Invalid counts line in SDF: {counts_line}")
        
        # Parse atoms
        coordinates = []
        atom_types = []
        
        atoms_parsed = 0
        i = 4
        
        while atoms_parsed < num_atoms and i < len(lines):
            line = lines[i].strip()
            i += 1
            
            # Skip empty lines and SDF terminators/metadata
            if (not line or line.startswith('M  END') or line.startswith('$$$$') or 
                line.startswith('M  ') or line.startswith('A  ') or line.startswith('V  ')):
                continue
                
            parts = line.split()
            if len(parts) < 4:
                continue  # Skip invalid atom lines
            
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                atom_type = parts[3]
                
                coordinates.append([x, y, z])
                atom_types.append(atom_type)
                atoms_parsed += 1
            except (ValueError, IndexError):
                # Skip lines that can't be parsed as atom coordinates
                continue
        
        # Parse bonds - start from current position after atoms
        bonds = []
        bonds_parsed = 0
        
        while bonds_parsed < num_bonds and i < len(lines):
            line = lines[i].strip()
            i += 1
            
            # Skip empty lines and SDF terminators/metadata
            if (not line or line.startswith('M  END') or line.startswith('$$$$') or 
                line.startswith('M  ') or line.startswith('A  ') or line.startswith('V  ')):
                continue
                
            parts = line.split()
            if len(parts) >= 3:
                try:
                    atom1 = int(parts[0]) - 1  # Convert to 0-based
                    atom2 = int(parts[1]) - 1
                    bond_type = parts[2]
                    bonds.append((atom1, atom2, bond_type))
                    bonds_parsed += 1
                except (ValueError, IndexError):
                    continue
        
        # Validate that we parsed the expected number of atoms
        if atoms_parsed != num_atoms:
            self.logger.warning(f"Expected {num_atoms} atoms but parsed {atoms_parsed} from {file_path}")
        
        ligand_data = {
            'name': name,
            'coordinates': np.array(coordinates),
            'atom_types': atom_types,
            'bonds': bonds,
            'partial_charges': self.calculate_partial_charges_from_structure(coordinates, atom_types),
            'rotatable_bonds': self.identify_rotatable_bonds(bonds),
            'num_atoms': num_atoms,
            'num_heavy_atoms': len([atom for atom in atom_types if atom != 'H']),
            'molecular_weight': self.calculate_molecular_weight(atom_types),
            'format': 'sdf'
        }
        
        return ligand_data
    
    def read_mol2_file(self, file_path: Path) -> Dict[str, Any]:
        """Read MOL2 file"""
        self.logger.debug(f"Reading MOL2 file: {file_path}")
        
        # Simplified MOL2 reader
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find sections
        atom_section = False
        bond_section = False
        
        coordinates = []
        atom_types = []
        bonds = []
        name = file_path.stem
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('@<TRIPOS>MOLECULE'):
                atom_section = False
                bond_section = False
                continue
            elif line.startswith('@<TRIPOS>ATOM'):
                atom_section = True
                bond_section = False
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                atom_section = False
                bond_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                atom_section = False
                bond_section = False
                continue
            
            if atom_section and line:
                parts = line.split()
                if len(parts) >= 6:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5].split('.')[0]  # Remove subtype
                    
                    coordinates.append([x, y, z])
                    atom_types.append(atom_type)
            
            elif bond_section and line:
                parts = line.split()
                if len(parts) >= 4:
                    atom1 = int(parts[1]) - 1  # Convert to 0-based
                    atom2 = int(parts[2]) - 1
                    bond_type = parts[3]
                    bonds.append((atom1, atom2, bond_type))
        
        ligand_data = {
            'name': name,
            'coordinates': np.array(coordinates),
            'atom_types': atom_types,
            'bonds': bonds,
            'partial_charges': self.calculate_partial_charges_from_structure(coordinates, atom_types),
            'rotatable_bonds': self.identify_rotatable_bonds(bonds),
            'num_atoms': len(coordinates),
            'num_heavy_atoms': len([atom for atom in atom_types if atom != 'H']),
            'molecular_weight': self.calculate_molecular_weight(atom_types),
            'format': 'mol2'
        }
        
        return ligand_data
    
    def read_pdb_file(self, file_path: Path) -> Dict[str, Any]:
        """Read PDB file (ligand only)"""
        self.logger.debug(f"Reading PDB file: {file_path}")
        
        coordinates = []
        atom_types = []
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    # Parse PDB atom line
                    atom_type = line[76:78].strip()
                    if not atom_type:
                        atom_type = line[12:16].strip()[0]  # Use first character of atom name
                    
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    coordinates.append([x, y, z])
                    atom_types.append(atom_type)
        
        if not coordinates:
            raise ValueError(f"No atoms found in PDB file: {file_path}")
        
        # Generate bonds (simplified)
        bonds = self.infer_bonds_from_coordinates(np.array(coordinates), atom_types)
        
        ligand_data = {
            'name': file_path.stem,
            'coordinates': np.array(coordinates),
            'atom_types': atom_types,
            'bonds': bonds,
            'partial_charges': self.calculate_partial_charges_from_structure(coordinates, atom_types),
            'rotatable_bonds': self.identify_rotatable_bonds(bonds),
            'num_atoms': len(coordinates),
            'num_heavy_atoms': len([atom for atom in atom_types if atom != 'H']),
            'molecular_weight': self.calculate_molecular_weight(atom_types),
            'format': 'pdb'
        }
        
        return ligand_data
    
    def infer_bonds_from_coordinates(self, coordinates: np.ndarray, atom_types: List[str]) -> List[Tuple[int, int, str]]:
        """Infer bonds from atomic coordinates"""
        bonds = []
        
        # Covalent radii (Angstroms)
        covalent_radii = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07,
            'S': 1.05, 'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33
        }
        
        for i in range(len(coordinates)):
            for j in range(i + 1, len(coordinates)):
                distance = np.linalg.norm(coordinates[i] - coordinates[j])
                
                # Get covalent radii
                radius1 = covalent_radii.get(atom_types[i], 1.0)
                radius2 = covalent_radii.get(atom_types[j], 1.0)
                
                # Check if atoms are bonded
                max_bond_distance = (radius1 + radius2) * 1.3  # 30% tolerance
                
                if distance < max_bond_distance:
                    bonds.append((i, j, 'single'))
        
        return bonds
    
    def calculate_partial_charges(self, smiles: str) -> List[float]:
        """Calculate partial charges from SMILES"""
        # Simplified charge calculation
        # In practice, would use methods like:
        # - Gasteiger charges
        # - AM1-BCC
        # - ESP charges from quantum calculations
        
        atom_types = self.extract_atom_types_from_smiles(smiles)
        charges = []
        
        for atom_type in atom_types:
            if atom_type == 'C':
                charges.append(0.0)
            elif atom_type == 'N':
                charges.append(-0.3)
            elif atom_type == 'O':
                charges.append(-0.4)
            elif atom_type == 'H':
                charges.append(0.1)
            else:
                charges.append(0.0)
        
        return charges
    
    def calculate_partial_charges_from_structure(self, coordinates: np.ndarray, atom_types: List[str]) -> List[float]:
        """Calculate partial charges from 3D structure"""
        # Simplified implementation
        charges = []
        
        for atom_type in atom_types:
            if atom_type == 'C':
                charges.append(0.0)
            elif atom_type == 'N':
                charges.append(-0.3)
            elif atom_type == 'O':
                charges.append(-0.4)
            elif atom_type == 'H':
                charges.append(0.1)
            else:
                charges.append(0.0)
        
        return charges
    
    def identify_rotatable_bonds(self, bonds: List[Tuple[int, int, str]]) -> List[int]:
        """Identify rotatable bonds"""
        # Simplified rotatable bond identification
        # In practice, would use proper chemical rules
        
        rotatable_bonds = []
        
        for i, (atom1, atom2, bond_type) in enumerate(bonds):
            # Single bonds are potentially rotatable
            if bond_type in ['single', '1']:
                rotatable_bonds.append(i)
        
        return rotatable_bonds
    
    def calculate_molecular_weight(self, atom_types: List[str]) -> float:
        """Calculate molecular weight"""
        atomic_weights = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'P': 30.974, 'S': 32.065, 'F': 18.998, 'Cl': 35.453,
            'Br': 79.904, 'I': 126.904
        }
        
        total_weight = 0.0
        for atom_type in atom_types:
            total_weight += atomic_weights.get(atom_type, 0.0)
        
        return total_weight
    
    def calculate_molecular_properties(self, ligand_data: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate molecular properties"""
        properties = {
            'molecular_weight': ligand_data['molecular_weight'],
            'num_atoms': ligand_data['num_atoms'],
            'num_heavy_atoms': ligand_data['num_heavy_atoms'],
            'num_rotatable_bonds': len(ligand_data['rotatable_bonds']),
            'formal_charge': sum(ligand_data['partial_charges']),
            'logp': self.estimate_logp(ligand_data['atom_types']),
            'hbd': self.count_hbd(ligand_data['atom_types']),
            'hba': self.count_hba(ligand_data['atom_types']),
            'tpsa': self.estimate_tpsa(ligand_data['atom_types']),
            'lipinski_violations': 0
        }
        
        # Check Lipinski's Rule of Five
        violations = 0
        if properties['molecular_weight'] > 500:
            violations += 1
        if properties['logp'] > 5:
            violations += 1
        if properties['hbd'] > 5:
            violations += 1
        if properties['hba'] > 10:
            violations += 1
        
        properties['lipinski_violations'] = violations
        
        return properties
    
    def estimate_logp(self, atom_types: List[str]) -> float:
        """Estimate logP using atomic contributions"""
        # Simplified logP estimation
        logp_contributions = {
            'C': 0.5, 'N': -0.8, 'O': -1.0, 'H': 0.0,
            'P': 0.0, 'S': 0.0, 'F': 0.0, 'Cl': 0.5, 'Br': 1.0
        }
        
        total_logp = 0.0
        for atom_type in atom_types:
            total_logp += logp_contributions.get(atom_type, 0.0)
        
        return total_logp
    
    def count_hbd(self, atom_types: List[str]) -> int:
        """Count hydrogen bond donors"""
        # Simplified HBD counting
        hbd_count = 0
        for atom_type in atom_types:
            if atom_type in ['N', 'O']:  # Simplified
                hbd_count += 1
        return hbd_count
    
    def count_hba(self, atom_types: List[str]) -> int:
        """Count hydrogen bond acceptors"""
        # Simplified HBA counting
        hba_count = 0
        for atom_type in atom_types:
            if atom_type in ['N', 'O']:  # Simplified
                hba_count += 1
        return hba_count
    
    def estimate_tpsa(self, atom_types: List[str]) -> float:
        """Estimate topological polar surface area"""
        # Simplified TPSA estimation
        tpsa_contributions = {
            'N': 23.79, 'O': 23.06, 'P': 13.59, 'S': 32.09
        }
        
        total_tpsa = 0.0
        for atom_type in atom_types:
            total_tpsa += tpsa_contributions.get(atom_type, 0.0)
        
        return total_tpsa
    
    def apply_geometric_constraints(self, coordinates: np.ndarray) -> np.ndarray:
        """Apply basic geometric constraints to coordinates"""
        # Ensure reasonable bond lengths
        # This is a simplified implementation
        
        if len(coordinates) < 2:
            return coordinates
        
        # Adjust coordinates to have reasonable distances
        for i in range(1, len(coordinates)):
            # Distance to previous atom
            distance = np.linalg.norm(coordinates[i] - coordinates[i-1])
            
            # Adjust if too close or too far
            if distance < 1.0:
                direction = coordinates[i] - coordinates[i-1]
                if np.linalg.norm(direction) > 0:
                    direction = direction / np.linalg.norm(direction)
                    coordinates[i] = coordinates[i-1] + direction * 1.4
            elif distance > 3.0:
                direction = coordinates[i] - coordinates[i-1]
                if np.linalg.norm(direction) > 0:
                    direction = direction / np.linalg.norm(direction)
                    coordinates[i] = coordinates[i-1] + direction * 1.8
        
        return coordinates
    
    def validate_ligand(self, ligand_data: Dict[str, Any]) -> bool:
        """Validate prepared ligand data"""
        # Check required fields
        required_fields = ['coordinates', 'atom_types', 'bonds', 'partial_charges']
        
        for field in required_fields:
            if field not in ligand_data:
                self.logger.error(f"Missing required field: {field}")
                return False
        
        # Check data consistency
        num_atoms = len(ligand_data['coordinates'])
        
        if len(ligand_data['atom_types']) != num_atoms:
            self.logger.error("Atom types count mismatch")
            return False
        
        if len(ligand_data['partial_charges']) != num_atoms:
            self.logger.error("Partial charges count mismatch")
            return False
        
        # Check for reasonable molecular weight
        if ligand_data['molecular_weight'] > 2000:
            self.logger.warning("Very high molecular weight")
        
        # Check for reasonable atom count
        if num_atoms > 300:
            self.logger.warning("Very large number of atoms")
        
        return True
    
    def save_ligand(self, ligand_data: Dict[str, Any], output_path: str, format: str = 'sdf'):
        """Save prepared ligand to file"""
        output_path = Path(output_path)
        
        if format == 'sdf':
            self.write_sdf(ligand_data, output_path)
        elif format == 'mol2':
            self.write_mol2(ligand_data, output_path)
        elif format == 'pdb':
            self.write_pdb(ligand_data, output_path)
        else:
            raise ValueError(f"Unsupported output format: {format}")
    
    def write_sdf(self, ligand_data: Dict[str, Any], output_path: Path):
        """Write ligand to SDF format"""
        with open(output_path, 'w') as f:
            # Header
            f.write(f"{ligand_data.get('name', 'ligand')}\n")
            f.write("  PandaDock\n")
            f.write("\n")
            
            # Counts line
            num_atoms = len(ligand_data['coordinates'])
            num_bonds = len(ligand_data['bonds'])
            f.write(f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000\n")
            
            # Atoms
            for i, (coord, atom_type) in enumerate(zip(ligand_data['coordinates'], ligand_data['atom_types'])):
                f.write(f"{coord[0]:10.4f}{coord[1]:10.4f}{coord[2]:10.4f} {atom_type:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
            
            # Bonds
            for atom1, atom2, bond_type in ligand_data['bonds']:
                bond_order = '1' if bond_type == 'single' else '2' if bond_type == 'double' else '1'
                f.write(f"{atom1+1:3d}{atom2+1:3d}  {bond_order}  0  0  0  0\n")
            
            f.write("M  END\n")
            f.write("$$$$\n")
    
    def write_mol2(self, ligand_data: Dict[str, Any], output_path: Path):
        """Write ligand to MOL2 format"""
        with open(output_path, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{ligand_data.get('name', 'ligand')}\n")
            f.write(f"{len(ligand_data['coordinates'])} {len(ligand_data['bonds'])} 0 0 0\n")
            f.write("SMALL\n")
            f.write("GASTEIGER\n")
            f.write("\n")
            
            f.write("@<TRIPOS>ATOM\n")
            for i, (coord, atom_type, charge) in enumerate(zip(
                ligand_data['coordinates'], 
                ligand_data['atom_types'], 
                ligand_data['partial_charges']
            )):
                f.write(f"{i+1:7d} {atom_type}{i+1:<8s} {coord[0]:9.4f} {coord[1]:9.4f} {coord[2]:9.4f} {atom_type}.3  1 LIG  {charge:8.4f}\n")
            
            f.write("@<TRIPOS>BOND\n")
            for i, (atom1, atom2, bond_type) in enumerate(ligand_data['bonds']):
                f.write(f"{i+1:6d} {atom1+1:5d} {atom2+1:5d} {bond_type}\n")
    
    def write_pdb(self, ligand_data: Dict[str, Any], output_path: Path):
        """Write ligand to PDB format"""
        with open(output_path, 'w') as f:
            f.write("HEADER    LIGAND\n")
            f.write("COMPND    LIGAND\n")
            
            for i, (coord, atom_type) in enumerate(zip(ligand_data['coordinates'], ligand_data['atom_types'])):
                f.write(f"HETATM{i+1:5d}  {atom_type:<4s} LIG A   1    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00 20.00          {atom_type:>2s}\n")
            
            f.write("END\n")
    
    def get_ligand_info(self) -> Dict[str, Any]:
        """Get information about the prepared ligand"""
        if not self.ligand_data:
            return {}
        
        return {
            'name': self.ligand_data.get('name', 'Unknown'),
            'format': self.ligand_data.get('format', 'Unknown'),
            'num_atoms': self.ligand_data['num_atoms'],
            'num_heavy_atoms': self.ligand_data['num_heavy_atoms'],
            'molecular_weight': self.ligand_data['molecular_weight'],
            'num_rotatable_bonds': len(self.ligand_data['rotatable_bonds']),
            'properties': self.ligand_data.get('properties', {})
        }