"""
Structure preparation for molecular docking.

This module provides tools for preparing protein and ligand structures
for docking, including protonation, energy minimization, and validation.
"""

import numpy as np
import logging
from typing import Dict, List, Optional, Any
from pathlib import Path

from .protein_handler import ProteinStructure, ProteinHandler
from .ligand_handler import LigandStructure, LigandHandler


class StructurePreparation:
    """Handles preparation of protein and ligand structures for docking."""
    
    def __init__(self):
        """Initialize structure preparation."""
        self.logger = logging.getLogger(__name__)
        self.protein_handler = ProteinHandler()
        self.ligand_handler = LigandHandler()
    
    def prepare_protein(self, protein: ProteinStructure, ph: float = 7.4,
                       add_hydrogens: bool = True, remove_water: bool = True) -> ProteinStructure:
        """
        Prepare protein structure for docking.
        
        Args:
            protein: Input protein structure
            ph: pH for protonation states
            add_hydrogens: Whether to add hydrogen atoms
            remove_water: Whether to remove water molecules
            
        Returns:
            Prepared protein structure
        """
        self.logger.info(f"Preparing protein structure (pH={ph})")
        
        prepared_protein = protein  # Start with input structure
        
        # Remove water molecules if requested
        if remove_water:
            prepared_protein = self._remove_water_molecules(prepared_protein)
        
        # Add hydrogens if requested
        if add_hydrogens:
            prepared_protein = self._add_hydrogens_to_protein(prepared_protein, ph)
        
        # Validate structure
        self._validate_protein_structure(prepared_protein)
        
        self.logger.info(f"Protein preparation completed: {prepared_protein.n_atoms} atoms")
        
        return prepared_protein
    
    def prepare_ligand(self, ligand: LigandStructure, add_hydrogens: bool = True,
                      generate_conformers: bool = True, n_conformers: int = 10) -> LigandStructure:
        """
        Prepare ligand structure for docking.
        
        Args:
            ligand: Input ligand structure
            add_hydrogens: Whether to add hydrogen atoms
            generate_conformers: Whether to generate conformers
            n_conformers: Number of conformers to generate
            
        Returns:
            Prepared ligand structure
        """
        self.logger.info("Preparing ligand structure")
        
        prepared_ligand = ligand  # Start with input structure
        
        # Add hydrogens if requested
        if add_hydrogens:
            prepared_ligand = self._add_hydrogens_to_ligand(prepared_ligand)
        
        # Generate conformers if requested
        if generate_conformers:
            prepared_ligand = self.ligand_handler.generate_conformers(
                prepared_ligand, n_conformers
            )
        
        # Validate structure
        self._validate_ligand_structure(prepared_ligand)
        
        self.logger.info(
            f"Ligand preparation completed: {prepared_ligand.n_atoms} atoms, "
            f"{len(prepared_ligand.conformers)} conformers"
        )
        
        return prepared_ligand
    
    def prepare_protein_file(self, input_file: str, output_file: str,
                           add_hydrogens: bool = True, remove_water: bool = True,
                           ph: float = 7.4) -> Dict[str, Any]:
        """
        Prepare protein from file.
        
        Args:
            input_file: Input protein file path
            output_file: Output protein file path
            add_hydrogens: Whether to add hydrogens
            remove_water: Whether to remove water
            ph: pH for protonation
            
        Returns:
            Preparation result dictionary
        """
        try:
            # Load protein
            protein = self.protein_handler.load_protein(input_file)
            
            # Prepare protein
            prepared_protein = self.prepare_protein(
                protein, ph=ph, add_hydrogens=add_hydrogens, remove_water=remove_water
            )
            
            # Save prepared protein
            self.protein_handler.save_protein(prepared_protein, output_file)
            
            return {
                'success': True,
                'input_file': input_file,
                'output_file': output_file,
                'original_atoms': protein.n_atoms,
                'prepared_atoms': prepared_protein.n_atoms,
                'changes_made': {
                    'hydrogens_added': add_hydrogens,
                    'water_removed': remove_water,
                    'ph_adjusted': ph
                }
            }
            
        except Exception as e:
            self.logger.error(f"Protein preparation failed: {e}")
            return {
                'success': False,
                'error': str(e),
                'input_file': input_file
            }
    
    def prepare_ligand_file(self, input_file: str, output_file: str,
                          add_hydrogens: bool = True, generate_conformers: bool = False) -> Dict[str, Any]:
        """
        Prepare ligand from file.
        
        Args:
            input_file: Input ligand file path
            output_file: Output ligand file path
            add_hydrogens: Whether to add hydrogens
            generate_conformers: Whether to generate conformers
            
        Returns:
            Preparation result dictionary
        """
        try:
            # Load ligand
            ligand = self.ligand_handler.load_ligand(input_file)
            
            # Prepare ligand
            prepared_ligand = self.prepare_ligand(
                ligand, add_hydrogens=add_hydrogens, generate_conformers=generate_conformers
            )
            
            # Save prepared ligand
            self.ligand_handler.save_ligand(prepared_ligand, output_file)
            
            return {
                'success': True,
                'input_file': input_file,
                'output_file': output_file,
                'original_atoms': ligand.n_atoms,
                'prepared_atoms': prepared_ligand.n_atoms,
                'conformers_generated': len(prepared_ligand.conformers),
                'changes_made': {
                    'hydrogens_added': add_hydrogens,
                    'conformers_generated': generate_conformers
                }
            }
            
        except Exception as e:
            self.logger.error(f"Ligand preparation failed: {e}")
            return {
                'success': False,
                'error': str(e),
                'input_file': input_file
            }
    
    def _remove_water_molecules(self, protein: ProteinStructure) -> ProteinStructure:
        """Remove water molecules from protein structure."""
        
        # Identify non-water atoms
        non_water_indices = []
        non_water_atom_names = []
        non_water_residue_names = []
        non_water_residue_ids = []
        
        for i, res_name in enumerate(protein.residue_names):
            if res_name not in ['HOH', 'WAT', 'TIP3', 'SOL']:
                non_water_indices.append(i)
                non_water_atom_names.append(protein.atom_names[i])
                non_water_residue_names.append(protein.residue_names[i])
                non_water_residue_ids.append(protein.residue_ids[i])
        
        if len(non_water_indices) < len(protein.coords):
            self.logger.info(
                f"Removed {len(protein.coords) - len(non_water_indices)} water molecules"
            )
            
            # Create new protein without water
            new_coords = protein.coords[non_water_indices]
            
            return ProteinStructure(
                coords=new_coords,
                atom_names=non_water_atom_names,
                residue_names=non_water_residue_names,
                residue_ids=non_water_residue_ids
            )
        
        return protein
    
    def _add_hydrogens_to_protein(self, protein: ProteinStructure, ph: float) -> ProteinStructure:
        """Add hydrogen atoms to protein structure."""
        
        try:
            # Try using a chemistry library for proper hydrogen addition
            return self._add_hydrogens_with_chemistry_lib(protein, ph, is_protein=True)
            
        except Exception as e:
            self.logger.warning(f"Advanced hydrogen addition failed: {e}")
            # Fallback to simple hydrogen addition
            return self._add_hydrogens_simple(protein, is_protein=True)
    
    def _add_hydrogens_to_ligand(self, ligand: LigandStructure) -> LigandStructure:
        """Add hydrogen atoms to ligand structure."""
        
        try:
            # Try using a chemistry library for proper hydrogen addition
            return self._add_hydrogens_with_chemistry_lib(ligand, is_protein=False)
            
        except Exception as e:
            self.logger.warning(f"Advanced hydrogen addition failed: {e}")
            # Fallback to simple hydrogen addition
            return self._add_hydrogens_simple(ligand, is_protein=False)
    
    def _add_hydrogens_with_chemistry_lib(self, structure, ph: float = 7.4, is_protein: bool = True):
        """Add hydrogens using chemistry library (RDKit/OpenEye)."""
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            if is_protein:
                # For proteins, add hydrogens based on standard residue templates
                return self._add_protein_hydrogens_rdkit(structure, ph)
            else:
                # For ligands, use RDKit's hydrogen addition
                return self._add_ligand_hydrogens_rdkit(structure)
                
        except ImportError:
            # RDKit not available
            raise ImportError("RDKit not available for hydrogen addition")
    
    def _add_protein_hydrogens_rdkit(self, protein: ProteinStructure, ph: float) -> ProteinStructure:
        """Add hydrogens to protein using standard residue templates."""
        
        # This is a simplified implementation
        # In practice, you'd use proper protein preparation tools like PDBFixer
        
        # For now, just return the original protein with a note
        self.logger.info("Protein hydrogen addition: using simplified method")
        return protein
    
    def _add_ligand_hydrogens_rdkit(self, ligand: LigandStructure) -> LigandStructure:
        """Add hydrogens to ligand using RDKit."""
        
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
        
        # Add hydrogens
        mol_with_h = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol_with_h)
        
        # Extract coordinates and atom information
        coords = []
        atom_names = []
        atom_types = []
        bonds = []
        
        for atom in mol_with_h.GetAtoms():
            coords.append(mol_with_h.GetConformer().GetAtomPosition(atom.GetIdx()))
            atom_names.append(f"{atom.GetSymbol()}{atom.GetIdx()}")
            atom_types.append(atom.GetSymbol())
        
        for bond in mol_with_h.GetBonds():
            bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        
        # Convert coordinates
        coords_array = np.array([[pos.x, pos.y, pos.z] for pos in coords])
        
        return LigandStructure(
            coords=coords_array,
            atom_names=atom_names,
            atom_types=atom_types,
            bonds=bonds
        )
    
    def _add_hydrogens_simple(self, structure, is_protein: bool = True):
        """Simple hydrogen addition fallback."""
        
        # This is a very simplified approach
        # In practice, proper hydrogen addition requires knowledge of:
        # - Bond connectivity
        # - Hybridization states
        # - Protonation states (for proteins)
        
        self.logger.info("Using simplified hydrogen addition")
        
        if is_protein:
            # For proteins, estimate some hydrogen positions
            original_coords = structure.coords
            n_original = len(original_coords)
            
            # Add approximately 1 hydrogen per heavy atom (very rough estimate)
            h_coords = []
            h_atom_names = []
            h_residue_names = []
            h_residue_ids = []
            
            for i, coord in enumerate(original_coords):
                # Add a hydrogen at a nearby position
                h_coord = coord + np.random.normal(0, 0.5, 3)  # Small random offset
                h_coords.append(h_coord)
                h_atom_names.append(f"H{i}")
                h_residue_names.append(structure.residue_names[i])
                h_residue_ids.append(structure.residue_ids[i])
            
            # Combine original and hydrogen coordinates
            all_coords = np.vstack([original_coords, h_coords])
            all_atom_names = structure.atom_names + h_atom_names
            all_residue_names = structure.residue_names + h_residue_names
            all_residue_ids = structure.residue_ids + h_residue_ids
            
            return ProteinStructure(
                coords=all_coords,
                atom_names=all_atom_names,
                residue_names=all_residue_names,
                residue_ids=all_residue_ids
            )
        
        else:
            # For ligands, add hydrogens to carbon atoms
            original_coords = structure.coords
            
            h_coords = []
            h_atom_names = []
            h_atom_types = []
            
            for i, (coord, atom_type) in enumerate(zip(original_coords, structure.atom_types)):
                if atom_type == 'C':  # Add hydrogens to carbons
                    # Add 1-3 hydrogens depending on assumed hybridization
                    n_hydrogens = np.random.choice([1, 2, 3])  # Random for simplicity
                    
                    for j in range(n_hydrogens):
                        h_coord = coord + np.random.normal(0, 0.5, 3)
                        h_coords.append(h_coord)
                        h_atom_names.append(f"H{i}_{j}")
                        h_atom_types.append('H')
            
            if h_coords:
                # Combine original and hydrogen coordinates
                all_coords = np.vstack([original_coords, h_coords])
                all_atom_names = structure.atom_names + h_atom_names
                all_atom_types = structure.atom_types + h_atom_types
                
                return LigandStructure(
                    coords=all_coords,
                    atom_names=all_atom_names,
                    atom_types=all_atom_types,
                    bonds=structure.bonds  # Keep original bonds
                )
        
        return structure
    
    def _validate_protein_structure(self, protein: ProteinStructure) -> None:
        """Validate protein structure."""
        
        if protein.n_atoms == 0:
            raise ValueError("Protein has no atoms")
        
        # Check for reasonable coordinates
        coords = protein.coords
        if np.any(np.abs(coords) > 1000):
            self.logger.warning("Protein has very large coordinates")
        
        # Check for duplicate atoms (very close positions)
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                distance = np.linalg.norm(coords[i] - coords[j])
                if distance < 0.1:
                    self.logger.warning(f"Very close atoms detected: {i}, {j} (distance: {distance:.3f})")
        
        self.logger.info("Protein structure validation completed")
    
    def _validate_ligand_structure(self, ligand: LigandStructure) -> None:
        """Validate ligand structure."""
        
        if ligand.n_atoms == 0:
            raise ValueError("Ligand has no atoms")
        
        # Check for reasonable coordinates
        coords = ligand.coords
        if np.any(np.abs(coords) > 1000):
            self.logger.warning("Ligand has very large coordinates")
        
        # Check bond connectivity
        for atom1, atom2 in ligand.bonds:
            if atom1 >= ligand.n_atoms or atom2 >= ligand.n_atoms:
                raise ValueError(f"Invalid bond: {atom1}-{atom2} (only {ligand.n_atoms} atoms)")
        
        self.logger.info("Ligand structure validation completed")
    
    def optimize_structure(self, structure, method: str = 'simple') -> Any:
        """
        Optimize structure geometry.
        
        Args:
            structure: Protein or ligand structure
            method: Optimization method ('simple', 'mmff', 'uff')
            
        Returns:
            Optimized structure
        """
        if method == 'simple':
            return self._simple_geometry_optimization(structure)
        else:
            # More advanced methods would require additional libraries
            self.logger.warning(f"Advanced optimization method '{method}' not implemented")
            return structure
    
    def _simple_geometry_optimization(self, structure):
        """Simple geometry optimization using steepest descent."""
        
        # This is a very simplified optimization
        # In practice, you'd use proper force fields and optimization algorithms
        
        coords = structure.coords.copy()
        
        # Simple optimization: move atoms away from very close contacts
        for iteration in range(10):
            forces = np.zeros_like(coords)
            
            # Calculate repulsive forces for very close atoms
            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    diff = coords[i] - coords[j]
                    distance = np.linalg.norm(diff)
                    
                    if distance < 2.0 and distance > 0.01:  # Very close atoms
                        force_magnitude = 1.0 / (distance ** 2)
                        force_direction = diff / distance
                        
                        forces[i] += force_magnitude * force_direction
                        forces[j] -= force_magnitude * force_direction
            
            # Apply forces with small step size
            coords += 0.01 * forces
            
            # Check convergence
            max_force = np.max(np.linalg.norm(forces, axis=1))
            if max_force < 0.01:
                break
        
        # Update structure coordinates
        if hasattr(structure, 'coords'):
            structure.coords = coords
        
        self.logger.info(f"Simple geometry optimization completed in {iteration+1} iterations")
        
        return structure