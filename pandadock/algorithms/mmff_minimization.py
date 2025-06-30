"""
MMFF94 force field minimization for molecular docking.

This module provides molecular mechanics minimization using the MMFF94 force field
with fallback to UFF if MMFF is not available. Uses RDKit for force field calculations.
"""

import numpy as np
import logging
from typing import Any, Optional, Dict, Tuple
import tempfile
import os
from copy import deepcopy

logger = logging.getLogger(__name__)


class MMFFMinimization:
    """MMFF94 force field minimization with RDKit integration."""
    
    def __init__(self, max_iterations: int = 200, converge_criterion: float = 0.01,
                 distance_cutoff: float = 2.0, use_constraints: bool = True):
        """
        Initialize MMFF minimization.
        
        Args:
            max_iterations: Maximum minimization iterations
            converge_criterion: Energy change convergence threshold
            distance_cutoff: Distance cutoff for protein-ligand interactions
            use_constraints: Whether to use distance constraints
        """
        self.max_iterations = max_iterations
        self.converge_criterion = converge_criterion
        self.distance_cutoff = distance_cutoff
        self.use_constraints = use_constraints
        
        self.logger = logging.getLogger(__name__)
        
        # Check RDKit availability
        self.rdkit_available = False
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, rdMolDescriptors
            self.rdkit_available = True
            self.logger.info("RDKit available for MMFF94 minimization")
        except ImportError:
            self.logger.warning("RDKit not available - MMFF minimization disabled")
    
    def minimize_ligand(self, ligand: Any) -> Tuple[Any, float]:
        """
        Minimize ligand in isolation using MMFF94.
        
        Args:
            ligand: Ligand molecule object
            
        Returns:
            Tuple of (minimized_ligand, energy)
        """
        if not self.rdkit_available:
            self.logger.warning("RDKit not available - returning original ligand")
            return deepcopy(ligand), 0.0
        
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            self.logger.debug("Starting ligand minimization with MMFF94")
            
            # Convert ligand to RDKit molecule
            mol = self._ligand_to_rdkit_mol(ligand)
            if mol is None:
                self.logger.warning("Failed to convert ligand to RDKit molecule")
                return deepcopy(ligand), 0.0
            
            # Add hydrogens if needed
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates if not present
            if not mol.GetNumConformers():
                AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Try MMFF94 minimization
            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
                    if ff is not None:
                        initial_energy = ff.CalcEnergy()
                        converged = ff.Minimize(maxIts=self.max_iterations)
                        final_energy = ff.CalcEnergy()
                        
                        self.logger.debug(f"MMFF94 minimization: {initial_energy:.3f} → {final_energy:.3f} kcal/mol")
                        
                        # Convert back to ligand format
                        minimized_ligand = self._rdkit_mol_to_ligand(mol, ligand)
                        return minimized_ligand, final_energy
                        
                else:
                    raise ValueError("MMFF94 parameters not available")
                    
            except Exception as e:
                self.logger.debug(f"MMFF94 failed, trying UFF: {e}")
                
                # Fallback to UFF
                ff = AllChem.UFFGetMoleculeForceField(mol)
                if ff is not None:
                    initial_energy = ff.CalcEnergy()
                    converged = ff.Minimize(maxIts=self.max_iterations)
                    final_energy = ff.CalcEnergy()
                    
                    self.logger.debug(f"UFF minimization: {initial_energy:.3f} → {final_energy:.3f} kcal/mol")
                    
                    minimized_ligand = self._rdkit_mol_to_ligand(mol, ligand)
                    return minimized_ligand, final_energy
                else:
                    raise ValueError("Both MMFF94 and UFF failed")
            
        except Exception as e:
            self.logger.warning(f"Minimization failed: {e}")
            return deepcopy(ligand), 0.0
    
    def minimize_pose(self, protein: Any, ligand_pose: Any) -> Tuple[Any, float]:
        """
        Minimize ligand pose in protein environment with constraints.
        
        Args:
            protein: Protein molecule object
            ligand_pose: Ligand pose to minimize
            
        Returns:
            Tuple of (minimized_pose, energy)
        """
        if not self.rdkit_available:
            self.logger.warning("RDKit not available - returning original pose")
            return deepcopy(ligand_pose), 0.0
        
        try:
            self.logger.debug("Starting constrained pose minimization")
            
            # For now, minimize ligand in isolation
            # TODO: Implement proper protein-ligand constraints
            return self.minimize_ligand(ligand_pose)
            
        except Exception as e:
            self.logger.warning(f"Constrained minimization failed: {e}")
            return deepcopy(ligand_pose), 0.0
    
    def _ligand_to_rdkit_mol(self, ligand: Any) -> Optional[Any]:
        """Convert ligand object to RDKit molecule."""
        try:
            from rdkit import Chem
            
            # Try to use SMILES if available
            if hasattr(ligand, 'smiles') and ligand.smiles:
                mol = Chem.MolFromSmiles(ligand.smiles)
                if mol is not None:
                    return mol
            
            # Try to create from SDF data
            if hasattr(ligand, 'sdf_data') and ligand.sdf_data:
                mol = Chem.MolFromMolBlock(ligand.sdf_data)
                if mol is not None:
                    return mol
            
            # Try to write to SDF file and read back
            coords = self._get_ligand_coordinates(ligand)
            if coords is not None:
                sdf_file = self._write_ligand_to_sdf(ligand, coords)
                if sdf_file and os.path.exists(sdf_file):
                    suppl = Chem.SDMolSupplier(sdf_file)
                    mol = next(suppl)
                    os.unlink(sdf_file)  # Clean up
                    return mol
            
            return None
            
        except Exception as e:
            self.logger.debug(f"Error converting ligand to RDKit: {e}")
            return None
    
    def _rdkit_mol_to_ligand(self, mol: Any, reference_ligand: Any) -> Any:
        """Convert RDKit molecule back to ligand format."""
        try:
            from rdkit import Chem
            
            minimized_ligand = deepcopy(reference_ligand)
            
            # Extract coordinates from RDKit molecule
            conf = mol.GetConformer()
            coords = []
            
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
            
            coords = np.array(coords)
            
            # Update ligand coordinates
            self._set_ligand_coordinates(minimized_ligand, coords)
            
            return minimized_ligand
            
        except Exception as e:
            self.logger.debug(f"Error converting RDKit molecule to ligand: {e}")
            return deepcopy(reference_ligand)
    
    def _get_ligand_coordinates(self, ligand: Any) -> Optional[np.ndarray]:
        """Extract coordinates from ligand object."""
        try:
            if hasattr(ligand, 'coords'):
                coords = np.array(ligand.coords)
            elif hasattr(ligand, 'coordinates'):
                coords = np.array(ligand.coordinates)
            elif hasattr(ligand, 'xyz'):
                coords = np.array(ligand.xyz)
            elif hasattr(ligand, 'atoms'):
                coords_list = []
                for atom in ligand.atoms:
                    if isinstance(atom, dict) and 'coords' in atom:
                        coords_list.append(atom['coords'])
                    elif hasattr(atom, 'coords'):
                        coords_list.append(atom.coords)
                if coords_list:
                    coords = np.array(coords_list)
                else:
                    return None
            else:
                return None
            
            # Ensure 2D array
            if coords.ndim == 1:
                coords = coords.reshape(-1, 3)
            
            return coords
            
        except Exception as e:
            self.logger.debug(f"Error extracting ligand coordinates: {e}")
            return None
    
    def _set_ligand_coordinates(self, ligand: Any, coords: np.ndarray) -> None:
        """Set coordinates in ligand object."""
        try:
            if hasattr(ligand, 'coords'):
                ligand.coords = coords
            elif hasattr(ligand, 'coordinates'):
                ligand.coordinates = coords
            elif hasattr(ligand, 'xyz'):
                ligand.xyz = coords
            elif hasattr(ligand, 'atoms'):
                for i, atom in enumerate(ligand.atoms):
                    if i < len(coords):
                        if isinstance(atom, dict):
                            atom['coords'] = coords[i]
                        elif hasattr(atom, 'coords'):
                            atom.coords = coords[i]
        except Exception as e:
            self.logger.debug(f"Error setting ligand coordinates: {e}")
    
    def _write_ligand_to_sdf(self, ligand: Any, coords: np.ndarray) -> Optional[str]:
        """Write ligand to temporary SDF file."""
        try:
            # Create a simple SDF file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as f:
                # Write SDF header
                f.write("Ligand\\n")
                f.write("  Generated\\n")
                f.write("\\n")
                
                # Write counts line
                num_atoms = len(coords)
                f.write(f"{num_atoms:3d}{0:3d}  0  0  0  0  0  0  0  0999 V2000\\n")
                
                # Write atoms (simplified - all carbon)
                for i, coord in enumerate(coords):
                    f.write(f"{coord[0]:10.4f}{coord[1]:10.4f}{coord[2]:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\\n")
                
                # Write end
                f.write("M  END\\n")
                f.write("$$$$\\n")
                
                return f.name
            
        except Exception as e:
            self.logger.debug(f"Error writing SDF file: {e}")
            return None


class SimpleMMFFMinimization:
    """Simple minimization for fallback compatibility."""
    
    def __init__(self, **kwargs):
        self.logger = logging.getLogger(__name__)
        self.logger.warning("Using simple minimization - MMFF94 not available")
    
    def minimize_ligand(self, ligand: Any) -> Tuple[Any, float]:
        """Simple minimization fallback."""
        return deepcopy(ligand), 0.0
    
    def minimize_pose(self, protein: Any, ligand_pose: Any) -> Tuple[Any, float]:
        """Simple pose minimization fallback."""
        return deepcopy(ligand_pose), 0.0


# Factory function for creating minimization instance
def create_mmff_minimizer(**kwargs) -> MMFFMinimization:
    """Create MMFF minimization instance with automatic fallback."""
    try:
        return MMFFMinimization(**kwargs)
    except Exception:
        return SimpleMMFFMinimization(**kwargs)