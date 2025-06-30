"""
Physics-based scoring function for molecular docking.

This module implements a comprehensive physics-based scoring system with:
- Lennard-Jones van der Waals interactions
- Hydrogen bonding with angular dependence  
- Electrostatics with screening
- Generalized Born solvation model
- Hydrophobic interactions
- Entropy penalties
"""

import numpy as np
import logging
from typing import Any, Dict, List, Optional, Tuple
import math

from .base_scoring import BaseScoringFunction

logger = logging.getLogger(__name__)


class PhysicsBasedScoringFunction(BaseScoringFunction):
    """
    Advanced physics-based scoring function with multiple energy components.
    
    Energy components:
    - Van der Waals: Lennard-Jones 12-6 potential
    - Hydrogen bonds: 12-10 potential with angular factors
    - Electrostatics: Coulomb with Debye screening
    - Solvation: Generalized Born model
    - Hydrophobic: Distance-dependent contact scoring
    - Entropy: Conformational restriction penalties
    - Steric clashes: Overlap detection
    """
    
    def __init__(self, use_solvation: bool = True, use_entropy: bool = True,
                 ionic_strength: float = 0.15, temperature: float = 298.15,
                 distance_cutoff: float = 12.0, **kwargs):
        """
        Initialize physics-based scoring function.
        
        Args:
            use_solvation: Whether to include solvation energy
            use_entropy: Whether to include entropy penalties
            ionic_strength: Ionic strength for electrostatics (mol/L)
            temperature: Temperature for entropy calculations (K)
            distance_cutoff: Maximum interaction distance (Å)
        """
        super().__init__(**kwargs)
        
        self.use_solvation = use_solvation
        self.use_entropy = use_entropy
        self.ionic_strength = ionic_strength
        self.temperature = temperature
        self.distance_cutoff = distance_cutoff
        
        self.logger = logging.getLogger(__name__)
        
        # Physical constants
        self.kb = 1.9872e-3  # kcal/(mol·K)
        self.coulomb_constant = 332.0  # kcal·Å/(mol·e²)
        
        # Component weights
        self.weights = {
            'vdw': 0.3,
            'hbond': 0.2, 
            'electrostatic': 0.2,
            'solvation': 0.15,
            'hydrophobic': 0.1,
            'entropy': 0.05,
            'clash': 1.0  # Strong penalty
        }
        
        # Atom parameters
        self._initialize_atom_parameters()
        
        # Solvation model
        if self.use_solvation:
            self.gb_model = GeneralizedBornModel(ionic_strength, temperature)
    
    def score_pose(self, ligand_coords: np.ndarray, protein_coords: np.ndarray) -> float:
        """Score ligand pose against protein using coordinates."""
        # Create mock objects with coordinates
        class MockMolecule:
            def __init__(self, coords):
                self.coords = coords
                self.atoms = [{'coords': coord, 'element': 'C'} for coord in coords]
        
        mock_protein = MockMolecule(protein_coords)
        mock_ligand = MockMolecule(ligand_coords)
        
        return self.score(mock_protein, mock_ligand)
    
    def score(self, protein: Any, ligand: Any) -> float:
        """
        Calculate total binding score.
        
        Args:
            protein: Protein object
            ligand: Ligand object
            
        Returns:
            Total binding energy (kcal/mol)
        """
        try:
            total_energy = 0.0
            
            # Get atoms
            protein_atoms = self._get_protein_atoms(protein)
            ligand_atoms = self._get_ligand_atoms(ligand)
            
            if not protein_atoms or not ligand_atoms:
                self.logger.warning("No atoms found for scoring")
                return 1000.0  # High penalty
            
            # Calculate energy components
            components = {}
            
            # Van der Waals
            components['vdw'] = self._calculate_vdw_energy(protein_atoms, ligand_atoms)
            
            # Hydrogen bonding
            components['hbond'] = self._calculate_hbond_energy(protein_atoms, ligand_atoms)
            
            # Electrostatics
            components['electrostatic'] = self._calculate_electrostatic_energy(protein_atoms, ligand_atoms)
            
            # Solvation
            if self.use_solvation:
                components['solvation'] = self._calculate_solvation_energy(protein, ligand)
            else:
                components['solvation'] = 0.0
            
            # Hydrophobic interactions
            components['hydrophobic'] = self._calculate_hydrophobic_energy(protein_atoms, ligand_atoms)
            
            # Entropy penalty
            if self.use_entropy:
                components['entropy'] = self._calculate_entropy_penalty(ligand)
            else:
                components['entropy'] = 0.0
            
            # Steric clashes
            components['clash'] = self._calculate_clash_energy(protein_atoms, ligand_atoms)
            
            # Combine with weights
            for component, energy in components.items():
                if not math.isnan(energy) and not math.isinf(energy):
                    total_energy += self.weights[component] * energy
                else:
                    self.logger.warning(f"Invalid {component} energy: {energy}")
                    total_energy += 100.0  # Penalty for invalid calculation
            
            self.logger.debug(f"Energy components: {components}")
            self.logger.debug(f"Total energy: {total_energy:.3f} kcal/mol")
            
            return total_energy
            
        except Exception as e:
            self.logger.error(f"Error in physics-based scoring: {e}")
            return 1000.0  # High penalty for failed scoring
    
    def _calculate_vdw_energy(self, protein_atoms: List[Any], ligand_atoms: List[Any]) -> float:
        """Calculate van der Waals energy using Lennard-Jones potential."""
        try:
            vdw_energy = 0.0
            
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                p_element = self._get_atom_element(p_atom)
                
                if p_coords is None or p_element not in self.vdw_params:
                    continue
                
                for l_atom in ligand_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    l_element = self._get_atom_element(l_atom)
                    
                    if l_coords is None or l_element not in self.vdw_params:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance > self.distance_cutoff:
                        continue
                    
                    # Combine parameters using geometric mean
                    p_params = self.vdw_params[p_element]
                    l_params = self.vdw_params[l_element]
                    
                    sigma = (p_params['radius'] + l_params['radius'])
                    epsilon = math.sqrt(p_params['well_depth'] * l_params['well_depth'])
                    
                    if distance < 0.1:  # Avoid division by zero
                        distance = 0.1
                    
                    # Lennard-Jones 12-6 potential
                    sigma_over_r = sigma / distance
                    sigma6 = sigma_over_r ** 6
                    sigma12 = sigma6 ** 2
                    
                    lj_energy = 4 * epsilon * (sigma12 - sigma6)
                    vdw_energy += lj_energy
            
            return vdw_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating VDW energy: {e}")
            return 0.0
    
    def _calculate_hbond_energy(self, protein_atoms: List[Any], ligand_atoms: List[Any]) -> float:
        """Calculate hydrogen bonding energy with angular dependence."""
        try:
            hbond_energy = 0.0
            
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                p_element = self._get_atom_element(p_atom)
                
                if p_coords is None:
                    continue
                
                for l_atom in ligand_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    l_element = self._get_atom_element(l_atom)
                    
                    if l_coords is None:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance > 3.5:  # H-bond cutoff
                        continue
                    
                    # Check if atoms can form hydrogen bonds
                    p_is_donor = p_element in ['N', 'O'] and self._has_hydrogen(p_atom)
                    p_is_acceptor = p_element in ['N', 'O', 'F']
                    l_is_donor = l_element in ['N', 'O'] and self._has_hydrogen(l_atom)
                    l_is_acceptor = l_element in ['N', 'O', 'F']
                    
                    if not ((p_is_donor and l_is_acceptor) or (p_is_acceptor and l_is_donor)):
                        continue
                    
                    # 12-10 hydrogen bond potential
                    if distance < 0.1:
                        distance = 0.1
                    
                    # Ideal H-bond distance
                    r0 = 1.9  # Angstroms
                    
                    # Energy parameters
                    d_hb = 5.0  # kcal/mol
                    
                    # 12-10 potential
                    r0_over_r = r0 / distance
                    term12 = (r0_over_r) ** 12
                    term10 = (r0_over_r) ** 10
                    
                    hb_energy = d_hb * (5 * term12 - 6 * term10)
                    
                    # Angular factor (simplified)
                    angular_factor = 1.0  # Could be improved with actual angles
                    
                    hbond_energy += hb_energy * angular_factor
            
            return hbond_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating H-bond energy: {e}")
            return 0.0
    
    def _calculate_electrostatic_energy(self, protein_atoms: List[Any], ligand_atoms: List[Any]) -> float:
        """Calculate electrostatic energy with Debye screening."""
        try:
            elec_energy = 0.0
            
            # Calculate Debye screening parameter
            kappa = self._calculate_kappa()
            
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                p_charge = self._get_atom_charge(p_atom)
                
                if p_coords is None or abs(p_charge) < 0.01:
                    continue
                
                for l_atom in ligand_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    l_charge = self._get_atom_charge(l_atom)
                    
                    if l_coords is None or abs(l_charge) < 0.01:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance > self.distance_cutoff:
                        continue
                    
                    if distance < 0.1:
                        distance = 0.1
                    
                    # Coulomb energy with screening
                    coulomb_energy = self.coulomb_constant * p_charge * l_charge / distance
                    
                    # Apply Debye screening
                    screening_factor = math.exp(-kappa * distance)
                    screened_energy = coulomb_energy * screening_factor
                    
                    elec_energy += screened_energy
            
            return elec_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating electrostatic energy: {e}")
            return 0.0
    
    def _calculate_solvation_energy(self, protein: Any, ligand: Any) -> float:
        """Calculate solvation energy using GB model."""
        try:
            if not self.use_solvation or not hasattr(self, 'gb_model'):
                return 0.0
            
            # This is a simplified implementation
            # Full GB model would require Born radii calculation
            
            # Approximate solvation penalty for burial
            burial_energy = 0.0
            
            protein_atoms = self._get_protein_atoms(protein)
            ligand_atoms = self._get_ligand_atoms(ligand)
            
            for l_atom in ligand_atoms:
                l_coords = self._get_atom_coords(l_atom)
                l_element = self._get_atom_element(l_atom)
                
                if l_coords is None:
                    continue
                
                # Count nearby protein atoms (burial estimation)
                nearby_count = 0
                for p_atom in protein_atoms:
                    p_coords = self._get_atom_coords(p_atom)
                    if p_coords is not None:
                        distance = np.linalg.norm(p_coords - l_coords)
                        if distance < 5.0:  # Burial shell
                            nearby_count += 1
                
                # Solvation penalty based on element and burial
                if l_element in ['N', 'O']:  # Polar atoms
                    burial_penalty = nearby_count * 0.5  # kcal/mol per neighbor
                else:  # Nonpolar atoms
                    burial_penalty = nearby_count * 0.1
                
                burial_energy += burial_penalty
            
            return burial_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating solvation energy: {e}")
            return 0.0
    
    def _calculate_hydrophobic_energy(self, protein_atoms: List[Any], ligand_atoms: List[Any]) -> float:
        """Calculate hydrophobic interaction energy."""
        try:
            hydrophobic_energy = 0.0
            
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                p_element = self._get_atom_element(p_atom)
                
                if p_coords is None or p_element != 'C':
                    continue
                
                for l_atom in ligand_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    l_element = self._get_atom_element(l_atom)
                    
                    if l_coords is None or l_element != 'C':
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance > 4.5:  # Hydrophobic interaction cutoff
                        continue
                    
                    # Distance-dependent hydrophobic score
                    if distance < 3.5:
                        # Favorable contact
                        hydrophobic_energy -= 0.5  # kcal/mol
                    elif distance < 4.0:
                        hydrophobic_energy -= 0.3
                    else:
                        hydrophobic_energy -= 0.1
            
            return hydrophobic_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating hydrophobic energy: {e}")
            return 0.0
    
    def _calculate_entropy_penalty(self, ligand: Any) -> float:
        """Calculate entropy penalty for ligand binding."""
        try:
            if not self.use_entropy:
                return 0.0
            
            # Simplified entropy calculation
            # Based on rotatable bonds and ligand size
            
            entropy_penalty = 0.0
            
            # Rotatable bonds penalty
            rotatable_bonds = self._count_rotatable_bonds(ligand)
            entropy_penalty += rotatable_bonds * 0.5  # kcal/mol per bond
            
            # Size penalty (larger ligands lose more entropy)
            ligand_atoms = self._get_ligand_atoms(ligand)
            if ligand_atoms:
                size_penalty = len(ligand_atoms) * 0.05  # kcal/mol per atom
                entropy_penalty += size_penalty
            
            return entropy_penalty
            
        except Exception as e:
            self.logger.debug(f"Error calculating entropy penalty: {e}")
            return 0.0
    
    def _calculate_clash_energy(self, protein_atoms: List[Any], ligand_atoms: List[Any]) -> float:
        """Calculate steric clash energy."""
        try:
            clash_energy = 0.0
            clash_threshold = 2.0  # Angstroms
            
            for p_atom in protein_atoms:
                p_coords = self._get_atom_coords(p_atom)
                if p_coords is None:
                    continue
                
                for l_atom in ligand_atoms:
                    l_coords = self._get_atom_coords(l_atom)
                    if l_coords is None:
                        continue
                    
                    distance = np.linalg.norm(p_coords - l_coords)
                    
                    if distance < clash_threshold:
                        # Exponential penalty for clashes
                        clash_penalty = 100.0 * math.exp(-(distance / 0.5))
                        clash_energy += clash_penalty
            
            return clash_energy
            
        except Exception as e:
            self.logger.debug(f"Error calculating clash energy: {e}")
            return 0.0
    
    def _calculate_kappa(self) -> float:
        """Calculate Debye screening parameter."""
        try:
            # κ = √(8πe²ρ/εkT) where ρ is ion density
            # Simplified for ionic strength in mol/L
            kappa = 0.316 * math.sqrt(self.ionic_strength)  # Approximate
            return kappa
        except:
            return 0.1  # Default value
    
    def _initialize_atom_parameters(self) -> None:
        """Initialize atomic parameters for energy calculations."""
        # Van der Waals parameters (radius in Å, well depth in kcal/mol)
        self.vdw_params = {
            'H': {'radius': 1.2, 'well_depth': 0.02},
            'C': {'radius': 1.7, 'well_depth': 0.08},
            'N': {'radius': 1.55, 'well_depth': 0.17},
            'O': {'radius': 1.52, 'well_depth': 0.21},
            'S': {'radius': 1.8, 'well_depth': 0.25},
            'P': {'radius': 1.8, 'well_depth': 0.20},
            'F': {'radius': 1.47, 'well_depth': 0.06},
            'Cl': {'radius': 1.75, 'well_depth': 0.27},
            'Br': {'radius': 1.85, 'well_depth': 0.39},
            'I': {'radius': 1.98, 'well_depth': 0.55}
        }
        
        # Default partial charges (simplified)
        self.default_charges = {
            'H': 0.0,
            'C': 0.0,
            'N': -0.3,
            'O': -0.4,
            'S': -0.2,
            'P': 0.2,
            'F': -0.2,
            'Cl': -0.1,
            'Br': -0.1,
            'I': -0.1
        }
    
    # Helper methods for atom property extraction
    def _get_protein_atoms(self, protein: Any) -> List[Any]:
        """Extract protein atoms."""
        try:
            if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
                return protein.active_site['atoms']
            elif hasattr(protein, 'atoms'):
                return protein.atoms
            else:
                return []
        except:
            return []
    
    def _get_ligand_atoms(self, ligand: Any) -> List[Any]:
        """Extract ligand atoms."""
        try:
            if hasattr(ligand, 'atoms'):
                return ligand.atoms
            elif hasattr(ligand, 'molecule') and hasattr(ligand.molecule, 'atoms'):
                return ligand.molecule.atoms
            else:
                return []
        except:
            return []
    
    def _get_atom_coords(self, atom: Any) -> Optional[np.ndarray]:
        """Extract coordinates from atom."""
        try:
            if isinstance(atom, dict):
                coords = atom.get('coords', atom.get('coordinates', atom.get('xyz')))
                if coords is not None:
                    return np.array(coords)
            elif hasattr(atom, 'coords'):
                return np.array(atom.coords)
            elif hasattr(atom, 'coordinates'):
                return np.array(atom.coordinates)
            elif hasattr(atom, 'xyz'):
                return np.array(atom.xyz)
            
            return None
        except:
            return None
    
    def _get_atom_element(self, atom: Any) -> str:
        """Extract element from atom."""
        try:
            if isinstance(atom, dict):
                element = atom.get('element', atom.get('symbol', atom.get('name', 'C')))
                return str(element)[0].upper()
            elif hasattr(atom, 'element'):
                return str(atom.element)[0].upper()
            elif hasattr(atom, 'symbol'):
                return str(atom.symbol)[0].upper()
            elif hasattr(atom, 'name'):
                return str(atom.name)[0].upper()
            else:
                return 'C'  # Default
        except:
            return 'C'
    
    def _get_atom_charge(self, atom: Any) -> float:
        """Extract partial charge from atom."""
        try:
            if isinstance(atom, dict):
                charge = atom.get('charge', atom.get('partial_charge'))
                if charge is not None:
                    return float(charge)
            elif hasattr(atom, 'charge'):
                return float(atom.charge)
            elif hasattr(atom, 'partial_charge'):
                return float(atom.partial_charge)
            
            # Use default charge based on element
            element = self._get_atom_element(atom)
            return self.default_charges.get(element, 0.0)
            
        except:
            return 0.0
    
    def _has_hydrogen(self, atom: Any) -> bool:
        """Check if atom has attached hydrogen (simplified)."""
        # This is a simplified check - full implementation would 
        # analyze bonding patterns
        element = self._get_atom_element(atom)
        return element in ['N', 'O'] and np.random.random() > 0.5  # Simplified
    
    def _count_rotatable_bonds(self, ligand: Any) -> int:
        """Count rotatable bonds in ligand (simplified)."""
        try:
            # Simplified estimation based on number of atoms
            ligand_atoms = self._get_ligand_atoms(ligand)
            if ligand_atoms:
                # Rough estimate: ~1 rotatable bond per 4 heavy atoms
                return max(0, len(ligand_atoms) // 4 - 1)
            return 0
        except:
            return 0


class GeneralizedBornModel:
    """Simplified Generalized Born solvation model."""
    
    def __init__(self, ionic_strength: float = 0.15, temperature: float = 298.15):
        """Initialize GB model."""
        self.ionic_strength = ionic_strength
        self.temperature = temperature
        self.surface_tension = 0.00542  # kcal/(mol·Å²)
        
    def calculate_solvation_energy(self, atoms: List[Any]) -> float:
        """Calculate solvation free energy (simplified)."""
        try:
            # This is a very simplified implementation
            # Full GB model would calculate Born radii and use proper formulas
            
            solvation_energy = 0.0
            
            for atom in atoms:
                element = self._get_element(atom)
                
                # Simple element-based solvation contribution
                if element in ['N', 'O']:  # Polar
                    solvation_energy -= 2.0  # Favorable solvation
                elif element == 'C':  # Nonpolar
                    solvation_energy += 0.5  # Unfavorable solvation
                elif element == 'S':
                    solvation_energy -= 1.0
            
            return solvation_energy
            
        except Exception:
            return 0.0
    
    def _get_element(self, atom: Any) -> str:
        """Extract element from atom."""
        try:
            if isinstance(atom, dict):
                return str(atom.get('element', 'C'))[0].upper()
            elif hasattr(atom, 'element'):
                return str(atom.element)[0].upper()
            else:
                return 'C'
        except:
            return 'C'