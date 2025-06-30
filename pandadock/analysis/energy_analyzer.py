"""
Energy decomposition and analysis for molecular docking results.

This module provides methods for decomposing binding energy into component
contributions and analyzing per-residue energy contributions.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Tuple, Union
from types import SimpleNamespace

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

logger = logging.getLogger(__name__)


class EnergyDecomposition:
    """Decompose binding energy into component contributions."""
    
    def __init__(self, scoring_function: Any):
        """
        Initialize energy decomposition analyzer.
        
        Args:
            scoring_function: Scoring function to use for decomposition
        """
        self.scoring_function = scoring_function
        self.logger = logging.getLogger(__name__)
    
    def decompose_energy(self, protein: Any, ligand: Any) -> Dict[str, float]:
        """
        Break down the binding energy into components.
        
        Args:
            protein: Protein object
            ligand: Ligand pose
            
        Returns:
            Energy components dictionary
        """
        try:
            self.logger.debug("Decomposing binding energy into components")
            
            # Initialize components dictionary
            components = {}
            
            # Calculate total energy first
            total_energy = self.scoring_function.score(protein, ligand)
            components['total'] = total_energy
            
            # Try to get individual energy components from scoring function
            component_methods = [
                ('vdw', ['_calculate_vdw', '_calculate_vdw_energy']),
                ('hbond', ['_calculate_hbond', '_calculate_hbond_energy']),
                ('electrostatic', ['_calculate_electrostatics', '_calculate_electrostatic_energy']),
                ('desolvation', ['_calculate_desolvation', '_calculate_desolvation_energy']),
                ('hydrophobic', ['_calculate_hydrophobic', '_calculate_hydrophobic_energy']),
                ('entropy', ['_calculate_entropy', '_calculate_entropy_energy']),
                ('clash', ['_calculate_clashes', '_calculate_clash_energy'])
            ]
            
            for comp_name, method_names in component_methods:
                energy = self._try_calculate_component(protein, ligand, method_names)
                if energy is not None:
                    components[comp_name] = energy
            
            # Compute remainder if some components are missing
            known_components = [v for k, v in components.items() if k != 'total']
            if known_components:
                known_sum = sum(known_components)
                if abs(known_sum - total_energy) > 0.1:
                    components['other'] = total_energy - known_sum
            
            self.logger.info(f"Energy decomposition: {len(components)-1} components identified")
            return components
            
        except Exception as e:
            self.logger.error(f"Error decomposing energy: {e}")
            return {'total': 0.0}
    
    def _try_calculate_component(self, protein: Any, ligand: Any, 
                                method_names: List[str]) -> Optional[float]:
        """Try to calculate an energy component using available methods."""
        for method_name in method_names:
            if hasattr(self.scoring_function, method_name):
                try:
                    method = getattr(self.scoring_function, method_name)
                    if method_name.endswith('_entropy'):
                        # Entropy methods typically take only ligand
                        return method(ligand)
                    else:
                        return method(protein, ligand)
                except Exception as e:
                    self.logger.debug(f"Method {method_name} failed: {e}")
                    continue
        return None
    
    def residue_contributions(self, protein: Any, ligand: Any, 
                            radius: float = 2.0) -> List[Tuple[str, float]]:
        """
        Calculate per-residue contributions to binding energy.
        
        Args:
            protein: Protein object
            ligand: Ligand pose
            radius: Radius around ligand to consider residues (Angstroms)
            
        Returns:
            List of (residue, energy) tuples, sorted by contribution
        """
        try:
            self.logger.debug(f"Calculating per-residue contributions within {radius} Å")
            
            # Extract appropriate protein atoms
            protein_atoms = self._get_protein_atoms(protein)
            if not protein_atoms:
                self.logger.warning("No protein atoms found")
                return []
            
            # Group atoms by residue
            residue_atoms = self._group_atoms_by_residue(protein_atoms)
            
            # Calculate energy contribution for each residue
            contributions = []
            ligand_atoms = self._get_ligand_atoms(ligand)
            
            for res_key, atoms in residue_atoms.items():
                # Check if any atom is near the ligand
                if not self._is_residue_near_ligand(atoms, ligand_atoms, radius):
                    continue
                
                # Create temporary protein with only this residue
                temp_protein = self._create_residue_protein(atoms)
                
                # Calculate energy
                try:
                    energy = self.scoring_function.score(temp_protein, ligand)
                    contributions.append((res_key, energy))
                except Exception as e:
                    self.logger.debug(f"Failed to calculate energy for {res_key}: {e}")
                    contributions.append((res_key, 0.0))
            
            # Sort by energy (largest contribution first)
            contributions.sort(key=lambda x: abs(x[1]), reverse=True)
            
            self.logger.info(f"Calculated contributions for {len(contributions)} residues")
            return contributions
            
        except Exception as e:
            self.logger.error(f"Error calculating residue contributions: {e}")
            return []
    
    def _get_protein_atoms(self, protein: Any) -> List[Any]:
        """Extract protein atoms from protein object."""
        if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
            return protein.active_site['atoms']
        elif hasattr(protein, 'atoms'):
            return protein.atoms
        else:
            return []
    
    def _get_ligand_atoms(self, ligand: Any) -> List[Any]:
        """Extract ligand atoms from ligand object."""
        if hasattr(ligand, 'atoms'):
            return ligand.atoms
        elif hasattr(ligand, 'molecule') and hasattr(ligand.molecule, 'atoms'):
            return ligand.molecule.atoms
        elif isinstance(ligand, dict) and 'atoms' in ligand:
            return ligand['atoms']
        else:
            return []
    
    def _group_atoms_by_residue(self, atoms: List[Any]) -> Dict[str, List[Any]]:
        """Group atoms by residue."""
        residue_atoms = {}
        
        for atom in atoms:
            # Extract residue information
            if isinstance(atom, dict):
                res_name = atom.get('residue_name', 'UNK')
                chain_id = atom.get('chain_id', 'A')
                res_id = atom.get('residue_id', 0)
            else:
                res_name = getattr(atom, 'residue_name', 'UNK')
                chain_id = getattr(atom, 'chain_id', 'A')
                res_id = getattr(atom, 'residue_id', 0)
            
            res_key = f"{res_name} {chain_id}:{res_id}"
            
            if res_key not in residue_atoms:
                residue_atoms[res_key] = []
            
            residue_atoms[res_key].append(atom)
        
        return residue_atoms
    
    def _is_residue_near_ligand(self, residue_atoms: List[Any], 
                               ligand_atoms: List[Any], radius: float) -> bool:
        """Check if residue has any atom near the ligand."""
        for atom in residue_atoms:
            atom_pos = self._get_atom_coords(atom)
            if atom_pos is None:
                continue
                
            for l_atom in ligand_atoms:
                l_pos = self._get_atom_coords(l_atom)
                if l_pos is None:
                    continue
                    
                if np.linalg.norm(atom_pos - l_pos) <= radius:
                    return True
        
        return False
    
    def _get_atom_coords(self, atom: Any) -> Optional[np.ndarray]:
        """Extract coordinates from atom object."""
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
    
    def _create_residue_protein(self, atoms: List[Any]) -> Any:
        """Create temporary protein object with only specified atoms."""
        temp_protein = SimpleNamespace()
        temp_protein.atoms = atoms
        temp_protein.active_site = None
        return temp_protein
    
    def visualize_decomposition(self, energy_components: Dict[str, float]) -> Optional[Any]:
        """
        Create visualization of energy component contributions.
        
        Args:
            energy_components: Energy components from decompose_energy method
            
        Returns:
            Matplotlib figure with visualization (if available)
        """
        if not MATPLOTLIB_AVAILABLE:
            self.logger.warning("matplotlib not available for energy visualization")
            return None
        
        try:
            # Extract components (excluding total)
            components = {k: v for k, v in energy_components.items() if k != 'total'}
            
            if not components:
                self.logger.warning("No energy components to visualize")
                return None
            
            # Create bar chart
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Sort components by absolute value
            sorted_components = sorted(components.items(), key=lambda x: abs(x[1]), reverse=True)
            labels = [item[0] for item in sorted_components]
            values = [item[1] for item in sorted_components]
            
            # Set bar colors based on sign
            colors = ['green' if v < 0 else 'red' for v in values]
            
            # Create the bar chart
            bars = ax.bar(labels, values, color=colors, alpha=0.7)
            
            # Add labels and title
            ax.set_xlabel('Energy Component')
            ax.set_ylabel('Energy (kcal/mol)')
            ax.set_title('Binding Energy Decomposition')
            
            # Add value labels on bars
            for bar in bars:
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width()/2.,
                    height + (0.1 if height > 0 else -0.1),
                    f'{height:.2f}',
                    ha='center',
                    va='bottom' if height > 0 else 'top',
                    fontsize=9
                )
            
            # Add total energy line
            total = energy_components.get('total', sum(values))
            ax.axhline(y=total, color='black', linestyle='-', alpha=0.7)
            ax.text(0, total, f'Total: {total:.2f}', va='bottom', ha='left', fontsize=9)
            
            # Adjust layout
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            
            return fig
            
        except Exception as e:
            self.logger.error(f"Error creating energy visualization: {e}")
            return None
    
    def visualize_residue_contributions(self, contributions: List[Tuple[str, float]], 
                                      top_n: int = 10) -> Optional[Any]:
        """
        Create visualization of per-residue energy contributions.
        
        Args:
            contributions: Residue contributions from residue_contributions method
            top_n: Number of top contributors to show
            
        Returns:
            Matplotlib figure with visualization (if available)
        """
        if not MATPLOTLIB_AVAILABLE:
            self.logger.warning("matplotlib not available for residue visualization")
            return None
        
        try:
            if not contributions:
                self.logger.warning("No residue contributions to visualize")
                return None
            
            # Get top contributors
            top_contributions = contributions[:top_n]
            
            # Create horizontal bar chart
            fig, ax = plt.subplots(figsize=(10, 8))
            
            residues = [item[0] for item in top_contributions]
            energies = [item[1] for item in top_contributions]
            
            # Set bar colors based on sign
            colors = ['green' if e < 0 else 'red' for e in energies]
            
            # Create horizontal bar chart
            bars = ax.barh(residues, energies, color=colors, alpha=0.7)
            
            # Add labels and title
            ax.set_xlabel('Energy Contribution (kcal/mol)')
            ax.set_ylabel('Residue')
            ax.set_title(f'Top {len(top_contributions)} Residue Energy Contributions')
            
            # Add value labels on bars
            for bar in bars:
                width = bar.get_width()
                ax.text(
                    width + (0.1 if width > 0 else -0.1),
                    bar.get_y() + bar.get_height()/2.,
                    f'{width:.2f}',
                    ha='left' if width > 0 else 'right',
                    va='center',
                    fontsize=9
                )
            
            # Adjust layout
            plt.tight_layout()
            
            return fig
            
        except Exception as e:
            self.logger.error(f"Error creating residue visualization: {e}")
            return None
    
    def export_analysis(self, energy_components: Dict[str, float], 
                       residue_contributions: List[Tuple[str, float]], 
                       output_file: str) -> None:
        """
        Export energy analysis to file.
        
        Args:
            energy_components: Energy components from decompose_energy
            residue_contributions: Residue contributions from residue_contributions  
            output_file: Output file path
        """
        try:
            import json
            
            export_data = {
                'energy_components': energy_components,
                'residue_contributions': [
                    {'residue': res, 'energy': energy} 
                    for res, energy in residue_contributions
                ],
                'summary': {
                    'total_energy': energy_components.get('total', 0.0),
                    'num_components': len(energy_components) - 1,  # Exclude total
                    'num_contributing_residues': len(residue_contributions)
                }
            }
            
            with open(output_file, 'w') as f:
                json.dump(export_data, f, indent=2)
            
            self.logger.info(f"Energy analysis exported to {output_file}")
            
        except Exception as e:
            self.logger.error(f"Error exporting energy analysis: {e}")