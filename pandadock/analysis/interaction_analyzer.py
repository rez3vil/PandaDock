"""
Protein-ligand interaction analysis for molecular docking results.

This module provides methods for generating interaction fingerprints and analyzing
protein-ligand binding patterns to understand docking results.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)


class InteractionFingerprinter:
    """Generate and analyze protein-ligand interaction fingerprints."""
    
    def __init__(self, interaction_types: Optional[List[str]] = None, 
                 distance_cutoff: float = 4.5, angle_cutoff: float = 60,
                 include_water: bool = False):
        """
        Initialize interaction fingerprinter.
        
        Args:
            interaction_types: Types of interactions to include in fingerprint
            distance_cutoff: Maximum distance for interactions (Angstroms)
            angle_cutoff: Angle cutoff for directional interactions (degrees)
            include_water: Whether to include water-mediated interactions
        """
        self.interaction_types = interaction_types or [
            'hbond', 'hydrophobic', 'ionic', 'aromatic', 'halogen'
        ]
        self.distance_cutoff = distance_cutoff
        self.angle_cutoff = angle_cutoff
        self.include_water = include_water
        self.logger = logging.getLogger(__name__)
        
        # Define atom types for different interactions
        self.hbond_donors = {'N', 'O', 'S'}
        self.hbond_acceptors = {'N', 'O', 'F', 'S'}
        self.hydrophobic_atoms = {'C'}
        self.halogen_atoms = {'F', 'Cl', 'Br', 'I'}
        self.charged_positive = {'N'}  # Simplified: N in Lys, Arg
        self.charged_negative = {'O'}  # Simplified: O in Asp, Glu
        self.aromatic_atoms = {'C'}    # Simplified: aromatic carbons
    
    def generate_fingerprint(self, protein: Any, ligand: Any) -> Dict[str, Any]:
        """
        Generate interaction fingerprint for a protein-ligand complex.
        
        Args:
            protein: Protein object
            ligand: Ligand pose
            
        Returns:
            Interaction fingerprint with counts of each interaction type
        """
        try:
            # Extract appropriate protein atoms
            if hasattr(protein, 'active_site') and protein.active_site and 'atoms' in protein.active_site:
                protein_atoms = protein.active_site['atoms']
            elif hasattr(protein, 'atoms'):
                protein_atoms = protein.atoms
            else:
                self.logger.warning("No atoms found in protein object")
                return self._empty_fingerprint()
            
            # Initialize fingerprint
            fingerprint = {
                'hbond': 0,
                'hydrophobic': 0,
                'ionic': 0,
                'aromatic': 0,
                'halogen': 0,
                'interactions': []  # Detailed interaction data
            }
            
            # Extract ligand atoms
            ligand_atoms = self._extract_ligand_atoms(ligand)
            if not ligand_atoms:
                self.logger.warning("No atoms found in ligand object")
                return fingerprint
            
            self.logger.debug(f"Analyzing interactions between {len(ligand_atoms)} ligand atoms and {len(protein_atoms)} protein atoms")
            
            # Check each ligand-protein atom pair for interactions
            for l_atom in ligand_atoms:
                l_symbol = self._get_atom_symbol(l_atom)
                l_coords = self._get_atom_coords(l_atom)
                
                if l_coords is None:
                    continue
                
                for p_atom in protein_atoms:
                    p_symbol = self._get_protein_atom_symbol(p_atom)
                    p_coords = self._get_atom_coords(p_atom)
                    
                    if p_coords is None:
                        continue
                    
                    # Calculate distance
                    distance = np.linalg.norm(l_coords - p_coords)
                    
                    # Skip if too far
                    if distance > self.distance_cutoff:
                        continue
                    
                    # Check for different types of interactions
                    self._check_hydrogen_bonds(fingerprint, l_atom, p_atom, l_symbol, p_symbol, distance)
                    self._check_hydrophobic_interactions(fingerprint, l_atom, p_atom, l_symbol, p_symbol, distance)
                    self._check_ionic_interactions(fingerprint, l_atom, p_atom, l_symbol, p_symbol, distance)
                    self._check_halogen_bonds(fingerprint, l_atom, p_atom, l_symbol, p_symbol, distance)
            
            # Check for aromatic interactions (requires special handling)
            if 'aromatic' in self.interaction_types:
                self._check_aromatic_interactions(fingerprint, ligand_atoms, protein_atoms)
            
            self.logger.info(f"Generated fingerprint: {sum(fingerprint[k] for k in fingerprint if k != 'interactions')} total interactions")
            
            return fingerprint
            
        except Exception as e:
            self.logger.error(f"Error generating fingerprint: {e}")
            return self._empty_fingerprint()
    
    def _check_hydrogen_bonds(self, fingerprint: Dict, l_atom: Any, p_atom: Any, 
                            l_symbol: str, p_symbol: str, distance: float) -> None:
        """Check for hydrogen bond interactions."""
        if 'hbond' not in self.interaction_types:
            return
            
        # Check donor-acceptor pairs
        if ((l_symbol in self.hbond_donors and p_symbol in self.hbond_acceptors) or
            (l_symbol in self.hbond_acceptors and p_symbol in self.hbond_donors)):
            # In a full implementation, we would check angles here
            if distance < 3.5:  # Typical H-bond distance
                fingerprint['hbond'] += 1
                fingerprint['interactions'].append({
                    'type': 'hbond',
                    'ligand_atom': l_atom,
                    'protein_atom': p_atom,
                    'distance': distance
                })
    
    def _check_hydrophobic_interactions(self, fingerprint: Dict, l_atom: Any, p_atom: Any,
                                      l_symbol: str, p_symbol: str, distance: float) -> None:
        """Check for hydrophobic interactions."""
        if 'hydrophobic' not in self.interaction_types:
            return
            
        if l_symbol in self.hydrophobic_atoms and p_symbol in self.hydrophobic_atoms:
            if distance < 4.0:  # Typical hydrophobic interaction distance
                fingerprint['hydrophobic'] += 1
                fingerprint['interactions'].append({
                    'type': 'hydrophobic',
                    'ligand_atom': l_atom,
                    'protein_atom': p_atom,
                    'distance': distance
                })
    
    def _check_ionic_interactions(self, fingerprint: Dict, l_atom: Any, p_atom: Any,
                                l_symbol: str, p_symbol: str, distance: float) -> None:
        """Check for ionic interactions."""
        if 'ionic' not in self.interaction_types:
            return
            
        if ((l_symbol in self.charged_positive and p_symbol in self.charged_negative) or
            (l_symbol in self.charged_negative and p_symbol in self.charged_positive)):
            if distance < 4.0:  # Typical ionic interaction distance
                fingerprint['ionic'] += 1
                fingerprint['interactions'].append({
                    'type': 'ionic',
                    'ligand_atom': l_atom,
                    'protein_atom': p_atom,
                    'distance': distance
                })
    
    def _check_halogen_bonds(self, fingerprint: Dict, l_atom: Any, p_atom: Any,
                           l_symbol: str, p_symbol: str, distance: float) -> None:
        """Check for halogen bond interactions."""
        if 'halogen' not in self.interaction_types:
            return
            
        if l_symbol in self.halogen_atoms and p_symbol in self.hbond_acceptors:
            if distance < 3.5:  # Typical halogen bond distance
                fingerprint['halogen'] += 1
                fingerprint['interactions'].append({
                    'type': 'halogen',
                    'ligand_atom': l_atom,
                    'protein_atom': p_atom,
                    'distance': distance
                })
    
    def _check_aromatic_interactions(self, fingerprint: Dict, ligand_atoms: List, protein_atoms: List) -> None:
        """Check for aromatic interactions (simplified implementation)."""
        # In a full implementation, we would identify aromatic rings
        # and check for π-π stacking and other aromatic interactions
        # For now, this is a placeholder
        pass
    
    def compare_fingerprints(self, fp1: Dict[str, Any], fp2: Dict[str, Any]) -> float:
        """
        Compare two interaction fingerprints and return similarity score.
        
        Args:
            fp1: First fingerprint
            fp2: Second fingerprint
            
        Returns:
            Similarity score (0.0 to 1.0)
        """
        # Calculate Tanimoto coefficient for interaction counts
        intersection = 0
        union = 0
        
        for interaction_type in self.interaction_types:
            count1 = fp1.get(interaction_type, 0)
            count2 = fp2.get(interaction_type, 0)
            
            intersection += min(count1, count2)
            union += max(count1, count2)
        
        if union == 0:
            return 0.0
        
        return intersection / union
    
    def analyze_key_interactions(self, protein: Any, ligand: Any) -> List[str]:
        """
        Identify key interactions in a protein-ligand complex.
        
        Args:
            protein: Protein object
            ligand: Ligand pose
            
        Returns:
            List of key interaction descriptions
        """
        # Generate fingerprint
        fingerprint = self.generate_fingerprint(protein, ligand)
        
        # Extract protein residue information
        residue_interactions = {}
        
        for interaction in fingerprint['interactions']:
            p_atom = interaction['protein_atom']
            
            # Get residue info
            res_info = self._get_residue_info(p_atom)
            res_key = f"{res_info['name']} {res_info['chain']}:{res_info['id']}"
            
            if res_key not in residue_interactions:
                residue_interactions[res_key] = []
            
            residue_interactions[res_key].append(interaction)
        
        # Create descriptive list of key interactions
        key_interactions = []
        
        for res_key, interactions in residue_interactions.items():
            # Count interaction types
            interaction_counts = {}
            for interaction in interactions:
                int_type = interaction['type']
                if int_type not in interaction_counts:
                    interaction_counts[int_type] = 0
                interaction_counts[int_type] += 1
            
            # Generate description
            description = f"{res_key}: "
            interaction_strs = []
            
            for int_type, count in interaction_counts.items():
                if count == 1:
                    interaction_strs.append(f"1 {int_type}")
                else:
                    interaction_strs.append(f"{count} {int_type}s")
            
            description += ", ".join(interaction_strs)
            key_interactions.append(description)
        
        # Sort by residue ID
        key_interactions.sort()
        
        return key_interactions
    
    def _extract_ligand_atoms(self, ligand: Any) -> List[Any]:
        """Extract atoms from ligand object."""
        if hasattr(ligand, 'atoms'):
            return ligand.atoms
        elif hasattr(ligand, 'molecule') and hasattr(ligand.molecule, 'atoms'):
            return ligand.molecule.atoms
        elif isinstance(ligand, dict) and 'atoms' in ligand:
            return ligand['atoms']
        else:
            return []
    
    def _get_atom_symbol(self, atom: Any) -> str:
        """Extract atom symbol from atom object."""
        if isinstance(atom, dict):
            return atom.get('symbol', atom.get('element', 'C'))
        elif hasattr(atom, 'symbol'):
            return atom.symbol
        elif hasattr(atom, 'element'):
            return atom.element
        else:
            return 'C'  # Default fallback
    
    def _get_protein_atom_symbol(self, atom: Any) -> str:
        """Extract atom symbol from protein atom object."""
        if isinstance(atom, dict):
            element = atom.get('element', atom.get('name', 'C'))
            return element[0] if isinstance(element, str) else str(element)
        elif hasattr(atom, 'element'):
            return str(atom.element)[0]
        elif hasattr(atom, 'name'):
            return str(atom.name)[0]
        else:
            return 'C'
    
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
    
    def _get_residue_info(self, atom: Any) -> Dict[str, Any]:
        """Extract residue information from protein atom."""
        if isinstance(atom, dict):
            return {
                'name': atom.get('residue_name', 'UNK'),
                'chain': atom.get('chain_id', 'A'),
                'id': atom.get('residue_id', 0)
            }
        else:
            return {
                'name': getattr(atom, 'residue_name', 'UNK'),
                'chain': getattr(atom, 'chain_id', 'A'),
                'id': getattr(atom, 'residue_id', 0)
            }
    
    def _empty_fingerprint(self) -> Dict[str, Any]:
        """Return empty fingerprint structure."""
        return {
            'hbond': 0,
            'hydrophobic': 0,
            'ionic': 0,
            'aromatic': 0,
            'halogen': 0,
            'interactions': []
        }