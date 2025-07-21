# -*- coding: utf-8 -*-
"""
Protein preparation module for PandaDock
"""

import logging
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path


class ProteinPreparer:
    """
    Handles protein preparation for docking
    
    Features:
    - PDB file parsing
    - Hydrogen addition
    - Charge assignment
    - Flexible residue handling
    - Binding site detection
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Standard amino acid residues
        self.standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # Flexible residue types
        self.flexible_residue_types = {
            'ARG', 'LYS', 'GLU', 'ASP', 'HIS', 'TYR', 'TRP', 'PHE', 'SER', 'THR', 'ASN', 'GLN'
        }
    
    def prepare_from_file(self, protein_file: str, flexible_residues: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Prepare protein from PDB file
        
        Args:
            protein_file: Path to protein PDB file
            flexible_residues: List of flexible residue identifiers (e.g., ['HIS57', 'SER195'])
            
        Returns:
            Dictionary containing prepared protein data
        """
        self.logger.info(f"Preparing protein from {protein_file}")
        
        protein_path = Path(protein_file)
        if not protein_path.exists():
            raise FileNotFoundError(f"Protein file not found: {protein_file}")
        
        # Parse PDB file
        protein_data = self._parse_pdb_file(protein_file)
        
        # Process flexible residues
        if flexible_residues:
            protein_data['flexible_residues'] = self._process_flexible_residues(
                protein_data, flexible_residues
            )
        else:
            protein_data['flexible_residues'] = []
        
        # Add hydrogens (placeholder - would use actual algorithm)
        protein_data = self._add_hydrogens(protein_data)
        
        # Assign charges (placeholder - would use actual algorithm)
        protein_data = self._assign_charges(protein_data)
        
        # Detect binding sites
        protein_data['binding_sites'] = self._detect_binding_sites(protein_data)
        
        self.logger.info(f"Prepared protein with {len(protein_data['atoms'])} atoms, "
                        f"{len(protein_data['residues'])} residues, "
                        f"{len(protein_data['flexible_residues'])} flexible residues")
        
        return protein_data
    
    def _parse_pdb_file(self, pdb_file: str) -> Dict[str, Any]:
        """Parse PDB file and extract structural information"""
        atoms = []
        residues = {}
        chains = set()
        
        with open(pdb_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        atom_data = self._parse_atom_line(line)
                        atoms.append(atom_data)
                        
                        # Track residues
                        res_id = f"{atom_data['chain']}{atom_data['res_num']}{atom_data['res_name']}"
                        if res_id not in residues:
                            residues[res_id] = {
                                'chain': atom_data['chain'],
                                'res_num': atom_data['res_num'],
                                'res_name': atom_data['res_name'],
                                'atoms': []
                            }
                        residues[res_id]['atoms'].append(len(atoms) - 1)  # Store atom index
                        
                        chains.add(atom_data['chain'])
                        
                    except Exception as e:
                        self.logger.warning(f"Error parsing line {line_num}: {e}")
        
        return {
            'atoms': atoms,
            'residues': residues,
            'chains': list(chains),
            'num_atoms': len(atoms),
            'num_residues': len(residues)
        }
    
    def _parse_atom_line(self, line: str) -> Dict[str, Any]:
        """Parse a single ATOM/HETATM line from PDB"""
        return {
            'record_type': line[0:6].strip(),
            'atom_num': int(line[6:11].strip()),
            'atom_name': line[12:16].strip(),
            'alt_loc': line[16:17].strip(),
            'res_name': line[17:20].strip(),
            'chain': line[21:22].strip(),
            'res_num': int(line[22:26].strip()),
            'insertion': line[26:27].strip(),
            'x': float(line[30:38].strip()),
            'y': float(line[38:46].strip()),
            'z': float(line[46:54].strip()),
            'occupancy': float(line[54:60].strip()) if line[54:60].strip() else 1.0,
            'b_factor': float(line[60:66].strip()) if line[60:66].strip() else 0.0,
            'element': line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0],
            'charge': line[78:80].strip() if len(line) > 78 else ''
        }
    
    def _process_flexible_residues(self, protein_data: Dict[str, Any], 
                                 flexible_residues: List[str]) -> List[Dict[str, Any]]:
        """Process flexible residue specifications"""
        flexible_res_data = []
        
        for flex_res in flexible_residues:
            # Parse residue identifier (e.g., 'HIS57', 'A:HIS57', 'SER195')
            if ':' in flex_res:
                chain, res_spec = flex_res.split(':')
            else:
                chain = None
                res_spec = flex_res
            
            # Extract residue name and number
            if res_spec[-1].isdigit():
                # Format: HIS57
                res_name = ''.join([c for c in res_spec if not c.isdigit()])
                res_num_str = ''.join([c for c in res_spec if c.isdigit()])
                res_num = int(res_num_str) if res_num_str else None
            else:
                # Format might be different
                res_name = res_spec
                res_num = None
            
            # Find matching residues
            for res_id, res_data in protein_data['residues'].items():
                if res_data['res_name'] == res_name:
                    if res_num is None or res_data['res_num'] == res_num:
                        if chain is None or res_data['chain'] == chain:
                            flexible_res_data.append({
                                'residue_id': res_id,
                                'chain': res_data['chain'],
                                'res_num': res_data['res_num'],
                                'res_name': res_data['res_name'],
                                'atom_indices': res_data['atoms']
                            })
                            self.logger.info(f"Added flexible residue: {res_id}")
                            break
            else:
                self.logger.warning(f"Flexible residue not found: {flex_res}")
        
        return flexible_res_data
    
    def _add_hydrogens(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:
        """Add hydrogen atoms to protein structure"""
        # Placeholder implementation
        self.logger.info("Adding hydrogens (placeholder implementation)")
        protein_data['hydrogens_added'] = True
        protein_data['ph'] = 7.0  # Default pH
        return protein_data
    
    def _assign_charges(self, protein_data: Dict[str, Any]) -> Dict[str, Any]:
        """Assign partial charges to atoms"""
        # Placeholder implementation
        self.logger.info("Assigning charges (placeholder implementation)")
        
        # Add default charges to atoms
        for atom in protein_data['atoms']:
            element = atom['element']
            if element == 'N':
                atom['partial_charge'] = -0.3
            elif element == 'O':
                atom['partial_charge'] = -0.4
            elif element == 'C':
                atom['partial_charge'] = 0.1
            else:
                atom['partial_charge'] = 0.0
        
        protein_data['charges_assigned'] = True
        return protein_data
    
    def _detect_binding_sites(self, protein_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Detect potential binding sites in the protein"""
        self.logger.info("Detecting binding sites (placeholder implementation)")
        
        # For now, create a default binding site at the center of the protein
        if protein_data['atoms']:
            coordinates = np.array([[atom['x'], atom['y'], atom['z']] for atom in protein_data['atoms']])
            center = np.mean(coordinates, axis=0)
            
            binding_site = {
                'site_id': 'site_1',
                'center': center.tolist(),
                'volume': 1000.0,
                'surface_area': 500.0,
                'druggability_score': 0.8,
                'residues': []
            }
            
            return [binding_site]
        
        return []
