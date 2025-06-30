"""
Flexible residue detection for PandaDock.

This module provides functionality to automatically detect residues that should
be made flexible during docking based on structural and energetic criteria.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple
import numpy as np


class FlexibleResidueDetector:
    """
    Detects protein residues that should be made flexible during docking.
    
    Uses various criteria including:
    - Proximity to binding site
    - Side chain flexibility
    - B-factor analysis
    - Secondary structure considerations
    """
    
    def __init__(self, distance_cutoff: float = 5.0):
        """
        Initialize the flexible residue detector.
        
        Args:
            distance_cutoff: Maximum distance from binding site for residue selection
        """
        self.distance_cutoff = distance_cutoff
        self.logger = logging.getLogger(__name__)
    
    def detect_flexible_residues(self, protein: Any, max_residues: int = 5,
                               binding_site_center: Optional[Tuple[float, float, float]] = None) -> List[str]:
        """
        Detect residues that should be made flexible during docking.
        
        Args:
            protein: Protein structure object
            max_residues: Maximum number of residues to select
            binding_site_center: Center of binding site (optional)
            
        Returns:
            List of residue identifiers for flexible residues
        """
        try:
            self.logger.info(f"Detecting flexible residues (max: {max_residues})")
            
            # Get binding site center if not provided
            if binding_site_center is None:
                binding_site_center = self._get_binding_site_center(protein)
            
            # Get candidate residues near binding site
            candidate_residues = self._get_nearby_residues(protein, binding_site_center)
            
            # Score residues based on flexibility criteria
            scored_residues = self._score_residues(protein, candidate_residues)
            
            # Select top residues
            flexible_residues = self._select_top_residues(scored_residues, max_residues)
            
            self.logger.info(f"Selected {len(flexible_residues)} flexible residues")
            
            return flexible_residues
            
        except Exception as e:
            self.logger.error(f"Flexible residue detection failed: {e}")
            return []
    
    def _get_binding_site_center(self, protein: Any) -> Tuple[float, float, float]:
        """Get the center of the binding site."""
        try:
            # Try to get active site center if defined
            if hasattr(protein, 'active_site_center'):
                return protein.active_site_center
            
            # Otherwise use geometric center of protein
            coords = protein.get_coordinates()
            center = np.mean(coords, axis=0)
            
            self.logger.info(f"Using protein center as binding site: {center}")
            return tuple(center)
            
        except Exception as e:
            self.logger.warning(f"Could not determine binding site center: {e}")
            # Return origin as fallback
            return (0.0, 0.0, 0.0)
    
    def _get_nearby_residues(self, protein: Any, center: Tuple[float, float, float]) -> List[Dict[str, Any]]:
        """Get residues within distance cutoff of binding site."""
        try:
            nearby_residues = []
            
            # Get all residues from protein
            residues = protein.get_residues() if hasattr(protein, 'get_residues') else []
            
            for residue in residues:
                # Calculate distance to binding site center
                residue_center = self._get_residue_center(residue)
                distance = np.linalg.norm(np.array(residue_center) - np.array(center))
                
                if distance <= self.distance_cutoff:
                    nearby_residues.append({
                        'residue': residue,
                        'distance': distance,
                        'residue_id': self._get_residue_id(residue)
                    })
            
            self.logger.info(f"Found {len(nearby_residues)} residues within {self.distance_cutoff}Å")
            
            return nearby_residues
            
        except Exception as e:
            self.logger.error(f"Error finding nearby residues: {e}")
            return []
    
    def _get_residue_center(self, residue: Any) -> Tuple[float, float, float]:
        """Get the geometric center of a residue."""
        try:
            if hasattr(residue, 'get_center'):
                return residue.get_center()
            elif hasattr(residue, 'get_coordinates'):
                coords = residue.get_coordinates()
                center = np.mean(coords, axis=0)
                return tuple(center)
            else:
                # Fallback - try to get CA atom position
                if hasattr(residue, 'get_atom'):
                    ca_atom = residue.get_atom('CA')
                    if ca_atom:
                        return ca_atom.get_coordinates()
                
                # Final fallback
                return (0.0, 0.0, 0.0)
                
        except Exception:
            return (0.0, 0.0, 0.0)
    
    def _get_residue_id(self, residue: Any) -> str:
        """Get a string identifier for the residue."""
        try:
            if hasattr(residue, 'get_id'):
                return str(residue.get_id())
            elif hasattr(residue, 'id'):
                return str(residue.id)
            elif hasattr(residue, 'resnum') and hasattr(residue, 'resname'):
                return f"{residue.resname}{residue.resnum}"
            else:
                return f"UNK{id(residue)}"
                
        except Exception:
            return f"UNK{id(residue)}"
    
    def _score_residues(self, protein: Any, candidate_residues: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Score residues based on flexibility criteria."""
        scored_residues = []
        
        for residue_info in candidate_residues:
            residue = residue_info['residue']
            
            # Calculate flexibility score
            score = 0.0
            
            # Distance score (closer is better)
            distance_score = max(0, 1.0 - residue_info['distance'] / self.distance_cutoff)
            score += distance_score * 0.4
            
            # Side chain flexibility score
            flexibility_score = self._calculate_side_chain_flexibility(residue)
            score += flexibility_score * 0.3
            
            # B-factor score (if available)
            bfactor_score = self._calculate_bfactor_score(residue)
            score += bfactor_score * 0.2
            
            # Secondary structure score
            ss_score = self._calculate_secondary_structure_score(residue)
            score += ss_score * 0.1
            
            scored_residues.append({
                'residue_id': residue_info['residue_id'],
                'residue': residue,
                'score': score,
                'distance': residue_info['distance']
            })
        
        # Sort by score (highest first)
        scored_residues.sort(key=lambda x: x['score'], reverse=True)
        
        return scored_residues
    
    def _calculate_side_chain_flexibility(self, residue: Any) -> float:
        """Calculate side chain flexibility based on residue type."""
        # Flexibility scores for different residue types
        flexibility_map = {
            'GLY': 0.9,  # Very flexible
            'PRO': 0.1,  # Very rigid
            'ALA': 0.3,  # Small, limited flexibility
            'VAL': 0.4, 'LEU': 0.5, 'ILE': 0.4,  # Branched aliphatic
            'MET': 0.7, 'CYS': 0.6,  # Sulfur-containing
            'PHE': 0.6, 'TYR': 0.7, 'TRP': 0.8,  # Aromatic
            'SER': 0.6, 'THR': 0.5,  # Small polar
            'ASN': 0.7, 'GLN': 0.8,  # Amide
            'ASP': 0.8, 'GLU': 0.9,  # Acidic
            'HIS': 0.7, 'ARG': 0.9, 'LYS': 0.8,  # Basic
        }
        
        try:
            # Get residue name
            resname = None
            if hasattr(residue, 'get_resname'):
                resname = residue.get_resname()
            elif hasattr(residue, 'resname'):
                resname = residue.resname
            
            if resname and resname in flexibility_map:
                return flexibility_map[resname]
            else:
                return 0.5  # Default moderate flexibility
                
        except Exception:
            return 0.5
    
    def _calculate_bfactor_score(self, residue: Any) -> float:
        """Calculate flexibility score based on B-factors."""
        try:
            # Try to get B-factors from residue atoms
            bfactors = []
            
            if hasattr(residue, 'get_atoms'):
                atoms = residue.get_atoms()
                for atom in atoms:
                    if hasattr(atom, 'get_bfactor'):
                        bfactors.append(atom.get_bfactor())
                    elif hasattr(atom, 'bfactor'):
                        bfactors.append(atom.bfactor)
            
            if bfactors:
                avg_bfactor = np.mean(bfactors)
                # Normalize B-factor (typical range 10-100)
                normalized_score = min(1.0, max(0.0, (avg_bfactor - 20) / 60))
                return normalized_score
            else:
                return 0.5  # Default if no B-factors available
                
        except Exception:
            return 0.5
    
    def _calculate_secondary_structure_score(self, residue: Any) -> float:
        """Calculate flexibility score based on secondary structure."""
        try:
            # Get secondary structure if available
            ss = None
            if hasattr(residue, 'get_secondary_structure'):
                ss = residue.get_secondary_structure()
            elif hasattr(residue, 'secondary_structure'):
                ss = residue.secondary_structure
            
            # Score based on secondary structure type
            if ss:
                if ss in ['C', 'L', '-']:  # Coil/loop
                    return 1.0
                elif ss in ['T', 'S']:  # Turn/bend
                    return 0.8
                elif ss in ['H', 'G', 'I']:  # Helix
                    return 0.3
                elif ss in ['E', 'B']:  # Sheet
                    return 0.2
            
            return 0.6  # Default for unknown
            
        except Exception:
            return 0.6
    
    def _select_top_residues(self, scored_residues: List[Dict[str, Any]], 
                           max_residues: int) -> List[str]:
        """Select the top-scoring residues."""
        selected = []
        
        for residue_info in scored_residues[:max_residues]:
            selected.append(residue_info['residue_id'])
            self.logger.debug(f"Selected flexible residue: {residue_info['residue_id']} "
                            f"(score: {residue_info['score']:.3f}, "
                            f"distance: {residue_info['distance']:.2f}Å)")
        
        return selected