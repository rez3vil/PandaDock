#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D Interaction Map Generator for PandaDock
==========================================

Creates Discovery Studio-style 2D interaction diagrams for protein-ligand complexes.
Analyzes hydrogen bonds, hydrophobic contacts, π-π stacking, and other interactions.

Author: PandaDock Team
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Wedge
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional, Any, Union
import pandas as pd
from dataclasses import dataclass
import warnings
from scipy.spatial.distance import cdist

warnings.filterwarnings('ignore')

@dataclass
class Atom:
    """Represents an atom with coordinates and properties"""
    name: str
    element: str
    residue: str
    residue_id: int
    chain: str
    x: float
    y: float
    z: float
    is_protein: bool = True

@dataclass
class Interaction:
    """Represents a protein-ligand interaction"""
    interaction_type: str
    protein_atom: Atom
    ligand_atom: Atom
    distance: float
    angle: float = None
    strength: str = "medium"  # weak, medium, strong

class InteractionAnalyzer:
    """
    Analyzes protein-ligand interactions and creates 2D interaction maps
    """
    
    def __init__(self, output_dir: str):
        """Initialize the interaction analyzer"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        
        # Interaction parameters
        self.interaction_params = {
            'hydrogen_bond': {
                'max_distance': 3.5,
                'min_angle': 120.0,
                'donor_elements': ['N', 'O', 'S'],
                'acceptor_elements': ['N', 'O', 'S', 'F']
            },
            'hydrophobic': {
                'max_distance': 4.0,
                'hydrophobic_atoms': ['C']
            },
            'pi_stacking': {
                'max_distance': 4.5,
                'aromatic_residues': ['PHE', 'TYR', 'TRP', 'HIS']
            },
            'electrostatic': {
                'max_distance': 6.0,
                'positive_residues': ['ARG', 'LYS', 'HIS'],
                'negative_residues': ['ASP', 'GLU']
            },
            'van_der_waals': {
                'min_distance': 2.5,
                'max_distance': 4.0
            }
        }
        
        # Color scheme for interactions
        self.interaction_colors = {
            'hydrogen_bond': '#2E86AB',
            'hydrophobic': '#F18F01',
            'pi_stacking': '#A23B72',
            'electrostatic': '#C73E1D',
            'van_der_waals': '#44AF69',
            'metal_coordination': '#8E44AD'
        }
        
        # Residue properties
        self.residue_properties = {
            'ALA': {'type': 'hydrophobic', 'color': '#FFA500'},
            'VAL': {'type': 'hydrophobic', 'color': '#FFA500'},
            'LEU': {'type': 'hydrophobic', 'color': '#FFA500'},
            'ILE': {'type': 'hydrophobic', 'color': '#FFA500'},
            'PHE': {'type': 'aromatic', 'color': '#FF69B4'},
            'TRP': {'type': 'aromatic', 'color': '#FF69B4'},
            'TYR': {'type': 'aromatic', 'color': '#FF69B4'},
            'PRO': {'type': 'hydrophobic', 'color': '#FFA500'},
            'SER': {'type': 'polar', 'color': '#00CED1'},
            'THR': {'type': 'polar', 'color': '#00CED1'},
            'CYS': {'type': 'polar', 'color': '#00CED1'},
            'MET': {'type': 'hydrophobic', 'color': '#FFA500'},
            'ASN': {'type': 'polar', 'color': '#00CED1'},
            'GLN': {'type': 'polar', 'color': '#00CED1'},
            'ASP': {'type': 'negative', 'color': '#FF4500'},
            'GLU': {'type': 'negative', 'color': '#FF4500'},
            'LYS': {'type': 'positive', 'color': '#0000FF'},
            'ARG': {'type': 'positive', 'color': '#0000FF'},
            'HIS': {'type': 'positive', 'color': '#0000FF'},
            'GLY': {'type': 'neutral', 'color': '#808080'}
        }
    
    def parse_pdb_file(self, pdb_file: str) -> Tuple[List[Atom], List[Atom]]:
        """
        Parse PDB file and extract protein and ligand atoms
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            Tuple of (protein_atoms, ligand_atoms)
        """
        protein_atoms = []
        ligand_atoms = []
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Parse PDB line
                        atom_name = line[12:16].strip()
                        residue = line[17:20].strip()
                        chain = line[21].strip()
                        residue_id = int(line[22:26].strip())
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        # Determine element from atom name
                        element = atom_name[0] if atom_name[0].isalpha() else atom_name[1]
                        element = element.upper()
                        
                        # Create atom object
                        atom = Atom(
                            name=atom_name,
                            element=element,
                            residue=residue,
                            residue_id=residue_id,
                            chain=chain,
                            x=x, y=y, z=z
                        )
                        
                        # Classify as protein or ligand
                        if residue in self.residue_properties or line.startswith('ATOM'):
                            atom.is_protein = True
                            protein_atoms.append(atom)
                        else:
                            atom.is_protein = False
                            ligand_atoms.append(atom)
                            
        except Exception as e:
            self.logger.warning(f"Error parsing PDB file {pdb_file}: {e}")
            
        self.logger.info(f"Parsed {len(protein_atoms)} protein atoms and {len(ligand_atoms)} ligand atoms")
        return protein_atoms, ligand_atoms
    
    def analyze_interactions(self, protein_atoms: List[Atom], ligand_atoms: List[Atom]) -> List[Interaction]:
        """
        Analyze all interactions between protein and ligand
        
        Args:
            protein_atoms: List of protein atoms
            ligand_atoms: List of ligand atoms
            
        Returns:
            List of interactions
        """
        interactions = []
        
        # Calculate distance matrix
        protein_coords = np.array([[atom.x, atom.y, atom.z] for atom in protein_atoms])
        ligand_coords = np.array([[atom.x, atom.y, atom.z] for atom in ligand_atoms])
        
        if len(protein_coords) == 0 or len(ligand_coords) == 0:
            return interactions
        
        distances = cdist(protein_coords, ligand_coords)
        
        # Analyze each protein-ligand atom pair
        for i, p_atom in enumerate(protein_atoms):
            for j, l_atom in enumerate(ligand_atoms):
                distance = distances[i, j]
                
                # Skip if too far
                if distance > 6.0:
                    continue
                
                # Check for hydrogen bonds
                if self._is_hydrogen_bond(p_atom, l_atom, distance):
                    interaction = Interaction(
                        interaction_type='hydrogen_bond',
                        protein_atom=p_atom,
                        ligand_atom=l_atom,
                        distance=distance,
                        strength=self._classify_strength('hydrogen_bond', distance)
                    )
                    interactions.append(interaction)
                
                # Check for hydrophobic contacts
                elif self._is_hydrophobic_contact(p_atom, l_atom, distance):
                    interaction = Interaction(
                        interaction_type='hydrophobic',
                        protein_atom=p_atom,
                        ligand_atom=l_atom,
                        distance=distance,
                        strength=self._classify_strength('hydrophobic', distance)
                    )
                    interactions.append(interaction)
                
                # Check for electrostatic interactions
                elif self._is_electrostatic(p_atom, l_atom, distance):
                    interaction = Interaction(
                        interaction_type='electrostatic',
                        protein_atom=p_atom,
                        ligand_atom=l_atom,
                        distance=distance,
                        strength=self._classify_strength('electrostatic', distance)
                    )
                    interactions.append(interaction)
                
                # Check for π-π stacking
                elif self._is_pi_stacking(p_atom, l_atom, distance):
                    interaction = Interaction(
                        interaction_type='pi_stacking',
                        protein_atom=p_atom,
                        ligand_atom=l_atom,
                        distance=distance,
                        strength=self._classify_strength('pi_stacking', distance)
                    )
                    interactions.append(interaction)
                
                # Check for van der Waals
                elif self._is_van_der_waals(p_atom, l_atom, distance):
                    interaction = Interaction(
                        interaction_type='van_der_waals',
                        protein_atom=p_atom,
                        ligand_atom=l_atom,
                        distance=distance,
                        strength=self._classify_strength('van_der_waals', distance)
                    )
                    interactions.append(interaction)
        
        self.logger.info(f"Found {len(interactions)} interactions")
        return interactions
    
    def _is_hydrogen_bond(self, p_atom: Atom, l_atom: Atom, distance: float) -> bool:
        """Check if atoms form a hydrogen bond"""
        params = self.interaction_params['hydrogen_bond']
        
        if distance > params['max_distance']:
            return False
        
        # Check if one is donor and other is acceptor
        donor_elements = params['donor_elements']
        acceptor_elements = params['acceptor_elements']
        
        is_hbond = (
            (p_atom.element in donor_elements and l_atom.element in acceptor_elements) or
            (p_atom.element in acceptor_elements and l_atom.element in donor_elements)
        )
        
        return is_hbond
    
    def _is_hydrophobic_contact(self, p_atom: Atom, l_atom: Atom, distance: float) -> bool:
        """Check if atoms form a hydrophobic contact"""
        params = self.interaction_params['hydrophobic']
        
        if distance > params['max_distance']:
            return False
        
        # Both atoms should be carbon
        is_hydrophobic = (
            p_atom.element == 'C' and l_atom.element == 'C' and
            p_atom.residue in ['ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'PRO', 'MET']
        )
        
        return is_hydrophobic
    
    def _is_electrostatic(self, p_atom: Atom, l_atom: Atom, distance: float) -> bool:
        """Check if atoms form an electrostatic interaction"""
        params = self.interaction_params['electrostatic']
        
        if distance > params['max_distance']:
            return False
        
        # Check for charged residues
        is_electrostatic = (
            p_atom.residue in params['positive_residues'] + params['negative_residues']
        )
        
        return is_electrostatic
    
    def _is_pi_stacking(self, p_atom: Atom, l_atom: Atom, distance: float) -> bool:
        """Check if atoms participate in π-π stacking"""
        params = self.interaction_params['pi_stacking']
        
        if distance > params['max_distance']:
            return False
        
        # Check for aromatic residues and atoms
        is_pi_stacking = (
            p_atom.residue in params['aromatic_residues'] and
            l_atom.element == 'C'  # Simplified check for aromatic carbon
        )
        
        return is_pi_stacking
    
    def _is_van_der_waals(self, p_atom: Atom, l_atom: Atom, distance: float) -> bool:
        """Check if atoms form van der Waals contact"""
        params = self.interaction_params['van_der_waals']
        
        return params['min_distance'] <= distance <= params['max_distance']
    
    def _classify_strength(self, interaction_type: str, distance: float) -> str:
        """Classify interaction strength based on distance"""
        if interaction_type == 'hydrogen_bond':
            if distance < 2.5:
                return 'strong'
            elif distance < 3.0:
                return 'medium'
            else:
                return 'weak'
        elif interaction_type in ['hydrophobic', 'van_der_waals']:
            if distance < 3.5:
                return 'strong'
            elif distance < 4.0:
                return 'medium'
            else:
                return 'weak'
        else:
            return 'medium'
    
    def create_2d_interaction_map(self, interactions: List[Interaction], 
                                 pose_id: str, binding_affinity: float) -> str:
        """
        Create a 2D interaction map in Discovery Studio style
        
        Args:
            interactions: List of protein-ligand interactions
            pose_id: Identifier for the pose
            binding_affinity: Binding affinity value
            
        Returns:
            Path to generated plot file
        """
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        # Group interactions by residue
        residue_interactions = {}
        for interaction in interactions:
            residue_key = f"{interaction.protein_atom.residue}{interaction.protein_atom.residue_id}"
            if residue_key not in residue_interactions:
                residue_interactions[residue_key] = []
            residue_interactions[residue_key].append(interaction)
        
        # Create circular layout for residues around central ligand
        center_x, center_y = 0.5, 0.5
        ligand_radius = 0.08
        residue_radius = 0.3
        
        # Draw central ligand
        ligand_circle = Circle((center_x, center_y), ligand_radius, 
                             facecolor='lightblue', edgecolor='navy', linewidth=2)
        ax.add_patch(ligand_circle)
        ax.text(center_x, center_y, 'LIGAND', ha='center', va='center', 
               fontweight='bold', fontsize=10)
        
        # Position residues around ligand
        num_residues = len(residue_interactions)
        if num_residues > 0:
            angle_step = 2 * np.pi / num_residues
            
            for i, (residue_key, res_interactions) in enumerate(residue_interactions.items()):
                angle = i * angle_step
                res_x = center_x + residue_radius * np.cos(angle)
                res_y = center_y + residue_radius * np.sin(angle)
                
                # Get residue info
                residue_name = res_interactions[0].protein_atom.residue
                residue_id = res_interactions[0].protein_atom.residue_id
                
                # Choose color based on residue type
                res_color = self.residue_properties.get(residue_name, {}).get('color', '#808080')
                
                # Draw residue circle
                residue_circle = Circle((res_x, res_y), 0.04, 
                                      facecolor=res_color, edgecolor='black', linewidth=1)
                ax.add_patch(residue_circle)
                
                # Add residue label
                ax.text(res_x, res_y - 0.08, f"{residue_name}{residue_id}", 
                       ha='center', va='center', fontsize=8, fontweight='bold')
                
                # Draw interaction lines
                for interaction in res_interactions:
                    color = self.interaction_colors.get(interaction.interaction_type, 'gray')
                    
                    # Line style based on strength
                    if interaction.strength == 'strong':
                        linewidth = 3
                        linestyle = '-'
                    elif interaction.strength == 'medium':
                        linewidth = 2
                        linestyle = '-'
                    else:
                        linewidth = 1
                        linestyle = '--'
                    
                    # Draw interaction line
                    ax.plot([center_x, res_x], [center_y, res_y], 
                           color=color, linewidth=linewidth, linestyle=linestyle, alpha=0.8)
                    
                    # Add interaction label
                    mid_x = (center_x + res_x) / 2
                    mid_y = (center_y + res_y) / 2
                    
                    # Create interaction symbol
                    if interaction.interaction_type == 'hydrogen_bond':
                        symbol = 'H'
                    elif interaction.interaction_type == 'hydrophobic':
                        symbol = '○'
                    elif interaction.interaction_type == 'electrostatic':
                        symbol = '±'
                    elif interaction.interaction_type == 'pi_stacking':
                        symbol = 'π'
                    else:
                        symbol = '•'
                    
                    ax.text(mid_x, mid_y, symbol, ha='center', va='center',
                           fontsize=8, fontweight='bold', 
                           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
        
        # Create legend
        legend_elements = []
        for interaction_type, color in self.interaction_colors.items():
            if any(int.interaction_type == interaction_type for int in interactions):
                count = sum(1 for int in interactions if int.interaction_type == interaction_type)
                legend_elements.append(
                    plt.Line2D([0], [0], color=color, linewidth=2, 
                              label=f"{interaction_type.replace('_', ' ').title()} ({count})")
                )
        
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
        
        # Add statistics box
        total_interactions = len(interactions)
        unique_residues = len(residue_interactions)
        
        stats_text = f"""Pose: {pose_id}
Binding Affinity: {binding_affinity:.2f} kcal/mol

Interaction Summary:
• Total Interactions: {total_interactions}
• Unique Residues: {unique_residues}
• H-bonds: {sum(1 for i in interactions if i.interaction_type == 'hydrogen_bond')}
• Hydrophobic: {sum(1 for i in interactions if i.interaction_type == 'hydrophobic')}
• Electrostatic: {sum(1 for i in interactions if i.interaction_type == 'electrostatic')}
• π-π Stacking: {sum(1 for i in interactions if i.interaction_type == 'pi_stacking')}
• van der Waals: {sum(1 for i in interactions if i.interaction_type == 'van_der_waals')}"""
        
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               verticalalignment='top', fontsize=9,
               bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        # Set plot properties
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(f'2D Interaction Map - {pose_id}', fontsize=14, fontweight='bold', pad=20)
        
        # Save plot
        output_file = self.output_dir / f'interaction_map_{pose_id}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Generated 2D interaction map: {output_file}")
        return str(output_file)
    
    def analyze_pose_interactions(self, pose_pdb_file: str, pose_id: str, 
                                binding_affinity: float) -> str:
        """
        Analyze interactions for a single pose and create 2D map
        
        Args:
            pose_pdb_file: Path to pose PDB file
            pose_id: Pose identifier
            binding_affinity: Binding affinity value
            
        Returns:
            Path to generated interaction map
        """
        try:
            # Parse PDB file
            protein_atoms, ligand_atoms = self.parse_pdb_file(pose_pdb_file)
            
            if not protein_atoms or not ligand_atoms:
                self.logger.warning(f"No protein or ligand atoms found in {pose_pdb_file}")
                return self._create_empty_interaction_map(pose_id, binding_affinity)
            
            # Analyze interactions
            interactions = self.analyze_interactions(protein_atoms, ligand_atoms)
            
            # Create 2D interaction map
            return self.create_2d_interaction_map(interactions, pose_id, binding_affinity)
            
        except Exception as e:
            self.logger.error(f"Error analyzing pose {pose_id}: {e}")
            return self._create_empty_interaction_map(pose_id, binding_affinity)
    
    def _create_empty_interaction_map(self, pose_id: str, binding_affinity: float) -> str:
        """Create placeholder interaction map when analysis fails"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        ax.text(0.5, 0.5, f'2D Interaction Map\n\nPose: {pose_id}\n'
                         f'Binding Affinity: {binding_affinity:.2f} kcal/mol\n\n'
                         f'[Unable to analyze interactions]\n'
                         f'Please check PDB file format',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=12, bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(f'Interaction Map - {pose_id}')
        ax.axis('off')
        
        output_file = self.output_dir / f'interaction_map_{pose_id}_placeholder.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(output_file)


def create_interaction_maps_for_poses(poses_df: pd.DataFrame, poses_dir: str, 
                                    output_dir: str, top_n: int = 3) -> List[str]:
    """
    Create interaction maps for top N poses
    
    Args:
        poses_df: DataFrame with pose information
        poses_dir: Directory containing pose PDB files
        output_dir: Output directory for maps
        top_n: Number of top poses to analyze
        
    Returns:
        List of generated interaction map file paths
    """
    analyzer = InteractionAnalyzer(output_dir)
    generated_files = []
    
    top_poses = poses_df.head(top_n)
    
    for _, pose in top_poses.iterrows():
        pose_pdb_file = Path(poses_dir) / f"{pose['Pose_ID']}.pdb"
        
        if pose_pdb_file.exists():
            interaction_map = analyzer.analyze_pose_interactions(
                str(pose_pdb_file), pose['Pose_ID'], pose['Binding_Affinity']
            )
            generated_files.append(interaction_map)
        else:
            # Create placeholder if PDB file doesn't exist
            placeholder_map = analyzer._create_empty_interaction_map(
                pose['Pose_ID'], pose['Binding_Affinity']
            )
            generated_files.append(placeholder_map)
    
    return generated_files


if __name__ == "__main__":
    print("PandaDock Interaction Analyzer")
    print("Creates 2D interaction maps for protein-ligand complexes")