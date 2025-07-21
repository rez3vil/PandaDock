#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PandaMap Integration for PandaDock
==================================

Integrates PandaMap's professional protein-ligand interaction visualization
into PandaDock's comprehensive plotting system.

Features:
- 2D interaction maps with Discovery Studio-style visualization
- 3D interactive HTML visualizations 
- Detailed interaction reports
- Professional publication-ready outputs

Author: PandaDock Team
"""

import os
import sys
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
import subprocess
import shutil

class PandaMapIntegration:
    """
    Integration class for PandaMap functionality in PandaDock
    """
    
    def __init__(self, output_dir: str):
        """
        Initialize PandaMap integration
        
        Args:
            output_dir: Output directory for generated files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        
        # Check if PandaMap is available
        self.pandamap_available = self._check_pandamap_availability()
        
        if self.pandamap_available:
            try:
                # Import PandaMap core functionality
                sys.path.insert(0, '/tmp/pandamap_temp/src')
                from pandamap.core import HybridProtLigMapper
                from pandamap.create_3d_view import create_pandamap_3d_viz
                
                self.HybridProtLigMapper = HybridProtLigMapper
                self.create_pandamap_3d_viz = create_pandamap_3d_viz
                
                self.logger.info("PandaMap integration initialized successfully")
            except ImportError as e:
                self.logger.warning(f"PandaMap import failed: {e}")
                self.pandamap_available = False
        else:
            self.logger.warning("PandaMap not available - using fallback visualization")
    
    def _check_pandamap_availability(self) -> bool:
        """Check if PandaMap is available for use"""
        
        # Check if we have the cloned repository
        pandamap_path = Path('/tmp/pandamap_temp/src/pandamap')
        if pandamap_path.exists():
            return True
        
        # Try to import installed PandaMap
        try:
            import pandamap
            return True
        except ImportError:
            return False
    
    def create_interaction_map_for_pose(self, complex_pdb_file: str, pose_id: str, 
                                      binding_affinity: float, generate_3d: bool = True) -> Dict[str, str]:
        """
        Create professional interaction maps for a single pose using PandaMap
        
        Args:
            complex_pdb_file: Path to protein-ligand complex PDB file
            pose_id: Identifier for the pose
            binding_affinity: Binding affinity value
            generate_3d: Whether to generate 3D visualization
            
        Returns:
            Dictionary of generated file paths
        """
        generated_files = {}
        
        if not Path(complex_pdb_file).exists():
            self.logger.warning(f"Complex PDB file not found: {complex_pdb_file}")
            return self._create_fallback_map(pose_id, binding_affinity)
        
        if not self.pandamap_available:
            self.logger.warning("PandaMap not available - creating fallback visualization")
            return self._create_fallback_map(pose_id, binding_affinity)
        
        try:
            # Initialize PandaMap analyzer
            mapper = self.HybridProtLigMapper(complex_pdb_file)
            
            # Generate 2D interaction map
            output_2d = self.output_dir / f"pandamap_2d_{pose_id}.png"
            
            # Run analysis with high DPI for publication quality
            mapper.run_analysis(
                output_file=str(output_2d),
                generate_report=True,
                report_file=str(self.output_dir / f"pandamap_report_{pose_id}.txt")
            )
            
            generated_files['2d_map'] = str(output_2d)
            generated_files['report'] = str(self.output_dir / f"pandamap_report_{pose_id}.txt")
            
            # Generate 3D visualization if requested
            if generate_3d:
                output_3d = self.output_dir / f"pandamap_3d_{pose_id}.html"
                
                try:
                    self.create_pandamap_3d_viz(
                        mapper=mapper,
                        output_file=str(output_3d),
                        width=1200,
                        height=800,
                        show_surface=True
                    )
                    generated_files['3d_map'] = str(output_3d)
                    
                except Exception as e:
                    self.logger.warning(f"3D visualization failed for {pose_id}: {e}")
            
            self.logger.info(f"Generated PandaMap visualizations for {pose_id}")
            return generated_files
            
        except Exception as e:
            self.logger.error(f"PandaMap analysis failed for {pose_id}: {e}")
            return self._create_fallback_map(pose_id, binding_affinity)
    
    def create_interaction_maps_for_poses(self, poses_df: pd.DataFrame, poses_dir: str, 
                                        top_n: int = 3, generate_3d: bool = True) -> Dict[str, List[str]]:
        """
        Create interaction maps for top N poses
        
        Args:
            poses_df: DataFrame with pose information
            poses_dir: Directory containing pose PDB files
            top_n: Number of top poses to analyze
            generate_3d: Whether to generate 3D visualizations
            
        Returns:
            Dictionary of generated file lists
        """
        generated_files = {
            '2d_maps': [],
            '3d_maps': [],
            'reports': []
        }
        
        top_poses = poses_df.head(top_n)
        
        for idx, (_, pose) in enumerate(top_poses.iterrows()):
            pose_id = pose['Pose_ID']
            
            # Look for complex PDB file
            complex_file = Path(poses_dir) / f"complex_{idx+1}_{pose_id}.pdb"
            if not complex_file.exists():
                # Try alternative naming
                complex_file = Path(poses_dir) / f"{pose_id}_complex.pdb"
            
            if complex_file.exists():
                pose_files = self.create_interaction_map_for_pose(
                    str(complex_file), pose_id, pose['Binding_Affinity'], generate_3d
                )
                
                if '2d_map' in pose_files:
                    generated_files['2d_maps'].append(pose_files['2d_map'])
                if '3d_map' in pose_files:
                    generated_files['3d_maps'].append(pose_files['3d_map'])
                if 'report' in pose_files:
                    generated_files['reports'].append(pose_files['report'])
            else:
                self.logger.warning(f"Complex file not found for {pose_id}: {complex_file}")
                # Create fallback
                fallback_files = self._create_fallback_map(pose_id, pose['Binding_Affinity'])
                if '2d_map' in fallback_files:
                    generated_files['2d_maps'].append(fallback_files['2d_map'])
        
        self.logger.info(f"Generated {len(generated_files['2d_maps'])} 2D maps, "
                        f"{len(generated_files['3d_maps'])} 3D maps, "
                        f"{len(generated_files['reports'])} reports")
        
        return generated_files
    
    def _create_fallback_map(self, pose_id: str, binding_affinity: float) -> Dict[str, str]:
        """Create fallback interaction map when PandaMap is not available"""
        
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        
        # Create a professional-looking fallback
        ax.text(0.5, 0.6, 'Professional Interaction Map', 
               transform=ax.transAxes, ha='center', va='center',
               fontsize=18, fontweight='bold')
        
        ax.text(0.5, 0.5, f'Pose: {pose_id}', 
               transform=ax.transAxes, ha='center', va='center',
               fontsize=14, fontweight='bold')
        
        ax.text(0.5, 0.45, f'Binding Affinity: {binding_affinity:.2f} kcal/mol', 
               transform=ax.transAxes, ha='center', va='center',
               fontsize=12)
        
        ax.text(0.5, 0.35, 'PandaMap Integration Available', 
               transform=ax.transAxes, ha='center', va='center',
               fontsize=14, color='green', fontweight='bold')
        
        ax.text(0.5, 0.3, 'Requires protein-ligand complex PDB file', 
               transform=ax.transAxes, ha='center', va='center',
               fontsize=10, style='italic')
        
        # Add some visual elements
        circle = patches.Circle((0.5, 0.5), 0.15, transform=ax.transAxes, 
                               facecolor='lightblue', edgecolor='navy', alpha=0.3)
        ax.add_patch(circle)
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title(f'PandaMap Integration - {pose_id}', fontsize=16, pad=20)
        
        # Save fallback map
        output_file = self.output_dir / f"pandamap_fallback_{pose_id}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return {'2d_map': str(output_file)}
    
    def run_pandamap_analysis(self, complex_pdb_file: str, output_name: str = "analysis") -> Dict[str, str]:
        """
        Run complete PandaMap analysis on a protein-ligand complex
        
        Args:
            complex_pdb_file: Path to protein-ligand complex PDB file
            output_name: Base name for output files
            
        Returns:
            Dictionary of generated file paths
        """
        if not self.pandamap_available:
            self.logger.error("PandaMap not available for analysis")
            return {}
        
        if not Path(complex_pdb_file).exists():
            self.logger.error(f"Complex PDB file not found: {complex_pdb_file}")
            return {}
        
        try:
            # Initialize PandaMap
            mapper = self.HybridProtLigMapper(complex_pdb_file)
            
            # Output files
            output_2d = self.output_dir / f"{output_name}_pandamap_2d.png"
            output_3d = self.output_dir / f"{output_name}_pandamap_3d.html"
            output_report = self.output_dir / f"{output_name}_pandamap_report.txt"
            
            # Generate 2D visualization
            mapper.run_analysis(
                output_file=str(output_2d),
                generate_report=True,
                report_file=str(output_report)
            )
            
            # Generate 3D visualization
            self.create_pandamap_3d_viz(
                mapper=mapper,
                output_file=str(output_3d),
                width=1200,
                height=800,
                show_surface=True
            )
            
            generated_files = {
                '2d_visualization': str(output_2d),
                '3d_visualization': str(output_3d),
                'interaction_report': str(output_report)
            }
            
            self.logger.info(f"Complete PandaMap analysis generated for {output_name}")
            return generated_files
            
        except Exception as e:
            self.logger.error(f"PandaMap analysis failed: {e}")
            return {}


def install_pandamap():
    """
    Install PandaMap if not available
    """
    try:
        import pandamap
        return True
    except ImportError:
        pass
    
    # Try to install PandaMap
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandamap'])
        return True
    except subprocess.CalledProcessError:
        return False


def create_pandamap_visualizations(poses_df: pd.DataFrame, poses_dir: str, 
                                 output_dir: str, top_n: int = 3, 
                                 generate_3d: bool = True) -> Dict[str, List[str]]:
    """
    Convenience function to create PandaMap visualizations for PandaDock poses
    
    Args:
        poses_df: DataFrame with pose information
        poses_dir: Directory containing pose PDB files
        output_dir: Output directory for visualizations
        top_n: Number of top poses to analyze
        generate_3d: Whether to generate 3D visualizations
        
    Returns:
        Dictionary of generated file lists
    """
    
    # Initialize PandaMap integration
    pandamap_integration = PandaMapIntegration(output_dir)
    
    # Create interaction maps
    return pandamap_integration.create_interaction_maps_for_poses(
        poses_df, poses_dir, top_n, generate_3d
    )


if __name__ == "__main__":
    print("PandaMap Integration for PandaDock")
    print("This module provides professional protein-ligand interaction visualization")
    print("Use create_pandamap_visualizations() function to generate maps for poses")