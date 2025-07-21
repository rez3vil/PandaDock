# -*- coding: utf-8 -*-
"""
Output writer for PandaDock results
"""

import os
import json
import logging
from typing import List, Dict, Any
from pathlib import Path


class OutputWriter:
    """
    Handles writing PandaDock results to various output formats
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def write_results(self, poses: List[Any], output_dir: str, format: str = 'html'):
        """
        Write docking results in specified format
        
        Args:
            poses: List of Pose objects
            output_dir: Output directory path
            format: Output format ('html', 'json', 'csv')
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        if format == 'json':
            self._write_json(poses, output_path / "pandadock_report.json")
        elif format == 'csv':
            self._write_csv(poses, output_path / "pandadock_report.csv")
        else:  # Default to HTML
            self._write_html(poses, output_path / "pandadock_report.html")
    
    def _write_json(self, poses: List[Any], filename: Path):
        """Write results as JSON"""
        self.logger.info(f"Writing JSON report to {filename}")
        
        results = {
            'metadata': {
                'num_poses': len(poses),
                'software': 'PandaDock',
                'version': '1.0.0'
            },
            'poses': [pose.to_dict() for pose in poses]
        }
        
        with open(filename, 'w') as f:
            json.dump(results, f, indent=2)
    
    def _write_csv(self, poses: List[Any], filename: Path):
        """Write results as CSV with IC50 and EC50 in scientific notation"""
        self.logger.info(f"Writing CSV report to {filename}")
        
        with open(filename, 'w') as f:
            f.write("Rank,Pose_ID,Score,Energy,Confidence,Binding_Affinity,IC50_uM,EC50_uM,Ligand_Efficiency,Clash_Score\n")
            for i, pose in enumerate(poses):
                ic50_um = pose.get_ic50(units='uM')
                ec50_um = pose.get_ec50(units='uM')
                f.write(f"{i+1},{pose.pose_id},{pose.score:.6f},{pose.energy:.6f},"
                       f"{pose.confidence:.6f},{pose.get_binding_affinity():.6f},"
                       f"{ic50_um:.2e},{ec50_um:.2e},{pose.ligand_efficiency:.6f},{pose.clash_score:.6f}\n")
    
    def _write_html(self, poses: List[Any], filename: Path):
        """Write results as HTML (placeholder for actual HTML report generation)"""
        self.logger.info(f"Writing HTML report to {filename}")
        
        # This would be implemented by the HTML report generator
        # For now, just create a simple HTML file
        with open(filename, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>PandaDock Results</title>
</head>
<body>
    <h1>PandaDock Docking Results</h1>
    <p>HTML report generation not implemented. Please use JSON or CSV format.</p>
</body>
</html>""")
    
    def save_complex_files(self, poses: List[Any], protein_file: str, output_dir: str):
        """
        Save protein-ligand complex files
        
        Args:
            poses: List of Pose objects
            protein_file: Path to protein file
            output_dir: Output directory
        """
        self.logger.info(f"Saving complex files to {output_dir}")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for i, pose in enumerate(poses):
            complex_filename = output_path / f"complex_{i+1}_{pose.pose_id}.pdb"
            self._save_complex_pdb(pose, protein_file, complex_filename)
    
    def _save_complex_pdb(self, pose: Any, protein_file: str, filename: Path):
        """Save a single protein-ligand complex as PDB"""
        with open(filename, 'w') as f:
            f.write("REMARK PandaDock protein-ligand complex\n")
            f.write(f"REMARK Pose ID: {pose.pose_id}\n")
            f.write(f"REMARK Score: {pose.score:.3f}\n")
            f.write(f"REMARK Binding Affinity: {pose.get_binding_affinity():.2f} kcal/mol\n")
            f.write(f"REMARK IC50: {pose.get_ic50():.1f} nM\n")
            f.write("REMARK\n")
            
            # Copy protein structure
            try:
                with open(protein_file, 'r') as protein_f:
                    for line in protein_f:
                        if line.startswith(('ATOM', 'HETATM')) and not line[17:20].strip() == 'LIG':
                            f.write(line)
            except FileNotFoundError:
                f.write("REMARK Warning: Could not read protein file\n")
            
            f.write("TER\n")
            
            # Write ligand
            for i, coord in enumerate(pose.coordinates):
                f.write(f"HETATM{i+1:5d}  C    LIG B   1    "
                       f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                       f"  1.00 20.00           C\n")
            
            f.write("END\n")
