# -*- coding: utf-8 -*-
"""
HTML report generation for PandaDock results
"""

import os
import json
import base64
from datetime import datetime
from typing import List, Dict, Any, Optional
import logging
from pathlib import Path

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.base_engine import Pose
from utils.ic50_calculator import IC50Calculator

# Import new plotting modules
try:
    from .plot_generator import create_plots_for_pandadock
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False


class HTMLReportGenerator:
    """
    HTML report generator for docking results
    
    Features:
    - Interactive pose visualization
    - Comprehensive scoring breakdown
    - Binding affinity analysis
    - Interaction summaries
    - Molecular property tables
    - Comparison charts
    """
    
    def __init__(self, config=None):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.ic50_calc = IC50Calculator()
        
        # Report template
        self.html_template = self._get_html_template()
        
        # Plotting configuration
        self.enable_plots = True
        self.generate_interaction_maps = True
        self.generate_txt_report = True
        
        self.logger.info("Initialized HTMLReportGenerator with plotting support: %s", PLOTTING_AVAILABLE)
    
    def generate_report(self, results: List[Pose], output_path: str = None) -> str:
        """
        Generate comprehensive HTML report
        
        Args:
            results: List of docked poses
            output_path: Output file path (optional)
            
        Returns:
            Path to generated HTML report
        """
        self.logger.info(f"Generating HTML report for {len(results)} poses")
        
        if not results:
            self.logger.warning("No results to report")
            return ""
        
        # Prepare data for report
        report_data = self._prepare_report_data(results)
        
        # Generate HTML content
        html_content = self._generate_html_content(report_data)
        
        # Determine output path
        if output_path is None:
            output_dir = self.config.io.output_dir if self.config else "output"
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, "pandadock_report.html")
        
        # Write HTML file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report generated: {output_path}")
        return output_path
    
    def generate_comprehensive_report(self, results: List[Pose], output_dir: str,
                                    protein_name: str = None, ligand_name: str = None,
                                    algorithm_info: Dict[str, Any] = None,
                                    command_info: Dict[str, Any] = None) -> Dict[str, str]:
        """
        Generate comprehensive report with plots, interaction maps, and TXT report
        
        Args:
            results: List of docked poses
            output_dir: Output directory for all files
            protein_name: Name of the protein target
            ligand_name: Name of the ligand
            algorithm_info: Algorithm information dictionary
            command_info: Command information dictionary
            
        Returns:
            Dictionary of generated file paths
        """
        self.logger.info(f"Generating comprehensive report for {len(results)} poses")
        
        if not results:
            self.logger.warning("No results to report")
            return {}
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        generated_files = {}
        
        # 1. Generate standard HTML report
        html_report_path = output_dir / "pandadock_report.html"
        html_path = self.generate_report(results, str(html_report_path))
        generated_files['html_report'] = html_path
        
        # 2. Generate CSV export
        csv_path = output_dir / "pandadock_report.csv"
        csv_export = self.export_data(results, format='csv', output_path=str(csv_path))
        generated_files['csv_report'] = csv_export
        
        # 3. Generate JSON export
        json_path = output_dir / "pandadock_report.json"
        json_export = self.export_data(results, format='json', output_path=str(json_path))
        generated_files['json_report'] = json_export
        
        # 4. Create poses summary CSV for plotting
        poses_csv_path = output_dir / "poses" / "poses_summary.csv"
        poses_csv_path.parent.mkdir(exist_ok=True)
        poses_summary = self._create_poses_summary_csv(results, str(poses_csv_path))
        generated_files['poses_csv'] = poses_summary
        
        # 5. Generate comprehensive plots if available
        if PLOTTING_AVAILABLE and self.enable_plots:
            try:
                self.logger.info("Generating comprehensive plots...")
                
                # Set default values
                if algorithm_info is None:
                    algorithm_info = {
                        'algorithm': 'PandaDock',
                        'version': 'Latest',
                        'scoring_function': self.config.scoring.scoring_function if self.config else 'Auto',
                        'engine': self.config.docking.mode if self.config else 'Auto',
                        'mode': self.config.docking.mode if self.config else 'Balanced'
                    }
                
                if command_info is None:
                    command_info = {
                        'command': 'python -m pandadock',
                        'protein': 'protein.pdb',
                        'ligand': 'ligand.pdb',
                        'center': 'Auto-detected',
                        'size': 'Auto-detected',
                        'exhaustiveness': self.config.docking.exhaustiveness if self.config else 'Default'
                    }
                
                # Generate all plots
                plot_files = create_plots_for_pandadock(
                    output_dir=str(output_dir),
                    poses_csv=poses_summary,
                    poses_dir=str(poses_csv_path.parent),
                    protein_name=protein_name,
                    ligand_name=ligand_name,
                    algorithm_info=algorithm_info,
                    command_info=command_info
                )
                
                generated_files.update(plot_files)
                self.logger.info(f"Generated {len(plot_files)} plot files")
                
            except Exception as e:
                self.logger.warning(f"Failed to generate plots: {e}")
        
        # 6. Generate interaction maps using PandaMap if available
        if PLOTTING_AVAILABLE and self.generate_interaction_maps:
            try:
                self.logger.info("Generating 2D interaction maps with PandaMap...")
                
                # Try to use PandaMap integration
                from .pandamap_integration import create_pandamap_visualizations
                
                poses_df = self._create_poses_dataframe(results)
                pandamap_files = create_pandamap_visualizations(
                    poses_df=poses_df,
                    poses_dir=str(poses_csv_path.parent),
                    output_dir=str(output_dir),
                    top_n=3,
                    generate_3d=False  # Only 2D for this function
                )
                
                # Add PandaMap files to generated files
                if pandamap_files and '2d_maps' in pandamap_files:
                    generated_files['interaction_maps'] = pandamap_files['2d_maps']
                    self.logger.info(f"Generated {len(pandamap_files['2d_maps'])} PandaMap interaction maps")
                else:
                    self.logger.info("PandaMap interaction maps not generated - skipping")
                
            except ImportError:
                self.logger.info("PandaMap not available - skipping interaction maps")
            except Exception as e:
                self.logger.warning(f"Failed to generate PandaMap interaction maps: {e}")
        
        self.logger.info(f"Comprehensive report generation completed. Generated {len(generated_files)} files.")
        return generated_files
    
    def _create_poses_summary_csv(self, results: List[Pose], output_path: str) -> str:
        """Create poses summary CSV for plotting"""
        import csv
        
        with open(output_path, 'w', newline='') as csvfile:
            fieldnames = ['Rank', 'Pose_ID', 'Score', 'Energy', 'Confidence', 
                         'Binding_Affinity', 'IC50_uM', 'EC50_uM', 'Ligand_Efficiency', 'Clash_Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for i, pose in enumerate(results):
                binding_affinity = pose.get_binding_affinity()
                ic50_um = pose.get_ic50(units='uM')
                ec50_um = pose.get_ec50(units='uM')
                num_heavy_atoms = len(pose.coordinates)
                le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
                
                writer.writerow({
                    'Rank': i + 1,
                    'Pose_ID': pose.pose_id,
                    'Score': pose.score,
                    'Energy': pose.energy,
                    'Confidence': pose.confidence,
                    'Binding_Affinity': binding_affinity,
                    'IC50_uM': f"{ic50_um:.2e}" if ic50_um != float('inf') else '',
                    'EC50_uM': f"{ec50_um:.2e}" if ec50_um != float('inf') else '',
                    'Ligand_Efficiency': le,
                    'Clash_Score': pose.clash_score
                })
        
        return output_path
    
    def _create_poses_dataframe(self, results: List[Pose]):
        """Create pandas DataFrame from poses for plotting"""
        import pandas as pd
        
        data = []
        for i, pose in enumerate(results):
            binding_affinity = pose.get_binding_affinity()
            ic50_um = pose.get_ic50(units='uM')
            ec50_um = pose.get_ec50(units='uM')
            num_heavy_atoms = len(pose.coordinates)
            le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
            
            data.append({
                'Rank': i + 1,
                'Pose_ID': pose.pose_id,
                'Score': pose.score,
                'Energy': pose.energy,
                'Confidence': pose.confidence,
                'Binding_Affinity': binding_affinity,
                'IC50_uM': f"{ic50_um:.2e}" if ic50_um != float('inf') else '',
                'EC50_uM': f"{ec50_um:.2e}" if ec50_um != float('inf') else '',
                'Ligand_Efficiency': le,
                'Clash_Score': pose.clash_score
            })
        
        return pd.DataFrame(data)
    
    def _prepare_report_data(self, results: List[Pose]) -> Dict[str, Any]:
        """Prepare data for HTML report"""
        report_data = {
            'metadata': self._get_metadata(),
            'summary': self._get_summary_stats(results),
            'poses': self._prepare_pose_data(results),
            'interactions': self._analyze_interactions(results),
            'energy_analysis': self._analyze_energies(results),
            'binding_analysis': self._analyze_binding_affinity(results),
            'molecular_properties': self._analyze_molecular_properties(results),
            'comparison_data': self._prepare_comparison_data(results)
        }
        
        return report_data
    
    def _get_metadata(self) -> Dict[str, Any]:
        """Get report metadata"""
        return {
            'timestamp': datetime.now().isoformat(),
            'pandadock_version': '1.0.0',
            'docking_mode': self.config.docking.mode if self.config else 'unknown',
            'scoring_function': self.config.scoring.scoring_function if self.config else 'unknown',
            'num_poses': self.config.docking.num_poses if self.config else 10,
            'flexible_residues': self.config.docking.flexible_residues if self.config else []
        }
    
    def _get_summary_stats(self, results: List[Pose]) -> Dict[str, Any]:
        """Get summary statistics"""
        if not results:
            return {}
        
        scores = [pose.score for pose in results]
        energies = [pose.energy for pose in results]
        
        summary = {
            'total_poses': len(results),
            'best_score': min(scores),
            'worst_score': max(scores),
            'mean_score': sum(scores) / len(scores),
            'score_range': max(scores) - min(scores),
            'best_energy': min(energies),
            'worst_energy': max(energies),
            'mean_energy': sum(energies) / len(energies),
            'energy_range': max(energies) - min(energies)
        }
        
        # Calculate binding affinities
        binding_affinities = []
        ic50_values_um = []
        ec50_values_um = []
        
        for pose in results:
            binding_affinity = pose.get_binding_affinity()
            ic50_um = pose.get_ic50(units='uM')
            ec50_um = pose.get_ec50(units='uM')
            
            if binding_affinity != float('inf'):
                binding_affinities.append(binding_affinity)
            if ic50_um != float('inf'):
                ic50_values_um.append(ic50_um)
            if ec50_um != float('inf'):
                ec50_values_um.append(ec50_um)
        
        if binding_affinities:
            summary['best_binding_affinity'] = min(binding_affinities)
            summary['mean_binding_affinity'] = sum(binding_affinities) / len(binding_affinities)
        
        if ic50_values_um:
            summary['best_ic50'] = min(ic50_values_um)  # For backward compatibility
            summary['best_ic50_um'] = min(ic50_values_um)
            summary['mean_ic50_um'] = sum(ic50_values_um) / len(ic50_values_um)
            
        if ec50_values_um:
            summary['best_ec50_um'] = min(ec50_values_um)
            summary['mean_ec50_um'] = sum(ec50_values_um) / len(ec50_values_um)
        
        return summary
    
    def _prepare_pose_data(self, results: List[Pose]) -> List[Dict[str, Any]]:
        """Prepare pose data for report"""
        pose_data = []
        
        for i, pose in enumerate(results):
            # Calculate additional metrics
            binding_affinity = pose.get_binding_affinity()
            ic50_um = pose.get_ic50(units='uM')
            ec50_um = pose.get_ec50(units='uM')
            
            # Estimate ligand efficiency
            num_heavy_atoms = len(pose.coordinates)  # Simplified
            le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
            
            data = {
                'rank': i + 1,
                'pose_id': pose.pose_id,
                'score': pose.score,
                'energy': pose.energy,
                'binding_affinity': binding_affinity,
                'ic50': ic50_um,  # For backward compatibility
                'ic50_nm': ic50_um * 1000 if ic50_um != float('inf') else float('inf'),  # For backward compatibility
                'ic50_um': ic50_um if ic50_um != float('inf') else float('inf'),
                'ec50_um': ec50_um if ec50_um != float('inf') else float('inf'),
                'ligand_efficiency': le,
                'rmsd': pose.rmsd,
                'confidence': pose.confidence,
                'clash_score': pose.clash_score,
                'binding_site_coverage': pose.binding_site_coverage,
                'energy_breakdown': {
                    'vdw': pose.vdw_energy,
                    'electrostatic': pose.electrostatic_energy,
                    'hbond': pose.hbond_energy,
                    'hydrophobic': pose.hydrophobic_energy,
                    'solvation': pose.solvation_energy,
                    'entropy': pose.entropy_energy
                },
                'interactions': {
                    'hbonds': len(pose.hbond_interactions),
                    'hydrophobic': len(pose.hydrophobic_interactions),
                    'salt_bridges': len(pose.salt_bridge_interactions)
                },
                'flexible_residues': pose.flexible_residues,
                'coordinates': pose.coordinates.tolist()
            }
            
            pose_data.append(data)
        
        return pose_data
    
    def _analyze_interactions(self, results: List[Pose]) -> Dict[str, Any]:
        """Analyze interactions across all poses"""
        interaction_analysis = {
            'hbond_stats': {'total': 0, 'average': 0, 'max': 0, 'min': float('inf')},
            'hydrophobic_stats': {'total': 0, 'average': 0, 'max': 0, 'min': float('inf')},
            'salt_bridge_stats': {'total': 0, 'average': 0, 'max': 0, 'min': float('inf')},
            'interaction_frequency': {},
            'best_interactions': []
        }
        
        if not results:
            return interaction_analysis
        
        # Collect interaction statistics
        hbond_counts = [len(pose.hbond_interactions) for pose in results]
        hydrophobic_counts = [len(pose.hydrophobic_interactions) for pose in results]
        salt_bridge_counts = [len(pose.salt_bridge_interactions) for pose in results]
        
        # H-bond statistics
        interaction_analysis['hbond_stats'] = {
            'total': sum(hbond_counts),
            'average': sum(hbond_counts) / len(hbond_counts),
            'max': max(hbond_counts),
            'min': min(hbond_counts)
        }
        
        # Hydrophobic statistics
        interaction_analysis['hydrophobic_stats'] = {
            'total': sum(hydrophobic_counts),
            'average': sum(hydrophobic_counts) / len(hydrophobic_counts),
            'max': max(hydrophobic_counts),
            'min': min(hydrophobic_counts)
        }
        
        # Salt bridge statistics
        interaction_analysis['salt_bridge_stats'] = {
            'total': sum(salt_bridge_counts),
            'average': sum(salt_bridge_counts) / len(salt_bridge_counts),
            'max': max(salt_bridge_counts),
            'min': min(salt_bridge_counts)
        }
        
        # Find best interactions (from top pose)
        if results:
            best_pose = results[0]
            interaction_analysis['best_interactions'] = {
                'hbonds': best_pose.hbond_interactions[:5],  # Top 5
                'hydrophobic': best_pose.hydrophobic_interactions[:5],
                'salt_bridges': best_pose.salt_bridge_interactions[:5]
            }
        
        return interaction_analysis
    
    def _analyze_energies(self, results: List[Pose]) -> Dict[str, Any]:
        """Analyze energy contributions"""
        energy_analysis = {
            'energy_contributions': {},
            'energy_correlations': {},
            'energy_distribution': {}
        }
        
        if not results:
            return energy_analysis
        
        # Calculate average energy contributions
        energy_terms = ['vdw_energy', 'electrostatic_energy', 'hbond_energy', 
                       'hydrophobic_energy', 'solvation_energy', 'entropy_energy']
        
        for term in energy_terms:
            values = [getattr(pose, term) for pose in results]
            energy_analysis['energy_contributions'][term] = {
                'average': sum(values) / len(values),
                'min': min(values),
                'max': max(values),
                'std': (sum((x - sum(values)/len(values))**2 for x in values) / len(values))**0.5
            }
        
        # Energy distribution for top poses
        top_poses = results[:min(10, len(results))]
        energy_analysis['energy_distribution'] = {
            'poses': [pose.pose_id for pose in top_poses],
            'total_energy': [pose.energy for pose in top_poses],
            'vdw': [pose.vdw_energy for pose in top_poses],
            'electrostatic': [pose.electrostatic_energy for pose in top_poses],
            'hbond': [pose.hbond_energy for pose in top_poses],
            'hydrophobic': [pose.hydrophobic_energy for pose in top_poses],
            'solvation': [pose.solvation_energy for pose in top_poses],
            'entropy': [pose.entropy_energy for pose in top_poses]
        }
        
        return energy_analysis
    
    def _analyze_binding_affinity(self, results: List[Pose]) -> Dict[str, Any]:
        """Analyze binding affinity and related metrics"""
        binding_analysis = {
            'affinity_distribution': [],
            'ic50_distribution': [],
            'ligand_efficiency_distribution': [],
            'drug_likeness_assessment': []
        }
        
        if not results:
            return binding_analysis
        
        for pose in results:
            binding_affinity = pose.get_binding_affinity()
            ic50 = pose.get_ic50()
            
            # Estimate molecular properties (simplified)
            num_heavy_atoms = len(pose.coordinates)
            molecular_weight = num_heavy_atoms * 15.0  # Rough estimate
            logp = 2.5  # Placeholder
            
            # Calculate ligand efficiency
            le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
            
            binding_analysis['affinity_distribution'].append(binding_affinity)
            binding_analysis['ic50_distribution'].append(ic50)
            binding_analysis['ligand_efficiency_distribution'].append(le)
            
            # Drug-likeness assessment
            drug_assessment = self.ic50_calc.assess_drug_likeness(
                binding_affinity, molecular_weight, logp, 2, 4, num_heavy_atoms
            )
            binding_analysis['drug_likeness_assessment'].append(drug_assessment)
        
        return binding_analysis
    
    def _analyze_molecular_properties(self, results: List[Pose]) -> Dict[str, Any]:
        """Analyze molecular properties"""
        molecular_analysis = {
            'size_distribution': [],
            'flexibility_distribution': [],
            'shape_descriptors': []
        }
        
        if not results:
            return molecular_analysis
        
        for pose in results:
            # Size metrics
            num_atoms = len(pose.coordinates)
            molecular_analysis['size_distribution'].append(num_atoms)
            
            # Flexibility (number of flexible residues)
            flexibility = len(pose.flexible_residues)
            molecular_analysis['flexibility_distribution'].append(flexibility)
            
            # Shape descriptors
            if len(pose.coordinates) > 0:
                # Calculate radius of gyration
                center = pose.coordinates.mean(axis=0)
                distances = ((pose.coordinates - center)**2).sum(axis=1)**0.5
                rg = (distances**2).mean()**0.5
                
                # Calculate dimensions
                min_coords = pose.coordinates.min(axis=0)
                max_coords = pose.coordinates.max(axis=0)
                dimensions = max_coords - min_coords
                
                molecular_analysis['shape_descriptors'].append({
                    'radius_of_gyration': rg,
                    'length': dimensions[0],
                    'width': dimensions[1],
                    'height': dimensions[2]
                })
        
        return molecular_analysis
    
    def _prepare_comparison_data(self, results: List[Pose]) -> Dict[str, Any]:
        """Prepare data for comparison charts"""
        comparison_data = {
            'score_vs_energy': [],
            'score_vs_interactions': [],
            'energy_vs_interactions': [],
            'affinity_vs_efficiency': []
        }
        
        if not results:
            return comparison_data
        
        for pose in results:
            binding_affinity = pose.get_binding_affinity()
            num_heavy_atoms = len(pose.coordinates)
            le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
            
            total_interactions = (len(pose.hbond_interactions) + 
                                len(pose.hydrophobic_interactions) + 
                                len(pose.salt_bridge_interactions))
            
            comparison_data['score_vs_energy'].append([pose.score, pose.energy])
            comparison_data['score_vs_interactions'].append([pose.score, total_interactions])
            comparison_data['energy_vs_interactions'].append([pose.energy, total_interactions])
            comparison_data['affinity_vs_efficiency'].append([binding_affinity, le])
        
        return comparison_data
    
    def _generate_html_content(self, report_data: Dict[str, Any]) -> str:
        """Generate HTML content from report data"""
        # Convert data to JSON for JavaScript
        json_data = json.dumps(report_data, indent=2, default=str)
        
        # Replace placeholders in template
        html_content = self.html_template.replace('{{REPORT_DATA}}', json_data)
        html_content = html_content.replace('{{TIMESTAMP}}', report_data['metadata']['timestamp'])
        html_content = html_content.replace('{{PANDADOCK_VERSION}}', report_data['metadata']['pandadock_version'])
        
        return html_content
    
    def _get_html_template(self) -> str:
        """Get HTML template for report"""
        return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PandaDock Docking Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            color: #333;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .header {
            text-align: center;
            margin-bottom: 30px;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 20px;
        }
        .header h1 {
            color: #4CAF50;
            margin: 0;
            font-size: 2.5em;
        }
        .subtitle {
            color: #666;
            margin-top: 10px;
        }
        .section {
            margin-bottom: 30px;
        }
        .section h2 {
            color: #4CAF50;
            border-bottom: 1px solid #ddd;
            padding-bottom: 10px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        .summary-card {
            background-color: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #4CAF50;
        }
        .summary-card h3 {
            margin-top: 0;
            color: #4CAF50;
        }
        .summary-value {
            font-size: 1.5em;
            font-weight: bold;
            color: #333;
        }
        .table-container {
            overflow-x: auto;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .pose-row {
            cursor: pointer;
        }
        .pose-row.selected {
            background-color: #e8f5e8;
        }
        .chart-container {
            margin: 20px 0;
            padding: 20px;
            background-color: #f8f9fa;
            border-radius: 8px;
        }
        .interaction-details {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin-top: 15px;
        }
        .energy-breakdown {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 10px;
            margin-top: 15px;
        }
        .energy-item {
            background-color: white;
            padding: 10px;
            border-radius: 5px;
            text-align: center;
            border: 1px solid #ddd;
        }
        .energy-label {
            font-size: 0.9em;
            color: #666;
        }
        .energy-value {
            font-size: 1.2em;
            font-weight: bold;
            color: #333;
        }
        .tabs {
            display: flex;
            background-color: #f0f0f0;
            border-radius: 8px 8px 0 0;
        }
        .tab {
            padding: 15px 25px;
            cursor: pointer;
            border: none;
            background-color: transparent;
            font-size: 16px;
            transition: background-color 0.3s;
        }
        .tab.active {
            background-color: #4CAF50;
            color: white;
        }
        .tab-content {
            display: none;
            padding: 20px;
            background-color: white;
            border-radius: 0 0 8px 8px;
            border: 1px solid #ddd;
        }
        .tab-content.active {
            display: block;
        }
        .footer {
            text-align: center;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
        }
        .logo {
            width: 60px;
            height: 60px;
            background-color: #4CAF50;
            border-radius: 50%;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-size: 24px;
            font-weight: bold;
            margin-right: 15px;
        }
        .header-content {
            display: flex;
            align-items: center;
            justify-content: center;
        }
        .pose-details {
            display: none;
            margin-top: 15px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 8px;
        }
        .pose-details.active {
            display: block;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <div class="logo">PD</div>
                <div>
                    <h1>PandaDock Docking Report</h1>
                    <div class="subtitle">
                        Generated on {{TIMESTAMP}} | Version {{PANDADOCK_VERSION}}
                    </div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Executive Summary</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Poses</h3>
                    <div class="summary-value" id="total-poses">-</div>
                </div>
                <div class="summary-card">
                    <h3>Best Score</h3>
                    <div class="summary-value" id="best-score">-</div>
                </div>
                <div class="summary-card">
                    <h3>Best Binding Affinity</h3>
                    <div class="summary-value" id="best-affinity">-</div>
                </div>
                <div class="summary-card">
                    <h3>Best IC50</h3>
                    <div class="summary-value" id="best-ic50">-</div>
                </div>
                <div class="summary-card">
                    <h3>Best EC50</h3>
                    <div class="summary-value" id="best-ec50">-</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Detailed Results</h2>
            <div class="tabs">
                <button class="tab active" onclick="showTab('poses')">Poses</button>
                <button class="tab" onclick="showTab('interactions')">Interactions</button>
                <button class="tab" onclick="showTab('energies')">Energies</button>
                <button class="tab" onclick="showTab('properties')">Properties</button>
            </div>

            <div id="poses-tab" class="tab-content active">
                <div class="table-container">
                    <table id="poses-table">
                        <thead>
                            <tr>
                                <th>Rank</th>
                                <th>Pose ID</th>
                                <th>Score</th>
                                <th>Energy (kcal/mol)</th>
                                <th>ΔG (kcal/mol)</th>
                                <th>IC50 (μM)</th>
                                <th>EC50 (μM)</th>
                                <th>LE (kcal/mol/atom)</th>
                                <th>Interactions</th>
                                <th>Clash Score</th>
                            </tr>
                        </thead>
                        <tbody id="poses-tbody">
                        </tbody>
                    </table>
                </div>
            </div>

            <div id="interactions-tab" class="tab-content">
                <div class="interaction-details">
                    <h3>Interaction Statistics</h3>
                    <div id="interaction-stats"></div>
                </div>
            </div>

            <div id="energies-tab" class="tab-content">
                <div class="chart-container">
                    <h3>Energy Breakdown</h3>
                    <div id="energy-chart"></div>
                </div>
            </div>

            <div id="properties-tab" class="tab-content">
                <div class="chart-container">
                    <h3>Molecular Properties</h3>
                    <div id="properties-chart"></div>
                </div>
            </div>
        </div>

        <div class="footer">
            <p>Generated by PandaDock v{{PANDADOCK_VERSION}} | 
               <a href="https://github.com/pandadock/pandadock" target="_blank">GitHub</a> | 
               <a href="https://pandadock.readthedocs.io" target="_blank">Documentation</a>
            </p>
        </div>
    </div>

    <script>
        // Report data
        const reportData = {{REPORT_DATA}};
        
        // Scientific notation formatter for HTML display
        function formatScientificNotation(value, decimals = 2) {
            if (value === null || value === undefined || value === Infinity || value === -Infinity) {
                return '-';
            }
            
            if (value === 0) {
                return '0';
            }
            
            const exponent = Math.floor(Math.log10(Math.abs(value)));
            const mantissa = value / Math.pow(10, exponent);
            
            if (Math.abs(exponent) <= 2) {
                // For small exponents, use regular notation
                return value.toFixed(decimals);
            }
            
            // Use scientific notation with × symbol
            const mantissaFormatted = mantissa.toFixed(decimals);
            return `${mantissaFormatted} × 10<sup>${exponent}</sup>`;
        }
        
        // Initialize report
        document.addEventListener('DOMContentLoaded', function() {
            initializeReport();
        });
        
        function initializeReport() {
            updateSummary();
            populatePosesTable();
            updateInteractionStats();
            updateEnergyChart();
            updatePropertiesChart();
        }
        
        function updateSummary() {
            const summary = reportData.summary;
            
            document.getElementById('total-poses').textContent = summary.total_poses || '-';
            document.getElementById('best-score').textContent = (summary.best_score || 0).toFixed(3);
            document.getElementById('best-affinity').textContent = (summary.best_binding_affinity || 0).toFixed(2) + ' kcal/mol';
            
            const bestIC50 = summary.best_ic50_um;
            if (bestIC50 && bestIC50 !== Infinity) {
                document.getElementById('best-ic50').innerHTML = formatScientificNotation(bestIC50, 2) + ' μM';
            } else {
                document.getElementById('best-ic50').textContent = '-';
            }
            
            const bestEC50 = summary.best_ec50_um;
            if (bestEC50 && bestEC50 !== Infinity) {
                document.getElementById('best-ec50').innerHTML = formatScientificNotation(bestEC50, 2) + ' μM';
            } else {
                document.getElementById('best-ec50').textContent = '-';
            }
        }
        
        function populatePosesTable() {
            const tbody = document.getElementById('poses-tbody');
            tbody.innerHTML = '';
            
            reportData.poses.forEach((pose, index) => {
                const row = document.createElement('tr');
                row.className = 'pose-row';
                row.onclick = () => showPoseDetails(index);
                
                const totalInteractions = pose.interactions.hbonds + 
                                        pose.interactions.hydrophobic + 
                                        pose.interactions.salt_bridges;
                
                row.innerHTML = `
                    <td>${pose.rank}</td>
                    <td>${pose.pose_id}</td>
                    <td>${pose.score.toFixed(3)}</td>
                    <td>${pose.energy.toFixed(2)}</td>
                    <td>${pose.binding_affinity.toFixed(2)}</td>
                    <td>${pose.ic50_um === Infinity ? '-' : formatScientificNotation(pose.ic50_um, 2)}</td>
                    <td>${pose.ec50_um === Infinity ? '-' : formatScientificNotation(pose.ec50_um, 2)}</td>
                    <td>${pose.ligand_efficiency.toFixed(3)}</td>
                    <td>${totalInteractions}</td>
                    <td>${pose.clash_score.toFixed(2)}</td>
                `;
                
                tbody.appendChild(row);
                
                // Add details row
                const detailRow = document.createElement('tr');
                detailRow.innerHTML = `
                    <td colspan="10">
                        <div class="pose-details" id="pose-details-${index}">
                            <h4>Energy Breakdown</h4>
                            <div class="energy-breakdown">
                                <div class="energy-item">
                                    <div class="energy-label">van der Waals</div>
                                    <div class="energy-value">${pose.energy_breakdown.vdw.toFixed(2)}</div>
                                </div>
                                <div class="energy-item">
                                    <div class="energy-label">Electrostatic</div>
                                    <div class="energy-value">${pose.energy_breakdown.electrostatic.toFixed(2)}</div>
                                </div>
                                <div class="energy-item">
                                    <div class="energy-label">H-bonds</div>
                                    <div class="energy-value">${pose.energy_breakdown.hbond.toFixed(2)}</div>
                                </div>
                                <div class="energy-item">
                                    <div class="energy-label">Hydrophobic</div>
                                    <div class="energy-value">${pose.energy_breakdown.hydrophobic.toFixed(2)}</div>
                                </div>
                                <div class="energy-item">
                                    <div class="energy-label">Solvation</div>
                                    <div class="energy-value">${pose.energy_breakdown.solvation.toFixed(2)}</div>
                                </div>
                                <div class="energy-item">
                                    <div class="energy-label">Entropy</div>
                                    <div class="energy-value">${pose.energy_breakdown.entropy.toFixed(2)}</div>
                                </div>
                            </div>
                            <h4>Interactions</h4>
                            <p>Hydrogen bonds: ${pose.interactions.hbonds} | 
                               Hydrophobic contacts: ${pose.interactions.hydrophobic} | 
                               Salt bridges: ${pose.interactions.salt_bridges}</p>
                            ${pose.flexible_residues.length > 0 ? 
                              `<h4>Flexible Residues</h4><p>${pose.flexible_residues.join(', ')}</p>` : ''}
                        </div>
                    </td>
                `;
                tbody.appendChild(detailRow);
            });
        }
        
        function showPoseDetails(index) {
            // Hide all details
            document.querySelectorAll('.pose-details').forEach(detail => {
                detail.classList.remove('active');
            });
            
            // Remove selection from all rows
            document.querySelectorAll('.pose-row').forEach(row => {
                row.classList.remove('selected');
            });
            
            // Show selected details
            const detail = document.getElementById(`pose-details-${index}`);
            const row = document.querySelectorAll('.pose-row')[index];
            
            if (detail.classList.contains('active')) {
                detail.classList.remove('active');
                row.classList.remove('selected');
            } else {
                detail.classList.add('active');
                row.classList.add('selected');
            }
        }
        
        function updateInteractionStats() {
            const interactions = reportData.interactions;
            const statsDiv = document.getElementById('interaction-stats');
            
            statsDiv.innerHTML = `
                <div class="energy-breakdown">
                    <div class="energy-item">
                        <div class="energy-label">Hydrogen Bonds</div>
                        <div class="energy-value">Avg: ${interactions.hbond_stats.average.toFixed(1)}</div>
                        <div class="energy-label">Max: ${interactions.hbond_stats.max}</div>
                    </div>
                    <div class="energy-item">
                        <div class="energy-label">Hydrophobic</div>
                        <div class="energy-value">Avg: ${interactions.hydrophobic_stats.average.toFixed(1)}</div>
                        <div class="energy-label">Max: ${interactions.hydrophobic_stats.max}</div>
                    </div>
                    <div class="energy-item">
                        <div class="energy-label">Salt Bridges</div>
                        <div class="energy-value">Avg: ${interactions.salt_bridge_stats.average.toFixed(1)}</div>
                        <div class="energy-label">Max: ${interactions.salt_bridge_stats.max}</div>
                    </div>
                </div>
            `;
        }
        
        function updateEnergyChart() {
            const energyDiv = document.getElementById('energy-chart');
            const energyDist = reportData.energy_analysis.energy_distribution;
            const energyContrib = reportData.energy_analysis.energy_contributions;
            
            if (energyDist.poses && energyDist.poses.length > 0) {
                let html = '<div class="energy-analysis">';
                
                // Energy contributions overview
                html += '<h4>Average Energy Contributions</h4>';
                html += '<div class="energy-breakdown">';
                
                Object.keys(energyContrib).forEach(component => {
                    const data = energyContrib[component];
                    const color = component.includes('hbond') || component.includes('hydrophobic') || component.includes('solvation') ? '#4CAF50' : 
                                 component.includes('vdw') || component.includes('entropy') ? '#f44336' : '#2196F3';
                    
                    html += `
                        <div class="energy-item" style="border-left: 4px solid ${color}">
                            <div class="energy-label">${component.replace('_', ' ').toUpperCase()}</div>
                            <div class="energy-value">${data.average.toFixed(2)}</div>
                            <div class="energy-label">Range: ${data.min.toFixed(2)} to ${data.max.toFixed(2)}</div>
                        </div>
                    `;
                });
                
                html += '</div>';
                
                // Energy distribution table
                html += '<h4>Energy Distribution Across Poses</h4>';
                html += '<div class="table-container">';
                html += '<table>';
                html += '<thead><tr><th>Pose</th><th>Total</th><th>vdW</th><th>H-bonds</th><th>Hydrophobic</th><th>Solvation</th><th>Entropy</th></tr></thead>';
                html += '<tbody>';
                
                energyDist.poses.forEach((pose, i) => {
                    html += `<tr>
                        <td>${pose}</td>
                        <td>${energyDist.total_energy[i].toFixed(2)}</td>
                        <td>${energyDist.vdw[i].toFixed(2)}</td>
                        <td>${energyDist.hbond[i].toFixed(2)}</td>
                        <td>${energyDist.hydrophobic[i].toFixed(2)}</td>
                        <td>${energyDist.solvation[i].toFixed(2)}</td>
                        <td>${energyDist.entropy[i].toFixed(2)}</td>
                    </tr>`;
                });
                
                html += '</tbody></table></div></div>';
                energyDiv.innerHTML = html;
            } else {
                energyDiv.innerHTML = '<p>No energy data available</p>';
            }
        }
        
        function updatePropertiesChart() {
            const propertiesDiv = document.getElementById('properties-chart');
            const molecular = reportData.molecular_properties;
            const binding = reportData.binding_analysis;
            
            let html = '<div class="properties-analysis">';
            
            // Molecular size and flexibility
            html += '<h4>Molecular Properties Overview</h4>';
            html += '<div class="energy-breakdown">';
            
            if (molecular.size_distribution && molecular.size_distribution.length > 0) {
                const avgSize = molecular.size_distribution.reduce((a, b) => a + b, 0) / molecular.size_distribution.length;
                const avgFlex = molecular.flexibility_distribution.reduce((a, b) => a + b, 0) / molecular.flexibility_distribution.length;
                
                html += `
                    <div class="energy-item">
                        <div class="energy-label">Average Molecular Size</div>
                        <div class="energy-value">${avgSize.toFixed(0)} atoms</div>
                    </div>
                    <div class="energy-item">
                        <div class="energy-label">Average Flexibility</div>
                        <div class="energy-value">${avgFlex.toFixed(1)} bonds</div>
                    </div>
                `;
            }
            
            // Shape descriptors
            if (molecular.shape_descriptors && molecular.shape_descriptors.length > 0) {
                const avgShape = molecular.shape_descriptors[0];
                html += `
                    <div class="energy-item">
                        <div class="energy-label">Radius of Gyration</div>
                        <div class="energy-value">${avgShape.radius_of_gyration.toFixed(2)} Å</div>
                    </div>
                    <div class="energy-item">
                        <div class="energy-label">Molecular Dimensions</div>
                        <div class="energy-value">${avgShape.length.toFixed(1)} × ${avgShape.width.toFixed(1)} × ${avgShape.height.toFixed(1)} Å</div>
                    </div>
                `;
            }
            
            html += '</div>';
            
            // Drug-likeness assessment
            if (binding.drug_likeness_assessment && binding.drug_likeness_assessment.length > 0) {
                html += '<h4>Drug-likeness Assessment</h4>';
                html += '<div class="table-container">';
                html += '<table>';
                html += '<thead><tr><th>Pose</th><th>LE</th><th>Lipinski</th><th>Veber</th><th>TPSA</th><th>Category</th></tr></thead>';
                html += '<tbody>';
                
                binding.drug_likeness_assessment.forEach((assessment, i) => {
                    const pose = reportData.poses[i];
                    html += `<tr>
                        <td>${pose.pose_id}</td>
                        <td>${assessment.ligand_efficiency.toFixed(3)}</td>
                        <td>${assessment.lipinski_compliant ? 'Pass' : 'Fail'}</td>
                        <td>${assessment.veber_compliant ? 'Pass' : 'Fail'}</td>
                        <td>${assessment.tpsa_estimate}</td>
                        <td style="color: ${assessment.drug_likeness_category === 'good' ? 'green' : 
                                              assessment.drug_likeness_category === 'moderate' ? 'orange' : 'red'}">
                            ${assessment.drug_likeness_category}
                        </td>
                    </tr>`;
                });
                
                html += '</tbody></table></div>';
            }
            
            html += '</div>';
            propertiesDiv.innerHTML = html;
        }
        
        function showTab(tabName) {
            // Hide all tabs
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Remove active class from all tab buttons
            document.querySelectorAll('.tab').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Show selected tab
            document.getElementById(tabName + '-tab').classList.add('active');
            event.target.classList.add('active');
        }
    </script>
</body>
</html>
        """
    
    def generate_summary_report(self, results: List[Pose], output_path: str = None) -> str:
        """Generate a simplified summary report"""
        if not results:
            return ""
        
        summary_data = {
            'total_poses': len(results),
            'best_score': min(pose.score for pose in results),
            'best_energy': min(pose.energy for pose in results),
            'top_poses': results[:5]  # Top 5 poses
        }
        
        # Generate simple text report
        report_lines = [
            "PandaDock Docking Summary Report",
            "=" * 40,
            f"Total poses: {summary_data['total_poses']}",
            f"Best score: {summary_data['best_score']:.3f}",
            f"Best energy: {summary_data['best_energy']:.2f} kcal/mol",
            "",
            "Top 5 Poses:",
            "-" * 20
        ]
        
        for i, pose in enumerate(summary_data['top_poses']):
            binding_affinity = pose.get_binding_affinity()
            ic50 = pose.get_ic50()
            
            report_lines.append(f"Pose {i+1}: {pose.pose_id}")
            report_lines.append(f"  Score: {pose.score:.3f}")
            report_lines.append(f"  Energy: {pose.energy:.2f} kcal/mol")
            report_lines.append(f"  ΔG: {binding_affinity:.2f} kcal/mol")
            report_lines.append(f"  IC50: {ic50:.1f} nM" if ic50 != float('inf') else "  IC50: -")
            report_lines.append("")
        
        report_content = "\n".join(report_lines)
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(report_content)
            return output_path
        
        return report_content
    
    def export_data(self, results: List[Pose], format: str = 'json', output_path: str = None) -> str:
        """Export docking results in various formats"""
        if not results:
            return ""
        
        if format == 'json':
            data = self._prepare_report_data(results)
            content = json.dumps(data, indent=2, default=str)
            extension = '.json'
        
        elif format == 'csv':
            content = self._export_csv(results)
            extension = '.csv'
        
        else:
            raise ValueError(f"Unsupported export format: {format}")
        
        if output_path is None:
            output_dir = self.config.io.output_dir if self.config else "output"
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"pandadock_results{extension}")
        
        with open(output_path, 'w') as f:
            f.write(content)
        
        return output_path
    
    def _export_csv(self, results: List[Pose]) -> str:
        """Export results to CSV format"""
        import csv
        import io
        
        output = io.StringIO()
        writer = csv.writer(output)
        
        # Header
        writer.writerow([
            'Rank', 'Pose_ID', 'Score', 'Energy', 'Binding_Affinity', 'IC50_uM', 'EC50_uM',
            'Ligand_Efficiency', 'VdW_Energy', 'Electrostatic_Energy', 'HBond_Energy',
            'Hydrophobic_Energy', 'Solvation_Energy', 'Entropy_Energy', 'Clash_Score',
            'HBond_Count', 'Hydrophobic_Count', 'Salt_Bridge_Count', 'Flexible_Residues'
        ])
        
        # Data rows
        for i, pose in enumerate(results):
            binding_affinity = pose.get_binding_affinity()
            ic50_um = pose.get_ic50(units='uM')
            ec50_um = pose.get_ec50(units='uM')
            num_heavy_atoms = len(pose.coordinates)
            le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
            
            writer.writerow([
                i + 1,
                pose.pose_id,
                pose.score,
                pose.energy,
                binding_affinity,
                f"{ic50_um:.2e}" if ic50_um != float('inf') else '',
                f"{ec50_um:.2e}" if ec50_um != float('inf') else '',
                le,
                pose.vdw_energy,
                pose.electrostatic_energy,
                pose.hbond_energy,
                pose.hydrophobic_energy,
                pose.solvation_energy,
                pose.entropy_energy,
                pose.clash_score,
                len(pose.hbond_interactions),
                len(pose.hydrophobic_interactions),
                len(pose.salt_bridge_interactions),
                ';'.join(pose.flexible_residues)
            ])
        
        return output.getvalue()