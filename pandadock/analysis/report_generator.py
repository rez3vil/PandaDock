"""
Comprehensive docking report generation for molecular docking results.

This module provides methods for generating detailed HTML, PDF, and text reports
that summarize docking results, clustering analysis, and interaction analysis.
"""

import numpy as np
import logging
from typing import List, Dict, Any, Optional, Union
import os
from datetime import datetime

logger = logging.getLogger(__name__)


class DockingReportGenerator:
    """Generate comprehensive docking analysis reports."""
    
    def __init__(self, report_format: str = 'html', 
                 include_sections: Optional[List[str]] = None):
        """
        Initialize report generator.
        
        Args:
            report_format: Report format ('html', 'pdf', or 'txt')
            include_sections: Sections to include in the report
        """
        self.report_format = report_format.lower()
        self.include_sections = include_sections or [
            'summary', 'top_poses', 'clusters', 'interactions', 'energetics'
        ]
        self.logger = logging.getLogger(__name__)
    
    def generate_report(self, protein: Any, poses: List[Any], scores: List[float], 
                       output_file: str, clustering_results: Optional[Dict] = None,
                       energy_decomposition: Optional[Dict] = None,
                       interaction_analysis: Optional[Dict] = None) -> str:
        """
        Generate a comprehensive docking report.
        
        Args:
            protein: Protein object
            poses: List of ligand poses
            scores: Corresponding scores for each pose
            output_file: Output file path
            clustering_results: Clustering results from PoseClusterer
            energy_decomposition: Energy decomposition from EnergyDecomposition
            interaction_analysis: Interaction analysis results
            
        Returns:
            Path to generated report
        """
        try:
            if not poses or not scores:
                raise ValueError("No poses or scores provided")
                
            if len(poses) != len(scores):
                raise ValueError("Number of poses and scores must match")
            
            self.logger.info(f"Generating {self.report_format.upper()} report for {len(poses)} poses")
            
            if self.report_format == 'html':
                return self._generate_html_report(
                    protein, poses, scores, output_file, 
                    clustering_results, energy_decomposition, interaction_analysis
                )
            elif self.report_format == 'pdf':
                return self._generate_pdf_report(
                    protein, poses, scores, output_file, 
                    clustering_results, energy_decomposition, interaction_analysis
                )
            elif self.report_format == 'txt':
                return self._generate_txt_report(
                    protein, poses, scores, output_file, 
                    clustering_results, energy_decomposition, interaction_analysis
                )
            else:
                raise ValueError(f"Unsupported report format: {self.report_format}")
                
        except Exception as e:
            self.logger.error(f"Error generating report: {e}")
            raise
    
    def _generate_html_report(self, protein: Any, poses: List[Any], scores: List[float],
                             output_file: str, clustering_results: Optional[Dict],
                             energy_decomposition: Optional[Dict],
                             interaction_analysis: Optional[Dict]) -> str:
        """Generate HTML report."""
        try:
            # Create HTML content
            html_content = self._create_html_header()
            
            # Add report sections
            if 'summary' in self.include_sections:
                html_content.extend(self._create_summary_section_html(protein, poses, scores))
            
            if 'top_poses' in self.include_sections:
                html_content.extend(self._create_top_poses_section_html(poses, scores))
            
            if 'clusters' in self.include_sections and clustering_results:
                html_content.extend(self._create_clusters_section_html(clustering_results))
            
            if 'interactions' in self.include_sections:
                html_content.extend(self._create_interactions_section_html(
                    protein, poses, scores, interaction_analysis
                ))
            
            if 'energetics' in self.include_sections and energy_decomposition:
                html_content.extend(self._create_energetics_section_html(energy_decomposition))
            
            # Add footer
            html_content.extend(self._create_html_footer())
            
            # Write to file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write('\\n'.join(html_content))
            
            self.logger.info(f"HTML report generated: {output_file}")
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error generating HTML report: {e}")
            raise
    
    def _generate_txt_report(self, protein: Any, poses: List[Any], scores: List[float],
                            output_file: str, clustering_results: Optional[Dict],
                            energy_decomposition: Optional[Dict],
                            interaction_analysis: Optional[Dict]) -> str:
        """Generate text report."""
        try:
            lines = []
            
            # Header
            lines.append("=" * 80)
            lines.append("PANDADOCK DOCKING REPORT")
            lines.append("=" * 80)
            lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            lines.append("")
            
            # Summary section
            if 'summary' in self.include_sections:
                lines.extend(self._create_summary_section_txt(protein, poses, scores))
            
            # Top poses section
            if 'top_poses' in self.include_sections:
                lines.extend(self._create_top_poses_section_txt(poses, scores))
            
            # Clustering section
            if 'clusters' in self.include_sections and clustering_results:
                lines.extend(self._create_clusters_section_txt(clustering_results))
            
            # Interactions section
            if 'interactions' in self.include_sections:
                lines.extend(self._create_interactions_section_txt(
                    protein, poses, scores, interaction_analysis
                ))
            
            # Energetics section
            if 'energetics' in self.include_sections and energy_decomposition:
                lines.extend(self._create_energetics_section_txt(energy_decomposition))
            
            # Footer
            lines.append("")
            lines.append("=" * 80)
            lines.append("End of Report")
            lines.append("=" * 80)
            
            # Write to file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write('\\n'.join(lines))
            
            self.logger.info(f"Text report generated: {output_file}")
            return output_file
            
        except Exception as e:
            self.logger.error(f"Error generating text report: {e}")
            raise
    
    def _generate_pdf_report(self, protein: Any, poses: List[Any], scores: List[float],
                            output_file: str, clustering_results: Optional[Dict],
                            energy_decomposition: Optional[Dict],
                            interaction_analysis: Optional[Dict]) -> str:
        """Generate PDF report."""
        try:
            # For PDF generation, we'll create HTML first then convert
            html_file = output_file.replace('.pdf', '.html')
            
            # Generate HTML report
            self._generate_html_report(
                protein, poses, scores, html_file, 
                clustering_results, energy_decomposition, interaction_analysis
            )
            
            # Try to convert to PDF using weasyprint or other libraries
            try:
                import weasyprint
                weasyprint.HTML(html_file).write_pdf(output_file)
                os.remove(html_file)  # Clean up temporary HTML file
                self.logger.info(f"PDF report generated: {output_file}")
                return output_file
            except ImportError:
                self.logger.warning("weasyprint not available for PDF generation")
                # Keep HTML file as fallback
                self.logger.info(f"HTML report generated instead: {html_file}")
                return html_file
                
        except Exception as e:
            self.logger.error(f"Error generating PDF report: {e}")
            raise
    
    def _create_html_header(self) -> List[str]:
        """Create HTML header with CSS styling."""
        return [
            '<!DOCTYPE html>',
            '<html lang="en">',
            '<head>',
            '    <meta charset="UTF-8">',
            '    <meta name="viewport" content="width=device-width, initial-scale=1.0">',
            '    <title>PandaDock Docking Report</title>',
            '    <style>',
            '        body {',
            '            font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;',
            '            line-height: 1.6;',
            '            padding: 20px;',
            '            max-width: 1200px;',
            '            margin: 0 auto;',
            '            background-color: #f9f9f9;',
            '        }',
            '        .container {',
            '            background-color: white;',
            '            padding: 30px;',
            '            border-radius: 8px;',
            '            box-shadow: 0 2px 10px rgba(0,0,0,0.1);',
            '        }',
            '        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }',
            '        h2 { color: #34495e; border-bottom: 2px solid #ecf0f1; padding-bottom: 5px; margin-top: 30px; }',
            '        h3 { color: #2c3e50; }',
            '        table {',
            '            border-collapse: collapse;',
            '            width: 100%;',
            '            margin: 20px 0;',
            '            box-shadow: 0 1px 3px rgba(0,0,0,0.2);',
            '        }',
            '        th, td {',
            '            border: 1px solid #ddd;',
            '            padding: 12px 8px;',
            '            text-align: left;',
            '        }',
            '        th {',
            '            background-color: #3498db;',
            '            color: white;',
            '            font-weight: bold;',
            '        }',
            '        tr:nth-child(even) { background-color: #f8f9fa; }',
            '        tr:hover { background-color: #e8f4f8; }',
            '        .highlight { background-color: #fff3cd; }',
            '        .good-score { color: #27ae60; font-weight: bold; }',
            '        .poor-score { color: #e74c3c; }',
            '        .chart-container { margin: 20px 0; padding: 20px; background-color: #ecf0f1; border-radius: 5px; }',
            '        .section { margin-bottom: 40px; }',
            '        .footnote { font-size: 0.9em; color: #7f8c8d; margin-top: 40px; }',
            '        .stats-grid {',
            '            display: grid;',
            '            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));',
            '            gap: 20px;',
            '            margin: 20px 0;',
            '        }',
            '        .stat-box {',
            '            background-color: #ecf0f1;',
            '            padding: 15px;',
            '            border-radius: 5px;',
            '            text-align: center;',
            '        }',
            '        .stat-value { font-size: 1.5em; font-weight: bold; color: #2c3e50; }',
            '        .stat-label { color: #7f8c8d; }',
            '    </style>',
            '</head>',
            '<body>',
            '    <div class="container">',
            f'        <h1>🐼 PandaDock Docking Report</h1>',
            f'        <p><strong>Generated:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>'
        ]
    
    def _create_html_footer(self) -> List[str]:
        """Create HTML footer."""
        return [
            '        <div class="footnote">',
            '            <p>This report was generated by PandaDock molecular docking software.</p>',
            '            <p>For questions or support, please refer to the documentation.</p>',
            '        </div>',
            '    </div>',
            '</body>',
            '</html>'
        ]
    
    def _create_summary_section_html(self, protein: Any, poses: List[Any], 
                                   scores: List[float]) -> List[str]:
        """Create summary section for HTML report."""
        lines = [
            '        <div class="section">',
            '            <h2>📊 Summary</h2>',
            '            <div class="stats-grid">'
        ]
        
        # Basic statistics
        lines.extend([
            '                <div class="stat-box">',
            f'                    <div class="stat-value">{len(poses)}</div>',
            '                    <div class="stat-label">Total Poses</div>',
            '                </div>',
            '                <div class="stat-box">',
            f'                    <div class="stat-value">{min(scores):.3f}</div>',
            '                    <div class="stat-label">Best Score</div>',
            '                </div>',
            '                <div class="stat-box">',
            f'                    <div class="stat-value">{np.mean(scores):.3f}</div>',
            '                    <div class="stat-label">Average Score</div>',
            '                </div>',
            '                <div class="stat-box">',
            f'                    <div class="stat-value">{np.std(scores):.3f}</div>',
            '                    <div class="stat-label">Score StdDev</div>',
            '                </div>'
        ])
        
        lines.extend([
            '            </div>',
            '        </div>'
        ])
        
        return lines
    
    def _create_top_poses_section_html(self, poses: List[Any], scores: List[float]) -> List[str]:
        """Create top poses section for HTML report."""
        lines = [
            '        <div class="section">',
            '            <h2>🏆 Top Poses</h2>',
            '            <table>',
            '                <tr><th>Rank</th><th>Score</th><th>Status</th></tr>'
        ]
        
        # Sort poses by score and show top 15
        sorted_indices = np.argsort(scores)
        for i, idx in enumerate(sorted_indices[:15]):
            score = scores[idx]
            score_class = "good-score" if score < -5.0 else "poor-score" if score > 0 else ""
            status = "Excellent" if score < -8.0 else "Good" if score < -5.0 else "Fair" if score < 0 else "Poor"
            
            lines.append(f'                <tr><td>{i+1}</td><td class="{score_class}">{score:.4f}</td><td>{status}</td></tr>')
        
        lines.extend([
            '            </table>',
            '        </div>'
        ])
        
        return lines
    
    def _create_clusters_section_html(self, clustering_results: Dict) -> List[str]:
        """Create clusters section for HTML report."""
        lines = [
            '        <div class="section">',
            '            <h2>🔗 Clustering Analysis</h2>'
        ]
        
        if 'clusters' in clustering_results and clustering_results['clusters']:
            rmsd_cutoff = clustering_results.get('rmsd_cutoff', 2.0)
            lines.append(f'            <p>Found <strong>{len(clustering_results["clusters"])}</strong> clusters using <strong>{rmsd_cutoff} Å</strong> RMSD cutoff.</p>')
            
            lines.extend([
                '            <table>',
                '                <tr><th>Cluster</th><th>Size</th><th>Best Score</th><th>Avg Score</th><th>Representative</th></tr>'
            ])
            
            for i, cluster in enumerate(clustering_results['clusters']):
                best_score = cluster.get('best_score', 0.0)
                avg_score = cluster.get('avg_score', 0.0)
                size = cluster.get('size', 0)
                representative = cluster.get('representative', 0)
                
                lines.append(f'                <tr><td>{i+1}</td><td>{size}</td><td>{best_score:.4f}</td><td>{avg_score:.4f}</td><td>Pose #{representative+1}</td></tr>')
            
            lines.append('            </table>')
        else:
            lines.append('            <p>No clustering results available.</p>')
        
        lines.append('        </div>')
        return lines
    
    def _create_interactions_section_html(self, protein: Any, poses: List[Any], 
                                        scores: List[float], interaction_analysis: Optional[Dict]) -> List[str]:
        """Create interactions section for HTML report."""
        lines = [
            '        <div class="section">',
            '            <h2>🔬 Protein-Ligand Interactions</h2>'
        ]
        
        if interaction_analysis:
            lines.append('            <p>Detailed interaction analysis provided.</p>')
            # Add interaction details here if available
        else:
            lines.append('            <p>Interaction analysis for best pose:</p>')
            lines.append('            <div class="chart-container">')
            lines.append('                <p>💡 <em>Interaction analysis would be displayed here with detailed binding mode information.</em></p>')
            lines.append('            </div>')
        
        lines.append('        </div>')
        return lines
    
    def _create_energetics_section_html(self, energy_decomposition: Dict) -> List[str]:
        """Create energetics section for HTML report."""
        lines = [
            '        <div class="section">',
            '            <h2>⚡ Energy Analysis</h2>'
        ]
        
        if energy_decomposition:
            total_energy = energy_decomposition.get('total', 0.0)
            lines.append(f'            <p><strong>Total Binding Energy:</strong> {total_energy:.3f} kcal/mol</p>')
            
            # Create energy components table
            components = {k: v for k, v in energy_decomposition.items() if k != 'total'}
            if components:
                lines.extend([
                    '            <table>',
                    '                <tr><th>Energy Component</th><th>Value (kcal/mol)</th><th>Contribution</th></tr>'
                ])
                
                for component, value in sorted(components.items(), key=lambda x: abs(x[1]), reverse=True):
                    contribution = "Favorable" if value < 0 else "Unfavorable"
                    color_class = "good-score" if value < 0 else "poor-score"
                    lines.append(f'                <tr><td>{component.title()}</td><td class="{color_class}">{value:.3f}</td><td>{contribution}</td></tr>')
                
                lines.append('            </table>')
        else:
            lines.append('            <p>No energy decomposition available.</p>')
        
        lines.append('        </div>')
        return lines
    
    def _create_summary_section_txt(self, protein: Any, poses: List[Any], scores: List[float]) -> List[str]:
        """Create summary section for text report."""
        return [
            "SUMMARY",
            "-" * 20,
            f"Number of Poses: {len(poses)}",
            f"Best Score: {min(scores):.4f}",
            f"Average Score: {np.mean(scores):.4f}",
            f"Score Std Dev: {np.std(scores):.4f}",
            f"Score Range: {min(scores):.4f} to {max(scores):.4f}",
            ""
        ]
    
    def _create_top_poses_section_txt(self, poses: List[Any], scores: List[float]) -> List[str]:
        """Create top poses section for text report."""
        lines = [
            "TOP POSES",
            "-" * 20,
            f"{'Rank':<6} {'Score':<12} {'Status':<10}",
            "-" * 30
        ]
        
        sorted_indices = np.argsort(scores)
        for i, idx in enumerate(sorted_indices[:15]):
            score = scores[idx]
            status = "Excellent" if score < -8.0 else "Good" if score < -5.0 else "Fair" if score < 0 else "Poor"
            lines.append(f"{i+1:<6} {score:<12.4f} {status:<10}")
        
        lines.append("")
        return lines
    
    def _create_clusters_section_txt(self, clustering_results: Dict) -> List[str]:
        """Create clusters section for text report."""
        lines = [
            "CLUSTERING ANALYSIS",
            "-" * 20
        ]
        
        if 'clusters' in clustering_results and clustering_results['clusters']:
            rmsd_cutoff = clustering_results.get('rmsd_cutoff', 2.0)
            lines.append(f"Found {len(clustering_results['clusters'])} clusters using {rmsd_cutoff} Å RMSD cutoff.")
            lines.append("")
            lines.append(f"{'Cluster':<8} {'Size':<6} {'Best Score':<12} {'Avg Score':<12}")
            lines.append("-" * 40)
            
            for i, cluster in enumerate(clustering_results['clusters']):
                best_score = cluster.get('best_score', 0.0)
                avg_score = cluster.get('avg_score', 0.0)
                size = cluster.get('size', 0)
                lines.append(f"{i+1:<8} {size:<6} {best_score:<12.4f} {avg_score:<12.4f}")
        else:
            lines.append("No clustering results available.")
        
        lines.append("")
        return lines
    
    def _create_interactions_section_txt(self, protein: Any, poses: List[Any], 
                                       scores: List[float], interaction_analysis: Optional[Dict]) -> List[str]:
        """Create interactions section for text report."""
        lines = [
            "PROTEIN-LIGAND INTERACTIONS",
            "-" * 30
        ]
        
        if interaction_analysis:
            lines.append("Detailed interaction analysis available.")
        else:
            lines.append("Basic interaction analysis for best pose.")
            lines.append("(Detailed analysis requires interaction fingerprinting)")
        
        lines.append("")
        return lines
    
    def _create_energetics_section_txt(self, energy_decomposition: Dict) -> List[str]:
        """Create energetics section for text report."""
        lines = [
            "ENERGY ANALYSIS",
            "-" * 20
        ]
        
        if energy_decomposition:
            total_energy = energy_decomposition.get('total', 0.0)
            lines.append(f"Total Binding Energy: {total_energy:.3f} kcal/mol")
            lines.append("")
            
            components = {k: v for k, v in energy_decomposition.items() if k != 'total'}
            if components:
                lines.append(f"{'Component':<15} {'Value':<12} {'Type':<12}")
                lines.append("-" * 40)
                
                for component, value in sorted(components.items(), key=lambda x: abs(x[1]), reverse=True):
                    contribution = "Favorable" if value < 0 else "Unfavorable"
                    lines.append(f"{component.title():<15} {value:<12.3f} {contribution:<12}")
        else:
            lines.append("No energy decomposition available.")
        
        lines.append("")
        return lines