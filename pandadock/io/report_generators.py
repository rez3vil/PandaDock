"""
Report generation utilities for PandaDock.
"""

import logging
from typing import Dict, Any
from pathlib import Path


class ReportGenerators:
    """Handles generation of various report formats."""
    
    def __init__(self, output_dir: str):
        """
        Initialize report generators.
        
        Args:
            output_dir: Output directory for reports
        """
        self.output_dir = Path(output_dir)
        self.logger = logging.getLogger(__name__)
    
    def generate_html_report(self, results: Dict[str, Any]) -> str:
        """
        Generate HTML report.
        
        Args:
            results: Docking results
            
        Returns:
            Path to generated HTML report
        """
        
        html_file = self.output_dir / "docking_report.html"
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PandaDock Results</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                .header {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
                .stats {{ background-color: #ecf0f1; padding: 15px; margin: 20px 0; border-radius: 5px; }}
                .pose {{ margin: 10px 0; padding: 10px; background-color: #f8f9fa; border-left: 4px solid #007bff; }}
            </style>
        </head>
        <body>
            <h1 class="header">🐼 PandaDock Results Report</h1>
            
            <div class="stats">
                <h2>Summary Statistics</h2>
                <p><strong>Total Poses:</strong> {results.get('statistics', {}).get('total_poses', 0)}</p>
                <p><strong>Best Score:</strong> {results.get('statistics', {}).get('best_score', 'N/A')}</p>
            </div>
            
            <h2>Top Poses</h2>
        """
        
        poses = results.get('poses', [])
        for pose in poses[:10]:
            html_content += f"""
            <div class="pose">
                <strong>Rank {pose['rank']}</strong> - Score: {pose['score']:.4f}
                <br>Atoms: {pose.get('n_atoms', 'N/A')}
            </div>
            """
        
        html_content += """
        </body>
        </html>
        """
        
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"Generated HTML report: {html_file}")
        return str(html_file)