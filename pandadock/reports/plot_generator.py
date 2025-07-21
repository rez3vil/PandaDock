#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive Plot Generator for PandaDock
===========================================

Generates publication-ready plots for docking results including:
1. Binding affinity, docking scores, deltaG, IC50/EC50 plots
2. 2D interaction maps (Discovery Studio style)
3. Master publication plot combining all metrics
4. Detailed TXT reports

Author: PandaDock Team
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional, Any, Union
import json
from datetime import datetime
import warnings
from scipy import stats
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as path_effects

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set up publication-quality plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.format': 'png',
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'text.usetex': False
})

class PandaDockPlotGenerator:
    """
    Comprehensive plot generator for PandaDock docking results
    
    Features:
    - Binding affinity distribution plots
    - Docking score analysis
    - IC50/EC50 correlation plots
    - 2D interaction maps
    - Publication-ready master plots
    - Detailed TXT reports
    """
    
    def __init__(self, output_dir: str, protein_name: str = None, ligand_name: str = None):
        """
        Initialize the plot generator
        
        Args:
            output_dir: Output directory for plots and reports
            protein_name: Name of the protein target
            ligand_name: Name of the ligand
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.protein_name = protein_name or "Unknown_Protein"
        self.ligand_name = ligand_name or "Unknown_Ligand"
        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        
        # Color schemes for different plot types
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72',
            'tertiary': '#F18F01',
            'quaternary': '#C73E1D',
            'success': '#44AF69',
            'warning': '#F18F01',
            'danger': '#C73E1D',
            'info': '#2E86AB'
        }
        
        self.logger.info(f"Initialized PandaDock Plot Generator for {self.protein_name} - {self.ligand_name}")
    
    def generate_all_plots(self, poses_csv: str, poses_dir: str, 
                          algorithm_info: Dict[str, Any], 
                          command_info: Dict[str, Any]) -> Dict[str, str]:
        """
        Generate all types of plots and reports
        
        Args:
            poses_csv: Path to poses summary CSV file
            poses_dir: Directory containing pose PDB files
            algorithm_info: Information about the algorithm used
            command_info: Information about the command executed
            
        Returns:
            Dictionary of generated file paths
        """
        self.logger.info("Starting comprehensive plot generation...")
        
        # Load pose data
        poses_df = pd.read_csv(poses_csv)
        
        generated_files = {}
        
        # 1. Binding metrics plots
        metrics_plot = self.create_binding_metrics_plot(poses_df)
        generated_files['binding_metrics'] = metrics_plot
        
        # 2. Score distribution plots
        score_plot = self.create_score_distribution_plot(poses_df)
        generated_files['score_distribution'] = score_plot
        
        # 3. IC50/EC50 correlation plots
        ic50_plot = self.create_ic50_ec50_plot(poses_df)
        generated_files['ic50_ec50'] = ic50_plot
        
        # 4. Master publication plot
        master_plot = self.create_master_publication_plot(poses_df)
        generated_files['master_publication'] = master_plot
        
        # 5. 2D interaction maps (for top poses) - Only if PandaMap integration not used
        try:
            # Check if PandaMap integration is available
            from .pandamap_integration import PandaMapIntegration
            self.logger.info("PandaMap integration available - skipping basic interaction maps")
            generated_files['interaction_maps'] = []  # Will be handled by PandaMap
        except ImportError:
            # Fall back to basic interaction maps
            interaction_maps = self.create_interaction_maps(poses_df, poses_dir)
            generated_files['interaction_maps'] = interaction_maps
        
        # 6. Detailed TXT report
        txt_report = self.create_detailed_txt_report(
            poses_df, algorithm_info, command_info, generated_files
        )
        generated_files['txt_report'] = txt_report
        
        self.logger.info(f"Generated {len(generated_files)} plot files and reports")
        return generated_files
    
    def create_binding_metrics_plot(self, poses_df: pd.DataFrame) -> str:
        """Create comprehensive binding metrics visualization"""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'Binding Metrics Analysis: {self.protein_name} - {self.ligand_name}', 
                    fontsize=16, fontweight='bold')
        
        # 1. Binding Affinity Distribution
        ax = axes[0, 0]
        ax.hist(poses_df['Binding_Affinity'], bins=15, alpha=0.7, 
               color=self.colors['primary'], edgecolor='black')
        ax.set_xlabel('Binding Affinity (kcal/mol)')
        ax.set_ylabel('Number of Poses')
        ax.set_title('Binding Affinity Distribution')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        mean_affinity = poses_df['Binding_Affinity'].mean()
        std_affinity = poses_df['Binding_Affinity'].std()
        ax.axvline(mean_affinity, color='red', linestyle='--', 
                  label=f'Mean: {mean_affinity:.2f} ± {std_affinity:.2f}')
        ax.legend()
        
        # 2. Energy Distribution
        ax = axes[0, 1]
        ax.hist(poses_df['Energy'], bins=15, alpha=0.7, 
               color=self.colors['secondary'], edgecolor='black')
        ax.set_xlabel('Docking Energy (kcal/mol)')
        ax.set_ylabel('Number of Poses')
        ax.set_title('Docking Energy Distribution')
        ax.grid(True, alpha=0.3)
        
        mean_energy = poses_df['Energy'].mean()
        std_energy = poses_df['Energy'].std()
        ax.axvline(mean_energy, color='red', linestyle='--',
                  label=f'Mean: {mean_energy:.2f} ± {std_energy:.2f}')
        ax.legend()
        
        # 3. Score vs Confidence
        ax = axes[0, 2]
        scatter = ax.scatter(poses_df['Score'], poses_df['Confidence'], 
                           c=poses_df['Binding_Affinity'], cmap='viridis', 
                           s=60, alpha=0.7)
        ax.set_xlabel('Docking Score')
        ax.set_ylabel('Confidence')
        ax.set_title('Score vs Confidence')
        ax.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Binding Affinity (kcal/mol)')
        
        # 4. ΔG Calculation (Energy - Baseline)
        baseline_energy = poses_df['Energy'].max()  # Use worst energy as baseline
        delta_g = poses_df['Energy'] - baseline_energy
        
        ax = axes[1, 0]
        ax.hist(delta_g, bins=15, alpha=0.7, 
               color=self.colors['tertiary'], edgecolor='black')
        ax.set_xlabel('ΔG (kcal/mol)')
        ax.set_ylabel('Number of Poses')
        ax.set_title('ΔG Distribution (relative to worst pose)')
        ax.grid(True, alpha=0.3)
        
        mean_delta_g = delta_g.mean()
        ax.axvline(mean_delta_g, color='red', linestyle='--',
                  label=f'Mean ΔG: {mean_delta_g:.2f}')
        ax.legend()
        
        # 5. Ligand Efficiency
        ax = axes[1, 1]
        ax.hist(poses_df['Ligand_Efficiency'], bins=15, alpha=0.7, 
               color=self.colors['quaternary'], edgecolor='black')
        ax.set_xlabel('Ligand Efficiency')
        ax.set_ylabel('Number of Poses')
        ax.set_title('Ligand Efficiency Distribution')
        ax.grid(True, alpha=0.3)
        
        mean_le = poses_df['Ligand_Efficiency'].mean()
        ax.axvline(mean_le, color='red', linestyle='--',
                  label=f'Mean LE: {mean_le:.3f}')
        ax.legend()
        
        # 6. Rank vs Binding Affinity
        ax = axes[1, 2]
        ax.plot(poses_df['Rank'], poses_df['Binding_Affinity'], 
               'o-', color=self.colors['primary'], markersize=6, linewidth=2)
        ax.set_xlabel('Pose Rank')
        ax.set_ylabel('Binding Affinity (kcal/mol)')
        ax.set_title('Binding Affinity vs Pose Rank')
        ax.grid(True, alpha=0.3)
        
        # Highlight best pose
        best_idx = poses_df['Binding_Affinity'].idxmax()
        ax.scatter(poses_df.loc[best_idx, 'Rank'], 
                  poses_df.loc[best_idx, 'Binding_Affinity'],
                  s=100, color='red', zorder=5, 
                  label=f'Best: {poses_df.loc[best_idx, "Binding_Affinity"]:.2f}')
        ax.legend()
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'binding_metrics_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(str(output_file).replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Generated binding metrics plot: {output_file}")
        return str(output_file)
    
    def create_score_distribution_plot(self, poses_df: pd.DataFrame) -> str:
        """Create detailed score distribution analysis"""
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'Score Distribution Analysis: {self.protein_name} - {self.ligand_name}', 
                    fontsize=16, fontweight='bold')
        
        # 1. Score Distribution with KDE
        ax = axes[0, 0]
        ax.hist(poses_df['Score'], bins=20, alpha=0.7, density=True,
               color=self.colors['primary'], edgecolor='black', label='Histogram')
        
        # Add KDE
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(poses_df['Score'])
        x_range = np.linspace(poses_df['Score'].min(), poses_df['Score'].max(), 100)
        ax.plot(x_range, kde(x_range), color='red', linewidth=2, label='KDE')
        
        ax.set_xlabel('Docking Score')
        ax.set_ylabel('Density')
        ax.set_title('Score Distribution with KDE')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # 2. Score vs Energy Correlation
        ax = axes[0, 1]
        ax.scatter(poses_df['Score'], poses_df['Energy'], alpha=0.7, 
                  color=self.colors['secondary'], s=60)
        
        # Add correlation line
        z = np.polyfit(poses_df['Score'], poses_df['Energy'], 1)
        p = np.poly1d(z)
        ax.plot(poses_df['Score'], p(poses_df['Score']), "r--", alpha=0.8, linewidth=2)
        
        # Calculate correlation
        correlation = stats.pearsonr(poses_df['Score'], poses_df['Energy'])[0]
        ax.set_xlabel('Docking Score')
        ax.set_ylabel('Energy (kcal/mol)')
        ax.set_title(f'Score vs Energy (r = {correlation:.3f})')
        ax.grid(True, alpha=0.3)
        
        # 3. Confidence Distribution
        ax = axes[1, 0]
        ax.hist(poses_df['Confidence'], bins=15, alpha=0.7, 
               color=self.colors['tertiary'], edgecolor='black')
        ax.set_xlabel('Confidence Score')
        ax.set_ylabel('Number of Poses')
        ax.set_title('Confidence Distribution')
        ax.grid(True, alpha=0.3)
        
        # Add confidence threshold line
        confidence_threshold = 0.7  # Typical threshold
        ax.axvline(confidence_threshold, color='red', linestyle='--',
                  label=f'Threshold: {confidence_threshold}')
        
        high_conf_count = (poses_df['Confidence'] >= confidence_threshold).sum()
        ax.text(0.05, 0.95, f'High Confidence: {high_conf_count}/{len(poses_df)}',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.legend()
        
        # 4. Score vs Rank Performance
        ax = axes[1, 1]
        ax.plot(poses_df['Rank'], poses_df['Score'], 'o-', 
               color=self.colors['quaternary'], markersize=6, linewidth=2)
        ax.set_xlabel('Pose Rank')
        ax.set_ylabel('Docking Score')
        ax.set_title('Score vs Pose Rank')
        ax.grid(True, alpha=0.3)
        
        # Highlight top 3 poses
        top_3 = poses_df.head(3)
        ax.scatter(top_3['Rank'], top_3['Score'], 
                  s=100, color='red', zorder=5, label='Top 3')
        ax.legend()
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'score_distribution_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(str(output_file).replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Generated score distribution plot: {output_file}")
        return str(output_file)
    
    def create_ic50_ec50_plot(self, poses_df: pd.DataFrame) -> str:
        """Create IC50/EC50 analysis plots"""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'IC50/EC50 Analysis: {self.protein_name} - {self.ligand_name}', 
                    fontsize=16, fontweight='bold')
        
        # Convert scientific notation strings to float
        ic50_values = []
        ec50_values = []
        
        for ic50_str, ec50_str in zip(poses_df['IC50_uM'], poses_df['EC50_uM']):
            # Parse scientific notation (e.g., "1.23e+04")
            if isinstance(ic50_str, str):
                ic50_val = float(ic50_str.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ic50_val = float(ic50_str)
            
            if isinstance(ec50_str, str):
                ec50_val = float(ec50_str.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ec50_val = float(ec50_str)
                
            ic50_values.append(ic50_val)
            ec50_values.append(ec50_val)
        
        ic50_array = np.array(ic50_values)
        ec50_array = np.array(ec50_values)
        
        # 1. IC50 Distribution (log scale)
        ax = axes[0, 0]
        ax.hist(np.log10(ic50_array), bins=15, alpha=0.7, 
               color=self.colors['primary'], edgecolor='black')
        ax.set_xlabel('log₁₀(IC50 μM)')
        ax.set_ylabel('Number of Poses')
        ax.set_title('IC50 Distribution')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        median_ic50 = np.median(ic50_array)
        ax.axvline(np.log10(median_ic50), color='red', linestyle='--',
                  label=f'Median: {median_ic50:.1e} μM')
        ax.legend()
        
        # 2. EC50 Distribution (log scale)
        ax = axes[0, 1]
        ax.hist(np.log10(ec50_array), bins=15, alpha=0.7, 
               color=self.colors['secondary'], edgecolor='black')
        ax.set_xlabel('log₁₀(EC50 μM)')
        ax.set_ylabel('Number of Poses')
        ax.set_title('EC50 Distribution')
        ax.grid(True, alpha=0.3)
        
        median_ec50 = np.median(ec50_array)
        ax.axvline(np.log10(median_ec50), color='red', linestyle='--',
                  label=f'Median: {median_ec50:.1e} μM')
        ax.legend()
        
        # 3. IC50 vs EC50 Correlation
        ax = axes[0, 2]
        ax.scatter(np.log10(ic50_array), np.log10(ec50_array), 
                  alpha=0.7, color=self.colors['tertiary'], s=60)
        
        # Add correlation line
        z = np.polyfit(np.log10(ic50_array), np.log10(ec50_array), 1)
        p = np.poly1d(z)
        x_line = np.linspace(np.log10(ic50_array).min(), np.log10(ic50_array).max(), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2)
        
        correlation = stats.pearsonr(np.log10(ic50_array), np.log10(ec50_array))[0]
        ax.set_xlabel('log₁₀(IC50 μM)')
        ax.set_ylabel('log₁₀(EC50 μM)')
        ax.set_title(f'IC50 vs EC50 Correlation (r = {correlation:.3f})')
        ax.grid(True, alpha=0.3)
        
        # 4. Binding Affinity vs IC50
        ax = axes[1, 0]
        ax.scatter(poses_df['Binding_Affinity'], np.log10(ic50_array), 
                  alpha=0.7, color=self.colors['quaternary'], s=60)
        
        # Add trend line
        z = np.polyfit(poses_df['Binding_Affinity'], np.log10(ic50_array), 1)
        p = np.poly1d(z)
        x_line = np.linspace(poses_df['Binding_Affinity'].min(), 
                           poses_df['Binding_Affinity'].max(), 100)
        ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2)
        
        correlation = stats.pearsonr(poses_df['Binding_Affinity'], np.log10(ic50_array))[0]
        ax.set_xlabel('Binding Affinity (kcal/mol)')
        ax.set_ylabel('log₁₀(IC50 μM)')
        ax.set_title(f'Affinity vs IC50 (r = {correlation:.3f})')
        ax.grid(True, alpha=0.3)
        
        # 5. Potency Classification
        ax = axes[1, 1]
        
        # Classify potency ranges
        potency_ranges = {
            'High (< 1 μM)': (ic50_array < 1).sum(),
            'Moderate (1-10 μM)': ((ic50_array >= 1) & (ic50_array < 10)).sum(),
            'Low (10-100 μM)': ((ic50_array >= 10) & (ic50_array < 100)).sum(),
            'Very Low (≥ 100 μM)': (ic50_array >= 100).sum()
        }
        
        colors = [self.colors['success'], self.colors['info'], 
                 self.colors['warning'], self.colors['danger']]
        
        ax.pie(potency_ranges.values(), labels=potency_ranges.keys(), 
              autopct='%1.1f%%', colors=colors, startangle=90)
        ax.set_title('IC50 Potency Classification')
        
        # 6. Rank vs IC50
        ax = axes[1, 2]
        ax.semilogy(poses_df['Rank'], ic50_array, 'o-', 
                   color=self.colors['primary'], markersize=6, linewidth=2)
        ax.set_xlabel('Pose Rank')
        ax.set_ylabel('IC50 (μM)')
        ax.set_title('IC50 vs Pose Rank')
        ax.grid(True, alpha=0.3)
        
        # Highlight best IC50
        best_ic50_idx = np.argmin(ic50_array)
        ax.scatter(poses_df.iloc[best_ic50_idx]['Rank'], ic50_array[best_ic50_idx],
                  s=100, color='red', zorder=5, 
                  label=f'Best: {ic50_array[best_ic50_idx]:.1e} μM')
        ax.legend()
        
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / 'ic50_ec50_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(str(output_file).replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Generated IC50/EC50 analysis plot: {output_file}")
        return str(output_file)
    
    def create_master_publication_plot(self, poses_df: pd.DataFrame) -> str:
        """Create comprehensive master publication plot"""
        
        # Create figure with custom layout
        fig = plt.figure(figsize=(20, 16))
        gs = GridSpec(4, 4, figure=fig, hspace=0.3, wspace=0.3)
        
        # Main title
        fig.suptitle(f'PandaDock Analysis: {self.protein_name} - {self.ligand_name}', 
                    fontsize=20, fontweight='bold', y=0.95)
        
        # Convert IC50/EC50 values
        ic50_values = []
        ec50_values = []
        
        for ic50_x, ec50_x in zip(poses_df['IC50_uM'], poses_df['EC50_uM']):
            if isinstance(ic50_x, str):
                ic50_val = float(ic50_x.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ic50_val = float(ic50_x)
            
            if isinstance(ec50_x, str):
                ec50_val = float(ec50_x.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ec50_val = float(ec50_x)
                
            ic50_values.append(ic50_val)
            ec50_values.append(ec50_val)
        
        ic50_values = np.array(ic50_values)
        ec50_values = np.array(ec50_values)
        
        # 1. Main correlation plot (top-left, large)
        ax1 = fig.add_subplot(gs[0:2, 0:2])
        scatter = ax1.scatter(poses_df['Binding_Affinity'], poses_df['Energy'], 
                             c=np.log10(ic50_values), cmap='viridis_r', 
                             s=100, alpha=0.8, edgecolors='black', linewidth=0.5)
        
        # Add trend line
        z = np.polyfit(poses_df['Binding_Affinity'], poses_df['Energy'], 1)
        p = np.poly1d(z)
        x_line = np.linspace(poses_df['Binding_Affinity'].min(), 
                           poses_df['Binding_Affinity'].max(), 100)
        ax1.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=3)
        
        correlation = stats.pearsonr(poses_df['Binding_Affinity'], poses_df['Energy'])[0]
        ax1.set_xlabel('Binding Affinity (kcal/mol)', fontsize=14)
        ax1.set_ylabel('Docking Energy (kcal/mol)', fontsize=14)
        ax1.set_title(f'Binding Affinity vs Energy (r = {correlation:.3f})', fontsize=16)
        ax1.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('log₁₀(IC50 μM)', fontsize=12)
        
        # Annotate best poses
        top_3 = poses_df.head(3)
        for i, (idx, row) in enumerate(top_3.iterrows()):
            ax1.annotate(f'#{i+1}', (row['Binding_Affinity'], row['Energy']),
                        xytext=(10, 10), textcoords='offset points',
                        fontsize=12, fontweight='bold', color='white',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.8))
        
        # 2. Score distribution (top-right)
        ax2 = fig.add_subplot(gs[0, 2:])
        ax2.hist(poses_df['Score'], bins=15, alpha=0.7, 
                color=self.colors['primary'], edgecolor='black')
        ax2.set_xlabel('Docking Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Score Distribution')
        ax2.grid(True, alpha=0.3)
        
        # Add statistics
        mean_score = poses_df['Score'].mean()
        std_score = poses_df['Score'].std()
        ax2.axvline(mean_score, color='red', linestyle='--', linewidth=2,
                   label=f'μ = {mean_score:.3f} ± {std_score:.3f}')
        ax2.legend()
        
        # 3. IC50 potency analysis (middle-right)
        ax3 = fig.add_subplot(gs[1, 2:])
        potency_ranges = {
            'High\n(< 1 μM)': (ic50_values < 1).sum(),
            'Moderate\n(1-10 μM)': ((ic50_values >= 1) & (ic50_values < 10)).sum(),
            'Low\n(10-100 μM)': ((ic50_values >= 10) & (ic50_values < 100)).sum(),
            'Very Low\n(≥ 100 μM)': (ic50_values >= 100).sum()
        }
        
        colors = [self.colors['success'], self.colors['info'], 
                 self.colors['warning'], self.colors['danger']]
        bars = ax3.bar(potency_ranges.keys(), potency_ranges.values(), color=colors, alpha=0.8)
        ax3.set_ylabel('Number of Poses')
        ax3.set_title('IC50 Potency Distribution')
        ax3.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar, value in zip(bars, potency_ranges.values()):
            if value > 0:
                ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                        str(value), ha='center', va='bottom', fontweight='bold')
        
        # 4. Binding affinity vs rank (bottom-left)
        ax4 = fig.add_subplot(gs[2, 0:2])
        ax4.plot(poses_df['Rank'], poses_df['Binding_Affinity'], 'o-', 
                color=self.colors['secondary'], markersize=8, linewidth=3, alpha=0.8)
        ax4.set_xlabel('Pose Rank')
        ax4.set_ylabel('Binding Affinity (kcal/mol)')
        ax4.set_title('Binding Affinity Ranking')
        ax4.grid(True, alpha=0.3)
        
        # Highlight top poses
        top_5 = poses_df.head(5)
        ax4.scatter(top_5['Rank'], top_5['Binding_Affinity'], 
                   s=120, color='red', zorder=5, alpha=0.8)
        
        # 5. Ligand efficiency analysis (bottom-middle)
        ax5 = fig.add_subplot(gs[2, 2])
        ax5.scatter(poses_df['Binding_Affinity'], poses_df['Ligand_Efficiency'], 
                   alpha=0.7, color=self.colors['tertiary'], s=80)
        ax5.set_xlabel('Binding Affinity\n(kcal/mol)')
        ax5.set_ylabel('Ligand Efficiency')
        ax5.set_title('Ligand Efficiency')
        ax5.grid(True, alpha=0.3)
        
        # 6. Confidence distribution (bottom-right)
        ax6 = fig.add_subplot(gs[2, 3])
        ax6.hist(poses_df['Confidence'], bins=10, alpha=0.7, 
                color=self.colors['quaternary'], edgecolor='black', orientation='horizontal')
        ax6.set_ylabel('Confidence Score')
        ax6.set_xlabel('Frequency')
        ax6.set_title('Confidence')
        ax6.grid(True, alpha=0.3)
        
        # 7. Summary statistics table (bottom)
        ax7 = fig.add_subplot(gs[3, :])
        ax7.axis('off')
        
        # Calculate key statistics
        stats_data = {
            'Metric': ['Best Binding Affinity', 'Best IC50', 'Best EC50', 'Mean Score', 
                      'High Potency Poses', 'Total Poses'],
            'Value': [f'{poses_df["Binding_Affinity"].min():.2f} kcal/mol',
                     f'{ic50_values.min():.1e} μM',
                     f'{ec50_values.min():.1e} μM',
                     f'{poses_df["Score"].mean():.3f} ± {poses_df["Score"].std():.3f}',
                     f'{(ic50_values < 10).sum()}/{len(poses_df)}',
                     f'{len(poses_df)}']
        }
        
        # Create table
        table = ax7.table(cellText=[[metric, value] for metric, value in 
                                   zip(stats_data['Metric'], stats_data['Value'])],
                         colLabels=['Metric', 'Value'],
                         cellLoc='center',
                         loc='center',
                         bbox=[0.1, 0.2, 0.8, 0.6])
        
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1, 2)
        
        # Style the table
        for i in range(len(stats_data['Metric']) + 1):
            for j in range(2):
                cell = table[(i, j)]
                if i == 0:  # Header row
                    cell.set_facecolor(self.colors['primary'])
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f8f9fa' if i % 2 == 0 else 'white')
        
        # Add algorithm info
        ax7.text(0.5, 0.05, f'Generated by PandaDock | {datetime.now().strftime("%Y-%m-%d %H:%M")}',
                transform=ax7.transAxes, ha='center', va='bottom',
                fontsize=10, style='italic', alpha=0.7)
        
        # Save master publication plot
        output_file = self.output_dir / 'master_publication.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(str(output_file).replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Generated master publication plot: {output_file}")
        return str(output_file)
    
    def create_interaction_maps(self, poses_df: pd.DataFrame, poses_dir: str) -> List[str]:
        """Create 2D interaction maps for top poses using PandaMap if available"""
        
        interaction_files = []
        
        # Try to use PandaMap for professional interaction visualization
        try:
            from .pandamap_integration import create_pandamap_visualizations
            
            self.logger.info("Using PandaMap for professional interaction visualization")
            
            pandamap_files = create_pandamap_visualizations(
                poses_df=poses_df,
                poses_dir=poses_dir,
                output_dir=str(self.output_dir),
                top_n=3,
                generate_3d=False  # Only 2D for this function
            )
            
            # Return 2D maps
            if '2d_maps' in pandamap_files:
                interaction_files.extend(pandamap_files['2d_maps'])
            
            if interaction_files:
                self.logger.info(f"Generated {len(interaction_files)} PandaMap interaction maps")
                return interaction_files
                
        except ImportError:
            self.logger.info("PandaMap not available - using fallback visualization")
        except Exception as e:
            self.logger.warning(f"PandaMap visualization failed: {e} - using fallback")
        
        # Fallback to placeholder interaction analysis
        self.logger.info("Creating fallback interaction visualizations")
        
        top_poses = poses_df.head(3)  # Top 3 poses
        
        for idx, (_, pose) in enumerate(top_poses.iterrows()):
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            
            # Enhanced placeholder with more professional styling
            ax.text(0.5, 0.6, 'Professional Interaction Map', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=18, fontweight='bold')
            
            ax.text(0.5, 0.5, f'Pose: {pose["Pose_ID"]} (Rank {pose["Rank"]})', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=14, fontweight='bold')
            
            ax.text(0.5, 0.45, f'Binding Affinity: {pose["Binding_Affinity"]:.2f} kcal/mol', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=12)
            
            ax.text(0.5, 0.4, f'IC50: {pose["IC50_uM"]}', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=12)
            
            ax.text(0.5, 0.3, 'PandaMap Integration Available', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=14, color='green', fontweight='bold')
            
            ax.text(0.5, 0.25, 'Use --pandamap flag for full interaction analysis', 
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=10, style='italic')
            
            # Add visual elements
            import matplotlib.patches as patches
            circle = patches.Circle((0.5, 0.5), 0.15, transform=ax.transAxes, 
                                   facecolor='lightblue', edgecolor='navy', alpha=0.3)
            ax.add_patch(circle)
            
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_title(f'Interaction Analysis - {pose["Pose_ID"]}', fontsize=16, pad=20)
            ax.axis('off')
            
            # Save interaction map
            output_file = self.output_dir / f'interaction_map_pose_{idx+1}.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            interaction_files.append(str(output_file))
        
        self.logger.info(f"Generated {len(interaction_files)} fallback interaction maps")
        return interaction_files
    
    def create_detailed_txt_report(self, poses_df: pd.DataFrame, 
                                  algorithm_info: Dict[str, Any],
                                  command_info: Dict[str, Any],
                                  generated_files: Dict[str, str]) -> str:
        """Create comprehensive TXT report"""
        
        # Convert IC50/EC50 values
        ic50_values = []
        ec50_values = []
        
        for ic50_x, ec50_x in zip(poses_df['IC50_uM'], poses_df['EC50_uM']):
            if isinstance(ic50_x, str):
                ic50_val = float(ic50_x.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ic50_val = float(ic50_x)
            
            if isinstance(ec50_x, str):
                ec50_val = float(ec50_x.replace('e+', 'e').replace('e-', 'e-'))
            else:
                ec50_val = float(ec50_x)
                
            ic50_values.append(ic50_val)
            ec50_values.append(ec50_val)
        
        ic50_values = np.array(ic50_values)
        ec50_values = np.array(ec50_values)
        
        report_content = f"""
PandaDock Comprehensive Analysis Report
=====================================

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
Target: {self.protein_name}
Ligand: {self.ligand_name}

ALGORITHM INFORMATION
====================
Algorithm: {algorithm_info.get('algorithm', 'PandaDock')}
Version: {algorithm_info.get('version', 'Latest')}
Scoring Function: {algorithm_info.get('scoring_function', 'Unknown')}
Engine: {algorithm_info.get('engine', 'Unknown')}
Mode: {algorithm_info.get('mode', 'Unknown')}

COMMAND EXECUTED
===============
Command: {command_info.get('command', 'Not specified')}
Protein File: {command_info.get('protein', 'Not specified')}
Ligand File: {command_info.get('ligand', 'Not specified')}
Search Center: {command_info.get('center', 'Not specified')}
Search Size: {command_info.get('size', 'Not specified')}
Number of Poses: {len(poses_df)}
Exhaustiveness: {command_info.get('exhaustiveness', 'Default')}

SUMMARY STATISTICS
=================
Total Poses Generated: {len(poses_df)}
Best Binding Affinity: {poses_df['Binding_Affinity'].max():.3f} kcal/mol
Worst Binding Affinity: {poses_df['Binding_Affinity'].min():.3f} kcal/mol
Mean Binding Affinity: {poses_df['Binding_Affinity'].mean():.3f} ± {poses_df['Binding_Affinity'].std():.3f} kcal/mol

Best Docking Score: {poses_df['Score'].min():.4f}
Mean Docking Score: {poses_df['Score'].mean():.4f} ± {poses_df['Score'].std():.4f}

Best IC50: {ic50_values.min():.2e} μM
Median IC50: {np.median(ic50_values):.2e} μM
Worst IC50: {ic50_values.max():.2e} μM

Best EC50: {ec50_values.min():.2e} μM
Median EC50: {np.median(ec50_values):.2e} μM
Worst EC50: {ec50_values.max():.2e} μM

High Potency Poses (IC50 < 10 μM): {(ic50_values < 10).sum()}/{len(poses_df)}
Moderate Potency Poses (10-100 μM): {((ic50_values >= 10) & (ic50_values < 100)).sum()}/{len(poses_df)}
Low Potency Poses (≥ 100 μM): {(ic50_values >= 100).sum()}/{len(poses_df)}

Mean Confidence: {poses_df['Confidence'].mean():.3f} ± {poses_df['Confidence'].std():.3f}
High Confidence Poses (≥ 0.7): {(poses_df['Confidence'] >= 0.7).sum()}/{len(poses_df)}

DETAILED POSE ANALYSIS
=====================
"""
        
        # Add detailed pose information
        for idx, (_, pose) in enumerate(poses_df.iterrows()):
            ic50_val = ic50_values[idx]
            ec50_val = ec50_values[idx]
            
            report_content += f"""
Pose {idx + 1}: {pose['Pose_ID']}
--------------------------------
Rank: {pose['Rank']}
Docking Score: {pose['Score']:.4f}
Energy: {pose['Energy']:.3f} kcal/mol
Confidence: {pose['Confidence']:.3f}
Binding Affinity: {pose['Binding_Affinity']:.3f} kcal/mol
IC50: {ic50_val:.2e} μM
EC50: {ec50_val:.2e} μM
Ligand Efficiency: {pose['Ligand_Efficiency']:.4f}
Clash Score: {pose['Clash_Score']:.4f}
"""
        
        # Add correlation analysis
        affinity_energy_corr = stats.pearsonr(poses_df['Binding_Affinity'], poses_df['Energy'])[0]
        score_energy_corr = stats.pearsonr(poses_df['Score'], poses_df['Energy'])[0]
        affinity_ic50_corr = stats.pearsonr(poses_df['Binding_Affinity'], np.log10(ic50_values))[0]
        
        report_content += f"""

CORRELATION ANALYSIS
===================
Binding Affinity vs Energy: r = {affinity_energy_corr:.3f}
Docking Score vs Energy: r = {score_energy_corr:.3f}
Binding Affinity vs log(IC50): r = {affinity_ic50_corr:.3f}

GENERATED FILES
==============
"""
        
        for file_type, file_path in generated_files.items():
            if file_path:
                if isinstance(file_path, list):
                    # Handle list of files (e.g., interaction maps)
                    report_content += f"{file_type.replace('_', ' ').title()}:\n"
                    for i, fp in enumerate(file_path, 1):
                        report_content += f"  - Map {i}: {Path(fp).name}\n"
                else:
                    report_content += f"{file_type.replace('_', ' ').title()}: {Path(file_path).name}\n"
        
        report_content += f"""

POSE DATA (CSV FORMAT)
=====================
{poses_df.to_csv(index=False)}

ANALYSIS NOTES
=============
• IC50 values represent half-maximal inhibitory concentration
• EC50 values represent half-maximal effective concentration  
• Binding affinity calculated using thermodynamic relationship
• Ligand efficiency normalizes binding affinity by molecular size
• Confidence scores reflect docking algorithm certainty
• All energy values in kcal/mol, concentrations in μM

This report was generated by PandaDock comprehensive analysis system.
For questions or support, please refer to the PandaDock documentation.
"""
        
        # Save report
        output_file = self.output_dir / 'detailed_analysis_report.txt'
        with open(output_file, 'w') as f:
            f.write(report_content)
        
        self.logger.info(f"Generated detailed TXT report: {output_file}")
        return str(output_file)


def create_plots_for_pandadock(output_dir: str, poses_csv: str, poses_dir: str,
                              protein_name: str = None, ligand_name: str = None,
                              algorithm_info: Dict[str, Any] = None,
                              command_info: Dict[str, Any] = None) -> Dict[str, str]:
    """
    Convenience function to generate all plots for PandaDock results
    
    Args:
        output_dir: Directory to save plots
        poses_csv: Path to poses summary CSV file
        poses_dir: Directory containing pose PDB files
        protein_name: Name of the protein
        ligand_name: Name of the ligand
        algorithm_info: Algorithm information dictionary
        command_info: Command information dictionary
        
    Returns:
        Dictionary of generated file paths
    """
    
    # Default values
    if algorithm_info is None:
        algorithm_info = {
            'algorithm': 'PandaDock',
            'version': 'Latest',
            'scoring_function': 'Auto',
            'engine': 'Auto',
            'mode': 'Balanced'
        }
    
    if command_info is None:
        command_info = {
            'command': 'python -m pandadock',
            'protein': 'protein.pdb',
            'ligand': 'ligand.pdb',
            'center': 'Auto-detected',
            'size': 'Auto-detected',
            'exhaustiveness': 'Default'
        }
    
    # Initialize plot generator
    plot_gen = PandaDockPlotGenerator(output_dir, protein_name, ligand_name)
    
    # Generate all plots
    generated_files = plot_gen.generate_all_plots(poses_csv, poses_dir, 
                                                 algorithm_info, command_info)
    
    return generated_files


if __name__ == "__main__":
    # Example usage
    print("PandaDock Plot Generator")
    print("This module provides comprehensive plotting for PandaDock results")
    print("Use create_plots_for_pandadock() function to generate all plots")