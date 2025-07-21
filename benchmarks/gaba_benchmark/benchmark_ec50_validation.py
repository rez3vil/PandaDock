#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PandaDock EC50 Benchmarking Script
==================================

Comprehensive benchmarking of PandaDock predictions against experimental GABA_A receptor EC50 data.
This script validates docking accuracy and binding affinity predictions for a series of phenolic compounds.

Usage:
    python benchmark_ec50_validation.py

Features:
- Automated docking of all ligands against beta2-alpha1 GABA_A receptor
- Correlation analysis between predicted and experimental EC50 values
- Comprehensive statistical evaluation (R², RMSE, MAE, Pearson correlation)
- Publication-quality plots and visualizations
- Detailed performance reports for different scoring functions
"""

import os
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class EC50BenchmarkSuite:
    """Comprehensive benchmarking suite for PandaDock EC50 predictions"""
    
    def __init__(self, 
                 protein_file="tests/beta-2_alpha-1.pdb",
                 sdf_directory="tests/sdf",
                 results_dir="ec50_benchmark_results"):
        
        self.protein_file = protein_file
        self.sdf_directory = Path(sdf_directory)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        
        # Grid parameters for GABA_A receptor binding site
        self.center = [-15.7, -17.7, 8.18]
        self.size = [40, 40, 40]
        
        # Experimental EC50 data (μM) and ln(EC50) values
        self.experimental_data = {
            '8': {'ec50': 8.7, 'ln_ec50': 2.163323026},
            '1a': {'ec50': 6.4, 'ln_ec50': 1.85629799},
            '1b': {'ec50': 4.2, 'ln_ec50': 1.435084525},
            '1c': {'ec50': 4.1, 'ln_ec50': 1.410986974},
            '1d': {'ec50': 9.1, 'ln_ec50': 2.208274414},
            '1e': {'ec50': 3.0, 'ln_ec50': 1.098612289},
            '1f': {'ec50': 1.5, 'ln_ec50': 0.405465108},
            '1g': {'ec50': 3.1, 'ln_ec50': 1.131402111},
            '1h': {'ec50': 2.3, 'ln_ec50': 0.832909123},
            '1i': {'ec50': 0.6, 'ln_ec50': -0.510825624},
            '1j': {'ec50': 1.8, 'ln_ec50': 0.587786665},
            '1k': {'ec50': 1.0, 'ln_ec50': 0.0},
            '1l': {'ec50': 1.9, 'ln_ec50': 0.641853886},
            '1m': {'ec50': 2.4, 'ln_ec50': 0.875468737},
            '1n': {'ec50': 0.19, 'ln_ec50': -1.660731207},
            '1o': {'ec50': 4.3, 'ln_ec50': 1.458615023},
            '1p': {'ec50': 0.45, 'ln_ec50': -0.798507696},
            '1t': {'ec50': 3.4, 'ln_ec50': 1.223775432},
            '1u': {'ec50': 1.7, 'ln_ec50': 0.530628251},
            '24disecbutphenol': {'ec50': 32.9, 'ln_ec50': 3.493472658},
            '35ditertbutphenol': {'ec50': 94.4, 'ln_ec50': 4.547541073},
            '4iododiisopropylphenol': {'ec50': 11.1, 'ln_ec50': 2.406945108},
            'carboetomR': {'ec50': 5.4, 'ln_ec50': 1.686398954},
            'diethylphenol': {'ec50': 13.6, 'ln_ec50': 2.610069793},
            'diisopropylcatechol': {'ec50': 36.5, 'ln_ec50': 3.597312261},
            'dimethylphenol': {'ec50': 135.0, 'ln_ec50': 4.905274778},
            'disecbutphenol': {'ec50': 2.8, 'ln_ec50': 1.029619417},
            'etomR': {'ec50': 3.5, 'ln_ec50': 1.252762968},
            'isopropylphenol': {'ec50': 39.7, 'ln_ec50': 3.681351188},
            'phenol': {'ec50': 920.0, 'ln_ec50': 6.82437367}
        }
        
        # SDF file mapping (handle naming differences)
        self.sdf_file_mapping = {
            '24disecbutphenol': '24disecbutylphenol.pdbqt.sdf',
            '35ditertbutphenol': '35Ditertbutylphenol.pdbqt.sdf',
            '4iododiisopropylphenol': '4Iodo26diisopropylphenol.pdbqt.sdf',
            'diisopropylcatechol': 'diisopropylphenol.pdbqt.sdf'  # Need to verify this mapping
        }
        
        # Results storage
        self.docking_results = {}
        self.benchmark_stats = {}
        
        print(f"🧪 EC50 Benchmark Suite Initialized")
        print(f"📁 Protein: {self.protein_file}")
        print(f"📂 SDF Directory: {self.sdf_directory}")
        print(f"📊 Results Directory: {self.results_dir}")
        print(f"🎯 Grid Center: {self.center}")
        print(f"📏 Grid Size: {self.size}")
        print(f"🧬 Total Ligands: {len(self.experimental_data)}")
        print("")

    def get_sdf_file_path(self, ligand_name):
        """Get the correct SDF file path for a ligand"""
        # Check if there's a specific mapping
        if ligand_name in self.sdf_file_mapping:
            sdf_filename = self.sdf_file_mapping[ligand_name]
        else:
            sdf_filename = f"{ligand_name}.pdbqt.sdf"
        
        sdf_path = self.sdf_directory / sdf_filename
        
        # If file doesn't exist, try alternative naming
        if not sdf_path.exists():
            print(f"⚠️  File not found: {sdf_path}")
            # List available files to help with mapping
            available_files = list(self.sdf_directory.glob("*.sdf"))
            print(f"Available SDF files: {[f.name for f in available_files]}")
            return None
        
        return sdf_path

    def run_docking(self, ligand_name, scoring_function='pandacore', num_poses=10):
        """Run PandaDock for a single ligand"""
        
        sdf_path = self.get_sdf_file_path(ligand_name)
        if not sdf_path:
            print(f"❌ Skipping {ligand_name}: SDF file not found")
            return None
        
        output_dir = self.results_dir / f"{ligand_name}_{scoring_function}"
        
        print(f"🔬 Docking {ligand_name} with {scoring_function}...")
        
        cmd = [
            "python", "-m", "pandadock",
            "--protein", self.protein_file,
            "--ligand", str(sdf_path),
            "--center"] + [str(c) for c in self.center] + [
            "--size"] + [str(s) for s in self.size] + [
            "--scoring", scoring_function,
            "--num-poses", str(num_poses),
            "--report-format", "json",
            "--out", str(output_dir)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print(f"✅ {ligand_name} completed successfully")
                return self.parse_docking_results(output_dir, ligand_name, scoring_function)
            else:
                print(f"❌ {ligand_name} failed:")
                print(f"   Error: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            print(f"⏰ {ligand_name} timed out")
            return None
        except Exception as e:
            print(f"💥 {ligand_name} crashed: {str(e)}")
            return None

    def parse_docking_results(self, output_dir, ligand_name, scoring_function):
        """Parse docking results from JSON output"""
        
        json_file = output_dir / "pandadock_report.json"
        csv_file = output_dir / "poses" / "poses_summary.csv"
        
        results = {
            'ligand_name': ligand_name,
            'scoring_function': scoring_function,
            'experimental_ec50': self.experimental_data[ligand_name]['ec50'],
            'experimental_ln_ec50': self.experimental_data[ligand_name]['ln_ec50']
        }
        
        try:
            # Parse JSON results
            if json_file.exists():
                with open(json_file, 'r') as f:
                    json_data = json.load(f)
                
                # Get best pose data
                if 'poses' in json_data and len(json_data['poses']) > 0:
                    best_pose = json_data['poses'][0]  # Best ranked pose
                    
                    results.update({
                        'best_score': best_pose.get('score', np.nan),
                        'best_energy': best_pose.get('energy', np.nan),
                        'best_binding_affinity': best_pose.get('binding_affinity', np.nan),
                        'predicted_ic50_um': best_pose.get('ic50_um', np.nan),
                        'predicted_ec50_um': best_pose.get('ec50_um', np.nan),
                        'ligand_efficiency': best_pose.get('ligand_efficiency', np.nan),
                        'confidence': best_pose.get('confidence', np.nan)
                    })
                
                # Get summary statistics
                if 'summary' in json_data:
                    summary = json_data['summary']
                    results.update({
                        'mean_score': summary.get('mean_score', np.nan),
                        'mean_binding_affinity': summary.get('mean_binding_affinity', np.nan),
                        'best_summary_ic50_um': summary.get('best_ic50_um', np.nan),
                        'best_summary_ec50_um': summary.get('best_ec50_um', np.nan)
                    })
            
            # Calculate ln(predicted_EC50) for correlation analysis
            if not np.isnan(results.get('predicted_ec50_um', np.nan)) and results['predicted_ec50_um'] > 0:
                results['predicted_ln_ec50'] = np.log(results['predicted_ec50_um'])
            else:
                results['predicted_ln_ec50'] = np.nan
            
            return results
            
        except Exception as e:
            print(f"⚠️  Error parsing results for {ligand_name}: {str(e)}")
            return results

    def run_comprehensive_benchmark(self):
        """Run docking for all ligands with different scoring functions"""
        
        scoring_functions = ['pandacore', 'pandaml', 'pandaphysics']
        
        print("🚀 Starting Comprehensive EC50 Benchmark")
        print("=" * 60)
        
        for scoring_func in scoring_functions:
            print(f"\n📊 Testing {scoring_func.upper()} scoring function")
            print("-" * 40)
            
            scoring_results = []
            
            for i, ligand_name in enumerate(self.experimental_data.keys(), 1):
                print(f"[{i:2d}/{len(self.experimental_data)}] ", end="")
                
                result = self.run_docking(ligand_name, scoring_func, num_poses=5)
                if result:
                    scoring_results.append(result)
                    
                    # Quick feedback
                    exp_ec50 = result['experimental_ec50']
                    pred_ec50 = result.get('predicted_ec50_um', 'N/A')
                    print(f"    Exp: {exp_ec50:.2f} μM, Pred: {pred_ec50:.1e} μM" if pred_ec50 != 'N/A' else f"    Exp: {exp_ec50:.2f} μM, Pred: N/A")
            
            self.docking_results[scoring_func] = scoring_results
            print(f"\n✅ {scoring_func.upper()} completed: {len(scoring_results)}/{len(self.experimental_data)} successful")

    def calculate_benchmark_statistics(self):
        """Calculate comprehensive benchmark statistics"""
        
        print("\n📈 Calculating Benchmark Statistics")
        print("=" * 50)
        
        for scoring_func, results in self.docking_results.items():
            
            # Convert to DataFrame for easier analysis
            df = pd.DataFrame(results)
            
            # Remove rows with missing predictions
            df_clean = df.dropna(subset=['predicted_ec50_um', 'predicted_ln_ec50'])
            
            if len(df_clean) == 0:
                print(f"⚠️  No valid predictions for {scoring_func}")
                continue
            
            print(f"\n🔬 {scoring_func.upper()} Statistics ({len(df_clean)}/{len(df)} valid predictions)")
            
            # Extract experimental and predicted values
            exp_ec50 = df_clean['experimental_ec50'].values
            pred_ec50 = df_clean['predicted_ec50_um'].values
            exp_ln_ec50 = df_clean['experimental_ln_ec50'].values
            pred_ln_ec50 = df_clean['predicted_ln_ec50'].values
            
            binding_affinity = df_clean['best_binding_affinity'].values
            scores = df_clean['best_score'].values
            
            stats_dict = {}
            
            # EC50 Correlations (log space is more meaningful for EC50)
            if len(exp_ln_ec50) > 1:
                # Pearson correlation
                r_ln_ec50, p_ln_ec50 = pearsonr(exp_ln_ec50, pred_ln_ec50)
                # Spearman correlation  
                rho_ln_ec50, p_rho_ln_ec50 = spearmanr(exp_ln_ec50, pred_ln_ec50)
                
                # R-squared
                r2_ln_ec50 = r_ln_ec50 ** 2
                
                # RMSE and MAE in log space
                rmse_ln_ec50 = np.sqrt(np.mean((exp_ln_ec50 - pred_ln_ec50) ** 2))
                mae_ln_ec50 = np.mean(np.abs(exp_ln_ec50 - pred_ln_ec50))
                
                stats_dict.update({
                    'ln_ec50_pearson_r': r_ln_ec50,
                    'ln_ec50_pearson_p': p_ln_ec50,
                    'ln_ec50_spearman_rho': rho_ln_ec50,
                    'ln_ec50_spearman_p': p_rho_ln_ec50,
                    'ln_ec50_r2': r2_ln_ec50,
                    'ln_ec50_rmse': rmse_ln_ec50,
                    'ln_ec50_mae': mae_ln_ec50,
                })
                
                print(f"   ln(EC50) Pearson r: {r_ln_ec50:.3f} (p={p_ln_ec50:.3e})")
                print(f"   ln(EC50) R²: {r2_ln_ec50:.3f}")
                print(f"   ln(EC50) RMSE: {rmse_ln_ec50:.3f}")
                print(f"   ln(EC50) MAE: {mae_ln_ec50:.3f}")
            
            # Binding Affinity Correlations
            if len(binding_affinity) > 1:
                binding_affinity_clean = binding_affinity[~np.isnan(binding_affinity)]
                exp_ln_ec50_clean = exp_ln_ec50[~np.isnan(binding_affinity)]
                
                if len(binding_affinity_clean) > 1:
                    r_affinity, p_affinity = pearsonr(exp_ln_ec50_clean, binding_affinity_clean)
                    rho_affinity, p_rho_affinity = spearmanr(exp_ln_ec50_clean, binding_affinity_clean)
                    
                    stats_dict.update({
                        'affinity_pearson_r': r_affinity,
                        'affinity_pearson_p': p_affinity,
                        'affinity_spearman_rho': rho_affinity,
                        'affinity_spearman_p': p_rho_affinity,
                        'affinity_r2': r_affinity ** 2
                    })
                    
                    print(f"   Binding Affinity vs ln(EC50) r: {r_affinity:.3f} (p={p_affinity:.3e})")
            
            # Score Correlations
            if len(scores) > 1:
                scores_clean = scores[~np.isnan(scores)]
                exp_ln_ec50_clean = exp_ln_ec50[~np.isnan(scores)]
                
                if len(scores_clean) > 1:
                    r_score, p_score = pearsonr(exp_ln_ec50_clean, scores_clean)
                    rho_score, p_rho_score = spearmanr(exp_ln_ec50_clean, scores_clean)
                    
                    stats_dict.update({
                        'score_pearson_r': r_score,
                        'score_pearson_p': p_score,
                        'score_spearman_rho': rho_score,
                        'score_spearman_p': p_rho_score,
                        'score_r2': r_score ** 2
                    })
                    
                    print(f"   Docking Score vs ln(EC50) r: {r_score:.3f} (p={p_score:.3e})")
            
            # Success rate
            stats_dict['success_rate'] = len(df_clean) / len(df)
            stats_dict['n_successful'] = len(df_clean)
            stats_dict['n_total'] = len(df)
            
            print(f"   Success Rate: {stats_dict['success_rate']:.1%} ({stats_dict['n_successful']}/{stats_dict['n_total']})")
            
            # Save detailed results
            self.benchmark_stats[scoring_func] = {
                'stats': stats_dict,
                'data': df_clean
            }

    def create_comprehensive_plots(self):
        """Generate comprehensive benchmark plots"""
        
        print("\n🎨 Generating Comprehensive Plots")
        print("=" * 40)
        
        # Set up the plotting style
        plt.rcParams.update({
            'figure.figsize': (12, 8),
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10
        })
        
        scoring_functions = list(self.docking_results.keys())
        n_scoring = len(scoring_functions)
        
        if n_scoring == 0:
            print("⚠️  No results to plot")
            return
        
        # Plot 1: ln(EC50) Correlation Plot
        self._plot_ln_ec50_correlations(scoring_functions)
        
        # Plot 2: Binding Affinity vs ln(EC50)
        self._plot_binding_affinity_correlations(scoring_functions)
        
        # Plot 3: Score vs ln(EC50)
        self._plot_score_correlations(scoring_functions)
        
        # Plot 4: Performance Summary
        self._plot_performance_summary(scoring_functions)
        
        # Plot 5: Prediction vs Experimental Scatter Matrix
        self._plot_prediction_matrix(scoring_functions)
        
        # Plot 6: Error Analysis
        self._plot_error_analysis(scoring_functions)
        
        # Plot 7: Individual Ligand Performance
        self._plot_individual_ligand_performance(scoring_functions)
        
        print("✅ All plots generated successfully")

    def _plot_ln_ec50_correlations(self, scoring_functions):
        """Plot ln(EC50) correlations for all scoring functions"""
        
        fig, axes = plt.subplots(1, len(scoring_functions), figsize=(6*len(scoring_functions), 5))
        if len(scoring_functions) == 1:
            axes = [axes]
        
        for i, scoring_func in enumerate(scoring_functions):
            ax = axes[i]
            
            if scoring_func not in self.benchmark_stats:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{scoring_func.upper()}\nNo Valid Data')
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            stats = self.benchmark_stats[scoring_func]['stats']
            
            # Scatter plot
            ax.scatter(df['experimental_ln_ec50'], df['predicted_ln_ec50'], 
                      alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
            
            # Perfect correlation line
            min_val = min(df['experimental_ln_ec50'].min(), df['predicted_ln_ec50'].min())
            max_val = max(df['experimental_ln_ec50'].max(), df['predicted_ln_ec50'].max())
            ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, linewidth=2, label='Perfect Correlation')
            
            # Regression line
            z = np.polyfit(df['experimental_ln_ec50'], df['predicted_ln_ec50'], 1)
            p = np.poly1d(z)
            ax.plot(df['experimental_ln_ec50'], p(df['experimental_ln_ec50']), 'b-', alpha=0.8, linewidth=2, label='Best Fit')
            
            # Statistics text
            r = stats.get('ln_ec50_pearson_r', np.nan)
            r2 = stats.get('ln_ec50_r2', np.nan)
            rmse = stats.get('ln_ec50_rmse', np.nan)
            
            stats_text = f'r = {r:.3f}\nR² = {r2:.3f}\nRMSE = {rmse:.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.set_xlabel('Experimental ln(EC50)')
            ax.set_ylabel('Predicted ln(EC50)')
            ax.set_title(f'{scoring_func.upper()}\nln(EC50) Correlation')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'ln_ec50_correlations.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'ln_ec50_correlations.pdf', bbox_inches='tight')
        plt.show()

    def _plot_binding_affinity_correlations(self, scoring_functions):
        """Plot binding affinity vs ln(EC50) correlations"""
        
        fig, axes = plt.subplots(1, len(scoring_functions), figsize=(6*len(scoring_functions), 5))
        if len(scoring_functions) == 1:
            axes = [axes]
        
        for i, scoring_func in enumerate(scoring_functions):
            ax = axes[i]
            
            if scoring_func not in self.benchmark_stats:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{scoring_func.upper()}\nNo Valid Data')
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            stats = self.benchmark_stats[scoring_func]['stats']
            
            # Remove NaN values
            valid_mask = ~(np.isnan(df['experimental_ln_ec50']) | np.isnan(df['best_binding_affinity']))
            if valid_mask.sum() == 0:
                ax.text(0.5, 0.5, 'No Valid Data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            x = df['experimental_ln_ec50'][valid_mask]
            y = df['best_binding_affinity'][valid_mask]
            
            # Scatter plot
            ax.scatter(x, y, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
            
            # Regression line
            if len(x) > 1:
                z = np.polyfit(x, y, 1)
                p = np.poly1d(z)
                ax.plot(x, p(x), 'r-', alpha=0.8, linewidth=2, label='Best Fit')
            
            # Statistics text
            r = stats.get('affinity_pearson_r', np.nan)
            r2 = stats.get('affinity_r2', np.nan)
            
            stats_text = f'r = {r:.3f}\nR² = {r2:.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.set_xlabel('Experimental ln(EC50)')
            ax.set_ylabel('Binding Affinity (kcal/mol)')
            ax.set_title(f'{scoring_func.upper()}\nBinding Affinity vs ln(EC50)')
            if len(x) > 1:
                ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'binding_affinity_correlations.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'binding_affinity_correlations.pdf', bbox_inches='tight')
        plt.show()

    def _plot_score_correlations(self, scoring_functions):
        """Plot docking score vs ln(EC50) correlations"""
        
        fig, axes = plt.subplots(1, len(scoring_functions), figsize=(6*len(scoring_functions), 5))
        if len(scoring_functions) == 1:
            axes = [axes]
        
        for i, scoring_func in enumerate(scoring_functions):
            ax = axes[i]
            
            if scoring_func not in self.benchmark_stats:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{scoring_func.upper()}\nNo Valid Data')
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            stats = self.benchmark_stats[scoring_func]['stats']
            
            # Remove NaN values
            valid_mask = ~(np.isnan(df['experimental_ln_ec50']) | np.isnan(df['best_score']))
            if valid_mask.sum() == 0:
                ax.text(0.5, 0.5, 'No Valid Data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            x = df['experimental_ln_ec50'][valid_mask]
            y = df['best_score'][valid_mask]
            
            # Scatter plot
            ax.scatter(x, y, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
            
            # Regression line
            if len(x) > 1:
                z = np.polyfit(x, y, 1)
                p = np.poly1d(z)
                ax.plot(x, p(x), 'r-', alpha=0.8, linewidth=2, label='Best Fit')
            
            # Statistics text
            r = stats.get('score_pearson_r', np.nan)
            r2 = stats.get('score_r2', np.nan)
            
            stats_text = f'r = {r:.3f}\nR² = {r2:.3f}'
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.set_xlabel('Experimental ln(EC50)')
            ax.set_ylabel('Docking Score')
            ax.set_title(f'{scoring_func.upper()}\nDocking Score vs ln(EC50)')
            if len(x) > 1:
                ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'score_correlations.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'score_correlations.pdf', bbox_inches='tight')
        plt.show()

    def _plot_performance_summary(self, scoring_functions):
        """Plot performance summary across all scoring functions"""
        
        # Collect statistics
        performance_data = []
        
        for scoring_func in scoring_functions:
            if scoring_func in self.benchmark_stats:
                stats = self.benchmark_stats[scoring_func]['stats']
                performance_data.append({
                    'Scoring Function': scoring_func.upper(),
                    'ln(EC50) R²': stats.get('ln_ec50_r2', 0),
                    'ln(EC50) RMSE': stats.get('ln_ec50_rmse', np.inf),
                    'Affinity R²': stats.get('affinity_r2', 0),
                    'Score R²': stats.get('score_r2', 0),
                    'Success Rate': stats.get('success_rate', 0)
                })
        
        if not performance_data:
            print("⚠️  No performance data to plot")
            return
        
        df_perf = pd.DataFrame(performance_data)
        
        # Create subplots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # R² comparison
        metrics = ['ln(EC50) R²', 'Affinity R²', 'Score R²']
        ax = axes[0, 0]
        x = np.arange(len(scoring_functions))
        width = 0.25
        
        for i, metric in enumerate(metrics):
            values = df_perf[metric].values
            ax.bar(x + i*width, values, width, label=metric)
        
        ax.set_xlabel('Scoring Function')
        ax.set_ylabel('R²')
        ax.set_title('R² Comparison Across Scoring Functions')
        ax.set_xticks(x + width)
        ax.set_xticklabels(df_perf['Scoring Function'])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # RMSE comparison
        ax = axes[0, 1]
        rmse_values = df_perf['ln(EC50) RMSE'].values
        ax.bar(df_perf['Scoring Function'], rmse_values)
        ax.set_xlabel('Scoring Function')
        ax.set_ylabel('ln(EC50) RMSE')
        ax.set_title('RMSE Comparison')
        ax.grid(True, alpha=0.3)
        
        # Success rate
        ax = axes[1, 0]
        success_values = df_perf['Success Rate'].values * 100
        bars = ax.bar(df_perf['Scoring Function'], success_values)
        ax.set_xlabel('Scoring Function')
        ax.set_ylabel('Success Rate (%)')
        ax.set_title('Success Rate Comparison')
        ax.set_ylim(0, 100)
        
        # Add value labels on bars
        for bar, value in zip(bars, success_values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{value:.1f}%', ha='center', va='bottom')
        ax.grid(True, alpha=0.3)
        
        # Overall performance radar chart
        ax = axes[1, 1]
        
        # Normalize metrics for radar chart (0-1 scale)
        normalized_data = []
        for _, row in df_perf.iterrows():
            normalized_data.append([
                row['ln(EC50) R²'],  # Already 0-1
                1 - min(row['ln(EC50) RMSE'] / 5, 1),  # Invert RMSE (lower is better)
                row['Success Rate']  # Already 0-1
            ])
        
        categories = ['ln(EC50) R²', 'Low RMSE', 'Success Rate']
        N = len(categories)
        
        # Angles for radar chart
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Complete the circle
        
        colors = ['blue', 'red', 'green']
        
        for i, (scoring_func, data) in enumerate(zip(df_perf['Scoring Function'], normalized_data)):
            data += data[:1]  # Complete the circle
            ax.plot(angles, data, 'o-', linewidth=2, label=scoring_func, color=colors[i % len(colors)])
            ax.fill(angles, data, alpha=0.25, color=colors[i % len(colors)])
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        ax.set_ylim(0, 1)
        ax.set_title('Overall Performance Comparison')
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
        ax.grid(True)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'performance_summary.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'performance_summary.pdf', bbox_inches='tight')
        plt.show()

    def _plot_prediction_matrix(self, scoring_functions):
        """Plot prediction vs experimental matrix for all metrics"""
        
        metrics = [
            ('predicted_ec50_um', 'experimental_ec50', 'EC50 (μM)', 'log'),
            ('predicted_ln_ec50', 'experimental_ln_ec50', 'ln(EC50)', 'linear'),
            ('best_binding_affinity', 'experimental_ln_ec50', 'Binding Affinity vs ln(EC50)', 'linear')
        ]
        
        fig, axes = plt.subplots(len(metrics), len(scoring_functions), 
                                figsize=(6*len(scoring_functions), 5*len(metrics)))
        
        if len(scoring_functions) == 1:
            axes = axes.reshape(-1, 1)
        if len(metrics) == 1:
            axes = axes.reshape(1, -1)
        
        for i, (pred_col, exp_col, title, scale) in enumerate(metrics):
            for j, scoring_func in enumerate(scoring_functions):
                ax = axes[i, j]
                
                if scoring_func not in self.benchmark_stats:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'{scoring_func.upper()}\n{title}')
                    continue
                
                df = self.benchmark_stats[scoring_func]['data']
                
                # Handle different column combinations
                if pred_col == 'best_binding_affinity' and exp_col == 'experimental_ln_ec50':
                    x_data = df['experimental_ln_ec50']
                    y_data = df['best_binding_affinity']
                    x_label = 'Experimental ln(EC50)'
                    y_label = 'Predicted Binding Affinity (kcal/mol)'
                else:
                    x_data = df[exp_col]
                    y_data = df[pred_col]
                    if 'experimental' in exp_col:
                        x_label = f'Experimental {title.split("vs")[0].strip() if "vs" in title else title}'
                    else:
                        x_label = f'Experimental {title}'
                    y_label = f'Predicted {title.split("vs")[0].strip() if "vs" in title else title}'
                
                # Remove NaN values
                valid_mask = ~(np.isnan(x_data) | np.isnan(y_data))
                if valid_mask.sum() == 0:
                    ax.text(0.5, 0.5, 'No Valid Data', ha='center', va='center', transform=ax.transAxes)
                    continue
                
                x_clean = x_data[valid_mask]
                y_clean = y_data[valid_mask]
                
                # Scatter plot
                ax.scatter(x_clean, y_clean, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
                
                # Perfect correlation line (only for pred vs exp of same metric)
                if pred_col.replace('predicted_', '') == exp_col.replace('experimental_', ''):
                    min_val = min(x_clean.min(), y_clean.min())
                    max_val = max(x_clean.max(), y_clean.max())
                    ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, linewidth=2, label='Perfect')
                
                # Regression line
                if len(x_clean) > 1:
                    z = np.polyfit(x_clean, y_clean, 1)
                    p = np.poly1d(z)
                    ax.plot(x_clean, p(x_clean), 'b-', alpha=0.8, linewidth=2, label='Best Fit')
                
                # Calculate correlation
                if len(x_clean) > 1:
                    r, p_val = pearsonr(x_clean, y_clean)
                    stats_text = f'r = {r:.3f}\n(p={p_val:.3e})'
                    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                ax.set_title(f'{scoring_func.upper()}\n{title}')
                
                if scale == 'log':
                    ax.set_xscale('log')
                    ax.set_yscale('log')
                
                if len(x_clean) > 1:
                    ax.legend()
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'prediction_matrix.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'prediction_matrix.pdf', bbox_inches='tight')
        plt.show()

    def _plot_error_analysis(self, scoring_functions):
        """Plot error analysis for predictions"""
        
        fig, axes = plt.subplots(2, len(scoring_functions), figsize=(6*len(scoring_functions), 10))
        if len(scoring_functions) == 1:
            axes = axes.reshape(-1, 1)
        
        for j, scoring_func in enumerate(scoring_functions):
            
            if scoring_func not in self.benchmark_stats:
                for i in range(2):
                    axes[i, j].text(0.5, 0.5, 'No Data', ha='center', va='center', transform=axes[i, j].transAxes)
                    axes[i, j].set_title(f'{scoring_func.upper()}\nNo Valid Data')
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            
            # Calculate errors
            ln_ec50_error = df['predicted_ln_ec50'] - df['experimental_ln_ec50']
            valid_errors = ln_ec50_error.dropna()
            
            if len(valid_errors) == 0:
                for i in range(2):
                    axes[i, j].text(0.5, 0.5, 'No Valid Data', ha='center', va='center', transform=axes[i, j].transAxes)
                continue
            
            # Error distribution histogram
            ax = axes[0, j]
            ax.hist(valid_errors, bins=10, alpha=0.7, edgecolor='black')
            ax.axvline(0, color='red', linestyle='--', alpha=0.8, linewidth=2, label='Perfect Prediction')
            ax.axvline(valid_errors.mean(), color='blue', linestyle='-', alpha=0.8, linewidth=2, 
                      label=f'Mean Error = {valid_errors.mean():.3f}')
            ax.set_xlabel('ln(EC50) Prediction Error')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{scoring_func.upper()}\nError Distribution')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Error vs experimental value
            ax = axes[1, j]
            exp_values = df['experimental_ln_ec50'][~ln_ec50_error.isna()]
            ax.scatter(exp_values, valid_errors, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
            ax.axhline(0, color='red', linestyle='--', alpha=0.8, linewidth=2, label='Perfect Prediction')
            ax.axhline(valid_errors.mean(), color='blue', linestyle='-', alpha=0.8, linewidth=2, 
                      label=f'Mean Error = {valid_errors.mean():.3f}')
            ax.set_xlabel('Experimental ln(EC50)')
            ax.set_ylabel('ln(EC50) Prediction Error')
            ax.set_title(f'{scoring_func.upper()}\nError vs Experimental Value')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'error_analysis.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'error_analysis.pdf', bbox_inches='tight')
        plt.show()

    def _plot_individual_ligand_performance(self, scoring_functions):
        """Plot performance for individual ligands"""
        
        # Collect all ligand data
        all_ligands = set()
        for scoring_func in scoring_functions:
            if scoring_func in self.benchmark_stats:
                all_ligands.update(self.benchmark_stats[scoring_func]['data']['ligand_name'])
        
        all_ligands = sorted(list(all_ligands))
        
        if len(all_ligands) == 0:
            print("⚠️  No ligand data to plot")
            return
        
        # Create figure
        fig, axes = plt.subplots(3, 1, figsize=(15, 12))
        
        # Prepare data for plotting
        exp_ec50_values = [self.experimental_data.get(lig, {}).get('ec50', np.nan) for lig in all_ligands]
        exp_ln_ec50_values = [self.experimental_data.get(lig, {}).get('ln_ec50', np.nan) for lig in all_ligands]
        
        # Plot 1: Experimental EC50 values
        ax = axes[0]
        bars = ax.bar(range(len(all_ligands)), exp_ec50_values, alpha=0.7)
        ax.set_xlabel('Ligands')
        ax.set_ylabel('Experimental EC50 (μM)')
        ax.set_title('Experimental EC50 Values for All Ligands')
        ax.set_xticks(range(len(all_ligands)))
        ax.set_xticklabels(all_ligands, rotation=45, ha='right')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, exp_ec50_values)):
            if not np.isnan(value):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{value:.2f}', ha='center', va='bottom', fontsize=8, rotation=45)
        
        # Plot 2: Predicted vs Experimental ln(EC50) for each scoring function
        ax = axes[1]
        x = np.arange(len(all_ligands))
        width = 0.8 / len(scoring_functions)
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        
        for i, scoring_func in enumerate(scoring_functions):
            if scoring_func not in self.benchmark_stats:
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            
            pred_values = []
            for lig in all_ligands:
                lig_data = df[df['ligand_name'] == lig]
                if len(lig_data) > 0:
                    pred_values.append(lig_data['predicted_ln_ec50'].iloc[0])
                else:
                    pred_values.append(np.nan)
            
            # Remove NaN for plotting
            valid_indices = [j for j, val in enumerate(pred_values) if not np.isnan(val)]
            valid_x = [x[j] + i*width for j in valid_indices]
            valid_pred = [pred_values[j] for j in valid_indices]
            
            ax.bar(valid_x, valid_pred, width, label=f'{scoring_func.upper()} Predicted', 
                  alpha=0.7, color=colors[i % len(colors)])
        
        # Add experimental values
        valid_exp_indices = [j for j, val in enumerate(exp_ln_ec50_values) if not np.isnan(val)]
        valid_exp_x = [x[j] + len(scoring_functions)*width for j in valid_exp_indices]
        valid_exp_values = [exp_ln_ec50_values[j] for j in valid_exp_indices]
        
        ax.bar(valid_exp_x, valid_exp_values, width, label='Experimental', 
               alpha=0.7, color='black', edgecolor='white', linewidth=1)
        
        ax.set_xlabel('Ligands')
        ax.set_ylabel('ln(EC50)')
        ax.set_title('Predicted vs Experimental ln(EC50) by Ligand')
        ax.set_xticks(x + width * len(scoring_functions) / 2)
        ax.set_xticklabels(all_ligands, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Prediction errors by ligand
        ax = axes[2]
        
        for i, scoring_func in enumerate(scoring_functions):
            if scoring_func not in self.benchmark_stats:
                continue
            
            df = self.benchmark_stats[scoring_func]['data']
            
            errors = []
            for lig in all_ligands:
                lig_data = df[df['ligand_name'] == lig]
                if len(lig_data) > 0:
                    pred = lig_data['predicted_ln_ec50'].iloc[0]
                    exp = self.experimental_data.get(lig, {}).get('ln_ec50', np.nan)
                    if not np.isnan(pred) and not np.isnan(exp):
                        errors.append(pred - exp)
                    else:
                        errors.append(np.nan)
                else:
                    errors.append(np.nan)
            
            # Remove NaN for plotting
            valid_indices = [j for j, val in enumerate(errors) if not np.isnan(val)]
            valid_x = [x[j] + i*width for j in valid_indices]
            valid_errors = [errors[j] for j in valid_indices]
            
            ax.bar(valid_x, valid_errors, width, label=f'{scoring_func.upper()}', 
                  alpha=0.7, color=colors[i % len(colors)])
        
        ax.axhline(0, color='black', linestyle='-', alpha=0.8, linewidth=2, label='Perfect Prediction')
        ax.set_xlabel('Ligands')
        ax.set_ylabel('ln(EC50) Prediction Error')
        ax.set_title('Prediction Errors by Ligand')
        ax.set_xticks(x + width * len(scoring_functions) / 2)
        ax.set_xticklabels(all_ligands, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'individual_ligand_performance.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.results_dir / 'individual_ligand_performance.pdf', bbox_inches='tight')
        plt.show()

    def generate_comprehensive_report(self):
        """Generate a comprehensive benchmark report"""
        
        report_file = self.results_dir / 'comprehensive_benchmark_report.md'
        
        print(f"\n📄 Generating Comprehensive Report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("# PandaDock EC50 Benchmarking Report\n\n")
            f.write(f"**Generated on:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Executive Summary\n\n")
            f.write("This report presents a comprehensive evaluation of PandaDock's ability to predict ")
            f.write("GABA_A receptor EC50 values for a series of phenolic compounds. The benchmark ")
            f.write("compares three scoring functions (PandaCore, PandaML, and PandaPhysics) against ")
            f.write("experimental data.\n\n")
            
            f.write("## Dataset Information\n\n")
            f.write(f"- **Total Ligands:** {len(self.experimental_data)}\n")
            f.write(f"- **Protein Target:** Beta2-Alpha1 GABA_A Receptor\n")
            f.write(f"- **Grid Center:** {self.center}\n")
            f.write(f"- **Grid Size:** {self.size}\n")
            f.write(f"- **EC50 Range:** {min(d['ec50'] for d in self.experimental_data.values()):.2f} - ")
            f.write(f"{max(d['ec50'] for d in self.experimental_data.values()):.1f} μM\n\n")
            
            f.write("## Scoring Function Performance\n\n")
            
            for scoring_func in self.docking_results.keys():
                if scoring_func not in self.benchmark_stats:
                    continue
                
                f.write(f"### {scoring_func.upper()}\n\n")
                
                stats = self.benchmark_stats[scoring_func]['stats']
                data = self.benchmark_stats[scoring_func]['data']
                
                f.write(f"- **Success Rate:** {stats['success_rate']:.1%} ({stats['n_successful']}/{stats['n_total']})\n")
                
                if 'ln_ec50_pearson_r' in stats:
                    f.write(f"- **ln(EC50) Correlation:**\n")
                    f.write(f"  - Pearson r: {stats['ln_ec50_pearson_r']:.3f} (p={stats['ln_ec50_pearson_p']:.3e})\n")
                    f.write(f"  - R²: {stats['ln_ec50_r2']:.3f}\n")
                    f.write(f"  - RMSE: {stats['ln_ec50_rmse']:.3f}\n")
                    f.write(f"  - MAE: {stats['ln_ec50_mae']:.3f}\n")
                
                if 'affinity_pearson_r' in stats:
                    f.write(f"- **Binding Affinity vs ln(EC50):**\n")
                    f.write(f"  - Pearson r: {stats['affinity_pearson_r']:.3f} (p={stats['affinity_pearson_p']:.3e})\n")
                    f.write(f"  - R²: {stats['affinity_r2']:.3f}\n")
                
                if 'score_pearson_r' in stats:
                    f.write(f"- **Docking Score vs ln(EC50):**\n")
                    f.write(f"  - Pearson r: {stats['score_pearson_r']:.3f} (p={stats['score_pearson_p']:.3e})\n")
                    f.write(f"  - R²: {stats['score_r2']:.3f}\n")
                
                f.write("\n")
            
            f.write("## Key Findings\n\n")
            
            # Find best performing scoring function
            best_r2 = 0
            best_scoring = None
            for scoring_func in self.benchmark_stats.keys():
                r2 = self.benchmark_stats[scoring_func]['stats'].get('ln_ec50_r2', 0)
                if r2 > best_r2:
                    best_r2 = r2
                    best_scoring = scoring_func
            
            if best_scoring:
                f.write(f"1. **Best Performing Scoring Function:** {best_scoring.upper()} (R² = {best_r2:.3f})\n")
            
            # Success rates
            success_rates = [(sf, stats['stats']['success_rate']) for sf, stats in self.benchmark_stats.items()]
            if success_rates:
                avg_success = np.mean([sr[1] for sr in success_rates])
                f.write(f"2. **Average Success Rate:** {avg_success:.1%}\n")
            
            # Correlation strength
            correlations = []
            for scoring_func in self.benchmark_stats.keys():
                r = self.benchmark_stats[scoring_func]['stats'].get('ln_ec50_pearson_r', 0)
                correlations.append(abs(r))
            
            if correlations:
                avg_correlation = np.mean(correlations)
                f.write(f"3. **Average |Correlation|:** {avg_correlation:.3f}\n")
            
            f.write("\n## Recommendations\n\n")
            
            if best_r2 < 0.5:
                f.write("1. **Model Improvement Needed:** Current R² values suggest moderate predictive power. ")
                f.write("Consider parameter optimization or alternative scoring approaches.\n")
            else:
                f.write("1. **Good Predictive Performance:** Current models show reasonable correlation with experimental data.\n")
            
            if success_rates and avg_success < 0.8:
                f.write("2. **Docking Reliability:** Some ligands failed to dock successfully. ")
                f.write("Consider reviewing grid settings or ligand preparation.\n")
            
            f.write("3. **Validation:** Consider expanding the dataset and performing cross-validation ")
            f.write("for more robust model evaluation.\n\n")
            
            f.write("## Files Generated\n\n")
            f.write("- `ln_ec50_correlations.png/pdf`: ln(EC50) correlation plots\n")
            f.write("- `binding_affinity_correlations.png/pdf`: Binding affinity correlation plots\n")
            f.write("- `score_correlations.png/pdf`: Docking score correlation plots\n")
            f.write("- `performance_summary.png/pdf`: Overall performance comparison\n")
            f.write("- `prediction_matrix.png/pdf`: Comprehensive prediction matrix\n")
            f.write("- `error_analysis.png/pdf`: Error distribution and analysis\n")
            f.write("- `individual_ligand_performance.png/pdf`: Per-ligand performance analysis\n")
            f.write("- `benchmark_results.csv`: Raw benchmark data\n")
            f.write("- `benchmark_statistics.json`: Statistical summary\n\n")
        
        print("✅ Comprehensive report generated successfully")

    def save_results(self):
        """Save all benchmark results to files"""
        
        print("\n💾 Saving Benchmark Results")
        print("=" * 30)
        
        # Save raw results as CSV
        all_results = []
        for scoring_func, results in self.docking_results.items():
            for result in results:
                result['scoring_function'] = scoring_func
                all_results.append(result)
        
        if all_results:
            results_df = pd.DataFrame(all_results)
            results_file = self.results_dir / 'benchmark_results.csv'
            results_df.to_csv(results_file, index=False)
            print(f"✅ Raw results saved: {results_file}")
        
        # Save statistics as JSON
        stats_file = self.results_dir / 'benchmark_statistics.json'
        with open(stats_file, 'w') as f:
            # Convert numpy types to Python types for JSON serialization
            json_stats = {}
            for scoring_func, data in self.benchmark_stats.items():
                json_stats[scoring_func] = {}
                for key, value in data['stats'].items():
                    if isinstance(value, (np.int64, np.int32)):
                        json_stats[scoring_func][key] = int(value)
                    elif isinstance(value, (np.float64, np.float32)):
                        json_stats[scoring_func][key] = float(value)
                    else:
                        json_stats[scoring_func][key] = value
            
            json.dump(json_stats, f, indent=2)
        print(f"✅ Statistics saved: {stats_file}")
        
        # Save experimental data reference
        exp_data_file = self.results_dir / 'experimental_data.csv'
        exp_df = pd.DataFrame([
            {'ligand_name': name, 'experimental_ec50_um': data['ec50'], 'experimental_ln_ec50': data['ln_ec50']}
            for name, data in self.experimental_data.items()
        ])
        exp_df.to_csv(exp_data_file, index=False)
        print(f"✅ Experimental data saved: {exp_data_file}")
        
        print("\n🎉 All results saved successfully!")

    def run_complete_benchmark(self):
        """Run the complete benchmarking pipeline"""
        
        print("🚀 Starting Complete EC50 Benchmark Pipeline")
        print("=" * 60)
        
        try:
            # Step 1: Run docking for all ligands and scoring functions
            self.run_comprehensive_benchmark()
            
            # Step 2: Calculate statistics
            self.calculate_benchmark_statistics()
            
            # Step 3: Generate plots
            self.create_comprehensive_plots()
            
            # Step 4: Generate report
            self.generate_comprehensive_report()
            
            # Step 5: Save results
            self.save_results()
            
            print("\n🎉 Complete Benchmark Pipeline Finished Successfully!")
            print(f"📂 All results saved in: {self.results_dir}")
            
        except Exception as e:
            print(f"\n💥 Benchmark pipeline failed: {str(e)}")
            import traceback
            traceback.print_exc()

def main():
    """Main function to run the benchmark"""
    
    print("🧪 PandaDock EC50 Benchmarking Suite")
    print("=" * 50)
    
    # Initialize benchmark suite
    benchmark = EC50BenchmarkSuite()
    
    # Run complete benchmark
    benchmark.run_complete_benchmark()

if __name__ == "__main__":
    main()