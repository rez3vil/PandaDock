#!/usr/bin/env python3
"""
Benchmark Results Analyzer for PandaDock

This script provides detailed analysis of benchmark results including:
- Performance statistics
- Binding affinity correlation analysis  
- Algorithm/device comparisons
- Visualization generation

Usage:
    python analyze_benchmark_results.py --results_dir benchmark_results_20231201_120000_comprehensive
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
import pandas as pd
import numpy as np
from datetime import datetime

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: matplotlib/seaborn not available. Plots will be skipped.")

try:
    from scipy import stats
    from sklearn.metrics import mean_squared_error, r2_score
    HAS_STATS = True
except ImportError:
    HAS_STATS = False
    print("Warning: scipy/sklearn not available. Advanced statistics will be limited.")


class BenchmarkAnalyzer:
    """Analyzer for PandaDock benchmark results."""
    
    def __init__(self, results_dir: str):
        """Initialize analyzer with results directory."""
        self.results_dir = Path(results_dir)
        self.results_file = self.results_dir / "benchmark_results.json"
        self.analysis_file = self.results_dir / "benchmark_analysis.json"
        
        if not self.results_dir.exists():
            raise FileNotFoundError(f"Results directory not found: {results_dir}")
        if not self.results_file.exists():
            raise FileNotFoundError(f"Results file not found: {self.results_file}")
        
        # Load results
        with open(self.results_file, 'r') as f:
            self.results = json.load(f)
        
        # Load analysis if available
        self.analysis = {}
        if self.analysis_file.exists():
            with open(self.analysis_file, 'r') as f:
                self.analysis = json.load(f)
    
    def print_summary(self):
        """Print a comprehensive summary of benchmark results."""
        print("=" * 80)
        print("PANDADOCK BENCHMARK RESULTS SUMMARY")
        print("=" * 80)
        print(f"Results directory: {self.results_dir}")
        print(f"Total benchmark runs: {len(self.results)}")
        print()
        
        # Basic statistics
        successful = [r for r in self.results if r['success']]
        failed = [r for r in self.results if not r['success']]
        
        print("BASIC STATISTICS")
        print("-" * 20)
        print(f"Successful runs: {len(successful)} ({len(successful)/len(self.results)*100:.1f}%)")
        print(f"Failed runs: {len(failed)} ({len(failed)/len(self.results)*100:.1f}%)")
        
        if successful:
            runtimes = [r['total_time'] for r in successful]
            print(f"Average runtime: {np.mean(runtimes):.1f} ± {np.std(runtimes):.1f} seconds")
            print(f"Median runtime: {np.median(runtimes):.1f} seconds")
            print(f"Runtime range: {min(runtimes):.1f} - {max(runtimes):.1f} seconds")
        print()
        
        # Device comparison
        self._print_device_comparison(successful)
        
        # Algorithm comparison
        self._print_algorithm_comparison(successful)
        
        # Scoring function comparison
        self._print_scoring_comparison(successful)
        
        # Binding affinity analysis
        self._print_binding_affinity_analysis(successful)
        
        # Top performers
        self._print_top_performers(successful)
        
        # Failure analysis
        if failed:
            self._print_failure_analysis(failed)
    
    def _print_device_comparison(self, successful_results):
        """Print device performance comparison."""
        devices = set(r['device'] for r in successful_results)
        
        if len(devices) < 2:
            return
        
        print("DEVICE PERFORMANCE COMPARISON")
        print("-" * 30)
        
        device_stats = {}
        for device in devices:
            device_results = [r for r in successful_results if r['device'] == device]
            runtimes = [r['total_time'] for r in device_results]
            
            device_stats[device] = {
                'count': len(device_results),
                'avg_runtime': np.mean(runtimes),
                'std_runtime': np.std(runtimes),
                'median_runtime': np.median(runtimes)
            }
            
            print(f"{device}:")
            print(f"  Runs: {len(device_results)}")
            print(f"  Avg runtime: {np.mean(runtimes):.1f} ± {np.std(runtimes):.1f}s")
            print(f"  Median runtime: {np.median(runtimes):.1f}s")
        
        # Calculate speedup if GPU and CPU both present
        if 'CPU' in device_stats and 'GPU' in device_stats:
            speedup = device_stats['CPU']['avg_runtime'] / device_stats['GPU']['avg_runtime']
            print(f"\nGPU speedup: {speedup:.1f}x")
        
        print()
    
    def _print_algorithm_comparison(self, successful_results):
        """Print algorithm performance comparison."""
        algorithms = set(r['algorithm'] for r in successful_results)
        
        if len(algorithms) < 2:
            return
        
        print("ALGORITHM PERFORMANCE COMPARISON")
        print("-" * 32)
        
        for algorithm in sorted(algorithms):
            alg_results = [r for r in successful_results if r['algorithm'] == algorithm]
            runtimes = [r['total_time'] for r in alg_results]
            scores = [r['best_score'] for r in alg_results if r['best_score'] != float('inf')]
            
            print(f"{algorithm.upper()}:")
            print(f"  Runs: {len(alg_results)}")
            print(f"  Avg runtime: {np.mean(runtimes):.1f} ± {np.std(runtimes):.1f}s")
            if scores:
                print(f"  Avg best score: {np.mean(scores):.3f} ± {np.std(scores):.3f}")
                print(f"  Best overall score: {min(scores):.3f}")
        
        print()
    
    def _print_scoring_comparison(self, successful_results):
        """Print scoring function comparison."""
        scoring_functions = set(r['scoring_function'] for r in successful_results)
        
        if len(scoring_functions) < 2:
            return
        
        print("SCORING FUNCTION COMPARISON")
        print("-" * 27)
        
        for scoring in sorted(scoring_functions):
            scoring_results = [r for r in successful_results if r['scoring_function'] == scoring]
            runtimes = [r['total_time'] for r in scoring_results]
            scores = [r['best_score'] for r in scoring_results if r['best_score'] != float('inf')]
            
            print(f"{scoring.upper()}:")
            print(f"  Runs: {len(scoring_results)}")
            print(f"  Avg runtime: {np.mean(runtimes):.1f} ± {np.std(runtimes):.1f}s")
            if scores:
                print(f"  Avg best score: {np.mean(scores):.3f} ± {np.std(scores):.3f}")
        
        print()
    
    def _print_binding_affinity_analysis(self, successful_results):
        """Print binding affinity correlation analysis."""
        # Filter results with both predicted and experimental values
        valid_results = [
            r for r in successful_results 
            if r.get('predicted_kd') is not None and r.get('experimental_value') is not None
        ]
        
        if len(valid_results) < 10:
            print("BINDING AFFINITY ANALYSIS")
            print("-" * 25)
            print(f"Insufficient data for correlation analysis ({len(valid_results)} valid results)")
            print("Need at least 10 results with both predicted and experimental values")
            print()
            return
        
        print("BINDING AFFINITY ANALYSIS")
        print("-" * 25)
        print(f"Valid results for correlation: {len(valid_results)}")
        
        # Extract values
        experimental = [r['experimental_value'] for r in valid_results]
        predicted_kd = [r['predicted_kd'] for r in valid_results]
        
        # Convert to log scale
        log_exp = np.log10(experimental)
        log_pred = np.log10(predicted_kd)
        
        # Calculate correlation
        if HAS_STATS:
            r_value, p_value = stats.pearsonr(log_exp, log_pred)
            rmse = np.sqrt(mean_squared_error(log_exp, log_pred))
            r2 = r2_score(log_exp, log_pred)
            
            print(f"Pearson correlation (R): {r_value:.3f}")
            print(f"P-value: {p_value:.2e}")
            print(f"R-squared: {r2:.3f}")
            print(f"RMSE (log scale): {rmse:.3f}")
            
            # Interpretation
            if abs(r_value) > 0.7:
                interpretation = "Strong correlation"
            elif abs(r_value) > 0.5:
                interpretation = "Moderate correlation"
            elif abs(r_value) > 0.3:
                interpretation = "Weak correlation"
            else:
                interpretation = "Very weak correlation"
            
            print(f"Interpretation: {interpretation}")
        
        # Experimental value distribution
        exp_types = {}
        for r in valid_results:
            exp_type = r.get('experimental_type', 'Unknown')
            if exp_type not in exp_types:
                exp_types[exp_type] = 0
            exp_types[exp_type] += 1
        
        print(f"Experimental data types: {dict(exp_types)}")
        print()
    
    def _print_top_performers(self, successful_results):
        """Print top performing complexes."""
        print("TOP PERFORMERS")
        print("-" * 14)
        
        # Sort by best score
        sorted_results = sorted(successful_results, key=lambda r: r['best_score'])
        
        print("Best docking scores:")
        for i, result in enumerate(sorted_results[:10], 1):
            print(f"  {i:2d}. {result['pdb_code']} ({result['algorithm']}/{result['device']}): "
                  f"{result['best_score']:.3f}")
        
        # Fastest runs
        sorted_by_time = sorted(successful_results, key=lambda r: r['total_time'])
        
        print("\nFastest runs:")
        for i, result in enumerate(sorted_by_time[:10], 1):
            print(f"  {i:2d}. {result['pdb_code']} ({result['algorithm']}/{result['device']}): "
                  f"{result['total_time']:.1f}s")
        
        print()
    
    def _print_failure_analysis(self, failed_results):
        """Print analysis of failed runs."""
        print("FAILURE ANALYSIS")
        print("-" * 16)
        print(f"Total failed runs: {len(failed_results)}")
        
        # Failure reasons
        error_types = {}
        for result in failed_results:
            error = result.get('error_message', 'Unknown error')
            # Categorize errors
            if 'timeout' in error.lower():
                category = 'Timeout'
            elif 'file' in error.lower() or 'not found' in error.lower():
                category = 'File not found'
            elif 'gpu' in error.lower() or 'cuda' in error.lower():
                category = 'GPU error'
            elif 'memory' in error.lower():
                category = 'Memory error'
            else:
                category = 'Other'
            
            if category not in error_types:
                error_types[category] = 0
            error_types[category] += 1
        
        print("Failure categories:")
        for category, count in sorted(error_types.items(), key=lambda x: x[1], reverse=True):
            print(f"  {category}: {count} ({count/len(failed_results)*100:.1f}%)")
        
        # Failures by device/algorithm
        failure_by_device = {}
        failure_by_algorithm = {}
        
        for result in failed_results:
            device = result['device']
            algorithm = result['algorithm']
            
            if device not in failure_by_device:
                failure_by_device[device] = 0
            failure_by_device[device] += 1
            
            if algorithm not in failure_by_algorithm:
                failure_by_algorithm[algorithm] = 0
            failure_by_algorithm[algorithm] += 1
        
        if failure_by_device:
            print(f"Failures by device: {dict(failure_by_device)}")
        if failure_by_algorithm:
            print(f"Failures by algorithm: {dict(failure_by_algorithm)}")
        
        print()
    
    def generate_detailed_report(self, output_file: Optional[str] = None):
        """Generate a detailed HTML report."""
        if output_file is None:
            output_file = self.results_dir / "detailed_report.html"
        
        # This would generate an HTML report with embedded plots and tables
        # For now, we'll create a markdown report
        self._generate_markdown_report(output_file.with_suffix('.md'))
    
    def _generate_markdown_report(self, output_file: Path):
        """Generate a markdown report."""
        successful = [r for r in self.results if r['success']]
        
        with open(output_file, 'w') as f:
            f.write("# PandaDock Benchmark Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Results directory: `{self.results_dir}`\n\n")
            
            # Summary table
            f.write("## Summary\n\n")
            f.write("| Metric | Value |\n")
            f.write("|--------|-------|\n")
            f.write(f"| Total runs | {len(self.results)} |\n")
            f.write(f"| Successful runs | {len(successful)} ({len(successful)/len(self.results)*100:.1f}%) |\n")
            
            if successful:
                runtimes = [r['total_time'] for r in successful]
                f.write(f"| Average runtime | {np.mean(runtimes):.1f} ± {np.std(runtimes):.1f} seconds |\n")
                f.write(f"| Median runtime | {np.median(runtimes):.1f} seconds |\n")
            
            f.write("\n")
            
            # Device comparison
            devices = set(r['device'] for r in successful)
            if len(devices) > 1:
                f.write("## Device Performance\n\n")
                f.write("| Device | Runs | Avg Runtime (s) | Std Runtime (s) |\n")
                f.write("|--------|------|-----------------|------------------|\n")
                
                for device in sorted(devices):
                    device_results = [r for r in successful if r['device'] == device]
                    runtimes = [r['total_time'] for r in device_results]
                    f.write(f"| {device} | {len(device_results)} | {np.mean(runtimes):.1f} | {np.std(runtimes):.1f} |\n")
                
                f.write("\n")
            
            # Algorithm comparison
            algorithms = set(r['algorithm'] for r in successful)
            if len(algorithms) > 1:
                f.write("## Algorithm Performance\n\n")
                f.write("| Algorithm | Runs | Avg Runtime (s) | Avg Best Score |\n")
                f.write("|-----------|------|-----------------|----------------|\n")
                
                for algorithm in sorted(algorithms):
                    alg_results = [r for r in successful if r['algorithm'] == algorithm]
                    runtimes = [r['total_time'] for r in alg_results]
                    scores = [r['best_score'] for r in alg_results if r['best_score'] != float('inf')]
                    avg_score = np.mean(scores) if scores else 'N/A'
                    f.write(f"| {algorithm} | {len(alg_results)} | {np.mean(runtimes):.1f} | {avg_score if avg_score == 'N/A' else f'{avg_score:.3f}'} |\n")
                
                f.write("\n")
            
            # Top performers
            f.write("## Top 10 Best Scores\n\n")
            sorted_results = sorted(successful, key=lambda r: r['best_score'])[:10]
            f.write("| Rank | PDB Code | Algorithm | Device | Score |\n")
            f.write("|------|----------|-----------|--------|---------|\n")
            
            for i, result in enumerate(sorted_results, 1):
                f.write(f"| {i} | {result['pdb_code']} | {result['algorithm']} | {result['device']} | {result['best_score']:.3f} |\n")
            
            f.write("\n")
            
            # Binding affinity analysis
            valid_affinity = [
                r for r in successful 
                if r.get('predicted_kd') is not None and r.get('experimental_value') is not None
            ]
            
            if len(valid_affinity) >= 10 and HAS_STATS:
                experimental = [r['experimental_value'] for r in valid_affinity]
                predicted = [r['predicted_kd'] for r in valid_affinity]
                
                log_exp = np.log10(experimental)
                log_pred = np.log10(predicted)
                
                r_value, p_value = stats.pearsonr(log_exp, log_pred)
                r2 = r2_score(log_exp, log_pred)
                
                f.write("## Binding Affinity Correlation\n\n")
                f.write(f"- Valid results: {len(valid_affinity)}\n")
                f.write(f"- Pearson correlation (R): {r_value:.3f}\n")
                f.write(f"- R-squared: {r2:.3f}\n")
                f.write(f"- P-value: {p_value:.2e}\n\n")
        
        print(f"Detailed report saved to: {output_file}")
    
    def export_csv_for_analysis(self, output_file: Optional[str] = None):
        """Export results to CSV for external analysis."""
        if output_file is None:
            output_file = self.results_dir / "results_for_analysis.csv"
        
        # Convert results to DataFrame
        df = pd.DataFrame(self.results)
        
        # Add derived columns
        df['success_binary'] = df['success'].astype(int)
        df['log_runtime'] = np.log10(df['total_time'].replace(0, 0.1))  # Avoid log(0)
        
        # Add experimental vs predicted comparison for successful runs
        if 'predicted_kd' in df.columns and 'experimental_value' in df.columns:
            valid_mask = (df['predicted_kd'].notna()) & (df['experimental_value'].notna()) & df['success']
            df.loc[valid_mask, 'log_experimental'] = np.log10(df.loc[valid_mask, 'experimental_value'])
            df.loc[valid_mask, 'log_predicted_kd'] = np.log10(df.loc[valid_mask, 'predicted_kd'])
            df.loc[valid_mask, 'affinity_error'] = (
                df.loc[valid_mask, 'log_predicted_kd'] - df.loc[valid_mask, 'log_experimental']
            )
        
        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"CSV exported to: {output_file}")
        
        return df


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Analyze PandaDock benchmark results",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'results_dir',
        help='Directory containing benchmark results'
    )
    
    parser.add_argument(
        '--summary', 
        action='store_true',
        help='Print summary to console (default)'
    )
    
    parser.add_argument(
        '--report', 
        help='Generate detailed report (HTML/Markdown)'
    )
    
    parser.add_argument(
        '--export_csv', 
        help='Export results to CSV for external analysis'
    )
    
    parser.add_argument(
        '--all', 
        action='store_true',
        help='Generate all outputs (summary, report, CSV)'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        analyzer = BenchmarkAnalyzer(args.results_dir)
        
        # Default to summary if no specific options given
        if not any([args.report, args.export_csv, args.all]):
            args.summary = True
        
        # Generate requested outputs
        if args.summary or args.all:
            analyzer.print_summary()
        
        if args.report or args.all:
            output_file = args.report if args.report else None
            analyzer.generate_detailed_report(output_file)
        
        if args.export_csv or args.all:
            output_file = args.export_csv if args.export_csv else None
            analyzer.export_csv_for_analysis(output_file)
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())