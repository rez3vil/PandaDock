#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scoring Function Performance Analysis
====================================

Analyzes the performance of different PandaDock scoring functions
to evaluate improvements and score quality.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def analyze_scoring_performance():
    """Analyze performance of all scoring functions"""
    
    print("🔬 PandaDock Scoring Function Performance Analysis")
    print("=" * 55)
    
    # Load results from different scoring functions
    results = {}
    scoring_functions = ['pandacore', 'pandaml', 'pandaphysics']
    
    for scoring_func in scoring_functions:
        csv_file = f"test_{scoring_func}/poses/poses_summary.csv"
        if Path(csv_file).exists():
            df = pd.read_csv(csv_file)
            results[scoring_func] = df
            print(f"✅ Loaded {scoring_func}: {len(df)} poses")
        else:
            print(f"❌ Missing {scoring_func} results")
    
    if not results:
        print("❌ No results found!")
        return
    
    # Performance analysis
    print(f"\n📊 Performance Analysis:")
    print("=" * 25)
    
    analysis_results = {}
    
    for name, df in results.items():
        print(f"\n🎯 {name.upper()} Performance:")
        print("-" * 25)
        
        # Extract key metrics
        scores = df['Score'].values
        energies = df['Energy'].values
        binding_affinities = df['Binding_Affinity'].values
        ic50s = df['IC50_uM'].values
        ec50s = df['EC50_uM'].values
        
        # Calculate statistics
        stats = {
            'score_range': f"{scores.min():.3f} to {scores.max():.3f}",
            'score_std': scores.std(),
            'energy_range': f"{energies.min():.2f} to {energies.max():.2f} kcal/mol",
            'energy_std': energies.std(),
            'affinity_range': f"{binding_affinities.min():.2f} to {binding_affinities.max():.2f} kcal/mol",
            'affinity_std': binding_affinities.std(),
            'unique_ic50s': len(np.unique(ic50s)),
            'unique_ec50s': len(np.unique(ec50s)),
            'best_binding_affinity': binding_affinities.max(),
            'best_ic50': ic50s.min(),
            'best_ec50': ec50s.min()
        }
        
        analysis_results[name] = stats
        
        print(f"Score range:      {stats['score_range']}")
        print(f"Score diversity:  σ = {stats['score_std']:.4f}")
        print(f"Energy range:     {stats['energy_range']}")
        print(f"Energy diversity: σ = {stats['energy_std']:.2f}")
        print(f"Affinity range:   {stats['affinity_range']}")
        print(f"Affinity diversity: σ = {stats['affinity_std']:.3f}")
        print(f"Unique IC50 values: {stats['unique_ic50s']}/{len(ic50s)}")
        print(f"Unique EC50 values: {stats['unique_ec50s']}/{len(ec50s)}")
        print(f"Best binding affinity: {stats['best_binding_affinity']:.3f} kcal/mol")
        print(f"Best IC50: {stats['best_ic50']:.1e} μM")
        print(f"Best EC50: {stats['best_ec50']:.1e} μM")
    
    # Comparison analysis
    print(f"\n🥇 Comparative Performance:")
    print("=" * 30)
    
    # Score diversity ranking
    score_diversity = {name: stats['score_std'] for name, stats in analysis_results.items()}
    best_score_diversity = max(score_diversity, key=score_diversity.get)
    
    # Energy diversity ranking
    energy_diversity = {name: stats['energy_std'] for name, stats in analysis_results.items()}
    best_energy_diversity = max(energy_diversity, key=energy_diversity.get)
    
    # IC50 diversity ranking
    ic50_diversity = {name: stats['unique_ic50s'] for name, stats in analysis_results.items()}
    best_ic50_diversity = max(ic50_diversity, key=ic50_diversity.get)
    
    # Best binding affinity
    best_affinity = {name: stats['best_binding_affinity'] for name, stats in analysis_results.items()}
    strongest_binder = max(best_affinity, key=best_affinity.get)
    
    print(f"🏆 Best Score Diversity:    {best_score_diversity.upper()} (σ = {score_diversity[best_score_diversity]:.4f})")
    print(f"🏆 Best Energy Diversity:   {best_energy_diversity.upper()} (σ = {energy_diversity[best_energy_diversity]:.2f})")
    print(f"🏆 Best IC50 Diversity:     {best_ic50_diversity.upper()} ({ic50_diversity[best_ic50_diversity]} unique)")
    print(f"🏆 Strongest Binding:       {strongest_binder.upper()} ({best_affinity[strongest_binder]:.3f} kcal/mol)")
    
    # Check for identical values issue
    print(f"\n⚠️  Issue Analysis:")
    print("=" * 20)
    
    for name, df in results.items():
        ic50s = df['IC50_uM'].values
        ec50s = df['EC50_uM'].values
        
        # Check for identical IC50/EC50 values
        unique_ic50_ratio = len(np.unique(ic50s)) / len(ic50s)
        unique_ec50_ratio = len(np.unique(ec50s)) / len(ec50s)
        
        if unique_ic50_ratio < 0.5:
            print(f"⚠️  {name.upper()}: Low IC50 diversity ({unique_ic50_ratio:.1%} unique)")
        else:
            print(f"✅ {name.upper()}: Good IC50 diversity ({unique_ic50_ratio:.1%} unique)")
        
        if unique_ec50_ratio < 0.5:
            print(f"⚠️  {name.upper()}: Low EC50 diversity ({unique_ec50_ratio:.1%} unique)")
        else:
            print(f"✅ {name.upper()}: Good EC50 diversity ({unique_ec50_ratio:.1%} unique)")
    
    # Create comparison visualization
    create_performance_visualization(results, analysis_results)
    
    # Overall assessment
    print(f"\n🎯 Overall Assessment:")
    print("=" * 22)
    
    # Calculate overall performance score
    performance_scores = {}
    for name in analysis_results.keys():
        score = (
            score_diversity[name] * 1000 +  # Score diversity (scaled)
            energy_diversity[name] * 10 +   # Energy diversity  
            ic50_diversity[name] +           # IC50 diversity
            best_affinity[name]              # Best binding strength
        )
        performance_scores[name] = score
    
    ranked_performance = sorted(performance_scores.items(), key=lambda x: x[1], reverse=True)
    
    print("Performance Ranking (higher = better):")
    for i, (name, score) in enumerate(ranked_performance, 1):
        print(f"{i}. {name.upper():12s}: {score:.2f}")
    
    return analysis_results

def create_performance_visualization(results, analysis_results):
    """Create visualization comparing scoring function performance"""
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Score distributions
    ax = axes[0, 0]
    for name, df in results.items():
        ax.hist(df['Score'], bins=20, alpha=0.6, label=name.upper(), density=True)
    ax.set_xlabel('Docking Score')
    ax.set_ylabel('Density')
    ax.set_title('Score Distributions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Energy distributions
    ax = axes[0, 1]
    for name, df in results.items():
        ax.hist(df['Energy'], bins=20, alpha=0.6, label=name.upper(), density=True)
    ax.set_xlabel('Energy (kcal/mol)')
    ax.set_ylabel('Density')
    ax.set_title('Energy Distributions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Binding affinity distributions
    ax = axes[0, 2]
    for name, df in results.items():
        ax.hist(df['Binding_Affinity'], bins=20, alpha=0.6, label=name.upper(), density=True)
    ax.set_xlabel('Binding Affinity (kcal/mol)')
    ax.set_ylabel('Density')
    ax.set_title('Binding Affinity Distributions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: IC50 distributions (log scale)
    ax = axes[1, 0]
    for name, df in results.items():
        ic50s = df['IC50_uM'].values
        ax.hist(np.log10(ic50s), bins=20, alpha=0.6, label=name.upper(), density=True)
    ax.set_xlabel('log10(IC50 μM)')
    ax.set_ylabel('Density')
    ax.set_title('IC50 Distributions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 5: Performance metrics comparison
    ax = axes[1, 1]
    metrics = ['score_std', 'energy_std', 'affinity_std']
    x_pos = np.arange(len(metrics))
    width = 0.25
    
    for i, (name, stats) in enumerate(analysis_results.items()):
        values = [stats['score_std']*1000, stats['energy_std'], stats['affinity_std']*10]
        ax.bar(x_pos + i*width, values, width, label=name.upper())
    
    ax.set_xlabel('Metric')
    ax.set_ylabel('Standard Deviation (scaled)')
    ax.set_title('Diversity Metrics Comparison')
    ax.set_xticks(x_pos + width)
    ax.set_xticklabels(['Score (×1000)', 'Energy', 'Affinity (×10)'])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 6: Unique value counts
    ax = axes[1, 2]
    scoring_funcs = list(analysis_results.keys())
    ic50_counts = [analysis_results[name]['unique_ic50s'] for name in scoring_funcs]
    ec50_counts = [analysis_results[name]['unique_ec50s'] for name in scoring_funcs]
    
    x_pos = np.arange(len(scoring_funcs))
    width = 0.35
    
    ax.bar(x_pos - width/2, ic50_counts, width, label='IC50', alpha=0.8)
    ax.bar(x_pos + width/2, ec50_counts, width, label='EC50', alpha=0.8)
    
    ax.set_xlabel('Scoring Function')
    ax.set_ylabel('Number of Unique Values')
    ax.set_title('Value Diversity Analysis')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([name.upper() for name in scoring_funcs])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    output_file = "scoring_function_performance_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    
    print(f"📊 Performance visualization saved: {output_file}")
    plt.show()

def test_ic50_calculation_fix():
    """Test if the IC50 calculation fix resolved the identical values issue"""
    
    print(f"\n🔧 Testing IC50 Calculation Fix:")
    print("=" * 35)
    
    # Load results and check for diversity
    scoring_functions = ['pandacore', 'pandaml', 'pandaphysics']
    
    for scoring_func in scoring_functions:
        csv_file = f"test_{scoring_func}/poses/poses_summary.csv"
        if Path(csv_file).exists():
            df = pd.read_csv(csv_file)
            
            ic50s = df['IC50_uM'].values
            ec50s = df['EC50_uM'].values
            
            unique_ic50s = len(np.unique(ic50s))
            unique_ec50s = len(np.unique(ec50s))
            
            print(f"\n{scoring_func.upper()}:")
            print(f"  IC50 diversity: {unique_ic50s}/{len(ic50s)} unique ({unique_ic50s/len(ic50s):.1%})")
            print(f"  EC50 diversity: {unique_ec50s}/{len(ec50s)} unique ({unique_ec50s/len(ec50s):.1%})")
            
            if unique_ic50s > len(ic50s) * 0.7:
                print(f"  ✅ IC50 diversity GOOD")
            else:
                print(f"  ⚠️  IC50 diversity LOW")
            
            if unique_ec50s > len(ec50s) * 0.7:
                print(f"  ✅ EC50 diversity GOOD")
            else:
                print(f"  ⚠️  EC50 diversity LOW")

if __name__ == "__main__":
    # Run comprehensive analysis
    analysis_results = analyze_scoring_performance()
    
    # Test IC50 fix
    test_ic50_calculation_fix()
    
    print(f"\n🎉 Scoring Function Analysis Complete!")
    print(f"📊 Check 'scoring_function_performance_analysis.png' for detailed visualizations")