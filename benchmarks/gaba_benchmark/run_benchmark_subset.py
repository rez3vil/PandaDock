#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run benchmark on a subset of ligands to validate the system
"""

import sys
sys.path.append('/Users/pritam/PandaDock')

from benchmark_ec50_validation import EC50BenchmarkSuite

def run_subset_benchmark():
    """Run benchmark on a manageable subset of ligands"""
    
    print("🧪 Running EC50 Benchmark on Subset")
    print("=" * 50)
    
    # Initialize benchmark suite
    benchmark = EC50BenchmarkSuite()
    
    # Select a representative subset of ligands across the EC50 range
    # High EC50: phenol (920), 35ditertbutphenol (94.4), isopropylphenol (39.7)
    # Medium EC50: 24disecbutphenol (32.9), 4iododiisopropylphenol (11.1), 1d (9.1), 8 (8.7)
    # Low EC50: 1a (6.4), 1b (4.2), 1e (3.0), 1j (1.8), 1k (1.0), 1n (0.19)
    
    subset_ligands = [
        'phenol',           # Very high EC50 (920)
        'isopropylphenol',  # High EC50 (39.7)
        '4iododiisopropylphenol',  # Medium-high EC50 (11.1)
        '8',               # Medium EC50 (8.7)
        '1a',              # Medium-low EC50 (6.4)
        '1e',              # Low EC50 (3.0)
        '1k',              # Low EC50 (1.0)
        '1n'               # Very low EC50 (0.19)
    ]
    
    print(f"🎯 Selected {len(subset_ligands)} representative ligands:")
    for ligand in subset_ligands:
        exp_data = benchmark.experimental_data[ligand]
        print(f"   {ligand}: {exp_data['ec50']:.2f} μM")
    
    # Temporarily override the experimental data to only include subset
    original_data = benchmark.experimental_data.copy()
    benchmark.experimental_data = {k: v for k, v in original_data.items() if k in subset_ligands}
    
    # Run complete benchmark
    benchmark.run_complete_benchmark()
    
    print(f"\n🎉 Subset Benchmark Complete!")
    print(f"📂 Results saved in: {benchmark.results_dir}")

if __name__ == "__main__":
    run_subset_benchmark()