#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick test of the benchmarking system with a few ligands
"""

import sys
sys.path.append('/Users/pritam/PandaDock')

from benchmark_ec50_validation import EC50BenchmarkSuite
import pandas as pd

def test_sample_ligands():
    """Test with just a few ligands to verify the system works"""
    
    print("🧪 Testing Benchmark System with Sample Ligands")
    print("=" * 50)
    
    # Initialize benchmark suite
    benchmark = EC50BenchmarkSuite()
    
    # Test with just a few ligands
    test_ligands = ['1a', '1b', '1c', '8', 'phenol']
    
    print(f"Testing with {len(test_ligands)} ligands: {test_ligands}")
    
    # Check SDF file availability
    print("\n📁 Checking SDF file availability:")
    for ligand in test_ligands:
        sdf_path = benchmark.get_sdf_file_path(ligand)
        if sdf_path:
            print(f"✅ {ligand}: {sdf_path.name}")
        else:
            print(f"❌ {ligand}: File not found")
    
    # Test docking for one ligand with each scoring function
    print("\n🔬 Testing docking for one ligand (1a):")
    
    for scoring_func in ['pandacore', 'pandaml']:  # Test just two scoring functions
        print(f"\nTesting {scoring_func}...")
        result = benchmark.run_docking('1a', scoring_func, num_poses=3)
        if result:
            print(f"✅ Success: Exp EC50={result['experimental_ec50']:.2f} μM, "
                  f"Pred EC50={result.get('predicted_ec50_um', 'N/A')}")
        else:
            print(f"❌ Failed")
    
    print("\n✅ Sample test completed!")

if __name__ == "__main__":
    test_sample_ligands()