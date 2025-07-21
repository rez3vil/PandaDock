#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the complete EC50 benchmark for all ligands

This script validates PandaDock's EC50 predictions against experimental GABA_A receptor data
for all 30 phenolic compounds. Due to the computational intensity, you can modify the 
script to run subsets or specific scoring functions.

Usage:
    python run_full_ec50_benchmark.py

Options:
    - Set QUICK_TEST = True for fast testing with fewer poses
    - Set SCORING_FUNCTIONS to test specific scoring algorithms
    - Set MAX_LIGANDS to limit the number of ligands tested
"""

import sys
sys.path.append('/Users/pritam/PandaDock')

from benchmark_ec50_validation import EC50BenchmarkSuite
import time

def run_full_benchmark():
    """Run the complete EC50 benchmark for all ligands"""
    
    print("🧪 PandaDock Complete EC50 Benchmark")
    print("=" * 60)
    print("🎯 Target: GABA_A Beta2-Alpha1 Receptor")
    print("🧬 Ligands: 30 Phenolic Compounds")
    print("⚖️  Scoring Functions: PandaCore, PandaML, PANDAPhysics")
    print("")
    
    # Configuration options
    QUICK_TEST = True  # Set to False for full benchmark with more poses
    SCORING_FUNCTIONS = ['pandacore', 'pandaml', 'pandaphysics']  # Can limit to ['pandaphysics'] for best performer
    MAX_LIGANDS = None  # Set to a number to limit ligands (e.g., 10 for testing)
    
    # Initialize benchmark suite
    benchmark = EC50BenchmarkSuite(results_dir="full_ec50_benchmark_results")
    
    print(f"⚙️  Configuration:")
    print(f"   Quick Test Mode: {QUICK_TEST}")
    print(f"   Scoring Functions: {SCORING_FUNCTIONS}")
    print(f"   Max Ligands: {MAX_LIGANDS if MAX_LIGANDS else 'All (30)'}")
    print("")
    
    # Optionally limit the dataset for testing
    if MAX_LIGANDS:
        original_data = benchmark.experimental_data.copy()
        ligand_names = list(original_data.keys())[:MAX_LIGANDS]
        benchmark.experimental_data = {k: v for k, v in original_data.items() if k in ligand_names}
        print(f"🔬 Limited to {len(benchmark.experimental_data)} ligands for testing")
    
    # Show experimental data range
    ec50_values = [data['ec50'] for data in benchmark.experimental_data.values()]
    print(f"📊 Experimental EC50 Range: {min(ec50_values):.2f} - {max(ec50_values):.1f} μM")
    print(f"📊 Log Range: {min(ec50_values)/max(ec50_values):.2e} (nearly 5000-fold difference)")
    print("")
    
    # Time estimation
    num_ligands = len(benchmark.experimental_data)
    num_scoring = len(SCORING_FUNCTIONS)
    poses_per_ligand = 3 if QUICK_TEST else 10
    estimated_time = num_ligands * num_scoring * 2  # ~2 minutes per ligand per scoring function
    
    print(f"⏱️  Estimated Runtime: {estimated_time:.0f} minutes ({estimated_time/60:.1f} hours)")
    print(f"   ({num_ligands} ligands × {num_scoring} scoring functions × {poses_per_ligand} poses)")
    print("")
    
    # Ask for confirmation for long runs
    if estimated_time > 30:
        response = input("⚠️  This will take a while. Continue? (y/N): ")
        if response.lower() != 'y':
            print("❌ Benchmark cancelled")
            return
    
    start_time = time.time()
    
    # Temporarily override scoring functions and poses
    original_run_method = benchmark.run_comprehensive_benchmark
    
    def custom_benchmark():
        """Modified benchmark method with custom settings"""
        print("🚀 Starting Custom EC50 Benchmark")
        print("=" * 60)
        
        for scoring_func in SCORING_FUNCTIONS:
            print(f"\n📊 Testing {scoring_func.upper()} scoring function")
            print("-" * 40)
            
            scoring_results = []
            
            for i, ligand_name in enumerate(benchmark.experimental_data.keys(), 1):
                print(f"[{i:2d}/{len(benchmark.experimental_data)}] ", end="")
                
                num_poses = poses_per_ligand
                result = benchmark.run_docking(ligand_name, scoring_func, num_poses=num_poses)
                if result:
                    scoring_results.append(result)
                    
                    # Quick feedback
                    exp_ec50 = result['experimental_ec50']
                    pred_ec50 = result.get('predicted_ec50_um', 'N/A')
                    if pred_ec50 != 'N/A':
                        error_fold = pred_ec50 / exp_ec50 if exp_ec50 > 0 else float('inf')
                        print(f"    Exp: {exp_ec50:6.2f} μM, Pred: {pred_ec50:8.1e} μM (×{error_fold:5.1f})")
                    else:
                        print(f"    Exp: {exp_ec50:6.2f} μM, Pred: N/A")
            
            benchmark.docking_results[scoring_func] = scoring_results
            success_rate = len(scoring_results) / len(benchmark.experimental_data) * 100
            print(f"\n✅ {scoring_func.upper()} completed: {len(scoring_results)}/{len(benchmark.experimental_data)} successful ({success_rate:.1f}%)")
    
    # Replace the benchmark method
    benchmark.run_comprehensive_benchmark = custom_benchmark
    
    try:
        # Run complete benchmark pipeline
        benchmark.run_complete_benchmark()
        
        end_time = time.time()
        runtime = end_time - start_time
        
        print(f"\n⏱️  Total Runtime: {runtime/60:.1f} minutes ({runtime/3600:.2f} hours)")
        print(f"📂 Results saved in: {benchmark.results_dir}")
        
        # Print key results
        print(f"\n📊 Key Results Summary:")
        print("=" * 40)
        
        for scoring_func in SCORING_FUNCTIONS:
            if scoring_func in benchmark.benchmark_stats:
                stats = benchmark.benchmark_stats[scoring_func]['stats']
                r2 = stats.get('ln_ec50_r2', 0)
                r = stats.get('ln_ec50_pearson_r', 0)
                p_val = stats.get('ln_ec50_pearson_p', 1)
                success = stats.get('success_rate', 0)
                
                significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
                
                print(f"{scoring_func.upper():12s}: R²={r2:.3f}, r={r:+.3f} ({significance}), Success={success:.1%}")
        
        print(f"\n🎉 Complete EC50 Benchmark Finished Successfully!")
        print(f"📈 Best performer: PANDAPhysics (based on subset results)")
        print(f"📁 View plots and detailed report in: {benchmark.results_dir}/")
        
    except KeyboardInterrupt:
        print(f"\n⚠️  Benchmark interrupted by user")
        print(f"📂 Partial results may be available in: {benchmark.results_dir}")
    except Exception as e:
        print(f"\n💥 Benchmark failed: {str(e)}")
        import traceback
        traceback.print_exc()

def show_experimental_data():
    """Display the experimental data for reference"""
    
    benchmark = EC50BenchmarkSuite()
    
    print("📋 Experimental GABA_A EC50 Data")
    print("=" * 50)
    print(f"{'Ligand Name':<25s} {'EC50 (μM)':<10s} {'ln(EC50)':<10s}")
    print("-" * 50)
    
    # Sort by EC50 value
    sorted_data = sorted(benchmark.experimental_data.items(), key=lambda x: x[1]['ec50'])
    
    for ligand_name, data in sorted_data:
        print(f"{ligand_name:<25s} {data['ec50']:>8.2f}  {data['ln_ec50']:>8.3f}")
    
    print("-" * 50)
    print(f"{'Total Ligands:':<25s} {len(benchmark.experimental_data):>8d}")
    
    ec50_values = [data['ec50'] for data in benchmark.experimental_data.values()]
    print(f"{'Range:':<25s} {min(ec50_values):>8.2f} - {max(ec50_values):8.1f} μM")
    print(f"{'Fold Difference:':<25s} {max(ec50_values)/min(ec50_values):>8.0f}×")

def main():
    """Main function with menu options"""
    
    print("🧪 PandaDock EC50 Benchmark Suite")
    print("=" * 50)
    print("1. Show experimental data")
    print("2. Run subset benchmark (8 ligands, ~20 min)")
    print("3. Run full benchmark (30 ligands, ~2-4 hours)")
    print("4. Run quick test (5 ligands, ~5 min)")
    print("0. Exit")
    
    while True:
        try:
            choice = input("\nSelect option (0-4): ").strip()
            
            if choice == '0':
                print("👋 Goodbye!")
                break
            elif choice == '1':
                show_experimental_data()
            elif choice == '2':
                # Run subset from previous test
                exec(open('run_benchmark_subset.py').read())
            elif choice == '3':
                run_full_benchmark()
            elif choice == '4':
                # Quick test with 5 ligands
                print("🚀 Quick Test Mode (5 ligands)")
                benchmark = EC50BenchmarkSuite(results_dir="quick_test_results")
                test_ligands = ['1n', '1k', '1a', '8', 'phenol']  # Representative range
                original_data = benchmark.experimental_data.copy()
                benchmark.experimental_data = {k: v for k, v in original_data.items() if k in test_ligands}
                benchmark.run_complete_benchmark()
            else:
                print("❌ Invalid option. Please choose 0-4.")
                
        except KeyboardInterrupt:
            print("\n👋 Goodbye!")
            break
        except Exception as e:
            print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    main()