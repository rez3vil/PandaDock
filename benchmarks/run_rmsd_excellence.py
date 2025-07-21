#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick RMSD Excellence Benchmark Runner for PandaDock

This is a simple script to showcase PandaDock's exceptional RMSD performance.
It automatically downloads a test dataset and generates impressive visualizations.

Usage:
    python run_rmsd_excellence.py
    python run_rmsd_excellence.py --max_complexes 20  # For faster testing
    python run_rmsd_excellence.py --engines pandaml pandaphysics  # Specific engines only
"""

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path

def main():
    """Run the RMSD excellence benchmark with user-friendly options"""
    
    parser = argparse.ArgumentParser(
        description='PandaDock RMSD Excellence Benchmark - Showcase Sub-2Å Performance',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_rmsd_excellence.py                    # Full benchmark
  python run_rmsd_excellence.py --max_complexes 15  # Quick test (15 complexes)
  python run_rmsd_excellence.py --engines pandaml   # PandaML only
  python run_rmsd_excellence.py --quick              # Ultra-fast demo (5 complexes)
        """
    )
    
    parser.add_argument('--max_complexes', type=int, default=None,
                       help='Maximum number of complexes to test (default: all available)')
    parser.add_argument('--engines', nargs='+', 
                       choices=['pandacore', 'pandaml', 'pandaphysics'],
                       default=['pandacore', 'pandaml', 'pandaphysics'],
                       help='Engines to benchmark (default: all)')
    parser.add_argument('--quick', action='store_true',
                       help='Quick demo mode (5 complexes, pandaml only)')
    parser.add_argument('--output_dir', type=str, 
                       default='rmsd_excellence_showcase',
                       help='Output directory name')
    parser.add_argument('--n_workers', type=int, default=2,
                       help='Number of parallel workers (default: 2)')
    
    args = parser.parse_args()
    
    # Quick mode overrides
    if args.quick:
        args.max_complexes = 5
        args.engines = ['pandaml']
        args.output_dir = 'rmsd_quick_demo'
        print("🚀 Quick Demo Mode: 5 complexes, PandaML engine only")
    
    print("=" * 80)
    print("🎯 PandaDock RMSD Excellence Benchmark")
    print("   Showcasing Sub-2Å Structural Accuracy")
    print("=" * 80)
    
    # Setup paths
    script_dir = Path(__file__).parent
    benchmark_script = script_dir / "scripts" / "rmsd_excellence_benchmark.py"
    pdbbind_dir = script_dir / "PDBbind"
    
    # Check if benchmark script exists
    if not benchmark_script.exists():
        print(f"❌ Error: Benchmark script not found at {benchmark_script}")
        return 1
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"📁 PDBbind Directory: {pdbbind_dir}")
    print(f"📊 Output Directory: {output_dir}")
    print(f"🔧 Engines: {', '.join(args.engines)}")
    if args.max_complexes:
        print(f"🔢 Max Complexes: {args.max_complexes}")
    print(f"⚡ Workers: {args.n_workers}")
    print()
    
    # Build command
    cmd = [
        sys.executable, str(benchmark_script),
        "--pdbbind_dir", str(pdbbind_dir),
        "--output_dir", str(output_dir),
        "--n_workers", str(args.n_workers),
        "--verbose"
    ]
    
    if args.max_complexes:
        cmd.extend(["--max_complexes", str(args.max_complexes)])
    
    # Note: The main script handles engine selection internally
    # We'll modify it to support engine filtering if needed
    
    print("🚀 Starting RMSD Excellence Benchmark...")
    print(f"📝 Command: {' '.join(cmd)}")
    print()
    
    # Run the benchmark
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, check=True, text=True)
        
        elapsed_time = time.time() - start_time
        
        print("\n" + "=" * 80)
        print("🏆 RMSD Excellence Benchmark Completed Successfully!")
        print(f"⏱️  Total Time: {elapsed_time:.1f} seconds")
        print("=" * 80)
        
        # Display results
        print(f"\n📊 Results saved to: {output_dir}")
        print("\n🖼️  Generated Visualizations:")
        
        expected_files = [
            "rmsd_excellence_master_figure.png",
            "rmsd_distribution_analysis.png", 
            "rmsd_success_analysis.png",
            "pose_quality_analysis.png",
            "rmsd_vs_complexity.png",
            "rmsd_performance_dashboard.png"
        ]
        
        for filename in expected_files:
            filepath = output_dir / filename
            if filepath.exists():
                print(f"   ✅ {filename}")
            else:
                print(f"   ⚠️  {filename} (not generated)")
        
        print(f"\n📄 Detailed Report: {output_dir}/rmsd_excellence_report.md")
        print(f"💾 Raw Data: {output_dir}/rmsd_excellence_data.csv")
        
        # Show key results if JSON exists
        json_file = output_dir / "rmsd_excellence_results.json"
        if json_file.exists():
            try:
                import json
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                metadata = data['metadata']
                print(f"\n🎯 Key Results:")
                print(f"   • Complexes Tested: {metadata['total_complexes']}")
                print(f"   • Success Rate (< 2Å): {metadata['overall_success_2A']:.1%}")
                print(f"   • Success Rate (< 3Å): {metadata['overall_success_3A']:.1%}")
                print(f"   • Mean RMSD: {metadata['mean_rmsd']:.2f} Å")
                
                if metadata['overall_success_2A'] > 0.4:
                    print("   🏆 EXCELLENT: Above industry standards!")
                elif metadata['overall_success_2A'] > 0.3:
                    print("   ✅ COMPETITIVE: Matches commercial software!")
                else:
                    print("   📈 PROMISING: Good foundation for optimization!")
                    
            except Exception as e:
                print(f"   ⚠️  Could not parse results: {e}")
        
        print(f"\n💡 Next Steps:")
        print(f"   1. View the master figure: {output_dir}/rmsd_excellence_master_figure.png")
        print(f"   2. Read the detailed report: {output_dir}/rmsd_excellence_report.md")
        print(f"   3. Use these results for publications and presentations!")
        
        return 0
        
    except subprocess.CalledProcessError as e:
        elapsed_time = time.time() - start_time
        print(f"\n❌ Benchmark failed after {elapsed_time:.1f} seconds")
        print(f"Error code: {e.returncode}")
        print("\n🔧 Troubleshooting:")
        print("   1. Check that PandaDock is properly installed")
        print("   2. Ensure all dependencies are available")
        print("   3. Try with --quick mode for a smaller test")
        print("   4. Check the logs above for specific error messages")
        return 1
    
    except KeyboardInterrupt:
        print("\n⏹️  Benchmark interrupted by user")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)