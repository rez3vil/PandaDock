#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Demo Benchmark Script for PandaDock

This script demonstrates the comprehensive benchmarking capabilities of PandaDock
by running a quick demo analysis with synthetic data.

Usage:
    python demo_benchmark.py --demo_type comprehensive
    python demo_benchmark.py --demo_type rmsd_analysis
    python demo_benchmark.py --demo_type full_workflow
"""

import os
import sys
import argparse
import tempfile
import shutil
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_comprehensive_demo(output_dir: str):
    """Run comprehensive benchmark demo with synthetic data"""
    logger.info("Running comprehensive benchmark demo...")
    
    # Create temporary PDBbind-like directory structure
    temp_dir = tempfile.mkdtemp(prefix="demo_pdbbind_")
    
    try:
        # Create mock index file
        index_file = Path(temp_dir) / "index.dat"
        with open(index_file, 'w') as f:
            f.write("# PDBbind demo index file\n")
            f.write("# code  resolution  year  -logKd/Ki  Kd/Ki  reference  ligand name\n")
            # Add some sample entries
            entries = [
                "1a1b  2.50  1998  5.22  6.03e-06  J.Med.Chem.(1998)41:1315  Demo_Ligand_1",
                "1c2d  1.80  2000  6.15  7.08e-07  Nature(2000)403:456  Demo_Ligand_2", 
                "1e3f  2.10  2002  4.89  1.29e-05  Science(2002)295:345  Demo_Ligand_3",
                "1g4h  1.95  2004  7.34  4.57e-08  Cell(2004)116:789  Demo_Ligand_4",
                "1i5j  2.30  2006  5.67  2.14e-06  PNAS(2006)103:123  Demo_Ligand_5"
            ]
            for entry in entries:
                f.write(f"{entry}\n")
        
        logger.info(f"Created demo index file: {index_file}")
        
        # Run comprehensive benchmark
        cmd = f"python pdbbind_comprehensive_benchmark.py --pdbbind_dir {temp_dir} --output_dir {output_dir} --max_complexes 5"
        logger.info(f"Running: {cmd}")
        
        import subprocess
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            logger.info("Comprehensive benchmark demo completed successfully!")
            logger.info(f"Results saved to: {output_dir}")
            
            # List generated files
            output_path = Path(output_dir)
            if output_path.exists():
                files = list(output_path.glob("*"))
                logger.info(f"Generated {len(files)} output files:")
                for file in files[:10]:  # Show first 10 files
                    logger.info(f"  - {file.name}")
                if len(files) > 10:
                    logger.info(f"  ... and {len(files) - 10} more files")
        else:
            logger.error(f"Benchmark failed: {result.stderr}")
            
    finally:
        # Clean up temp directory
        shutil.rmtree(temp_dir)
        logger.info("Cleaned up temporary files")

def run_rmsd_demo(output_dir: str):
    """Run RMSD analysis demo"""
    logger.info("Running RMSD analysis demo...")
    
    # Run RMSD analysis (will generate synthetic data if no input found)
    cmd = f"python rmsd_detailed_analysis.py --input_dir nonexistent --output_dir {output_dir}"
    logger.info(f"Running: {cmd}")
    
    import subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        logger.info("RMSD analysis demo completed successfully!")
        logger.info(f"Results saved to: {output_dir}")
        
        # List generated files
        output_path = Path(output_dir)
        if output_path.exists():
            files = list(output_path.glob("*"))
            logger.info(f"Generated {len(files)} output files:")
            for file in files:
                logger.info(f"  - {file.name}")
    else:
        logger.error(f"RMSD analysis failed: {result.stderr}")

def run_full_workflow_demo(output_dir: str):
    """Run full workflow demo combining both analyses"""
    logger.info("Running full workflow demo...")
    
    # Step 1: Run comprehensive benchmark
    benchmark_dir = Path(output_dir) / "benchmark_results"
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Step 1: Running comprehensive benchmark...")
    run_comprehensive_demo(str(benchmark_dir))
    
    # Step 2: Run RMSD analysis on benchmark results
    rmsd_dir = Path(output_dir) / "rmsd_analysis"
    rmsd_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Step 2: Running RMSD analysis...")
    cmd = f"python rmsd_detailed_analysis.py --input_dir {benchmark_dir} --output_dir {rmsd_dir}"
    
    import subprocess
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode == 0:
        logger.info("Full workflow demo completed successfully!")
        
        # Generate summary report
        summary_file = Path(output_dir) / "demo_summary.md"
        with open(summary_file, 'w') as f:
            f.write("# PandaDock Benchmark Demo Summary\n\n")
            f.write("## Workflow Completed\n\n")
            f.write("1. **Comprehensive Benchmark**: Generated synthetic PDBbind-style analysis\n")
            f.write("2. **RMSD Analysis**: Performed detailed structural accuracy evaluation\n\n")
            f.write("## Generated Outputs\n\n")
            
            # List all files
            for subdir in ["benchmark_results", "rmsd_analysis"]:
                subdir_path = Path(output_dir) / subdir
                if subdir_path.exists():
                    f.write(f"### {subdir.replace('_', ' ').title()}\n")
                    files = list(subdir_path.glob("*"))
                    for file in files:
                        f.write(f"- `{file.name}`\n")
                    f.write("\n")
            
            f.write("## Key Features Demonstrated\n\n")
            f.write("- Publication-quality visualizations\n")
            f.write("- Statistical analysis and validation\n")
            f.write("- Algorithm performance comparison\n")
            f.write("- Metal complex specialized analysis\n")
            f.write("- RMSD excellence evaluation\n")
            f.write("- Comprehensive reporting\n")
        
        logger.info(f"Demo summary saved to: {summary_file}")
        
    else:
        logger.error(f"RMSD analysis step failed: {result.stderr}")

def main():
    parser = argparse.ArgumentParser(description='PandaDock Benchmark Demo')
    parser.add_argument('--demo_type', default='full_workflow',
                       choices=['comprehensive', 'rmsd_analysis', 'full_workflow'],
                       help='Type of demo to run')
    parser.add_argument('--output_dir', default='demo_results',
                       help='Output directory for demo results')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Starting PandaDock benchmark demo: {args.demo_type}")
    logger.info(f"Output directory: {output_dir.absolute()}")
    
    try:
        if args.demo_type == 'comprehensive':
            run_comprehensive_demo(str(output_dir))
        elif args.demo_type == 'rmsd_analysis':
            run_rmsd_demo(str(output_dir))
        elif args.demo_type == 'full_workflow':
            run_full_workflow_demo(str(output_dir))
        
        logger.info("Demo completed successfully!")
        logger.info(f"\nTo view results, check: {output_dir.absolute()}")
        
    except Exception as e:
        logger.error(f"Demo failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()