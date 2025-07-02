#!/usr/bin/env python3
"""
Test script for PandaDock benchmarking system.

This script creates mock data to test the benchmarking system without requiring
the full PDBbind dataset. Useful for development and testing.
"""

import os
import sys
import json
import tempfile
import shutil
from pathlib import Path
from typing import List

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from pdbbind_benchmark import PDBBindEntry, BenchmarkResult, PandaDockBenchmark
from analyze_benchmark_results import BenchmarkAnalyzer


def create_mock_pdbbind_data(base_dir: Path) -> Path:
    """Create mock PDBbind dataset for testing."""
    pdbbind_dir = base_dir / "mock_pdbbind"
    pdbbind_dir.mkdir(exist_ok=True)
    
    # Create mock index file
    index_content = """# ==============================================================================
# List of protein-ligand complexes in the PDBbind refined set v.2020
# 5316 protein-ligand complexes in total, sorted by their release year.
# Latest update: July 2021
# PDB code, resolution, release year, binding data, reference, ligand name
# ==============================================================================
1abc  2.10  1995  Kd=49uM       // 1abc.pdf (ABC)
2def  1.80  1998  Ki=190nM      // 2def.pdf (DEF)
3ghi  2.30  2001  IC50=2.5uM    // 3ghi.pdf (GHI)
4jkl  1.90  2003  Kd=150nM      // 4jkl.pdf (JKL)
5mno  2.50  2005  Ki=8.7nM      // 5mno.pdf (MNO)
"""
    
    with open(pdbbind_dir / "INDEX_refined_set.2020", "w") as f:
        f.write(index_content)
    
    # Create mock PDB directories and files
    pdb_codes = ["1abc", "2def", "3ghi", "4jkl", "5mno"]
    
    for pdb_code in pdb_codes:
        pdb_dir = pdbbind_dir / pdb_code
        pdb_dir.mkdir(exist_ok=True)
        
        # Create mock protein file
        protein_content = """ATOM      1  N   ALA A   1      20.154  10.000  15.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1      21.618  10.000  15.000  1.00 20.00           C  
ATOM      3  C   ALA A   1      22.154  11.421  15.000  1.00 20.00           C  
ATOM      4  O   ALA A   1      21.618  12.421  15.000  1.00 20.00           O  
END
"""
        with open(pdb_dir / f"{pdb_code}_protein.pdb", "w") as f:
            f.write(protein_content)
        
        # Create mock ligand file (SDF format)
        ligand_content = """
  Marvin  01230100002D          

  3  2  0  0  0  0            999 V2000
   -0.7145    0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.2063    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.2063    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
        with open(pdb_dir / f"{pdb_code}_ligand.sdf", "w") as f:
            f.write(ligand_content)
        
        # Create mock pocket file
        pocket_content = """ATOM      1  C   LIG A   1      20.000  10.000  15.000  1.00 20.00           C  
ATOM      2  C   LIG A   1      21.000  10.000  15.000  1.00 20.00           C  
ATOM      3  N   LIG A   1      22.000  10.000  15.000  1.00 20.00           N  
END
"""
        with open(pdb_dir / f"{pdb_code}_pocket.pdb", "w") as f:
            f.write(pocket_content)
    
    print(f"Created mock PDBbind dataset at: {pdbbind_dir}")
    return pdbbind_dir


def create_mock_benchmark_results(output_dir: Path) -> List[BenchmarkResult]:
    """Create mock benchmark results for testing analysis."""
    import random
    random.seed(42)  # For reproducible results
    
    # Mock experimental data
    pdb_codes = ["1abc", "2def", "3ghi", "4jkl", "5mno"]
    experimental_values = [49e-6, 190e-9, 2.5e-6, 150e-9, 8.7e-9]  # Convert to M
    exp_types = ["KD", "KI", "IC50", "KD", "KI"]
    
    algorithms = ["genetic", "random"]
    devices = ["CPU", "GPU"]
    scoring_functions = ["standard", "enhanced"]
    
    results = []
    
    for i, pdb_code in enumerate(pdb_codes):
        for algorithm in algorithms:
            for device in devices:
                for scoring in scoring_functions:
                    # Generate realistic mock data
                    success = random.random() > 0.1  # 90% success rate
                    
                    if success:
                        # Generate correlated predicted values with some noise
                        exp_val = experimental_values[i]
                        log_exp = -np.log10(exp_val)
                        
                        # Add noise to correlation
                        noise = random.gauss(0, 1.5)
                        predicted_log = log_exp + noise
                        predicted_kd = 10**(-predicted_log)
                        
                        # Generate other metrics
                        best_score = -predicted_log + random.gauss(0, 0.5)
                        runtime = random.uniform(10, 120)
                        if device == "GPU":
                            runtime *= random.uniform(0.2, 0.6)  # GPU speedup
                        
                        poses = random.randint(1, 10)
                        
                        result = BenchmarkResult(
                            pdb_code=pdb_code,
                            algorithm=algorithm,
                            device=device,
                            scoring_function=scoring,
                            total_time=runtime,
                            setup_time=runtime * 0.1,
                            docking_time=runtime * 0.9,
                            best_score=best_score,
                            poses_generated=poses,
                            success=True,
                            predicted_delta_g=best_score,
                            predicted_kd=predicted_kd,
                            predicted_ic50=predicted_kd * 2.0,
                            predicted_ki=predicted_kd * 0.5,
                            experimental_value=exp_val,
                            experimental_type=exp_types[i]
                        )
                    else:
                        # Failed run
                        result = BenchmarkResult(
                            pdb_code=pdb_code,
                            algorithm=algorithm,
                            device=device,
                            scoring_function=scoring,
                            total_time=random.uniform(5, 30),
                            setup_time=0,
                            docking_time=0,
                            best_score=float('inf'),
                            poses_generated=0,
                            success=False,
                            error_message="Mock failure for testing",
                            experimental_value=exp_val,
                            experimental_type=exp_types[i]
                        )
                    
                    results.append(result)
    
    return results


def test_benchmark_system():
    """Test the complete benchmarking system."""
    import numpy as np
    
    print("Testing PandaDock Benchmark System")
    print("=" * 40)
    
    # Create temporary directory for test
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        try:
            # Test 1: Create mock PDBbind data
            print("1. Creating mock PDBbind dataset...")
            pdbbind_dir = create_mock_pdbbind_data(temp_path)
            
            # Test 2: Parse PDBbind data
            print("2. Testing PDBbind parser...")
            from pdbbind_benchmark import PDBBindParser
            
            parser = PDBBindParser(str(pdbbind_dir), "INDEX_refined_set.2020")
            entries = parser.parse_index_file()
            
            assert len(entries) == 5, f"Expected 5 entries, got {len(entries)}"
            assert entries[0].pdb_code == "1abc", f"Expected 1abc, got {entries[0].pdb_code}"
            assert entries[0].binding_type == "KD", f"Expected KD, got {entries[0].binding_type}"
            assert abs(entries[0].binding_value - 49e-6) < 1e-9, "Binding value conversion failed"
            
            print("   ✓ PDBbind parser working correctly")
            
            # Test 3: Create mock benchmark results
            print("3. Creating mock benchmark results...")
            output_dir = temp_path / "test_results"
            output_dir.mkdir()
            
            results = create_mock_benchmark_results(output_dir)
            
            # Save results
            results_data = [result.__dict__ for result in results]
            with open(output_dir / "benchmark_results.json", "w") as f:
                json.dump(results_data, f, indent=2)
            
            print(f"   ✓ Created {len(results)} mock results")
            
            # Test 4: Test benchmark analysis
            print("4. Testing benchmark analyzer...")
            analyzer = BenchmarkAnalyzer(str(output_dir))
            
            # Test summary
            print("   Testing summary generation...")
            analyzer.print_summary()
            
            print("   ✓ Analysis completed successfully")
            
            # Test 5: Test exports
            print("5. Testing export functions...")
            
            # Test CSV export
            df = analyzer.export_csv_for_analysis()
            assert len(df) == len(results), "CSV export length mismatch"
            
            # Test markdown report
            analyzer.generate_detailed_report()
            
            report_file = output_dir / "detailed_report.md"
            assert report_file.exists(), "Markdown report not created"
            
            print("   ✓ Export functions working correctly")
            
            # Test 6: Validate correlation analysis
            print("6. Testing correlation analysis...")
            
            successful_results = [r for r in results if r.success]
            if len(successful_results) > 0:
                experimental = [r.experimental_value for r in successful_results]
                predicted = [r.predicted_kd for r in successful_results]
                
                if len(experimental) > 3:
                    import scipy.stats as stats
                    r_value, p_value = stats.pearsonr(np.log10(experimental), np.log10(predicted))
                    print(f"   Mock correlation R = {r_value:.3f}, p = {p_value:.3f}")
                    
                    # Should have some correlation since we generated correlated data
                    assert abs(r_value) > 0.1, "Correlation too weak for mock data"
            
            print("   ✓ Correlation analysis working correctly")
            
            print("\n" + "=" * 40)
            print("✅ ALL TESTS PASSED!")
            print("The benchmarking system is working correctly.")
            print("\nTest output saved to:", output_dir)
            
            # Show some sample results
            print("\nSample Results:")
            print("-" * 20)
            for i, result in enumerate(results[:3]):
                status = "✓" if result.success else "✗"
                print(f"{status} {result.pdb_code} ({result.algorithm}/{result.device}): "
                      f"Score={result.best_score:.3f}, Time={result.total_time:.1f}s")
            
            return True
            
        except Exception as e:
            print(f"\n❌ TEST FAILED: {e}")
            import traceback
            traceback.print_exc()
            return False


def test_pdbbind_parser_edge_cases():
    """Test edge cases for PDBbind parser."""
    print("\nTesting PDBbind Parser Edge Cases")
    print("-" * 35)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        pdbbind_dir = temp_path / "edge_case_pdbbind"
        pdbbind_dir.mkdir()
        
        # Test various binding data formats
        index_content = """# Test index file
1abc  2.10  1995  Kd=49uM       // 1abc.pdf (ABC)
2def  1.80  1998  Ki=190nM      // 2def.pdf (DEF)
3ghi  2.30  2001  IC50=2.5mM    // 3ghi.pdf (GHI)
4jkl  1.90  2003  Kd=0.15uM     // 4jkl.pdf (JKL)
5mno  2.50  2005  Ki=8.7pM      // 5mno.pdf (MNO)
# This line should be ignored
"""
        
        with open(pdbbind_dir / "INDEX_test.2020", "w") as f:
            f.write(index_content)
        
        # Don't create directories - test missing files handling
        
        from pdbbind_benchmark import PDBBindParser
        
        parser = PDBBindParser(str(pdbbind_dir), "INDEX_test.2020")
        entries = parser.parse_index_file()
        
        # Should get 0 entries because files don't exist
        assert len(entries) == 0, f"Expected 0 entries due to missing files, got {len(entries)}"
        
        print("✓ Edge case testing passed")


if __name__ == "__main__":
    # Import required modules
    try:
        import numpy as np
        import pandas as pd
        
        # Run tests
        success = test_benchmark_system()
        test_pdbbind_parser_edge_cases()
        
        if success:
            print("\n🎉 Test completed successfully!")
            print("\nNext steps:")
            print("1. Download the PDBbind dataset")
            print("2. Edit run_benchmark_examples.sh with your PDBbind path")
            print("3. Run: ./run_benchmark_examples.sh")
        else:
            sys.exit(1)
            
    except ImportError as e:
        print(f"Missing required package: {e}")
        print("Please install required packages:")
        print("pip install pandas numpy scipy sklearn matplotlib seaborn")
        sys.exit(1)