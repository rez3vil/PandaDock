#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test Integrated Plotting in PandaDock Command Line
================================================

Tests the integrated plotting functionality in the main PandaDock command line interface.

Usage:
    python test_integrated_plots.py
"""

import subprocess
import sys
import os
from pathlib import Path
import time

def test_pandadock_with_plots():
    """Test PandaDock command line with plotting features"""
    
    print("🧪 Testing PandaDock Command Line Integration with Plots")
    print("=" * 65)
    
    # Check if test files exist
    protein_file = "tests/beta-2_alpha-1.pdb"
    ligand_file = "tests/propofol.pdb"
    
    if not Path(protein_file).exists():
        print(f"❌ Protein file not found: {protein_file}")
        return False
    
    if not Path(ligand_file).exists():
        print(f"❌ Ligand file not found: {ligand_file}")
        return False
    
    print(f"✅ Found test files:")
    print(f"   📁 Protein: {protein_file}")
    print(f"   💊 Ligand: {ligand_file}")
    
    # Test different command line combinations
    tests = [
        {
            'name': 'Basic docking with plots',
            'cmd': [
                'python', '-m', 'pandadock',
                '--protein', protein_file,
                '--ligand', ligand_file,
                '--mode', 'balanced',
                '--scoring', 'pandacore',
                '--center', '-15.7', '-17.7', '8.18',
                '--size', '40', '40', '40',
                '--plots',
                '--out', 'test_plots_basic'
            ],
            'expected_files': [
                'pandadock_report.html',
                'binding_metrics_analysis.png',
                'master_publication.png'
            ]
        },
        {
            'name': 'Comprehensive outputs',
            'cmd': [
                'python', '-m', 'pandadock',
                '--protein', protein_file,
                '--ligand', ligand_file,
                '--mode', 'balanced',
                '--scoring', 'pandacore',
                '--flexible-residues', 'ASN265',
                '--center', '-15.7', '-17.7', '8.18',
                '--size', '40', '40', '40',
                '--all-outputs',
                '--protein-name', 'GABA_A Receptor',
                '--ligand-name', 'Propofol',
                '--out', 'test_plots_comprehensive'
            ],
            'expected_files': [
                'pandadock_report.html',
                'pandadock_report.csv',
                'pandadock_report.json',
                'binding_metrics_analysis.png',
                'ic50_ec50_analysis.png',
                'master_publication.png',
                'detailed_analysis_report.txt'
            ]
        },
        {
            'name': 'Just interaction maps',
            'cmd': [
                'python', '-m', 'pandadock',
                '--protein', protein_file,
                '--ligand', ligand_file,
                '--mode', 'balanced',
                '--scoring', 'pandacore',
                '--center', '-15.7', '-17.7', '8.18',
                '--size', '40', '40', '40',
                '--interaction-maps',
                '--master-plot',
                '--out', 'test_plots_maps'
            ],
            'expected_files': [
                'pandadock_report.html',
                'master_publication.png'
            ]
        }
    ]
    
    success_count = 0
    
    for i, test in enumerate(tests, 1):
        print(f"\n🔬 Test {i}: {test['name']}")
        print("-" * 50)
        
        # Clean up previous test results
        output_dir = Path(test['cmd'][-1])  # Last argument is output directory
        if output_dir.exists():
            import shutil
            shutil.rmtree(output_dir)
        
        print(f"Running: {' '.join(test['cmd'])}")
        
        try:
            # Run the command
            start_time = time.time()
            result = subprocess.run(
                test['cmd'],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            end_time = time.time()
            
            print(f"⏱️  Execution time: {end_time - start_time:.1f} seconds")
            
            if result.returncode == 0:
                print("✅ Command executed successfully")
                
                # Check for expected files
                missing_files = []
                found_files = []
                
                for expected_file in test['expected_files']:
                    file_path = output_dir / expected_file
                    if file_path.exists():
                        file_size = file_path.stat().st_size / 1024
                        found_files.append(f"{expected_file} ({file_size:.1f} KB)")
                    else:
                        missing_files.append(expected_file)
                
                if found_files:
                    print("✅ Generated files:")
                    for file_info in found_files:
                        print(f"   📄 {file_info}")
                
                if missing_files:
                    print("⚠️  Missing expected files:")
                    for missing_file in missing_files:
                        print(f"   ❌ {missing_file}")
                
                # Check for any additional plot files
                if output_dir.exists():
                    all_files = list(output_dir.glob("**/*"))
                    plot_files = [f for f in all_files if f.suffix in ['.png', '.pdf']]
                    
                    if plot_files:
                        print(f"📊 Total plot files generated: {len(plot_files)}")
                        for plot_file in plot_files:
                            if plot_file.name not in [f.split(' (')[0] for f in found_files]:
                                file_size = plot_file.stat().st_size / 1024
                                print(f"   📊 {plot_file.name} ({file_size:.1f} KB)")
                
                if not missing_files:
                    success_count += 1
                    print("🎉 Test PASSED")
                else:
                    print("⚠️  Test PARTIAL - some files missing")
                
            else:
                print("❌ Command failed")
                print(f"Return code: {result.returncode}")
                if result.stderr:
                    print("Error output:")
                    print(result.stderr[:500])  # First 500 characters
                
        except subprocess.TimeoutExpired:
            print("❌ Command timed out after 5 minutes")
        except Exception as e:
            print(f"❌ Test failed with error: {e}")
    
    print(f"\n📊 Test Summary")
    print("=" * 30)
    print(f"Tests passed: {success_count}/{len(tests)}")
    print(f"Success rate: {success_count/len(tests)*100:.1f}%")
    
    if success_count == len(tests):
        print("🎉 All tests passed! Plotting integration is working correctly.")
        return True
    else:
        print("⚠️  Some tests failed. Check the output above for details.")
        return False

def test_help_output():
    """Test that help output includes new plotting options"""
    
    print(f"\n📖 Testing Help Output")
    print("-" * 25)
    
    try:
        result = subprocess.run(
            ['python', '-m', 'pandadock', '--help'],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            help_text = result.stdout
            
            # Check for plotting-related options
            plotting_options = [
                '--plots',
                '--interaction-maps',
                '--master-plot',
                '--txt-report',
                '--all-outputs',
                '--protein-name',
                '--ligand-name'
            ]
            
            found_options = []
            missing_options = []
            
            for option in plotting_options:
                if option in help_text:
                    found_options.append(option)
                else:
                    missing_options.append(option)
            
            print(f"✅ Found {len(found_options)}/{len(plotting_options)} plotting options in help")
            
            if missing_options:
                print("❌ Missing options:")
                for option in missing_options:
                    print(f"   {option}")
            
            # Check for examples
            if '--all-outputs' in help_text and 'GABA_A Receptor' in help_text:
                print("✅ Help includes comprehensive plotting examples")
            else:
                print("⚠️  Help missing comprehensive plotting examples")
            
            return len(missing_options) == 0
            
        else:
            print("❌ Help command failed")
            return False
            
    except Exception as e:
        print(f"❌ Help test failed: {e}")
        return False

def main():
    """Main test function"""
    
    print("Starting PandaDock plotting integration tests...\n")
    
    # Test help output first
    help_success = test_help_output()
    
    # Test actual functionality
    functionality_success = test_pandadock_with_plots()
    
    print("\n" + "=" * 70)
    if help_success and functionality_success:
        print("🎉 ALL TESTS PASSED!")
        print("\n✨ PandaDock plotting integration is working correctly!")
        print("\nUsers can now run:")
        print("python -m pandadock --protein target.pdb --ligand compound.sdf --all-outputs")
        print("\nTo get comprehensive plots, interaction maps, and publication-ready figures!")
    else:
        print("❌ Some tests failed.")
        if not help_success:
            print("   - Help integration needs work")
        if not functionality_success:
            print("   - Functionality integration needs work")
    
    return help_success and functionality_success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)