"""
Test script to verify PDB file formatting fixes.

This script tests that the corrected PDB writing functions produce
properly formatted files that can be opened in molecular viewers.
"""

import sys
import numpy as np
import tempfile
from pathlib import Path

# Add pandadock to path
sys.path.insert(0, '/Users/pritam/PandaDock')

def test_pdb_writing():
    """Test that PDB files are properly formatted."""
    print("Testing PDB file formatting fixes...")
    
    try:
        from pandadock.io.pdb_writer import PDBWriter, validate_pdb
        
        # Create test coordinates
        ligand_coords = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0], 
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
            [13.0, 14.0, 15.0]
        ])
        
        protein_coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 1.5, 1.5],
            [3.0, 3.0, 3.0]
        ])
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Test 1: Write ligand pose
            print("  Testing ligand pose writing...")
            ligand_file = temp_path / "test_ligand.pdb"
            
            PDBWriter.write_ligand_pose(
                filename=ligand_file,
                coordinates=ligand_coords,
                ligand_name="LIG",
                score=-8.5,
                pose_rank=1
            )
            
            # Validate ligand file
            is_valid, message = validate_pdb(ligand_file)
            if is_valid:
                print("    ✓ Ligand PDB file is properly formatted")
            else:
                print(f"    ✗ Ligand PDB validation failed: {message}")
                return False
            
            # Check file content
            with open(ligand_file, 'r') as f:
                content = f.read()
                print("    Sample ligand PDB content:")
                lines = content.split('\n')
                for line in lines[:10]:  # Show first 10 lines
                    if line.strip():
                        print(f"      {line}")
            
            # Test 2: Write protein-ligand complex
            print("  Testing protein-ligand complex writing...")
            complex_file = temp_path / "test_complex.pdb"
            
            PDBWriter.write_protein_ligand_complex(
                filename=complex_file,
                protein_coords=protein_coords,
                ligand_coords=ligand_coords,
                ligand_name="LIG",
                score=-8.5
            )
            
            # Validate complex file
            is_valid, message = validate_pdb(complex_file)
            if is_valid:
                print("    ✓ Complex PDB file is properly formatted")
            else:
                print(f"    ✗ Complex PDB validation failed: {message}")
                return False
            
            # Check complex content
            with open(complex_file, 'r') as f:
                content = f.read()
                print("    Sample complex PDB content:")
                lines = content.split('\n')
                for line in lines[:8]:  # Show first 8 lines
                    if line.strip():
                        print(f"      {line}")
            
            # Test 3: Check atom naming format
            print("  Testing atom naming format...")
            with open(ligand_file, 'r') as f:
                hetatm_lines = [line for line in f.readlines() if line.startswith('HETATM')]
            
            if hetatm_lines:
                # Check first HETATM line format
                first_line = hetatm_lines[0]
                atom_name = first_line[12:16]  # Atom name field
                print(f"    First atom name: '{atom_name}' (should be properly formatted)")
                
                # Check if atom name is properly formatted (no extra spaces)
                if atom_name.strip() and not atom_name.startswith(' C '):
                    print("    ✓ Atom names are properly formatted")
                else:
                    print("    ✗ Atom names still have formatting issues")
                    return False
            
            print("✓ All PDB formatting tests passed!")
            print("\nPDB files should now open correctly in:")
            print("  - ChimeraX")
            print("  - PyMOL") 
            print("  - VMD")
            print("  - Other molecular visualization software")
            
            return True
            
    except Exception as e:
        print(f"✗ PDB testing failed: {e}")
        return False


def test_result_writers_integration():
    """Test that result writers use the fixed PDB formatting."""
    print("\nTesting ResultWriters integration...")
    
    try:
        from pandadock.io.result_writers import ResultWriters
        
        # Create mock results
        mock_results = {
            'poses': [
                {
                    'rank': 1,
                    'score': -8.5,
                    'coordinates': [
                        [1.0, 2.0, 3.0],
                        [4.0, 5.0, 6.0],
                        [7.0, 8.0, 9.0]
                    ]
                }
            ]
        }
        
        mock_config = {
            'algorithm': 'genetic',
            'scoring_function': 'physics_based'
        }
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test result writers
            writer = ResultWriters(temp_dir)
            saved_files = writer.save_all_results(mock_results, mock_config)
            
            # Check if pose files were created
            poses_dir = Path(temp_dir) / "poses"
            if poses_dir.exists():
                pose_files = list(poses_dir.glob("*.pdb"))
                if pose_files:
                    print("  ✓ ResultWriters created PDB files")
                    
                    # Validate the created file
                    from pandadock.io.pdb_writer import validate_pdb
                    is_valid, message = validate_pdb(pose_files[0])
                    if is_valid:
                        print("  ✓ ResultWriters PDB files are properly formatted")
                        return True
                    else:
                        print(f"  ✗ ResultWriters PDB validation failed: {message}")
                        return False
                else:
                    print("  ✗ No PDB files created by ResultWriters")
                    return False
            else:
                print("  ✗ Poses directory not created")
                return False
                
    except Exception as e:
        print(f"✗ ResultWriters integration test failed: {e}")
        return False


def main():
    """Run all PDB formatting tests."""
    print("🧪 Testing PDB File Formatting Fixes")
    print("=" * 50)
    
    tests_passed = 0
    total_tests = 2
    
    # Test 1: Basic PDB writing
    if test_pdb_writing():
        tests_passed += 1
    
    # Test 2: Result writers integration
    if test_result_writers_integration():
        tests_passed += 1
    
    print("\n" + "=" * 50)
    print(f"📊 Test Results: {tests_passed}/{total_tests} passed")
    
    if tests_passed == total_tests:
        print("🎉 ALL PDB FORMATTING TESTS PASSED!")
        print("✅ PDB files should now open correctly in ChimeraX and other molecular viewers")
        return True
    else:
        print("❌ Some PDB formatting tests failed")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)