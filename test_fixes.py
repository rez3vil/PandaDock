#!/usr/bin/env python3
"""
Test script to verify all fixes work correctly
"""

import numpy as np
import sys
import os

# Add current directory to path to import modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from docking.base_engine import Pose
from scoring.scoring_functions import ScoringFunctions
from utils.ic50_calculator import IC50Calculator

def test_ic50_calculation():
    """Test IC50 calculation fix"""
    print("Testing IC50 calculation...")
    
    # Create a pose with sample data
    coordinates = np.random.randn(13, 3) * 3.0
    pose = Pose(coordinates=coordinates, score=0.165337, energy=-8.615635)
    
    # Test that binding affinity is now negative (favorable)
    binding_affinity = pose.get_binding_affinity()
    print(f"  Binding affinity: {binding_affinity:.6f} kcal/mol")
    
    # Test that IC50 is not infinite
    ic50 = pose.get_ic50()
    print(f"  IC50: {ic50:.6e} nM")
    
    # Test IC50Calculator directly
    calc = IC50Calculator()
    ic50_calc = calc.delta_g_to_ic50(binding_affinity)
    print(f"  IC50 (calculator): {ic50_calc:.6e} M")
    
    assert binding_affinity < 0, "Binding affinity should be negative (favorable)"
    assert ic50 != float('inf'), "IC50 should not be infinite"
    assert ic50_calc != float('inf'), "IC50 from calculator should not be infinite"
    
    print("  ✅ IC50 calculation fix verified!")
    return True

def test_interaction_detection():
    """Test interaction detection fix"""
    print("Testing interaction detection...")
    
    scoring = ScoringFunctions()
    coordinates = np.random.randn(13, 3) * 3.0
    
    # Test H-bond detection
    hbonds = scoring.find_hbond_interactions(coordinates, None)
    print(f"  H-bonds detected: {len(hbonds)}")
    
    # Test hydrophobic detection
    hydrophobic = scoring.find_hydrophobic_interactions(coordinates, None)
    print(f"  Hydrophobic interactions: {len(hydrophobic)}")
    
    # Test salt bridge detection
    salt_bridges = scoring.find_salt_bridge_interactions(coordinates, None)
    print(f"  Salt bridges: {len(salt_bridges)}")
    
    assert len(hbonds) > 0, "Should detect some H-bonds"
    assert len(hydrophobic) > 0, "Should detect some hydrophobic interactions"
    
    print("  ✅ Interaction detection fix verified!")
    return True

def test_sdf_bond_generation():
    """Test SDF bond generation fix by checking bond inclusion"""
    print("Testing SDF bond generation...")
    
    from docking.base_engine import DockingEngine
    import tempfile
    import os
    
    # Create a simple mock engine class that doesn't require config
    class MockEngine(DockingEngine):
        def dock(self, ligand, receptor):
            pass
        
        def score(self, pose):
            return 0.0
    
    # Test with a simple molecule structure
    test_coordinates = np.array([
        [0.0, 0.0, 0.0],
        [1.4, 0.0, 0.0], 
        [2.1, 1.2, 0.0],
        [3.5, 1.2, 0.0]
    ])
    test_atom_types = ['C', 'C', 'N', 'O']
    test_bonds = [(0, 1, 'single'), (1, 2, 'single'), (2, 3, 'double')]
    
    ligand_data = {
        'coordinates': test_coordinates,
        'atom_types': test_atom_types,
        'bonds': test_bonds,
        'num_atoms': len(test_coordinates),
        'molecular_weight': 58.0,
        'partial_charges': [0.0] * len(test_coordinates)
    }
    
    # Create a mock config
    class MockConfig:
        pass
    
    # Create a mock engine and set ligand data
    engine = MockEngine(MockConfig())
    engine.ligand = ligand_data
    
    # Create a test pose
    pose = Pose(coordinates=test_coordinates, score=0.2, energy=-5.0)
    pose.pose_id = "test_pose"
    
    # Test SDF generation
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.sdf', delete=False) as f:
        temp_filename = f.name
    
    try:
        engine._save_pose_as_sdf(pose, temp_filename)
        
        # Read the generated SDF file
        with open(temp_filename, 'r') as f:
            content = f.read()
        
        print(f"  Generated SDF content preview:")
        lines = content.split('\n')
        for i, line in enumerate(lines[:15]):  # Show first 15 lines
            print(f"    {i+1:2d}: {line}")
        
        # Check that bonds are included
        counts_line = lines[3]
        num_atoms = int(counts_line[:3])
        num_bonds = int(counts_line[3:6])
        
        print(f"  Atoms: {num_atoms}, Bonds: {num_bonds}")
        
        assert num_atoms == 4, f"Expected 4 atoms, got {num_atoms}"
        assert num_bonds == 3, f"Expected 3 bonds, got {num_bonds}"
        
        # Check that bond block exists by looking for the bond lines
        bond_block_found = False
        for i, line in enumerate(lines):
            if 'M  END' in line:
                # Bond lines should be before M  END
                bond_lines = lines[4+num_atoms:i]  # Skip header + atom lines
                valid_bond_lines = 0
                for bond_line in bond_lines:
                    parts = bond_line.strip().split()
                    if len(parts) >= 3:
                        try:
                            int(parts[0])  # atom 1
                            int(parts[1])  # atom 2  
                            int(parts[2])  # bond order
                            valid_bond_lines += 1
                        except ValueError:
                            pass
                
                print(f"  Valid bond lines found: {valid_bond_lines}")
                assert valid_bond_lines >= 3, f"Expected at least 3 bond lines, got {valid_bond_lines}"
                bond_block_found = True
                break
        
        assert bond_block_found, "Bond block not found in SDF file"
        
        print("  ✅ SDF bond generation fix verified!")
        return True
        
    finally:
        # Clean up
        if os.path.exists(temp_filename):
            os.unlink(temp_filename)

def test_filename_fix():
    """Test SDF filename fix"""
    print("Testing SDF filename fix...")
    
    from docking.base_engine import DockingEngine
    import tempfile
    import os
    
    # Create a simple mock engine class that doesn't require config
    class MockEngine(DockingEngine):
        def dock(self, ligand, receptor):
            pass
        
        def score(self, pose):
            return 0.0
    
    # Create a mock config
    class MockConfig:
        pass
    
    engine = MockEngine(MockConfig())
    
    # Create test poses
    test_coordinates = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])
    poses = [
        Pose(coordinates=test_coordinates, score=0.2, energy=-5.0, pose_id="test_pose_1"),
        Pose(coordinates=test_coordinates, score=0.3, energy=-4.0, pose_id="test_pose_2")
    ]
    
    # Test saving poses
    with tempfile.TemporaryDirectory() as temp_dir:
        engine.save_poses(poses, temp_dir)
        
        # Check that all.ligands.sdf was created (not all_poses.sdf)
        expected_file = os.path.join(temp_dir, "all.ligands.sdf")
        old_file = os.path.join(temp_dir, "all_poses.sdf")
        
        assert os.path.exists(expected_file), f"Expected file {expected_file} not found"
        assert not os.path.exists(old_file), f"Old filename {old_file} should not exist"
        
        print(f"  ✅ Filename changed to 'all.ligands.sdf'")
        
        # Check file content
        with open(expected_file, 'r') as f:
            content = f.read()
        
        print(f"  File contains {content.count('$$$$')} structures")
        assert content.count('$$$$') == 2, "Should contain 2 structures"
        
        print("  ✅ SDF filename fix verified!")
        return True

def main():
    """Run all tests"""
    print("🧪 Testing PandaDock fixes...\n")
    
    tests = [
        test_ic50_calculation,
        test_interaction_detection,
        test_sdf_bond_generation,
        test_filename_fix
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            print(f"\n{'='*50}")
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ❌ Test failed with error: {e}")
            failed += 1
    
    print(f"\n{'='*50}")
    print(f"📊 Test Results:")
    print(f"  ✅ Passed: {passed}")
    print(f"  ❌ Failed: {failed}")
    print(f"  📝 Total:  {passed + failed}")
    
    if failed == 0:
        print("\n🎉 All fixes verified successfully!")
        return True
    else:
        print(f"\n⚠️  {failed} test(s) failed. Please check the issues above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)