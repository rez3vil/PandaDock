#!/usr/bin/env python3
"""
Trace where atoms are being lost in the docking process.
"""

import sys
import os
sys.path.append('.')

from pandadock.molecules.ligand_handler import LigandHandler
from pandadock.core.pose_generator import PoseGenerator

def test_ligand_loading():
    """Test ligand loading to see where atoms are lost."""
    print("Testing ligand loading process...")
    
    handler = LigandHandler()
    
    # Test loading the 46-atom ligand
    try:
        ligand = handler.load_ligand("run/ligand.sdf")
        print(f"✅ Loaded ligand: {ligand.n_atoms} atoms")
        print(f"   Atom types: {len(set(ligand.atom_types))} unique types")
        print(f"   Bonds: {len(ligand.bonds)}")
        print(f"   Coordinates shape: {ligand.coords.shape}")
        
        # Check coordinate ranges
        import numpy as np
        coord_ranges = np.max(ligand.coords, axis=0) - np.min(ligand.coords, axis=0)
        print(f"   Coordinate ranges: {coord_ranges}")
        
        # Check first few atoms
        print(f"   First 5 atoms:")
        for i in range(min(5, ligand.n_atoms)):
            print(f"     {i+1}: {ligand.atom_types[i]} at {ligand.coords[i]}")
        
        return ligand
        
    except Exception as e:
        print(f"❌ Error loading ligand: {e}")
        return None

def test_pose_generation(ligand):
    """Test pose generation to see if atoms are lost there."""
    if ligand is None:
        return
        
    print(f"\n{'='*60}")
    print("Testing pose generation...")
    
    try:
        # Create a mock protein for testing
        class MockProtein:
            def __init__(self):
                import numpy as np
                # Create some protein coordinates
                self.coords = np.random.randn(100, 3) * 10
                self.active_site = {
                    'center': np.array([0.0, 0.0, 0.0]),
                    'radius': 10.0
                }
        
        protein = MockProtein()
        
        # Create pose generator
        generator = PoseGenerator()
        
        print(f"Original ligand: {ligand.n_atoms} atoms")
        print(f"Generating poses...")
        
        # Generate poses
        poses = generator.generate_initial_poses(ligand, protein, n_poses=3)
        
        print(f"Generated {len(poses)} poses")
        
        for i, pose in enumerate(poses):
            print(f"  Pose {i+1}: shape {pose.shape}, atoms: {len(pose)}")
        
        return poses
        
    except Exception as e:
        print(f"❌ Error in pose generation: {e}")
        import traceback
        traceback.print_exc()
        return None

def test_genetic_algorithm_coordinate_handling(ligand):
    """Test genetic algorithm coordinate handling."""
    print(f"\n{'='*60}")
    print("Testing genetic algorithm coordinate handling...")
    
    try:
        from pandadock.algorithms.genetic_algorithm_enhanced import GeneticAlgorithm
        
        # Mock scoring function
        class MockScoring:
            def score(self, protein, pose):
                return -5.0  # Random score
        
        # Create GA
        ga = GeneticAlgorithm(MockScoring(), population_size=5, max_generations=1)
        
        # Test coordinate extraction
        coords = ga._get_ligand_coordinates(ligand)
        print(f"GA extracted coordinates: {coords.shape if coords is not None else 'None'}")
        
        if coords is not None:
            print(f"  Expected: {ligand.n_atoms} atoms")
            print(f"  Got: {len(coords)} atoms")
            if len(coords) != ligand.n_atoms:
                print("❌ ATOM COUNT MISMATCH in genetic algorithm!")
            else:
                print("✅ Atom count preserved in genetic algorithm")
        
        return coords
        
    except Exception as e:
        print(f"❌ Error in genetic algorithm testing: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    # Test the loading and processing pipeline
    ligand = test_ligand_loading()
    
    if ligand:
        poses = test_pose_generation(ligand)
        coords = test_genetic_algorithm_coordinate_handling(ligand)
    
    print(f"\n{'='*60}")
    print("SUMMARY:")
    if ligand:
        print(f"✅ Ligand loading: {ligand.n_atoms} atoms preserved")
        if poses:
            print(f"✅ Pose generation: Works, generates poses with {len(poses[0])} atoms each")
        if coords is not None:
            print(f"{'✅' if len(coords) == ligand.n_atoms else '❌'} GA coordinate extraction: {len(coords)}/{ligand.n_atoms} atoms")