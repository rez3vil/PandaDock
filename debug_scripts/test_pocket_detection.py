#!/usr/bin/env python3
"""
Test script to validate pocket detection and sphere.pdb generation.
"""

import numpy as np
import tempfile
import os
from pathlib import Path

# Add the project to Python path
import sys
sys.path.insert(0, '/mnt/7b616197-a2a7-4736-bd58-c500d1a8c523/Panda-Software/PandaDock')

from pandadock.molecules import ProteinHandler, ProteinStructure

def create_test_protein():
    """Create a simple test protein structure."""
    # Create a simple protein with a cavity
    coords = []
    
    # Create a ring-like structure to have a pocket
    for i in range(20):
        angle = 2 * np.pi * i / 20
        x = 10 * np.cos(angle)
        y = 10 * np.sin(angle)
        z = 0
        coords.append([x, y, z])
    
    # Add some atoms above and below
    for i in range(10):
        angle = 2 * np.pi * i / 10
        x = 8 * np.cos(angle)
        y = 8 * np.sin(angle)
        z = 5
        coords.append([x, y, z])
        coords.append([x, y, -5])
    
    protein = ProteinStructure(
        coords=np.array(coords),
        atom_names=[f"CA{i}" for i in range(len(coords))],
        residue_names=["ALA"] * len(coords),
        residue_ids=[f"A_{i}" for i in range(len(coords))]
    )
    
    return protein

def test_pocket_detection():
    """Test pocket detection functionality."""
    print("Testing pocket detection...")
    
    # Create test protein
    protein = create_test_protein()
    print(f"Created test protein with {protein.n_atoms} atoms")
    
    # Test pocket detection
    pockets = protein.detect_pockets(min_volume=10.0)
    print(f"Detected {len(pockets)} pockets")
    
    for i, pocket in enumerate(pockets):
        print(f"  Pocket {i+1}:")
        print(f"    Center: {pocket['center']}")
        print(f"    Radius: {pocket['radius']:.2f} Å")
        print(f"    Volume: {pocket['volume']:.2f} Å³")
    
    # Define active site based on first pocket
    if pockets:
        pocket = pockets[0]
        protein.define_active_site(pocket['center'], min(pocket['radius'], 10.0))
        print(f"Defined active site at {pocket['center']} with radius {pocket['radius']:.2f}")
    
    return protein, pockets

def test_sphere_generation():
    """Test sphere.pdb generation."""
    print("\nTesting sphere.pdb generation...")
    
    protein, pockets = test_pocket_detection()
    
    if not protein.active_site:
        print("No active site defined, cannot test sphere generation")
        return
    
    # Create temporary file for sphere
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as temp_file:
        sphere_file = temp_file.name
    
    try:
        handler = ProteinHandler()
        handler.save_binding_site_sphere(protein, sphere_file)
        
        # Check if file was created
        if os.path.exists(sphere_file):
            file_size = os.path.getsize(sphere_file)
            print(f"✓ sphere.pdb generated successfully ({file_size} bytes)")
            
            # Read first few lines to verify content
            with open(sphere_file, 'r') as f:
                lines = f.readlines()[:10]
                print("First few lines of sphere.pdb:")
                for line in lines:
                    print(f"  {line.strip()}")
        else:
            print("✗ sphere.pdb was not created")
    
    except Exception as e:
        print(f"✗ Error generating sphere.pdb: {e}")
    
    finally:
        # Clean up
        if os.path.exists(sphere_file):
            os.unlink(sphere_file)

def test_plotting():
    """Test plot generation functionality."""
    print("\nTesting plot generation...")
    
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        print("✓ Matplotlib is available")
        
        # Create some test data
        scores = np.random.normal(-5, 2, 50)
        ranks = list(range(1, len(scores) + 1))
        
        # Test basic plotting
        plt.figure(figsize=(8, 6))
        plt.plot(ranks, scores, 'o-')
        plt.title('Test Plot')
        plt.xlabel('Rank')
        plt.ylabel('Score')
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as temp_file:
            plot_file = temp_file.name
        
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        if os.path.exists(plot_file):
            file_size = os.path.getsize(plot_file)
            print(f"✓ Test plot generated successfully ({file_size} bytes)")
            os.unlink(plot_file)
        else:
            print("✗ Plot file was not created")
    
    except ImportError:
        print("✗ Matplotlib not available - plots will not be generated")
    except Exception as e:
        print(f"✗ Error during plotting test: {e}")

if __name__ == "__main__":
    print("=" * 60)
    print("PandaDock Pocket Detection and Visualization Test")
    print("=" * 60)
    
    test_pocket_detection()
    test_sphere_generation()
    test_plotting()
    
    print("\n" + "=" * 60)
    print("Test completed")
    print("=" * 60)