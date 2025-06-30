"""
Create a sample PDB file using the fixed PDB writer to demonstrate proper formatting.
"""

import sys
import numpy as np
from pathlib import Path

# Add pandadock to path
sys.path.insert(0, '/Users/pritam/PandaDock')

def create_sample_pdb():
    """Create a sample PDB file that can be opened in ChimeraX."""
    
    from pandadock.io.pdb_writer import PDBWriter
    
    # Create sample ligand coordinates (a simple molecule)
    ligand_coords = np.array([
        [0.000,  0.000,  0.000],   # C1
        [1.540,  0.000,  0.000],   # C2
        [2.310,  1.335,  0.000],   # C3
        [1.540,  2.670,  0.000],   # C4
        [0.000,  2.670,  0.000],   # C5
        [-0.770, 1.335,  0.000],   # C6
    ])
    
    # Save the ligand
    output_file = Path('/Users/pritam/PandaDock/sample_ligand.pdb')
    
    PDBWriter.write_ligand_pose(
        filename=output_file,
        coordinates=ligand_coords,
        ligand_name="BZN",  # Benzene
        score=-12.34,
        pose_rank=1
    )
    
    print(f"✅ Created sample PDB file: {output_file}")
    print("\nTo test in ChimeraX:")
    print(f"1. Open ChimeraX")
    print(f"2. File > Open > {output_file}")
    print("3. The molecule should display correctly")
    
    # Show file content
    print("\n📄 File content:")
    with open(output_file, 'r') as f:
        content = f.read()
        print(content)
    
    return output_file

if __name__ == "__main__":
    create_sample_pdb()