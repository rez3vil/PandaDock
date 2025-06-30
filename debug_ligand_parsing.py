#!/usr/bin/env python3
"""
Debug ligand parsing to identify the coordinate loading issue.
"""

import numpy as np
from pathlib import Path

def debug_simple_mol_parser(file_path: str):
    """Debug the simple MOL file parser."""
    print(f"Debugging parser for: {file_path}")
    
    coords = []
    atom_names = []
    atom_types = []
    bonds = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    print(f"Total lines in file: {len(lines)}")
    
    # Skip header lines
    if len(lines) < 4:
        print("ERROR: Invalid MOL file format")
        return
    
    # Parse counts line
    counts_line = lines[3].strip()
    print(f"Counts line: '{counts_line}'")
    
    n_atoms = int(counts_line[:3])
    n_bonds = int(counts_line[3:6]) if len(counts_line) >= 6 else 0
    
    print(f"Expected atoms: {n_atoms}")
    print(f"Expected bonds: {n_bonds}")
    
    # Parse atoms
    print(f"\nParsing atoms from lines 4 to {4 + n_atoms - 1}:")
    
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            print(f"  Line {i}: MISSING (file too short)")
            break
        
        line = lines[i]
        print(f"  Line {i}: '{line.strip()}'")
        
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            atom_type = line[31:34].strip()
            
            coords.append([x, y, z])
            atom_names.append(f"{atom_type}{i-3}")
            atom_types.append(atom_type)
            
            print(f"    Parsed: x={x:.3f}, y={y:.3f}, z={z:.3f}, type='{atom_type}'")
            
        except (ValueError, IndexError) as e:
            print(f"    ERROR parsing line {i}: {e}")
            break
    
    print(f"\nSuccessfully parsed {len(coords)} atoms")
    
    # Check bond section
    print(f"\nBond section starts at line {4 + n_atoms}:")
    bond_start = 4 + n_atoms
    
    for i in range(bond_start, min(bond_start + 5, len(lines))):  # Show first 5 bond lines
        if i < len(lines):
            print(f"  Line {i}: '{lines[i].strip()}'")
    
    return len(coords), n_atoms

if __name__ == "__main__":
    # Test both ligand files
    test_files = [
        "run/ligand.sdf",
        "tests/ligand.sdf"
    ]
    
    for file_path in test_files:
        if Path(file_path).exists():
            print("=" * 80)
            parsed_atoms, expected_atoms = debug_simple_mol_parser(file_path)
            print(f"\nRESULT: Parsed {parsed_atoms}/{expected_atoms} atoms")
            if parsed_atoms != expected_atoms:
                print("❌ ATOM COUNT MISMATCH!")
            else:
                print("✅ All atoms parsed correctly")
            print()
        else:
            print(f"File not found: {file_path}")