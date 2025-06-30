#!/usr/bin/env python3
"""
Compare original ligand coordinates with generated poses.
"""

import numpy as np
from pathlib import Path

def parse_sdf_coordinates(sdf_file: str) -> np.ndarray:
    """Parse coordinates from SDF file."""
    coords = []
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines
    if len(lines) < 4:
        return np.array([])
    
    # Parse counts line
    counts_line = lines[3].strip()
    n_atoms = int(counts_line[:3])
    
    # Parse atoms
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            coords.append([x, y, z])
        except ValueError:
            continue
    
    return np.array(coords)

def parse_pdb_coordinates(pdb_file: str) -> np.ndarray:
    """Parse coordinates from PDB file."""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
    return np.array(coords)

def calculate_geometry_stats(coords: np.ndarray) -> dict:
    """Calculate basic geometry statistics."""
    if len(coords) == 0:
        return {}
    
    # Pairwise distances
    distances = []
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = np.linalg.norm(coords[i] - coords[j])
            distances.append(dist)
    
    distances = np.array(distances)
    
    # Center and spread
    center = np.mean(coords, axis=0)
    coord_ranges = np.max(coords, axis=0) - np.min(coords, axis=0)
    
    return {
        'center': center,
        'coord_ranges': coord_ranges,
        'mean_distance': np.mean(distances),
        'min_distance': np.min(distances),
        'max_distance': np.max(distances),
        'total_span': np.linalg.norm(coord_ranges)
    }

def analyze_transformation_preservation(original_coords: np.ndarray, 
                                      generated_coords: np.ndarray) -> dict:
    """Analyze how well molecular geometry is preserved."""
    if len(original_coords) != len(generated_coords):
        return {'error': 'Different number of atoms'}
    
    # Calculate internal distances for both
    orig_distances = []
    gen_distances = []
    
    n_atoms = len(original_coords)
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            orig_dist = np.linalg.norm(original_coords[i] - original_coords[j])
            gen_dist = np.linalg.norm(generated_coords[i] - generated_coords[j])
            orig_distances.append(orig_dist)
            gen_distances.append(gen_dist)
    
    orig_distances = np.array(orig_distances)
    gen_distances = np.array(gen_distances)
    
    # Calculate preservation metrics
    distance_rmsd = np.sqrt(np.mean((orig_distances - gen_distances)**2))
    distance_changes = np.abs(orig_distances - gen_distances)
    
    return {
        'distance_rmsd': distance_rmsd,
        'max_distance_change': np.max(distance_changes),
        'mean_distance_change': np.mean(distance_changes),
        'geometry_preserved': distance_rmsd < 0.1  # Threshold for good preservation
    }

if __name__ == "__main__":
    # Load original ligand
    original_sdf = Path("run/ligand.sdf")
    test_sdf = Path("tests/ligand.sdf")
    
    original_coords = None
    if original_sdf.exists():
        original_coords = parse_sdf_coordinates(original_sdf)
        print(f"Loaded original ligand from {original_sdf}: {len(original_coords)} atoms")
    elif test_sdf.exists():
        original_coords = parse_sdf_coordinates(test_sdf)
        print(f"Loaded original ligand from {test_sdf}: {len(original_coords)} atoms")
    else:
        print("No original ligand file found")
        exit(1)
    
    if len(original_coords) == 0:
        print("Could not parse original coordinates")
        exit(1)
    
    # Analyze original geometry
    orig_stats = calculate_geometry_stats(original_coords)
    print(f"\nOriginal ligand geometry:")
    print(f"  Center: [{orig_stats['center'][0]:.1f}, {orig_stats['center'][1]:.1f}, {orig_stats['center'][2]:.1f}]")
    print(f"  Coord ranges: {orig_stats['coord_ranges']}")
    print(f"  Distance stats: mean={orig_stats['mean_distance']:.2f}, min={orig_stats['min_distance']:.2f}, max={orig_stats['max_distance']:.2f}")
    print(f"  Total span: {orig_stats['total_span']:.2f}")
    
    # Analyze generated poses
    test_dirs = ["tests/new/poses", "tests/new1/poses"]
    
    for test_dir in test_dirs:
        pose_dir = Path(test_dir)
        if not pose_dir.exists():
            continue
            
        print(f"\n{'='*60}")
        print(f"ANALYZING: {pose_dir}")
        print(f"{'='*60}")
        
        pose_files = list(pose_dir.glob("pose_*.pdb"))[:5]  # First 5 poses
        
        geometry_preserved_count = 0
        total_valid_poses = 0
        
        for pose_file in pose_files:
            try:
                gen_coords = parse_pdb_coordinates(pose_file)
                
                if len(gen_coords) != len(original_coords):
                    print(f"{pose_file.name}: Atom count mismatch ({len(gen_coords)} vs {len(original_coords)})")
                    continue
                
                # Analyze geometry preservation
                preservation = analyze_transformation_preservation(original_coords, gen_coords)
                
                if 'error' in preservation:
                    print(f"{pose_file.name}: {preservation['error']}")
                    continue
                
                total_valid_poses += 1
                if preservation['geometry_preserved']:
                    geometry_preserved_count += 1
                
                gen_stats = calculate_geometry_stats(gen_coords)
                
                print(f"\n{pose_file.name}:")
                print(f"  Center: [{gen_stats['center'][0]:.1f}, {gen_stats['center'][1]:.1f}, {gen_stats['center'][2]:.1f}]")
                print(f"  Coord ranges: {gen_stats['coord_ranges']}")
                print(f"  Distance preservation RMSD: {preservation['distance_rmsd']:.4f}")
                print(f"  Max distance change: {preservation['max_distance_change']:.4f}")
                print(f"  Geometry preserved: {'YES' if preservation['geometry_preserved'] else 'NO'}")
                
            except Exception as e:
                print(f"Error analyzing {pose_file}: {e}")
        
        print(f"\n{test_dir} Summary:")
        print(f"  Valid poses analyzed: {total_valid_poses}")
        print(f"  Geometry preserved: {geometry_preserved_count}/{total_valid_poses}")
        if total_valid_poses > 0:
            print(f"  Preservation rate: {geometry_preserved_count/total_valid_poses:.1%}")