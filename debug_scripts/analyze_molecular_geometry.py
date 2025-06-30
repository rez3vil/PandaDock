#!/usr/bin/env python3
"""
Analyze molecular geometry preservation in generated poses.
"""

import numpy as np
from pathlib import Path

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

def parse_sdf_coordinates(sdf_file: str) -> np.ndarray:
    """Parse coordinates from SDF file."""
    coords = []
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 4:
        return np.array([])
    
    counts_line = lines[3].strip()
    n_atoms = int(counts_line[:3])
    
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

def calculate_internal_distances(coords: np.ndarray) -> np.ndarray:
    """Calculate all pairwise distances within molecule."""
    n_atoms = len(coords)
    distances = np.zeros((n_atoms, n_atoms))
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            distances[i, j] = dist
            distances[j, i] = dist
    
    return distances

def analyze_bond_preservation(original_coords: np.ndarray, generated_coords: np.ndarray) -> dict:
    """Analyze how well bond lengths are preserved."""
    if len(original_coords) != len(generated_coords):
        return {'error': 'Atom count mismatch'}
    
    orig_distances = calculate_internal_distances(original_coords)
    gen_distances = calculate_internal_distances(generated_coords)
    
    # Extract upper triangle (unique pairs)
    n_atoms = len(original_coords)
    orig_upper = []
    gen_upper = []
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            orig_upper.append(orig_distances[i, j])
            gen_upper.append(gen_distances[i, j])
    
    orig_upper = np.array(orig_upper)
    gen_upper = np.array(gen_upper)
    
    # Calculate differences
    distance_diffs = np.abs(orig_upper - gen_upper)
    
    # Identify likely bonds (distances < 2.0 Å in original)
    bond_mask = orig_upper < 2.0
    bond_orig = orig_upper[bond_mask]
    bond_gen = gen_upper[bond_mask]
    bond_diffs = distance_diffs[bond_mask]
    
    return {
        'n_atoms': n_atoms,
        'n_distances': len(orig_upper),
        'n_likely_bonds': len(bond_orig),
        'bond_rmsd': np.sqrt(np.mean(bond_diffs**2)) if len(bond_diffs) > 0 else 0,
        'max_bond_change': np.max(bond_diffs) if len(bond_diffs) > 0 else 0,
        'mean_bond_change': np.mean(bond_diffs) if len(bond_diffs) > 0 else 0,
        'distance_rmsd': np.sqrt(np.mean(distance_diffs**2)),
        'max_distance_change': np.max(distance_diffs),
        'severely_distorted_bonds': np.sum(bond_diffs > 0.5) if len(bond_diffs) > 0 else 0,
        'geometry_preserved': np.mean(bond_diffs) < 0.1 if len(bond_diffs) > 0 else False
    }

def check_for_triangular_clustering(coords: np.ndarray) -> dict:
    """Check for unnatural triangular clustering patterns."""
    n_atoms = len(coords)
    
    if n_atoms < 3:
        return {'triangular_clustering': False, 'reason': 'Too few atoms'}
    
    # Calculate center of mass
    center = np.mean(coords, axis=0)
    
    # Calculate distances from center
    center_distances = [np.linalg.norm(coord - center) for coord in coords]
    center_distances = np.array(center_distances)
    
    # Check for atoms clustered in small triangles
    triangular_patterns = 0
    
    # Look for sets of 3 atoms that form very small triangles
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            for k in range(j+1, n_atoms):
                # Calculate triangle side lengths
                d_ij = np.linalg.norm(coords[i] - coords[j])
                d_jk = np.linalg.norm(coords[j] - coords[k])
                d_ik = np.linalg.norm(coords[i] - coords[k])
                
                # Check if all sides are very short (< 1.5 Å)
                if d_ij < 1.5 and d_jk < 1.5 and d_ik < 1.5:
                    # Calculate triangle area using Heron's formula
                    s = (d_ij + d_jk + d_ik) / 2
                    area = np.sqrt(s * (s - d_ij) * (s - d_jk) * (s - d_ik))
                    
                    if area < 0.5:  # Very small triangle
                        triangular_patterns += 1
    
    # Check coordinate spreads
    coord_ranges = np.max(coords, axis=0) - np.min(coords, axis=0)
    min_range = np.min(coord_ranges)
    max_range = np.max(coord_ranges)
    aspect_ratio = max_range / min_range if min_range > 0 else float('inf')
    
    return {
        'triangular_clustering': triangular_patterns > 3,
        'n_small_triangles': triangular_patterns,
        'coord_ranges': coord_ranges,
        'aspect_ratio': aspect_ratio,
        'extreme_flatness': aspect_ratio > 20,
        'center_distances_std': np.std(center_distances),
        'center_distances_mean': np.mean(center_distances),
    }

def main():
    """Main analysis function."""
    print("Analyzing molecular geometry preservation...")
    
    # Load original ligand
    original_coords = parse_sdf_coordinates("tests/ligand.sdf")
    print(f"Original ligand: {len(original_coords)} atoms")
    
    if len(original_coords) == 0:
        print("Could not load original ligand!")
        return
    
    # Analyze original geometry
    orig_clustering = check_for_triangular_clustering(original_coords)
    print(f"Original triangular clustering: {orig_clustering['triangular_clustering']}")
    print(f"Original coord ranges: {orig_clustering['coord_ranges']}")
    print(f"Original aspect ratio: {orig_clustering['aspect_ratio']:.2f}")
    
    # Analyze generated poses
    poses_dir = Path("tests/new/poses")
    pose_files = list(poses_dir.glob("pose_*.pdb"))[:10]  # First 10 poses
    
    print(f"\nAnalyzing {len(pose_files)} generated poses...")
    print("=" * 80)
    
    distorted_poses = 0
    triangular_poses = 0
    
    for pose_file in pose_files:
        gen_coords = parse_pdb_coordinates(pose_file)
        
        if len(gen_coords) != len(original_coords):
            print(f"{pose_file.name}: Atom count mismatch ({len(gen_coords)} vs {len(original_coords)})")
            continue
        
        # Analyze bond preservation
        bond_analysis = analyze_bond_preservation(original_coords, gen_coords)
        
        # Check for triangular clustering
        clustering_analysis = check_for_triangular_clustering(gen_coords)
        
        is_distorted = not bond_analysis.get('geometry_preserved', False)
        has_triangular = clustering_analysis['triangular_clustering']
        
        if is_distorted:
            distorted_poses += 1
        if has_triangular:
            triangular_poses += 1
        
        print(f"\n{pose_file.name}:")
        print(f"  Bond RMSD: {bond_analysis.get('bond_rmsd', 0):.4f}")
        print(f"  Max bond change: {bond_analysis.get('max_bond_change', 0):.4f}")
        print(f"  Severely distorted bonds: {bond_analysis.get('severely_distorted_bonds', 0)}")
        print(f"  Geometry preserved: {'✅' if bond_analysis.get('geometry_preserved', False) else '❌'}")
        print(f"  Triangular clustering: {'❌' if has_triangular else '✅'}")
        print(f"  Small triangles: {clustering_analysis['n_small_triangles']}")
        print(f"  Aspect ratio: {clustering_analysis['aspect_ratio']:.2f}")
        print(f"  Extreme flatness: {'❌' if clustering_analysis['extreme_flatness'] else '✅'}")
    
    print(f"\n{'='*80}")
    print("SUMMARY:")
    print(f"Total poses analyzed: {len(pose_files)}")
    print(f"Poses with distorted geometry: {distorted_poses}")
    print(f"Poses with triangular clustering: {triangular_poses}")
    print(f"Geometry preservation rate: {(len(pose_files) - distorted_poses) / len(pose_files) * 100:.1f}%")

if __name__ == "__main__":
    main()