#!/usr/bin/env python3
"""
Simple script to analyze pose coordinates and identify issues.
"""

import numpy as np
import os
from pathlib import Path

def parse_pdb_coordinates(pdb_file: str) -> np.ndarray:
    """Parse coordinates from PDB file."""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46]) 
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

def analyze_pose_geometry(coords: np.ndarray) -> dict:
    """Analyze geometric properties of pose coordinates."""
    n_atoms = len(coords)
    
    # Calculate pairwise distances
    distances = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            distances.append(dist)
    
    distances = np.array(distances)
    
    # Calculate center of mass
    center = np.mean(coords, axis=0)
    
    # Calculate spread in each dimension
    coord_ranges = np.max(coords, axis=0) - np.min(coords, axis=0)
    
    return {
        'n_atoms': n_atoms,
        'center_of_mass': center,
        'coord_ranges': coord_ranges,
        'mean_pairwise_distance': np.mean(distances),
        'min_pairwise_distance': np.min(distances),
        'max_pairwise_distance': np.max(distances),
        'aspect_ratios': coord_ranges / np.min(coord_ranges) if np.min(coord_ranges) > 0 else [0, 0, 0]
    }

def check_triangular_clustering(coords: np.ndarray) -> tuple:
    """Check if coordinates show triangular/unnatural clustering."""
    analysis = analyze_pose_geometry(coords)
    
    issues = []
    
    # Check for extreme aspect ratios (flat/linear structures)
    max_aspect_ratio = np.max(analysis['aspect_ratios'])
    if max_aspect_ratio > 10.0:
        issues.append(f"High aspect ratio: {max_aspect_ratio:.1f}")
    
    # Check for unusually small minimum distances
    if analysis['min_pairwise_distance'] < 0.5:
        issues.append(f"Very small min distance: {analysis['min_pairwise_distance']:.2f}")
    
    # Check for unusually large distances
    if analysis['max_pairwise_distance'] > 50.0:
        issues.append(f"Very large max distance: {analysis['max_pairwise_distance']:.2f}")
    
    # Check coordinate ranges
    min_range = np.min(analysis['coord_ranges'])
    if min_range < 0.1:
        issues.append(f"Very small coordinate range: {min_range:.2f}")
        
    return len(issues) > 0, issues, analysis

if __name__ == "__main__":
    # Analyze pose directories
    base_dir = Path("tests")
    
    for pose_dir in ["new/poses", "new1/poses", "results/poses"]:
        full_path = base_dir / pose_dir
        if not full_path.exists():
            continue
            
        print(f"\n{'='*80}")
        print(f"ANALYZING: {full_path}")
        print(f"{'='*80}")
        
        pose_files = list(full_path.glob("pose_*.pdb"))
        pose_files.sort()
        
        problematic_poses = 0
        total_poses = len(pose_files)
        
        for i, pose_file in enumerate(pose_files[:10]):  # Analyze first 10
            try:
                coords = parse_pdb_coordinates(pose_file)
                has_issues, issues, analysis = check_triangular_clustering(coords)
                
                if has_issues:
                    problematic_poses += 1
                
                print(f"\n{pose_file.name}:")
                print(f"  Atoms: {analysis['n_atoms']}")
                print(f"  Center: [{analysis['center_of_mass'][0]:.1f}, {analysis['center_of_mass'][1]:.1f}, {analysis['center_of_mass'][2]:.1f}]")
                print(f"  Coord ranges: {analysis['coord_ranges']}")
                print(f"  Aspect ratios: {analysis['aspect_ratios']}")
                print(f"  Distance range: {analysis['min_pairwise_distance']:.2f} - {analysis['max_pairwise_distance']:.2f}")
                print(f"  Issues: {', '.join(issues) if issues else 'None'}")
            
            except Exception as e:
                print(f"Error analyzing {pose_file}: {e}")
        
        print(f"\nSUMMARY for {pose_dir}:")
        print(f"  Analyzed: {min(10, total_poses)}/{total_poses} poses")
        print(f"  Problematic: {problematic_poses}")