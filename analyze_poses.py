#!/usr/bin/env python3
"""
Script to analyze pose coordinates and identify triangular/clustering issues.
"""

import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple

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

def analyze_pose_geometry(coords: np.ndarray) -> Dict:
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
    
    # Calculate distances from center
    center_distances = [np.linalg.norm(coord - center) for coord in coords]
    
    # Calculate spread in each dimension
    coord_ranges = np.max(coords, axis=0) - np.min(coords, axis=0)
    
    return {
        'n_atoms': n_atoms,
        'center_of_mass': center,
        'coord_ranges': coord_ranges,
        'mean_pairwise_distance': np.mean(distances),
        'min_pairwise_distance': np.min(distances),
        'max_pairwise_distance': np.max(distances),
        'std_pairwise_distance': np.std(distances),
        'mean_center_distance': np.mean(center_distances),
        'max_center_distance': np.max(center_distances),
        'aspect_ratios': coord_ranges / np.min(coord_ranges) if np.min(coord_ranges) > 0 else [0, 0, 0]
    }

def check_triangular_clustering(coords: np.ndarray, threshold_ratio: float = 10.0) -> bool:
    """Check if coordinates show triangular/unnatural clustering."""
    analysis = analyze_pose_geometry(coords)
    
    # Check for extreme aspect ratios (flat/linear structures)
    max_aspect_ratio = np.max(analysis['aspect_ratios'])
    if max_aspect_ratio > threshold_ratio:
        return True
    
    # Check for unusually small minimum distances
    if analysis['min_pairwise_distance'] < 0.5:  # Less than typical bond length
        return True
    
    # Check for unusually large distances (scattered atoms)
    if analysis['max_pairwise_distance'] > 50.0:
        return True
        
    return False

def analyze_all_poses(poses_dir: str) -> None:
    """Analyze all poses in directory."""
    poses_dir = Path(poses_dir)
    
    if not poses_dir.exists():
        print(f"Directory not found: {poses_dir}")
        return
    
    pose_files = list(poses_dir.glob("pose_*.pdb"))
    pose_files.sort()
    
    print(f"Analyzing {len(pose_files)} poses in {poses_dir}")
    print("=" * 80)
    
    triangular_poses = []
    normal_poses = []
    
    for pose_file in pose_files:
        try:
            coords = parse_pdb_coordinates(pose_file)
            analysis = analyze_pose_geometry(coords)
            is_triangular = check_triangular_clustering(coords)
            
            if is_triangular:
                triangular_poses.append((pose_file.name, analysis))
            else:
                normal_poses.append((pose_file.name, analysis))
            
            # Print summary for first few poses
            if len(triangular_poses) + len(normal_poses) <= 5:
                print(f"\n{pose_file.name}:")
                print(f"  Center: [{analysis['center_of_mass'][0]:.2f}, {analysis['center_of_mass'][1]:.2f}, {analysis['center_of_mass'][2]:.2f}]")
                print(f"  Coord ranges: [{analysis['coord_ranges'][0]:.2f}, {analysis['coord_ranges'][1]:.2f}, {analysis['coord_ranges'][2]:.2f}]")
                print(f"  Aspect ratios: [{analysis['aspect_ratios'][0]:.2f}, {analysis['aspect_ratios'][1]:.2f}, {analysis['aspect_ratios'][2]:.2f}]")
                print(f"  Distance stats: mean={analysis['mean_pairwise_distance']:.2f}, min={analysis['min_pairwise_distance']:.2f}, max={analysis['max_pairwise_distance']:.2f}")
                print(f"  Triangular clustering: {'YES' if is_triangular else 'NO'}")
        
        except Exception as e:
            print(f"Error analyzing {pose_file}: {e}")
    
    print(f"\n{'='*80}")
    print(f"SUMMARY:")
    print(f"Total poses analyzed: {len(pose_files)}")
    print(f"Poses with triangular/distorted geometry: {len(triangular_poses)}")
    print(f"Poses with normal geometry: {len(normal_poses)}")
    
    if triangular_poses:
        print(f"\nDistorted poses:")
        for name, analysis in triangular_poses[:10]:  # Show first 10
            print(f"  {name}: aspect_ratios=[{analysis['aspect_ratios'][0]:.1f}, {analysis['aspect_ratios'][1]:.1f}, {analysis['aspect_ratios'][2]:.1f}], min_dist={analysis['min_pairwise_distance']:.2f}")
    
    # Visualize coordinate distributions
    create_visualization(triangular_poses, normal_poses, poses_dir)

def create_visualization(triangular_poses: List, normal_poses: List, poses_dir: Path) -> None:
    """Create visualization of coordinate patterns."""
    if not triangular_poses and not normal_poses:
        return
        
    plt.figure(figsize=(12, 8))
    
    # Plot aspect ratios
    plt.subplot(2, 2, 1)
    if triangular_poses:
        tri_aspects = [max(analysis['aspect_ratios']) for _, analysis in triangular_poses]
        plt.hist(tri_aspects, bins=20, alpha=0.7, label='Distorted poses', color='red')
    
    if normal_poses:
        norm_aspects = [max(analysis['aspect_ratios']) for _, analysis in normal_poses]
        plt.hist(norm_aspects, bins=20, alpha=0.7, label='Normal poses', color='blue')
    
    plt.xlabel('Max Aspect Ratio')
    plt.ylabel('Count')
    plt.title('Aspect Ratio Distribution')
    plt.legend()
    
    # Plot minimum distances
    plt.subplot(2, 2, 2)
    if triangular_poses:
        tri_min_dist = [analysis['min_pairwise_distance'] for _, analysis in triangular_poses]
        plt.hist(tri_min_dist, bins=20, alpha=0.7, label='Distorted poses', color='red')
    
    if normal_poses:
        norm_min_dist = [analysis['min_pairwise_distance'] for _, analysis in normal_poses]
        plt.hist(norm_min_dist, bins=20, alpha=0.7, label='Normal poses', color='blue')
    
    plt.xlabel('Minimum Pairwise Distance (Å)')
    plt.ylabel('Count')
    plt.title('Minimum Distance Distribution')
    plt.legend()
    
    # Plot coordinate ranges
    plt.subplot(2, 2, 3)
    if triangular_poses:
        tri_ranges = [np.mean(analysis['coord_ranges']) for _, analysis in triangular_poses]
        plt.hist(tri_ranges, bins=20, alpha=0.7, label='Distorted poses', color='red')
    
    if normal_poses:
        norm_ranges = [np.mean(analysis['coord_ranges']) for _, analysis in normal_poses]
        plt.hist(norm_ranges, bins=20, alpha=0.7, label='Normal poses', color='blue')
    
    plt.xlabel('Mean Coordinate Range (Å)')
    plt.ylabel('Count') 
    plt.title('Coordinate Range Distribution')
    plt.legend()
    
    # Plot center distances
    plt.subplot(2, 2, 4)
    if triangular_poses:
        tri_center = [analysis['max_center_distance'] for _, analysis in triangular_poses]
        plt.hist(tri_center, bins=20, alpha=0.7, label='Distorted poses', color='red')
    
    if normal_poses:
        norm_center = [analysis['max_center_distance'] for _, analysis in normal_poses]
        plt.hist(norm_center, bins=20, alpha=0.7, label='Normal poses', color='blue')
    
    plt.xlabel('Max Distance from Center (Å)')
    plt.ylabel('Count')
    plt.title('Center Distance Distribution')
    plt.legend()
    
    plt.tight_layout()
    
    # Save plot
    output_file = poses_dir / 'pose_analysis.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nVisualization saved to: {output_file}")
    
    plt.show()

if __name__ == "__main__":
    # Analyze different pose directories
    base_dir = Path("tests")
    
    for pose_dir in ["new/poses", "new1/poses", "results/poses"]:
        full_path = base_dir / pose_dir
        if full_path.exists():
            print(f"\n{'='*100}")
            print(f"ANALYZING: {full_path}")
            print(f"{'='*100}")
            analyze_all_poses(full_path)