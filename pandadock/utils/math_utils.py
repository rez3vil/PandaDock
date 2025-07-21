# -*- coding: utf-8 -*-
"""
Mathematical utility functions for PandaDock
"""

import numpy as np
from typing import Tuple, Union, List
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize


def distance_matrix(coords1: np.ndarray, coords2: np.ndarray = None) -> np.ndarray:
    """
    Calculate distance matrix between two sets of coordinates
    
    Args:
        coords1: First set of coordinates (N x 3)
        coords2: Second set of coordinates (M x 3), optional
        
    Returns:
        Distance matrix (N x M) or (N x N) if coords2 is None
    """
    if coords2 is None:
        # Calculate pairwise distances within coords1
        distances = pdist(coords1, metric='euclidean')
        return squareform(distances)
    else:
        # Calculate distances between coords1 and coords2
        distances = np.zeros((len(coords1), len(coords2)))
        for i, coord1 in enumerate(coords1):
            for j, coord2 in enumerate(coords2):
                distances[i, j] = np.linalg.norm(coord1 - coord2)
        return distances


def rotation_matrix(angles: np.ndarray) -> np.ndarray:
    """
    Create rotation matrix from Euler angles (ZYX convention)
    
    Args:
        angles: Euler angles [alpha, beta, gamma] in radians
        
    Returns:
        3x3 rotation matrix
    """
    alpha, beta, gamma = angles
    
    # Rotation matrices for each axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(alpha), -np.sin(alpha)],
        [0, np.sin(alpha), np.cos(alpha)]
    ])
    
    Ry = np.array([
        [np.cos(beta), 0, np.sin(beta)],
        [0, 1, 0],
        [-np.sin(beta), 0, np.cos(beta)]
    ])
    
    Rz = np.array([
        [np.cos(gamma), -np.sin(gamma), 0],
        [np.sin(gamma), np.cos(gamma), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation matrix (ZYX order)
    R = Rz @ Ry @ Rx
    
    return R


def quaternion_to_matrix(quaternion: np.ndarray) -> np.ndarray:
    """
    Convert quaternion to rotation matrix
    
    Args:
        quaternion: Quaternion [w, x, y, z] (scalar first)
        
    Returns:
        3x3 rotation matrix
    """
    w, x, y, z = quaternion
    
    # Normalize quaternion
    norm = np.sqrt(w*w + x*x + y*y + z*z)
    if norm == 0:
        return np.eye(3)
    
    w, x, y, z = w/norm, x/norm, y/norm, z/norm
    
    # Rotation matrix from quaternion
    R = np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)],
        [2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)],
        [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)]
    ])
    
    return R


def matrix_to_quaternion(matrix: np.ndarray) -> np.ndarray:
    """
    Convert rotation matrix to quaternion
    
    Args:
        matrix: 3x3 rotation matrix
        
    Returns:
        Quaternion [w, x, y, z] (scalar first)
    """
    trace = np.trace(matrix)
    
    if trace > 0:
        s = np.sqrt(trace + 1.0) * 2  # s = 4 * w
        w = 0.25 * s
        x = (matrix[2, 1] - matrix[1, 2]) / s
        y = (matrix[0, 2] - matrix[2, 0]) / s
        z = (matrix[1, 0] - matrix[0, 1]) / s
    elif matrix[0, 0] > matrix[1, 1] and matrix[0, 0] > matrix[2, 2]:
        s = np.sqrt(1.0 + matrix[0, 0] - matrix[1, 1] - matrix[2, 2]) * 2  # s = 4 * x
        w = (matrix[2, 1] - matrix[1, 2]) / s
        x = 0.25 * s
        y = (matrix[0, 1] + matrix[1, 0]) / s
        z = (matrix[0, 2] + matrix[2, 0]) / s
    elif matrix[1, 1] > matrix[2, 2]:
        s = np.sqrt(1.0 + matrix[1, 1] - matrix[0, 0] - matrix[2, 2]) * 2  # s = 4 * y
        w = (matrix[0, 2] - matrix[2, 0]) / s
        x = (matrix[0, 1] + matrix[1, 0]) / s
        y = 0.25 * s
        z = (matrix[1, 2] + matrix[2, 1]) / s
    else:
        s = np.sqrt(1.0 + matrix[2, 2] - matrix[0, 0] - matrix[1, 1]) * 2  # s = 4 * z
        w = (matrix[1, 0] - matrix[0, 1]) / s
        x = (matrix[0, 2] + matrix[2, 0]) / s
        y = (matrix[1, 2] + matrix[2, 1]) / s
        z = 0.25 * s
    
    return np.array([w, x, y, z])


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray, align: bool = True) -> float:
    """
    Calculate Root Mean Square Deviation (RMSD) between two coordinate sets
    
    Args:
        coords1: First set of coordinates (N x 3)
        coords2: Second set of coordinates (N x 3)
        align: Whether to align coordinates before calculating RMSD
        
    Returns:
        RMSD value in Angstroms
    """
    if coords1.shape != coords2.shape:
        # Handle coordinate shape mismatch by using the minimum number of atoms
        min_atoms = min(len(coords1), len(coords2))
        if min_atoms == 0:
            return float('inf')  # No atoms to compare
        
        # Use only the first min_atoms from each coordinate set
        coords1 = coords1[:min_atoms]
        coords2 = coords2[:min_atoms]
        
        # Log warning if shapes were different
        import logging
        logger = logging.getLogger(__name__)
        logger.warning(f"Coordinate shape mismatch. Using first {min_atoms} atoms for RMSD calculation.")
    
    if len(coords1) == 0:
        return 0.0
    
    # Center coordinates
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)
    
    if align:
        # Align using Kabsch algorithm
        coords2_aligned = kabsch_align(coords1_centered, coords2_centered)
    else:
        coords2_aligned = coords2_centered
    
    # Calculate RMSD
    diff = coords1_centered - coords2_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd


def kabsch_align(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    """
    Align coords2 to coords1 using Kabsch algorithm
    
    Args:
        coords1: Reference coordinates (N x 3)
        coords2: Coordinates to align (N x 3)
        
    Returns:
        Aligned coords2
    """
    # Calculate cross-covariance matrix
    H = coords2.T @ coords1
    
    # SVD decomposition
    U, S, Vt = np.linalg.svd(H)
    
    # Calculate rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Apply rotation
    aligned_coords = coords2 @ R.T
    
    return aligned_coords


def angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
    """
    Calculate angle between two vectors in radians
    
    Args:
        v1: First vector
        v2: Second vector
        
    Returns:
        Angle in radians
    """
    # Normalize vectors
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    
    # Calculate angle
    cos_angle = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
    angle = np.arccos(cos_angle)
    
    return angle


def dihedral_angle(coords: np.ndarray) -> float:
    """
    Calculate dihedral angle for four atoms
    
    Args:
        coords: Coordinates of 4 atoms (4 x 3)
        
    Returns:
        Dihedral angle in radians
    """
    if coords.shape != (4, 3):
        raise ValueError("Dihedral angle requires exactly 4 atoms")
    
    # Vectors between atoms
    v1 = coords[1] - coords[0]
    v2 = coords[2] - coords[1]
    v3 = coords[3] - coords[2]
    
    # Normal vectors to planes
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    
    # Normalize
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    
    # Calculate angle
    cos_angle = np.clip(np.dot(n1, n2), -1.0, 1.0)
    angle = np.arccos(cos_angle)
    
    # Determine sign
    if np.dot(np.cross(n1, n2), v2) < 0:
        angle = -angle
    
    return angle


def center_of_mass(coords: np.ndarray, masses: np.ndarray = None) -> np.ndarray:
    """
    Calculate center of mass
    
    Args:
        coords: Coordinates (N x 3)
        masses: Atomic masses (N,), optional
        
    Returns:
        Center of mass coordinates
    """
    if masses is None:
        # Use equal masses (geometric center)
        return np.mean(coords, axis=0)
    else:
        # Weighted center of mass
        total_mass = np.sum(masses)
        if total_mass == 0:
            return np.mean(coords, axis=0)
        
        weighted_coords = coords * masses[:, np.newaxis]
        return np.sum(weighted_coords, axis=0) / total_mass


def radius_of_gyration(coords: np.ndarray, masses: np.ndarray = None) -> float:
    """
    Calculate radius of gyration
    
    Args:
        coords: Coordinates (N x 3)
        masses: Atomic masses (N,), optional
        
    Returns:
        Radius of gyration
    """
    com = center_of_mass(coords, masses)
    
    if masses is None:
        masses = np.ones(len(coords))
    
    # Calculate radius of gyration
    distances_squared = np.sum((coords - com)**2, axis=1)
    rg_squared = np.sum(masses * distances_squared) / np.sum(masses)
    
    return np.sqrt(rg_squared)


def principal_moments(coords: np.ndarray, masses: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate principal moments of inertia
    
    Args:
        coords: Coordinates (N x 3)
        masses: Atomic masses (N,), optional
        
    Returns:
        Tuple of (eigenvalues, eigenvectors) of inertia tensor
    """
    if masses is None:
        masses = np.ones(len(coords))
    
    # Center coordinates
    com = center_of_mass(coords, masses)
    centered_coords = coords - com
    
    # Calculate inertia tensor
    inertia_tensor = np.zeros((3, 3))
    
    for i, (coord, mass) in enumerate(zip(centered_coords, masses)):
        r_squared = np.dot(coord, coord)
        inertia_tensor += mass * (r_squared * np.eye(3) - np.outer(coord, coord))
    
    # Calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)
    
    # Sort by eigenvalue
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    return eigenvalues, eigenvectors


def transform_coordinates(coords: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    """
    Apply rotation and translation to coordinates
    
    Args:
        coords: Coordinates (N x 3)
        rotation: 3x3 rotation matrix
        translation: Translation vector (3,)
        
    Returns:
        Transformed coordinates
    """
    # Apply rotation
    rotated = coords @ rotation.T
    
    # Apply translation
    transformed = rotated + translation
    
    return transformed


def random_rotation_matrix() -> np.ndarray:
    """
    Generate random rotation matrix
    
    Returns:
        Random 3x3 rotation matrix
    """
    # Generate random quaternion
    q = np.random.randn(4)
    q = q / np.linalg.norm(q)
    
    # Convert to rotation matrix
    return quaternion_to_matrix(q)


def fibonacci_sphere(n_points: int) -> np.ndarray:
    """
    Generate points on sphere using Fibonacci spiral
    
    Args:
        n_points: Number of points to generate
        
    Returns:
        Points on unit sphere (n_points x 3)
    """
    points = np.zeros((n_points, 3))
    
    golden_ratio = (1 + 5**0.5) / 2
    
    for i in range(n_points):
        # Latitude
        lat = np.arcsin(-1 + 2 * i / (n_points - 1))
        
        # Longitude
        lon = 2 * np.pi * i / golden_ratio
        
        # Convert to Cartesian coordinates
        points[i, 0] = np.cos(lat) * np.cos(lon)
        points[i, 1] = np.cos(lat) * np.sin(lon)
        points[i, 2] = np.sin(lat)
    
    return points


def volume_of_convex_hull(coords: np.ndarray) -> float:
    """
    Calculate volume of convex hull
    
    Args:
        coords: Coordinates (N x 3)
        
    Returns:
        Volume of convex hull
    """
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(coords)
        return hull.volume
    except ImportError:
        # Fallback: approximate using bounding box
        min_coords = np.min(coords, axis=0)
        max_coords = np.max(coords, axis=0)
        dimensions = max_coords - min_coords
        return np.prod(dimensions)


def surface_area_of_convex_hull(coords: np.ndarray) -> float:
    """
    Calculate surface area of convex hull
    
    Args:
        coords: Coordinates (N x 3)
        
    Returns:
        Surface area of convex hull
    """
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(coords)
        return hull.area
    except ImportError:
        # Fallback: approximate using sphere
        radius = np.max(np.linalg.norm(coords - np.mean(coords, axis=0), axis=1))
        return 4 * np.pi * radius**2


def optimize_rigid_body_transformation(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find optimal rigid body transformation to align coords2 to coords1
    
    Args:
        coords1: Target coordinates (N x 3)
        coords2: Coordinates to transform (N x 3)
        
    Returns:
        Tuple of (rotation_matrix, translation_vector)
    """
    # Center coordinates
    center1 = np.mean(coords1, axis=0)
    center2 = np.mean(coords2, axis=0)
    
    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2
    
    # Find optimal rotation using Kabsch algorithm
    rotation = kabsch_align(coords1_centered, coords2_centered)
    
    # Calculate translation
    translation = center1 - center2 @ rotation.T
    
    return rotation, translation


def gaussian_overlap(center1: np.ndarray, sigma1: float, center2: np.ndarray, sigma2: float) -> float:
    """
    Calculate overlap between two Gaussian distributions
    
    Args:
        center1: Center of first Gaussian
        sigma1: Standard deviation of first Gaussian
        center2: Center of second Gaussian
        sigma2: Standard deviation of second Gaussian
        
    Returns:
        Overlap integral
    """
    # Distance between centers
    distance = np.linalg.norm(center1 - center2)
    
    # Combined variance
    combined_sigma = np.sqrt(sigma1**2 + sigma2**2)
    
    # Overlap integral
    overlap = np.exp(-0.5 * (distance / combined_sigma)**2)
    
    return overlap


def molecular_surface_area(coords: np.ndarray, radii: np.ndarray, probe_radius: float = 1.4) -> float:
    """
    Calculate molecular surface area (simplified)
    
    Args:
        coords: Atomic coordinates (N x 3)
        radii: Atomic radii (N,)
        probe_radius: Probe radius (water = 1.4 Å)
        
    Returns:
        Molecular surface area
    """
    # This is a simplified implementation
    # In practice, would use algorithms like MSMS or NACCESS
    
    total_area = 0.0
    
    for i, (coord, radius) in enumerate(zip(coords, radii)):
        # Surface area of sphere
        sphere_area = 4 * np.pi * (radius + probe_radius)**2
        
        # Approximate burial by nearby atoms
        buried_fraction = 0.0
        for j, (other_coord, other_radius) in enumerate(zip(coords, radii)):
            if i != j:
                distance = np.linalg.norm(coord - other_coord)
                overlap_distance = (radius + other_radius + 2 * probe_radius)
                
                if distance < overlap_distance:
                    buried_fraction += 0.1  # Simplified
        
        # Add exposed area
        exposed_area = sphere_area * max(0, 1 - buried_fraction)
        total_area += exposed_area
    
    return total_area


def grid_search_optimization(objective_function, bounds: List[Tuple[float, float]], n_points: int = 10) -> Tuple[np.ndarray, float]:
    """
    Simple grid search optimization
    
    Args:
        objective_function: Function to minimize
        bounds: List of (min, max) bounds for each dimension
        n_points: Number of points per dimension
        
    Returns:
        Tuple of (best_params, best_value)
    """
    best_params = None
    best_value = float('inf')
    
    # Generate grid points
    grid_points = []
    for bound in bounds:
        points = np.linspace(bound[0], bound[1], n_points)
        grid_points.append(points)
    
    # Evaluate all combinations
    import itertools
    
    for params in itertools.product(*grid_points):
        try:
            value = objective_function(np.array(params))
            if value < best_value:
                best_value = value
                best_params = np.array(params)
        except:
            continue
    
    return best_params, best_value