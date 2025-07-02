"""
Mathematical utility functions for PandaDock.
"""

import numpy as np
from typing import Union, List, Tuple


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculate Root Mean Square Deviation between two coordinate sets.
    
    Args:
        coords1: First coordinate set (N, 3)
        coords2: Second coordinate set (N, 3)
        
    Returns:
        RMSD value
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    diff = coords1 - coords2
    squared_diff = np.sum(diff**2)
    msd = squared_diff / len(coords1)
    return np.sqrt(msd)


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    """
    Normalize a vector to unit length.
    
    Args:
        vector: Input vector
        
    Returns:
        Normalized vector
    """
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector
    return vector / norm


def rotation_matrix_from_euler(angles: np.ndarray) -> np.ndarray:
    """
    Create rotation matrix from Euler angles (ZYX convention).
    
    Args:
        angles: Euler angles [rx, ry, rz] in radians
        
    Returns:
        3x3 rotation matrix
    """
    rx, ry, rz = angles
    
    # Rotation around X-axis
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rx), -np.sin(rx)],
        [0, np.sin(rx), np.cos(rx)]
    ])
    
    # Rotation around Y-axis
    Ry = np.array([
        [np.cos(ry), 0, np.sin(ry)],
        [0, 1, 0],
        [-np.sin(ry), 0, np.cos(ry)]
    ])
    
    # Rotation around Z-axis
    Rz = np.array([
        [np.cos(rz), -np.sin(rz), 0],
        [np.sin(rz), np.cos(rz), 0],
        [0, 0, 1]
    ])
    
    # Combined rotation (ZYX order)
    return Rz @ Ry @ Rx


def apply_rotation_translation(coords: np.ndarray, rotation: np.ndarray, 
                             translation: np.ndarray) -> np.ndarray:
    """
    Apply rotation and translation to coordinates.
    
    Args:
        coords: Original coordinates (N, 3)
        rotation: 3x3 rotation matrix
        translation: Translation vector (3,)
        
    Returns:
        Transformed coordinates (N, 3)
    """
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("Coordinates must be (N, 3) array")
    if rotation.shape != (3, 3):
        raise ValueError("Rotation must be (3, 3) matrix")
    if translation.shape != (3,):
        raise ValueError("Translation must be (3,) vector")
    
    # Apply rotation then translation
    rotated_coords = coords @ rotation.T  # Apply rotation
    transformed_coords = rotated_coords + translation  # Apply translation
    
    return transformed_coords


def distance_matrix(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    """
    Calculate distance matrix between two coordinate sets.
    
    Args:
        coords1: First coordinate set (N, 3)
        coords2: Second coordinate set (M, 3)
        
    Returns:
        Distance matrix (N, M)
    """
    # Expand dimensions for broadcasting
    coords1_expanded = coords1[:, np.newaxis, :]  # (N, 1, 3)
    coords2_expanded = coords2[np.newaxis, :, :]  # (1, M, 3)
    
    # Calculate squared differences
    diff = coords1_expanded - coords2_expanded  # (N, M, 3)
    squared_distances = np.sum(diff**2, axis=2)  # (N, M)
    
    return np.sqrt(squared_distances)


def center_of_mass(coords: np.ndarray, masses: np.ndarray = None) -> np.ndarray:
    """
    Calculate center of mass for a set of coordinates.
    
    Args:
        coords: Coordinates (N, 3)
        masses: Atomic masses (N,). If None, assumes equal masses
        
    Returns:
        Center of mass coordinates (3,)
    """
    if masses is None:
        return np.mean(coords, axis=0)
    
    if len(masses) != len(coords):
        raise ValueError("Number of masses must match number of coordinates")
    
    total_mass = np.sum(masses)
    weighted_coords = coords * masses[:, np.newaxis]
    return np.sum(weighted_coords, axis=0) / total_mass


def align_structures(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Align two structures using Kabsch algorithm.
    
    Args:
        coords1: Reference coordinates (N, 3)
        coords2: Coordinates to align (N, 3)
        
    Returns:
        Tuple of (rotation_matrix, translation_vector)
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    # Center both structures
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    
    centered1 = coords1 - centroid1
    centered2 = coords2 - centroid2
    
    # Calculate covariance matrix
    H = centered2.T @ centered1
    
    # SVD
    U, S, Vt = np.linalg.svd(H)
    
    # Calculate rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Calculate translation
    t = centroid1 - R @ centroid2
    
    return R, t