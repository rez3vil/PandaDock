"""
Mathematical utility functions.
"""

import numpy as np
from typing import Union, Tuple

def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculate RMSD between two coordinate sets.
    
    Args:
        coords1: First coordinate set (N, 3)
        coords2: Second coordinate set (N, 3)
        
    Returns:
        RMSD value
    """
    
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
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

def rotation_matrix_from_euler(angles: Tuple[float, float, float]) -> np.ndarray:
    """
    Create rotation matrix from Euler angles.
    
    Args:
        angles: Euler angles (rx, ry, rz) in radians
        
    Returns:
        3x3 rotation matrix
    """
    
    rx, ry, rz = angles
    
    # Rotation matrices for each axis
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(rx), -np.sin(rx)],
                   [0, np.sin(rx), np.cos(rx)]])
    
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)],
                   [0, 1, 0],
                   [-np.sin(ry), 0, np.cos(ry)]])
    
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0],
                   [np.sin(rz), np.cos(rz), 0],
                   [0, 0, 1]])
    
    return np.dot(Rz, np.dot(Ry, Rx))

def apply_rotation_translation(coords: np.ndarray, rotation: np.ndarray, 
                             translation: np.ndarray) -> np.ndarray:
    """
    Apply rotation and translation to coordinates.
    
    Args:
        coords: Coordinate array (N, 3)
        rotation: 3x3 rotation matrix
        translation: Translation vector (3,)
        
    Returns:
        Transformed coordinates
    """
    if coords.ndim == 1:
        coords = coords.reshape(-1, 3)
    
    # Apply rotation then translation
    rotated = np.dot(coords, rotation.T)
    translated = rotated + translation
    
    return translated