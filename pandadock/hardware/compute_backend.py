"""
Unified compute backend for CPU/GPU operations.

This module provides a consistent interface for mathematical operations
across different compute devices, with automatic error handling and fallback.
"""

import numpy as np
import logging
from typing import Union, Optional, Any, Dict
from abc import ABC, abstractmethod

from .device_manager import DeviceManager, DeviceInfo


class ComputeBackendError(Exception):
    """Exception raised by compute backend operations."""
    pass


class BaseComputeBackend(ABC):
    """Abstract base class for compute backends."""
    
    @abstractmethod
    def compute_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """Compute pairwise distances between coordinate sets."""
        pass
    
    @abstractmethod
    def compute_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """Compute energies for given coordinates."""
        pass
    
    @abstractmethod
    def optimize_pose(self, initial_coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """Optimize molecular pose using local optimization."""
        pass


class CPUComputeBackend(BaseComputeBackend):
    """CPU-based compute backend using NumPy."""
    
    def __init__(self, device_info: DeviceInfo):
        """Initialize CPU compute backend."""
        self.device_info = device_info
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Initialized CPU backend: {device_info.name}")
    
    def compute_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """
        Compute pairwise distances using NumPy.
        
        Args:
            coords1: First set of coordinates (N, 3)
            coords2: Second set of coordinates (M, 3)
            
        Returns:
            Distance matrix (N, M)
        """
        try:
            # Efficient vectorized distance computation
            diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
            distances = np.sqrt(np.sum(diff**2, axis=2))
            return distances
            
        except Exception as e:
            raise ComputeBackendError(f"CPU distance computation failed: {e}")
    
    def compute_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Compute energies using CPU.
        
        Args:
            coords: Molecular coordinates
            params: Energy computation parameters
            
        Returns:
            Energy values
        """
        try:
            # Simple Lennard-Jones potential as example
            # In practice, this would use the actual scoring function
            
            # Extract parameters
            epsilon = params.get('epsilon', 1.0)
            sigma = params.get('sigma', 3.5)
            
            # Compute pairwise distances
            n_atoms = coords.shape[0]
            energies = np.zeros(n_atoms)
            
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    r = np.linalg.norm(coords[i] - coords[j])
                    if r > 0:
                        # Lennard-Jones potential
                        r6 = (sigma / r) ** 6
                        lj_energy = 4 * epsilon * (r6**2 - r6)
                        energies[i] += lj_energy
                        energies[j] += lj_energy
            
            return energies
            
        except Exception as e:
            raise ComputeBackendError(f"CPU energy computation failed: {e}")
    
    def optimize_pose(self, initial_coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Optimize pose using CPU-based optimization.
        
        Args:
            initial_coords: Initial coordinates
            params: Optimization parameters
            
        Returns:
            Optimized coordinates
        """
        try:
            from scipy.optimize import minimize
            
            def objective(coords_flat):
                coords = coords_flat.reshape(initial_coords.shape)
                energies = self.compute_energies(coords, params)
                return np.sum(energies)
            
            # Flatten coordinates for optimization
            x0 = initial_coords.flatten()
            
            # Run optimization
            result = minimize(objective, x0, method='L-BFGS-B')
            
            if result.success:
                return result.x.reshape(initial_coords.shape)
            else:
                self.logger.warning(f"CPU optimization did not converge: {result.message}")
                return initial_coords
                
        except Exception as e:
            raise ComputeBackendError(f"CPU pose optimization failed: {e}")


class GPUComputeBackend(BaseComputeBackend):
    """GPU-based compute backend using PyTorch or similar."""
    
    def __init__(self, device_info: DeviceInfo):
        """Initialize GPU compute backend."""
        self.device_info = device_info
        self.logger = logging.getLogger(__name__)
        
        # Initialize GPU backend
        self._initialize_gpu()
        
        self.logger.info(f"Initialized GPU backend: {device_info.name}")
    
    def _initialize_gpu(self):
        """Initialize GPU computation libraries."""
        try:
            if 'CUDA' in self.device_info.name:
                import torch
                self.torch = torch
                self.device = torch.device(f'cuda:{self.device_info.device_id}')
                
                # Test GPU availability
                test_tensor = torch.randn(10, 10, device=self.device)
                _ = torch.sum(test_tensor)
                
                self.backend_type = 'pytorch'
                
            else:
                # OpenCL fallback
                import pyopencl as cl
                self.cl = cl
                self.backend_type = 'opencl'
                
                # Initialize OpenCL context
                platforms = cl.get_platforms()
                devices = platforms[0].get_devices(cl.device_type.GPU)
                self.context = cl.Context([devices[self.device_info.device_id]])
                self.queue = cl.CommandQueue(self.context)
                
        except Exception as e:
            raise ComputeBackendError(f"GPU initialization failed: {e}")
    
    def compute_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """
        Compute pairwise distances using GPU.
        
        Args:
            coords1: First set of coordinates (N, 3)
            coords2: Second set of coordinates (M, 3)
            
        Returns:
            Distance matrix (N, M)
        """
        try:
            if self.backend_type == 'pytorch':
                return self._pytorch_distances(coords1, coords2)
            else:
                return self._opencl_distances(coords1, coords2)
                
        except Exception as e:
            raise ComputeBackendError(f"GPU distance computation failed: {e}")
    
    def _pytorch_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """Compute distances using PyTorch."""
        # Convert to PyTorch tensors
        tensor1 = self.torch.from_numpy(coords1.astype(np.float32)).to(self.device)
        tensor2 = self.torch.from_numpy(coords2.astype(np.float32)).to(self.device)
        
        # Compute distances
        diff = tensor1.unsqueeze(1) - tensor2.unsqueeze(0)
        distances = self.torch.sqrt(self.torch.sum(diff**2, dim=2))
        
        # Convert back to NumPy
        return distances.cpu().numpy()
    
    def _opencl_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """Compute distances using OpenCL."""
        # Simple OpenCL implementation
        # In practice, this would use optimized kernels
        
        n1, n2 = coords1.shape[0], coords2.shape[0]
        distances = np.zeros((n1, n2), dtype=np.float32)
        
        # Fallback to CPU for now
        diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
        distances = np.sqrt(np.sum(diff**2, axis=2)).astype(np.float32)
        
        return distances
    
    def compute_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Compute energies using GPU.
        
        Args:
            coords: Molecular coordinates
            params: Energy computation parameters
            
        Returns:
            Energy values
        """
        try:
            if self.backend_type == 'pytorch':
                return self._pytorch_energies(coords, params)
            else:
                return self._opencl_energies(coords, params)
                
        except Exception as e:
            raise ComputeBackendError(f"GPU energy computation failed: {e}")
    
    def _pytorch_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """Compute energies using PyTorch."""
        # Convert to PyTorch tensor
        coords_tensor = self.torch.from_numpy(coords.astype(np.float32)).to(self.device)
        
        # Extract parameters
        epsilon = params.get('epsilon', 1.0)
        sigma = params.get('sigma', 3.5)
        
        n_atoms = coords_tensor.shape[0]
        energies = self.torch.zeros(n_atoms, device=self.device)
        
        # Compute pairwise interactions
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r = self.torch.norm(coords_tensor[i] - coords_tensor[j])
                if r > 0:
                    r6 = (sigma / r) ** 6
                    lj_energy = 4 * epsilon * (r6**2 - r6)
                    energies[i] += lj_energy
                    energies[j] += lj_energy
        
        return energies.cpu().numpy()
    
    def _opencl_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """Compute energies using OpenCL."""
        # Fallback to CPU implementation for now
        # In practice, this would use optimized OpenCL kernels
        
        epsilon = params.get('epsilon', 1.0)
        sigma = params.get('sigma', 3.5)
        
        n_atoms = coords.shape[0]
        energies = np.zeros(n_atoms, dtype=np.float32)
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r = np.linalg.norm(coords[i] - coords[j])
                if r > 0:
                    r6 = (sigma / r) ** 6
                    lj_energy = 4 * epsilon * (r6**2 - r6)
                    energies[i] += lj_energy
                    energies[j] += lj_energy
        
        return energies
    
    def optimize_pose(self, initial_coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Optimize pose using GPU acceleration.
        
        Args:
            initial_coords: Initial coordinates
            params: Optimization parameters
            
        Returns:
            Optimized coordinates
        """
        try:
            if self.backend_type == 'pytorch':
                return self._pytorch_optimization(initial_coords, params)
            else:
                # Fallback to CPU for OpenCL
                cpu_backend = CPUComputeBackend(self.device_info)
                return cpu_backend.optimize_pose(initial_coords, params)
                
        except Exception as e:
            raise ComputeBackendError(f"GPU pose optimization failed: {e}")
    
    def _pytorch_optimization(self, initial_coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """Optimize pose using PyTorch."""
        # Convert to PyTorch tensor with gradient tracking
        coords = self.torch.from_numpy(initial_coords.astype(np.float32)).to(self.device)
        coords.requires_grad_(True)
        
        # Set up optimizer
        optimizer = self.torch.optim.Adam([coords], lr=0.01)
        
        # Optimization loop
        for _ in range(100):  # Max iterations
            optimizer.zero_grad()
            
            # Compute energy
            energy = self._compute_torch_energy(coords, params)
            
            # Backward pass
            energy.backward()
            
            # Optimization step
            optimizer.step()
            
            # Check convergence
            if energy.item() < 1e-6:
                break
        
        return coords.detach().cpu().numpy()
    
    def _compute_torch_energy(self, coords: 'torch.Tensor', params: Dict[str, Any]) -> 'torch.Tensor':
        """Compute energy using PyTorch tensors."""
        epsilon = params.get('epsilon', 1.0)
        sigma = params.get('sigma', 3.5)
        
        n_atoms = coords.shape[0]
        total_energy = self.torch.tensor(0.0, device=self.device)
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r = self.torch.norm(coords[i] - coords[j])
                if r > 0:
                    r6 = (sigma / r) ** 6
                    lj_energy = 4 * epsilon * (r6**2 - r6)
                    total_energy += lj_energy
        
        return total_energy


class ComputeBackend:
    """
    Main compute backend that automatically selects CPU or GPU implementation.
    
    Provides a unified interface for all compute operations with automatic
    error handling and fallback capabilities.
    """
    
    def __init__(self, device_manager: DeviceManager):
        """
        Initialize compute backend.
        
        Args:
            device_manager: Configured device manager instance
        """
        self.device_manager = device_manager
        self.logger = logging.getLogger(__name__)
        
        # Initialize appropriate backend
        self._initialize_backend()
    
    def _initialize_backend(self):
        """Initialize the appropriate compute backend."""
        device_info = self.device_manager.selected_device
        
        try:
            if device_info.device_type == 'gpu':
                self.backend = GPUComputeBackend(device_info)
                self.logger.info("Using GPU compute backend")
            else:
                self.backend = CPUComputeBackend(device_info)
                self.logger.info("Using CPU compute backend")
                
        except ComputeBackendError as e:
            self.logger.warning(f"Failed to initialize primary backend: {e}")
            
            # Force fallback to CPU
            self.device_manager.force_cpu_fallback("Backend initialization failed")
            cpu_device = self.device_manager.selected_device
            self.backend = CPUComputeBackend(cpu_device)
            self.logger.info("Using CPU fallback backend")
    
    def compute_distances(self, coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
        """
        Compute pairwise distances with automatic fallback.
        
        Args:
            coords1: First set of coordinates
            coords2: Second set of coordinates
            
        Returns:
            Distance matrix
        """
        try:
            return self.backend.compute_distances(coords1, coords2)
            
        except ComputeBackendError as e:
            if self.device_manager.is_gpu_selected:
                self.logger.warning(f"GPU distance computation failed: {e}")
                self._fallback_to_cpu()
                return self.backend.compute_distances(coords1, coords2)
            else:
                raise
    
    def compute_energies(self, coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Compute energies with automatic fallback.
        
        Args:
            coords: Molecular coordinates
            params: Energy computation parameters
            
        Returns:
            Energy values
        """
        try:
            return self.backend.compute_energies(coords, params)
            
        except ComputeBackendError as e:
            if self.device_manager.is_gpu_selected:
                self.logger.warning(f"GPU energy computation failed: {e}")
                self._fallback_to_cpu()
                return self.backend.compute_energies(coords, params)
            else:
                raise
    
    def optimize_pose(self, initial_coords: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
        """
        Optimize pose with automatic fallback.
        
        Args:
            initial_coords: Initial coordinates
            params: Optimization parameters
            
        Returns:
            Optimized coordinates
        """
        try:
            return self.backend.optimize_pose(initial_coords, params)
            
        except ComputeBackendError as e:
            if self.device_manager.is_gpu_selected:
                self.logger.warning(f"GPU pose optimization failed: {e}")
                self._fallback_to_cpu()
                return self.backend.optimize_pose(initial_coords, params)
            else:
                raise
    
    def _fallback_to_cpu(self):
        """Fall back to CPU backend."""
        self.device_manager.force_cpu_fallback("Compute operation failed")
        cpu_device = self.device_manager.selected_device
        self.backend = CPUComputeBackend(cpu_device)
        self.logger.info("Switched to CPU fallback backend")
    
    @property
    def device_info(self) -> DeviceInfo:
        """Get information about the current compute device."""
        return self.device_manager.selected_device
    
    def cleanup(self):
        """Clean up compute backend resources."""
        self.device_manager.cleanup()