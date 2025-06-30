"""
Device management for CPU/GPU resources.

This module handles automatic detection and configuration of available
compute devices, providing transparent CPU/GPU switching.
"""

import os
import logging
from typing import Optional, Dict, Any, List
from dataclasses import dataclass


@dataclass
class DeviceInfo:
    """Information about a compute device."""
    device_type: str  # 'cpu' or 'gpu'
    device_id: int
    name: str
    memory_gb: float
    is_available: bool = True


class DeviceManager:
    """
    Manages CPU and GPU resources for molecular docking computations.
    
    Provides automatic device detection, configuration, and resource allocation
    with graceful fallback from GPU to CPU when needed.
    """
    
    def __init__(self, prefer_gpu: bool = True, gpu_id: Optional[int] = None):
        """
        Initialize device manager.
        
        Args:
            prefer_gpu: Whether to prefer GPU over CPU if available
            gpu_id: Specific GPU device ID to use (None for auto-select)
        """
        self.logger = logging.getLogger(__name__)
        self.prefer_gpu = prefer_gpu
        self.requested_gpu_id = gpu_id
        
        # Device information
        self._cpu_info: Optional[DeviceInfo] = None
        self._gpu_devices: List[DeviceInfo] = []
        self._selected_device: Optional[DeviceInfo] = None
        
        # Initialize devices
        self._detect_devices()
        self._select_optimal_device()
        
        # Initialize compute backend
        self._compute_backend = None
    
    def _detect_devices(self) -> None:
        """Detect available CPU and GPU devices."""
        # Always detect CPU
        self._detect_cpu()
        
        # Try to detect GPU devices
        if self.prefer_gpu:
            self._detect_gpu_devices()
    
    def _detect_cpu(self) -> None:
        """Detect CPU information."""
        try:
            import psutil
            
            # Get CPU information
            cpu_count = psutil.cpu_count(logical=True)
            memory_gb = psutil.virtual_memory().total / (1024**3)
            
            self._cpu_info = DeviceInfo(
                device_type='cpu',
                device_id=0,
                name=f'CPU ({cpu_count} cores)',
                memory_gb=memory_gb,
                is_available=True
            )
            
            self.logger.info(f"Detected CPU: {self._cpu_info.name}, "
                           f"Memory: {memory_gb:.1f} GB")
            
        except ImportError:
            # Fallback without psutil
            self._cpu_info = DeviceInfo(
                device_type='cpu',
                device_id=0,
                name='CPU',
                memory_gb=8.0,  # Conservative estimate
                is_available=True
            )
            self.logger.warning("psutil not available, using CPU fallback")
    
    def _detect_gpu_devices(self) -> None:
        """Detect available GPU devices."""
        # Try CUDA first
        if self._detect_cuda_devices():
            return
        
        # Try OpenCL as fallback
        self._detect_opencl_devices()
    
    def _detect_cuda_devices(self) -> bool:
        """Detect CUDA-compatible GPUs."""
        try:
            import torch
            
            if not torch.cuda.is_available():
                self.logger.info("CUDA not available")
                return False
            
            device_count = torch.cuda.device_count()
            self.logger.info(f"Found {device_count} CUDA device(s)")
            
            for i in range(device_count):
                props = torch.cuda.get_device_properties(i)
                memory_gb = props.total_memory / (1024**3)
                
                device_info = DeviceInfo(
                    device_type='gpu',
                    device_id=i,
                    name=f'CUDA:{i} ({props.name})',
                    memory_gb=memory_gb,
                    is_available=True
                )
                
                self._gpu_devices.append(device_info)
                self.logger.info(f"Detected GPU {i}: {props.name}, "
                               f"Memory: {memory_gb:.1f} GB")
            
            return True
            
        except ImportError:
            self.logger.info("PyTorch not available for CUDA detection")
            return False
        except Exception as e:
            self.logger.warning(f"Error detecting CUDA devices: {e}")
            return False
    
    def _detect_opencl_devices(self) -> None:
        """Detect OpenCL-compatible devices."""
        try:
            import pyopencl as cl
            
            platforms = cl.get_platforms()
            gpu_count = 0
            
            for platform in platforms:
                devices = platform.get_devices(cl.device_type.GPU)
                for device in devices:
                    memory_gb = device.global_mem_size / (1024**3)
                    
                    device_info = DeviceInfo(
                        device_type='gpu',
                        device_id=gpu_count,
                        name=f'OpenCL:{gpu_count} ({device.name.strip()})',
                        memory_gb=memory_gb,
                        is_available=True
                    )
                    
                    self._gpu_devices.append(device_info)
                    self.logger.info(f"Detected OpenCL GPU {gpu_count}: "
                                   f"{device.name.strip()}, "
                                   f"Memory: {memory_gb:.1f} GB")
                    gpu_count += 1
            
            if gpu_count == 0:
                self.logger.info("No OpenCL GPU devices found")
                
        except ImportError:
            self.logger.info("PyOpenCL not available")
        except Exception as e:
            self.logger.warning(f"Error detecting OpenCL devices: {e}")
    
    def _select_optimal_device(self) -> None:
        """Select the best available device for computation."""
        if self.prefer_gpu and self._gpu_devices:
            # Select requested GPU or best available GPU
            if self.requested_gpu_id is not None:
                # Use specific GPU if available
                for device in self._gpu_devices:
                    if device.device_id == self.requested_gpu_id:
                        self._selected_device = device
                        self.logger.info(f"Selected requested GPU: {device.name}")
                        return
                
                self.logger.warning(f"Requested GPU {self.requested_gpu_id} not found, "
                                  "selecting best available GPU")
            
            # Select GPU with most memory
            best_gpu = max(self._gpu_devices, key=lambda d: d.memory_gb)
            self._selected_device = best_gpu
            self.logger.info(f"Selected best GPU: {best_gpu.name}")
            
        else:
            # Use CPU
            self._selected_device = self._cpu_info
            self.logger.info(f"Selected CPU: {self._cpu_info.name}")
    
    @property
    def selected_device(self) -> DeviceInfo:
        """Get the currently selected compute device."""
        return self._selected_device
    
    @property
    def is_gpu_selected(self) -> bool:
        """Check if a GPU is currently selected."""
        return self._selected_device.device_type == 'gpu'
    
    @property
    def is_cpu_selected(self) -> bool:
        """Check if CPU is currently selected."""
        return self._selected_device.device_type == 'cpu'
    
    def get_device_config(self) -> Dict[str, Any]:
        """
        Get device configuration for compute backends.
        
        Returns:
            Dictionary with device configuration parameters
        """
        device = self._selected_device
        
        return {
            'device_type': device.device_type,
            'device_id': device.device_id,
            'device_name': device.name,
            'memory_gb': device.memory_gb,
            'use_gpu': device.device_type == 'gpu'
        }
    
    def force_cpu_fallback(self, reason: str = "GPU computation failed") -> None:
        """
        Force fallback to CPU computation.
        
        Args:
            reason: Reason for the fallback
        """
        if self._cpu_info:
            self.logger.warning(f"Falling back to CPU: {reason}")
            self._selected_device = self._cpu_info
        else:
            raise RuntimeError("No CPU device available for fallback")
    
    def get_memory_info(self) -> Dict[str, float]:
        """
        Get memory information for the selected device.
        
        Returns:
            Dictionary with memory information in GB
        """
        if self.is_gpu_selected:
            try:
                if 'CUDA' in self._selected_device.name:
                    import torch
                    torch.cuda.empty_cache()
                    
                    device_id = self._selected_device.device_id
                    total_memory = torch.cuda.get_device_properties(device_id).total_memory
                    allocated_memory = torch.cuda.memory_allocated(device_id)
                    cached_memory = torch.cuda.memory_reserved(device_id)
                    
                    return {
                        'total_gb': total_memory / (1024**3),
                        'allocated_gb': allocated_memory / (1024**3),
                        'cached_gb': cached_memory / (1024**3),
                        'free_gb': (total_memory - allocated_memory) / (1024**3)
                    }
                    
            except Exception as e:
                self.logger.warning(f"Could not get GPU memory info: {e}")
        
        # CPU memory info
        try:
            import psutil
            memory = psutil.virtual_memory()
            return {
                'total_gb': memory.total / (1024**3),
                'available_gb': memory.available / (1024**3),
                'used_gb': memory.used / (1024**3),
                'free_gb': memory.free / (1024**3)
            }
        except ImportError:
            return {'total_gb': self._selected_device.memory_gb}
    
    def cleanup(self) -> None:
        """Clean up device resources."""
        if self.is_gpu_selected:
            try:
                if 'CUDA' in self._selected_device.name:
                    import torch
                    torch.cuda.empty_cache()
                    self.logger.info("GPU memory cache cleared")
            except Exception as e:
                self.logger.warning(f"Error during GPU cleanup: {e}")
    
    def get_device_info(self) -> Dict[str, Any]:
        """
        Get comprehensive device information.
        
        Returns:
            Dictionary with device information
        """
        device = self._selected_device
        
        info = {
            'selected_device': device.name,
            'device_type': device.device_type,
            'device_id': device.device_id,
            'memory_gb': device.memory_gb,
            'is_available': device.is_available,
            'total_gpus': len(self._gpu_devices),
            'gpu_devices': [gpu.name for gpu in self._gpu_devices],
            'cpu_info': self._cpu_info.name if self._cpu_info else "Unknown"
        }
        
        # Add memory info
        try:
            memory_info = self.get_memory_info()
            info['memory_info'] = memory_info
        except Exception as e:
            self.logger.debug(f"Could not get memory info: {e}")
            info['memory_info'] = {}
        
        return info
    
    @property
    def compute_backend(self):
        """Get or create compute backend."""
        if self._compute_backend is None:
            from .compute_backend import ComputeBackend
            self._compute_backend = ComputeBackend(self)
        return self._compute_backend