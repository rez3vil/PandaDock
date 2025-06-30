"""
Hardware detection and management module for PandaDock.
This module provides a unified interface for CPU and GPU operations.
"""

import os
import multiprocessing
import logging
from typing import Dict, Any, Optional, List, Union
from dataclasses import dataclass


@dataclass
class HardwareInfo:
    """Hardware information container."""
    cpu_cores: int
    cpu_threads: int
    memory_gb: float
    has_gpu: bool
    gpu_count: int = 0
    gpu_names: List[str] = None
    gpu_memory_gb: List[float] = None
    
    def __post_init__(self):
        if self.gpu_names is None:
            self.gpu_names = []
        if self.gpu_memory_gb is None:
            self.gpu_memory_gb = []


class HardwareDetector:
    """Detects available hardware resources."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def detect_hardware(self) -> HardwareInfo:
        """Detect all available hardware."""
        return HardwareInfo(
            cpu_cores=self._detect_cpu_cores(),
            cpu_threads=self._detect_cpu_threads(),
            memory_gb=self._detect_memory(),
            has_gpu=self._detect_gpu_availability(),
            gpu_count=self._detect_gpu_count(),
            gpu_names=self._detect_gpu_names(),
            gpu_memory_gb=self._detect_gpu_memory()
        )
    
    def _detect_cpu_cores(self) -> int:
        """Detect number of physical CPU cores."""
        try:
            import psutil
            return psutil.cpu_count(logical=False) or multiprocessing.cpu_count()
        except ImportError:
            return multiprocessing.cpu_count()
    
    def _detect_cpu_threads(self) -> int:
        """Detect number of CPU threads."""
        try:
            import psutil
            return psutil.cpu_count(logical=True) or multiprocessing.cpu_count()
        except ImportError:
            return multiprocessing.cpu_count()
    
    def _detect_memory(self) -> float:
        """Detect available system memory in GB."""
        try:
            import psutil
            return psutil.virtual_memory().total / (1024**3)
        except ImportError:
            return 8.0  # Default assumption
    
    def _detect_gpu_availability(self) -> bool:
        """Check if GPU is available."""
        # Check PyTorch
        try:
            import torch
            return torch.cuda.is_available()
        except ImportError:
            pass
        
        # Check TensorFlow
        try:
            import tensorflow as tf
            return len(tf.config.list_physical_devices('GPU')) > 0
        except ImportError:
            pass
        
        # Check CUDA directly
        try:
            import pynvml
            pynvml.nvmlInit()
            return pynvml.nvmlDeviceGetCount() > 0
        except ImportError:
            pass
        
        return False
    
    def _detect_gpu_count(self) -> int:
        """Detect number of available GPUs."""
        try:
            import torch
            return torch.cuda.device_count()
        except ImportError:
            pass
        
        try:
            import pynvml
            pynvml.nvmlInit()
            return pynvml.nvmlDeviceGetCount()
        except ImportError:
            pass
        
        return 0
    
    def _detect_gpu_names(self) -> List[str]:
        """Detect GPU names."""
        names = []
        try:
            import torch
            for i in range(torch.cuda.device_count()):
                names.append(torch.cuda.get_device_name(i))
            return names
        except ImportError:
            pass
        
        try:
            import pynvml
            pynvml.nvmlInit()
            for i in range(pynvml.nvmlDeviceGetCount()):
                handle = pynvml.nvmlDeviceGetHandleByIndex(i)
                name = pynvml.nvmlDeviceGetName(handle).decode('utf-8')
                names.append(name)
            return names
        except ImportError:
            pass
        
        return []
    
    def _detect_gpu_memory(self) -> List[float]:
        """Detect GPU memory in GB."""
        memory = []
        try:
            import torch
            for i in range(torch.cuda.device_count()):
                props = torch.cuda.get_device_properties(i)
                memory.append(props.total_memory / (1024**3))
            return memory
        except ImportError:
            pass
        
        try:
            import pynvml
            pynvml.nvmlInit()
            for i in range(pynvml.nvmlDeviceGetCount()):
                handle = pynvml.nvmlDeviceGetHandleByIndex(i)
                info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                memory.append(info.total / (1024**3))
            return memory
        except ImportError:
            pass
        
        return []


class HardwareManager:
    """Manages hardware resources and provides unified interface."""
    
    def __init__(self, prefer_gpu: bool = False, gpu_id: int = 0, cpu_workers: Optional[int] = None):
        """
        Initialize hardware manager.
        
        Parameters:
        -----------
        prefer_gpu : bool
            Whether to prefer GPU over CPU
        gpu_id : int
            GPU device ID to use
        cpu_workers : int, optional
            Number of CPU workers to use
        """
        self.logger = logging.getLogger(__name__)
        self.detector = HardwareDetector()
        self.hardware_info = self.detector.detect_hardware()
        
        self.prefer_gpu = prefer_gpu
        self.gpu_id = gpu_id
        self.cpu_workers = cpu_workers or self.hardware_info.cpu_cores
        
        self._setup_hardware()
    
    def _setup_hardware(self):
        """Setup hardware based on preferences and availability."""
        if self.prefer_gpu and self.hardware_info.has_gpu:
            self._setup_gpu()
        else:
            self._setup_cpu()
    
    def _setup_gpu(self):
        """Setup GPU for computations."""
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.set_device(self.gpu_id)
                self.device = f'cuda:{self.gpu_id}'
                self.logger.info(f"GPU setup complete: {torch.cuda.get_device_name(self.gpu_id)}")
            else:
                self.logger.warning("GPU requested but not available. Falling back to CPU.")
                self._setup_cpu()
        except ImportError:
            self.logger.warning("PyTorch not available. Falling back to CPU.")
            self._setup_cpu()
    
    def _setup_cpu(self):
        """Setup CPU for computations."""
        self.device = 'cpu'
        
        # Set CPU affinity if possible
        try:
            import psutil
            process = psutil.Process()
            available_cores = list(range(self.cpu_workers))
            process.cpu_affinity(available_cores)
            self.logger.info(f"CPU setup complete: {self.cpu_workers} workers")
        except ImportError:
            self.logger.info(f"CPU setup complete: {self.cpu_workers} workers (psutil not available)")
    
    def get_device(self) -> str:
        """Get the current device string."""
        return getattr(self, 'device', 'cpu')
    
    def is_gpu_available(self) -> bool:
        """Check if GPU is available and being used."""
        return hasattr(self, 'device') and 'cuda' in self.device
    
    def get_optimal_batch_size(self, problem_size: int, base_batch_size: int = 32) -> int:
        """
        Calculate optimal batch size based on hardware and problem size.
        
        Parameters:
        -----------
        problem_size : int
            Size of the problem (e.g., number of poses)
        base_batch_size : int
            Base batch size to scale from
            
        Returns:
        --------
        int
            Optimal batch size
        """
        if self.is_gpu_available():
            # GPU can handle larger batches
            gpu_mem_gb = self.hardware_info.gpu_memory_gb[self.gpu_id]
            if gpu_mem_gb > 8:
                return min(base_batch_size * 4, problem_size)
            elif gpu_mem_gb > 4:
                return min(base_batch_size * 2, problem_size)
            else:
                return min(base_batch_size, problem_size)
        else:
            # CPU batch size based on cores and memory
            cpu_factor = max(1, self.cpu_workers // 4)
            mem_factor = max(1, int(self.hardware_info.memory_gb // 8))
            return min(base_batch_size * min(cpu_factor, mem_factor), problem_size)
    
    def get_hardware_summary(self) -> Dict[str, Any]:
        """Get summary of hardware configuration."""
        return {
            'device': self.get_device(),
            'cpu_cores': self.hardware_info.cpu_cores,
            'cpu_threads': self.hardware_info.cpu_threads,
            'cpu_workers': self.cpu_workers,
            'memory_gb': self.hardware_info.memory_gb,
            'has_gpu': self.hardware_info.has_gpu,
            'gpu_count': self.hardware_info.gpu_count,
            'gpu_names': self.hardware_info.gpu_names,
            'gpu_memory_gb': self.hardware_info.gpu_memory_gb,
            'using_gpu': self.is_gpu_available()
        }
    
    def print_hardware_info(self):
        """Print detailed hardware information."""
        print("=" * 60)
        print("🔧 Hardware Configuration")
        print("=" * 60)
        
        # CPU Info
        print(f"🖥️  CPU: {self.hardware_info.cpu_cores} cores, {self.hardware_info.cpu_threads} threads")
        print(f"💾 Memory: {self.hardware_info.memory_gb:.1f} GB")
        print(f"⚙️  Workers: {self.cpu_workers}")
        
        # GPU Info
        if self.hardware_info.has_gpu:
            print(f"🎮 GPU: {self.hardware_info.gpu_count} device(s) available")
            for i, (name, mem) in enumerate(zip(self.hardware_info.gpu_names, self.hardware_info.gpu_memory_gb)):
                status = "✅ ACTIVE" if (self.is_gpu_available() and i == self.gpu_id) else "⏸️  IDLE"
                print(f"     Device {i}: {name} ({mem:.1f} GB) {status}")
        else:
            print("🎮 GPU: Not available")
        
        # Current configuration
        print(f"🚀 Active Device: {self.get_device().upper()}")
        print("=" * 60)
    
    def cleanup(self):
        """Cleanup hardware resources."""
        if self.is_gpu_available():
            try:
                import torch
                torch.cuda.empty_cache()
                self.logger.info("GPU cache cleared")
            except ImportError:
                pass
        
        self.logger.info("Hardware cleanup complete")


class AdaptiveProcessor:
    """Processor that adapts workload between CPU and GPU."""
    
    def __init__(self, hardware_manager: HardwareManager, workload_balance: float = 0.8):
        """
        Initialize adaptive processor.
        
        Parameters:
        -----------
        hardware_manager : HardwareManager
            Hardware manager instance
        workload_balance : float
            Balance between GPU and CPU (0.0 = all CPU, 1.0 = all GPU)
        """
        self.hardware_manager = hardware_manager
        self.workload_balance = workload_balance
        self.logger = logging.getLogger(__name__)
    
    def distribute_workload(self, total_work: int) -> Dict[str, int]:
        """
        Distribute workload between CPU and GPU.
        
        Parameters:
        -----------
        total_work : int
            Total amount of work to distribute
            
        Returns:
        --------
        Dict[str, int]
            Dictionary with 'gpu' and 'cpu' work amounts
        """
        if not self.hardware_manager.is_gpu_available():
            return {'gpu': 0, 'cpu': total_work}
        
        gpu_work = int(total_work * self.workload_balance)
        cpu_work = total_work - gpu_work
        
        self.logger.info(f"Workload distribution: GPU={gpu_work}, CPU={cpu_work}")
        
        return {'gpu': gpu_work, 'cpu': cpu_work}
    
    def process_batch(self, batch_data, processor_func, **kwargs):
        """
        Process a batch of data using optimal hardware.
        
        Parameters:
        -----------
        batch_data : list
            Data to process
        processor_func : callable
            Function to process the data
        **kwargs
            Additional arguments for processor function
            
        Returns:
        --------
        list
            Processed results
        """
        device = self.hardware_manager.get_device()
        batch_size = self.hardware_manager.get_optimal_batch_size(len(batch_data))
        
        results = []
        for i in range(0, len(batch_data), batch_size):
            batch = batch_data[i:i+batch_size]
            batch_results = processor_func(batch, device=device, **kwargs)
            results.extend(batch_results)
        
        return results


def detect_and_setup_hardware(prefer_gpu: bool = False, gpu_id: int = 0, 
                             cpu_workers: Optional[int] = None, 
                             verbose: bool = True) -> HardwareManager:
    """
    Convenience function to detect and setup hardware.
    
    Parameters:
    -----------
    prefer_gpu : bool
        Whether to prefer GPU over CPU
    gpu_id : int
        GPU device ID to use
    cpu_workers : int, optional
        Number of CPU workers
    verbose : bool
        Whether to print hardware info
        
    Returns:
    --------
    HardwareManager
        Configured hardware manager
    """
    manager = HardwareManager(prefer_gpu=prefer_gpu, gpu_id=gpu_id, cpu_workers=cpu_workers)
    
    if verbose:
        manager.print_hardware_info()
    
    return manager


def get_recommended_settings(target_performance: str = 'balanced') -> Dict[str, Any]:
    """
    Get recommended hardware settings based on target performance.
    
    Parameters:
    -----------
    target_performance : str
        Target performance level ('fast', 'balanced', 'accurate')
        
    Returns:
    --------
    Dict[str, Any]
        Recommended settings
    """
    detector = HardwareDetector()
    hardware_info = detector.detect_hardware()
    
    settings = {
        'use_gpu': hardware_info.has_gpu,
        'cpu_workers': hardware_info.cpu_cores,
        'workload_balance': 0.8
    }
    
    if target_performance == 'fast':
        settings.update({
            'use_gpu': hardware_info.has_gpu,
            'workload_balance': 0.9 if hardware_info.has_gpu else 0.0,
            'cpu_workers': max(1, hardware_info.cpu_cores // 2)
        })
    elif target_performance == 'accurate':
        settings.update({
            'use_gpu': False,  # CPU often more accurate for scientific computing
            'workload_balance': 0.0,
            'cpu_workers': hardware_info.cpu_cores
        })
    # 'balanced' uses default settings
    
    return settings
