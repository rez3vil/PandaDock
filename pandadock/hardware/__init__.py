"""
Hardware abstraction layer for PandaDock.

This module provides unified CPU/GPU compute interfaces and automatic
hardware detection and configuration.
"""

from .device_manager import DeviceManager
from .compute_backend import ComputeBackend
from .performance_monitor import PerformanceMonitor

__all__ = ['DeviceManager', 'ComputeBackend', 'PerformanceMonitor']