"""
Performance monitoring for compute operations.

This module provides tools to monitor and optimize performance across
different compute devices and configurations.
"""

import time
import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from collections import defaultdict
import threading


@dataclass
class PerformanceMetric:
    """Single performance measurement."""
    operation: str
    duration: float
    device_type: str
    memory_used: float = 0.0
    success: bool = True
    error_message: Optional[str] = None
    timestamp: float = field(default_factory=time.time)


class PerformanceMonitor:
    """
    Monitor and analyze performance of compute operations.
    
    Tracks execution times, memory usage, and success rates across
    different devices and operation types.
    """
    
    def __init__(self, device_manager=None):
        """
        Initialize performance monitor.
        
        Args:
            device_manager: Device manager instance for memory monitoring (optional)
        """
        self.device_manager = device_manager
        self.logger = logging.getLogger(__name__)
        
        # Performance data
        self._metrics: List[PerformanceMetric] = []
        self._operation_stats: Dict[str, Dict[str, List[float]]] = defaultdict(
            lambda: defaultdict(list)
        )
        
        # Thread safety
        self._lock = threading.Lock()
        
        # Monitoring state
        self._monitoring_enabled = True
    
    def start_operation(self, operation_name: str) -> 'OperationTimer':
        """
        Start timing an operation.
        
        Args:
            operation_name: Name of the operation being timed
            
        Returns:
            Context manager for timing the operation
        """
        return OperationTimer(self, operation_name)
    
    def record_metric(self, metric: PerformanceMetric) -> None:
        """
        Record a performance metric.
        
        Args:
            metric: Performance metric to record
        """
        if not self._monitoring_enabled:
            return
        
        with self._lock:
            self._metrics.append(metric)
            
            # Update operation statistics
            key = f"{metric.operation}_{metric.device_type}"
            self._operation_stats[key]['durations'].append(metric.duration)
            self._operation_stats[key]['memory'].append(metric.memory_used)
            self._operation_stats[key]['success'].append(1.0 if metric.success else 0.0)
    
    def get_operation_stats(self, operation: str, device_type: Optional[str] = None) -> Dict[str, Any]:
        """
        Get statistics for a specific operation.
        
        Args:
            operation: Operation name
            device_type: Specific device type (None for all devices)
            
        Returns:
            Dictionary with operation statistics
        """
        with self._lock:
            if device_type:
                key = f"{operation}_{device_type}"
                stats_data = self._operation_stats.get(key, {})
            else:
                # Aggregate across all device types
                stats_data = defaultdict(list)
                for key, data in self._operation_stats.items():
                    if key.startswith(f"{operation}_"):
                        for metric_type, values in data.items():
                            stats_data[metric_type].extend(values)
            
            if not stats_data or not stats_data.get('durations'):
                return {'count': 0}
            
            durations = stats_data['durations']
            memory_usage = stats_data.get('memory', [])
            success_rates = stats_data.get('success', [])
            
            return {
                'count': len(durations),
                'avg_duration': sum(durations) / len(durations),
                'min_duration': min(durations),
                'max_duration': max(durations),
                'total_duration': sum(durations),
                'avg_memory_mb': sum(memory_usage) / len(memory_usage) if memory_usage else 0.0,
                'success_rate': sum(success_rates) / len(success_rates) if success_rates else 1.0,
                'operations_per_second': len(durations) / sum(durations) if sum(durations) > 0 else 0.0
            }
    
    def get_device_performance_comparison(self) -> Dict[str, Dict[str, Any]]:
        """
        Compare performance across different devices.
        
        Returns:
            Dictionary comparing performance metrics by device type
        """
        device_performance = {}
        
        # Get unique device types
        device_types = set()
        for key in self._operation_stats.keys():
            device_type = key.split('_')[-1]
            device_types.add(device_type)
        
        for device_type in device_types:
            device_stats = {
                'total_operations': 0,
                'total_duration': 0.0,
                'avg_duration': 0.0,
                'success_rate': 0.0,
                'operations_per_second': 0.0
            }
            
            all_durations = []
            all_success = []
            
            for key, data in self._operation_stats.items():
                if key.endswith(f"_{device_type}"):
                    durations = data.get('durations', [])
                    success = data.get('success', [])
                    
                    all_durations.extend(durations)
                    all_success.extend(success)
            
            if all_durations:
                device_stats['total_operations'] = len(all_durations)
                device_stats['total_duration'] = sum(all_durations)
                device_stats['avg_duration'] = sum(all_durations) / len(all_durations)
                device_stats['success_rate'] = sum(all_success) / len(all_success) if all_success else 1.0
                device_stats['operations_per_second'] = len(all_durations) / sum(all_durations)
            
            device_performance[device_type] = device_stats
        
        return device_performance
    
    def get_recent_metrics(self, minutes: int = 10) -> List[PerformanceMetric]:
        """
        Get metrics from the recent time window.
        
        Args:
            minutes: Number of minutes to look back
            
        Returns:
            List of recent performance metrics
        """
        cutoff_time = time.time() - (minutes * 60)
        
        with self._lock:
            return [metric for metric in self._metrics if metric.timestamp >= cutoff_time]
    
    def detect_performance_issues(self) -> List[str]:
        """
        Detect potential performance issues.
        
        Returns:
            List of performance issue descriptions
        """
        issues = []
        
        # Check for high failure rates
        device_performance = self.get_device_performance_comparison()
        
        for device_type, stats in device_performance.items():
            if stats['success_rate'] < 0.8:
                issues.append(f"High failure rate on {device_type}: {stats['success_rate']:.1%}")
            
            if stats['operations_per_second'] < 1.0 and stats['total_operations'] > 10:
                issues.append(f"Low performance on {device_type}: {stats['operations_per_second']:.2f} ops/sec")
        
        # Check for memory issues
        recent_metrics = self.get_recent_metrics(5)  # Last 5 minutes
        memory_failures = [m for m in recent_metrics if not m.success and 'memory' in (m.error_message or '').lower()]
        
        if len(memory_failures) > 3:
            issues.append(f"Multiple memory-related failures detected: {len(memory_failures)} in last 5 minutes")
        
        # Check for GPU fallbacks
        gpu_failures = [m for m in recent_metrics if not m.success and m.device_type == 'gpu']
        if len(gpu_failures) > 2:
            issues.append(f"Frequent GPU failures detected: {len(gpu_failures)} in last 5 minutes")
        
        return issues
    
    def optimize_device_selection(self) -> str:
        """
        Suggest optimal device based on performance history.
        
        Returns:
            Recommended device type ('cpu' or 'gpu')
        """
        device_performance = self.get_device_performance_comparison()
        
        if not device_performance:
            return 'cpu'  # Default fallback
        
        # Score devices based on success rate and performance
        device_scores = {}
        
        for device_type, stats in device_performance.items():
            if stats['total_operations'] < 5:
                continue  # Not enough data
            
            # Calculate composite score
            success_weight = 0.6
            performance_weight = 0.4
            
            success_score = stats['success_rate']
            performance_score = min(stats['operations_per_second'] / 10.0, 1.0)  # Normalize to 0-1
            
            composite_score = (success_weight * success_score + 
                             performance_weight * performance_score)
            
            device_scores[device_type] = composite_score
        
        if device_scores:
            best_device = max(device_scores.items(), key=lambda x: x[1])[0]
            return best_device
        
        return 'cpu'  # Safe default
    
    def generate_performance_report(self) -> str:
        """
        Generate a comprehensive performance report.
        
        Returns:
            Formatted performance report string
        """
        report_lines = [
            "=== PandaDock Performance Report ===",
            ""
        ]
        
        # Overall statistics
        with self._lock:
            total_operations = len(self._metrics)
            if total_operations > 0:
                total_duration = sum(m.duration for m in self._metrics)
                avg_duration = total_duration / total_operations
                success_rate = sum(1 for m in self._metrics if m.success) / total_operations
                
                report_lines.extend([
                    f"Total Operations: {total_operations}",
                    f"Total Duration: {total_duration:.2f} seconds",
                    f"Average Duration: {avg_duration:.3f} seconds",
                    f"Overall Success Rate: {success_rate:.1%}",
                    ""
                ])
        
        # Device comparison
        device_performance = self.get_device_performance_comparison()
        if device_performance:
            report_lines.append("Device Performance Comparison:")
            report_lines.append("-" * 40)
            
            for device_type, stats in device_performance.items():
                report_lines.extend([
                    f"{device_type.upper()}:",
                    f"  Operations: {stats['total_operations']}",
                    f"  Success Rate: {stats['success_rate']:.1%}",
                    f"  Avg Duration: {stats['avg_duration']:.3f}s",
                    f"  Ops/Second: {stats['operations_per_second']:.2f}",
                    ""
                ])
        
        # Performance issues
        issues = self.detect_performance_issues()
        if issues:
            report_lines.append("Performance Issues Detected:")
            report_lines.append("-" * 30)
            for issue in issues:
                report_lines.append(f"• {issue}")
            report_lines.append("")
        
        # Recommendations
        optimal_device = self.optimize_device_selection()
        report_lines.extend([
            "Recommendations:",
            f"• Optimal device type: {optimal_device.upper()}",
            ""
        ])
        
        return "\n".join(report_lines)
    
    def reset_metrics(self) -> None:
        """Clear all recorded metrics."""
        with self._lock:
            self._metrics.clear()
            self._operation_stats.clear()
        
        self.logger.info("Performance metrics reset")
    
    def enable_monitoring(self) -> None:
        """Enable performance monitoring."""
        self._monitoring_enabled = True
    
    def disable_monitoring(self) -> None:
        """Disable performance monitoring."""
        self._monitoring_enabled = False


class OperationTimer:
    """Context manager for timing operations."""
    
    def __init__(self, monitor: PerformanceMonitor, operation_name: str):
        """
        Initialize operation timer.
        
        Args:
            monitor: Performance monitor instance
            operation_name: Name of the operation being timed
        """
        self.monitor = monitor
        self.operation_name = operation_name
        self.start_time = None
        self.memory_before = 0.0
    
    def __enter__(self) -> 'OperationTimer':
        """Start timing the operation."""
        self.start_time = time.time()
        
        # Record initial memory usage
        try:
            memory_info = self.monitor.device_manager.get_memory_info()
            if 'allocated_gb' in memory_info:
                self.memory_before = memory_info['allocated_gb'] * 1024  # Convert to MB
            elif 'used_gb' in memory_info:
                self.memory_before = memory_info['used_gb'] * 1024
        except Exception:
            self.memory_before = 0.0
        
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop timing and record the metric."""
        if self.start_time is None:
            return
        
        duration = time.time() - self.start_time
        success = exc_type is None
        error_message = str(exc_val) if exc_val else None
        
        # Calculate memory usage change
        memory_used = 0.0
        try:
            memory_info = self.monitor.device_manager.get_memory_info()
            if 'allocated_gb' in memory_info:
                memory_after = memory_info['allocated_gb'] * 1024
                memory_used = max(0, memory_after - self.memory_before)
            elif 'used_gb' in memory_info:
                memory_after = memory_info['used_gb'] * 1024
                memory_used = max(0, memory_after - self.memory_before)
        except Exception:
            pass
        
        # Create and record metric
        metric = PerformanceMetric(
            operation=self.operation_name,
            duration=duration,
            device_type=self.monitor.device_manager.selected_device.device_type,
            memory_used=memory_used,
            success=success,
            error_message=error_message
        )
        
        self.monitor.record_metric(metric)