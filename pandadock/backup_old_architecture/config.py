"""
Configuration management module for PandaDock.
This module handles all configuration settings in a centralized manner.
"""

import os
import json
from pathlib import Path
from typing import Dict, Any, Optional, List


class DockingConfig:
    """Centralized configuration manager for PandaDock."""
    
    def __init__(self):
        """Initialize configuration with default values."""
        self.reset_to_defaults()
    
    def reset_to_defaults(self):
        """Reset all configuration to default values."""
        self.hardware = HardwareConfig()
        self.algorithm = AlgorithmConfig()
        self.scoring = ScoringConfig()
        self.output = OutputConfig()
        self.analysis = AnalysisConfig()
    
    def from_args(self, args):
        """Update configuration from command-line arguments."""
        self.hardware.from_args(args)
        self.algorithm.from_args(args)
        self.scoring.from_args(args)
        self.output.from_args(args)
        self.analysis.from_args(args)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'hardware': self.hardware.to_dict(),
            'algorithm': self.algorithm.to_dict(),
            'scoring': self.scoring.to_dict(),
            'output': self.output.to_dict(),
            'analysis': self.analysis.to_dict()
        }
    
    def save(self, filepath: str):
        """Save configuration to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    def load(self, filepath: str):
        """Load configuration from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        self.hardware.from_dict(data.get('hardware', {}))
        self.algorithm.from_dict(data.get('algorithm', {}))
        self.scoring.from_dict(data.get('scoring', {}))
        self.output.from_dict(data.get('output', {}))
        self.analysis.from_dict(data.get('analysis', {}))


class HardwareConfig:
    """Hardware acceleration configuration."""
    
    def __init__(self):
        self.use_gpu = False
        self.gpu_id = 0
        self.gpu_precision = 'float32'
        self.cpu_workers = None
        self.cpu_affinity = False
        self.workload_balance = 0.8
        self.auto_tune = False
    
    def from_args(self, args):
        """Update from command-line arguments."""
        self.use_gpu = getattr(args, 'use_gpu', self.use_gpu)
        self.gpu_id = getattr(args, 'gpu_id', self.gpu_id)
        self.gpu_precision = getattr(args, 'gpu_precision', self.gpu_precision)
        self.cpu_workers = getattr(args, 'cpu_workers', self.cpu_workers)
        self.cpu_affinity = getattr(args, 'cpu_affinity', self.cpu_affinity)
        self.workload_balance = getattr(args, 'workload_balance', self.workload_balance)
        self.auto_tune = getattr(args, 'auto_tune', self.auto_tune)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'use_gpu': self.use_gpu,
            'gpu_id': self.gpu_id,
            'gpu_precision': self.gpu_precision,
            'cpu_workers': self.cpu_workers,
            'cpu_affinity': self.cpu_affinity,
            'workload_balance': self.workload_balance,
            'auto_tune': self.auto_tune
        }
    
    def from_dict(self, data: Dict[str, Any]):
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)
    
    def detect_optimal_settings(self):
        """Automatically detect optimal hardware settings."""
        import multiprocessing
        
        # Detect CPU cores
        if self.cpu_workers is None:
            self.cpu_workers = multiprocessing.cpu_count()
        
        # Try to detect GPU
        if not self.use_gpu:
            try:
                import torch
                if torch.cuda.is_available():
                    self.use_gpu = True
                    print(f"GPU detected: {torch.cuda.get_device_name()}")
            except ImportError:
                pass
        
        return self


class AlgorithmConfig:
    """Docking algorithm configuration."""
    
    def __init__(self):
        self.algorithm_type = 'genetic'
        self.max_iterations = 100
        self.population_size = 100
        self.mutation_rate = 0.2
        self.crossover_rate = 0.8
        self.exhaustiveness = 1
        self.local_optimization = False
        
        # Monte Carlo settings
        self.mc_steps = 1000
        self.temperature = 300.0
        self.cooling_factor = 0.95
        
        # PandaDock settings
        self.high_temp = 1000.0
        self.target_temp = 300.0
        self.num_conformers = 10
        self.num_orientations = 10
        self.md_steps = 1000
        self.minimize_steps = 200
        self.use_grid = False
        
        # Grid settings
        self.grid_spacing = 0.375
        self.grid_radius = 10.0
        self.grid_center = None
    
    def from_args(self, args):
        """Update from command-line arguments."""
        self.algorithm_type = getattr(args, 'algorithm', self.algorithm_type)
        self.max_iterations = getattr(args, 'iterations', self.max_iterations)
        self.population_size = getattr(args, 'population_size', self.population_size)
        self.exhaustiveness = getattr(args, 'exhaustiveness', self.exhaustiveness)
        self.local_optimization = getattr(args, 'local_opt', self.local_optimization)
        
        # Monte Carlo
        self.mc_steps = getattr(args, 'mc_steps', self.mc_steps)
        self.temperature = getattr(args, 'temperature', self.temperature)
        self.cooling_factor = getattr(args, 'cooling_factor', self.cooling_factor)
        
        # PandaDock
        self.high_temp = getattr(args, 'high_temp', self.high_temp)
        self.target_temp = getattr(args, 'target_temp', self.target_temp)
        self.num_conformers = getattr(args, 'num_conformers', self.num_conformers)
        self.num_orientations = getattr(args, 'num_orientations', self.num_orientations)
        self.md_steps = getattr(args, 'md_steps', self.md_steps)
        self.minimize_steps = getattr(args, 'minimize_steps', self.minimize_steps)
        self.use_grid = getattr(args, 'use_grid', self.use_grid)
        
        # Grid settings
        self.grid_spacing = getattr(args, 'grid_spacing', self.grid_spacing)
        self.grid_radius = getattr(args, 'grid_radius', self.grid_radius)
        if hasattr(args, 'site') and args.site:
            self.grid_center = args.site
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'algorithm_type': self.algorithm_type,
            'max_iterations': self.max_iterations,
            'population_size': self.population_size,
            'mutation_rate': self.mutation_rate,
            'crossover_rate': self.crossover_rate,
            'exhaustiveness': self.exhaustiveness,
            'local_optimization': self.local_optimization,
            'mc_steps': self.mc_steps,
            'temperature': self.temperature,
            'cooling_factor': self.cooling_factor,
            'high_temp': self.high_temp,
            'target_temp': self.target_temp,
            'num_conformers': self.num_conformers,
            'num_orientations': self.num_orientations,
            'md_steps': self.md_steps,
            'minimize_steps': self.minimize_steps,
            'use_grid': self.use_grid,
            'grid_spacing': self.grid_spacing,
            'grid_radius': self.grid_radius,
            'grid_center': self.grid_center
        }
    
    def from_dict(self, data: Dict[str, Any]):
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)


class ScoringConfig:
    """Scoring function configuration."""
    
    def __init__(self):
        self.scoring_type = 'standard'
        self.enhanced_scoring = False
        self.physics_based = False
        self.mmff_minimization = False
        self.tethered_docking = False
        self.tether_weight = 10.0
        self.weights = None
        self.ph = 7.4
    
    def from_args(self, args):
        """Update from command-line arguments."""
        if getattr(args, 'physics_based', False):
            self.scoring_type = 'physics'
            self.physics_based = True
        elif getattr(args, 'enhanced_scoring', False):
            self.scoring_type = 'enhanced'
            self.enhanced_scoring = True
        else:
            self.scoring_type = 'standard'
        
        self.mmff_minimization = getattr(args, 'mmff_minimization', self.mmff_minimization)
        self.tethered_docking = getattr(args, 'tethered_docking', self.tethered_docking)
        self.tether_weight = getattr(args, 'tether_weight', self.tether_weight)
        self.ph = getattr(args, 'ph', self.ph)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'scoring_type': self.scoring_type,
            'enhanced_scoring': self.enhanced_scoring,
            'physics_based': self.physics_based,
            'mmff_minimization': self.mmff_minimization,
            'tethered_docking': self.tethered_docking,
            'tether_weight': self.tether_weight,
            'weights': self.weights,
            'ph': self.ph
        }
    
    def from_dict(self, data: Dict[str, Any]):
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)


class OutputConfig:
    """Output configuration."""
    
    def __init__(self):
        self.output_dir = 'docking_results'
        self.report_format = 'all'
        self.report_name = None
        self.detailed_energy = False
        self.skip_plots = False
        self.prepare_molecules = False
    
    def from_args(self, args):
        """Update from command-line arguments."""
        self.output_dir = getattr(args, 'output', self.output_dir)
        self.report_format = getattr(args, 'report_format', self.report_format)
        self.report_name = getattr(args, 'report_name', self.report_name)
        self.detailed_energy = getattr(args, 'detailed_energy', self.detailed_energy)
        self.skip_plots = getattr(args, 'skip_plots', self.skip_plots)
        self.prepare_molecules = getattr(args, 'prepare_molecules', self.prepare_molecules)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'output_dir': self.output_dir,
            'report_format': self.report_format,
            'report_name': self.report_name,
            'detailed_energy': self.detailed_energy,
            'skip_plots': self.skip_plots,
            'prepare_molecules': self.prepare_molecules
        }
    
    def from_dict(self, data: Dict[str, Any]):
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)


class AnalysisConfig:
    """Analysis configuration."""
    
    def __init__(self):
        self.cluster_poses = False
        self.clustering_method = 'hierarchical'
        self.rmsd_cutoff = 10.0
        self.analyze_interactions = False
        self.interaction_types = ['hbond', 'hydrophobic', 'ionic']
        self.classify_modes = False
        self.discover_modes = False
        self.n_modes = 5
        self.energy_decomposition = False
        self.per_residue_energy = False
        self.generate_analysis_report = False
        self.analysis_report_format = 'html'
        self.analysis_report_sections = ['summary', 'clusters', 'interactions', 'energetics']
    
    def from_args(self, args):
        """Update from command-line arguments."""
        self.cluster_poses = getattr(args, 'cluster_poses', self.cluster_poses)
        self.clustering_method = getattr(args, 'clustering_method', self.clustering_method)
        self.rmsd_cutoff = getattr(args, 'rmsd_cutoff', self.rmsd_cutoff)
        self.analyze_interactions = getattr(args, 'analyze_interactions', self.analyze_interactions)
        self.interaction_types = getattr(args, 'interaction_types', self.interaction_types)
        self.classify_modes = getattr(args, 'classify_modes', self.classify_modes)
        self.discover_modes = getattr(args, 'discover_modes', self.discover_modes)
        self.n_modes = getattr(args, 'n_modes', self.n_modes)
        self.energy_decomposition = getattr(args, 'energy_decomposition', self.energy_decomposition)
        self.per_residue_energy = getattr(args, 'per_residue_energy', self.per_residue_energy)
        self.generate_analysis_report = getattr(args, 'generate_analysis_report', self.generate_analysis_report)
        self.analysis_report_format = getattr(args, 'analysis_report_format', self.analysis_report_format)
        self.analysis_report_sections = getattr(args, 'analysis_report_sections', self.analysis_report_sections)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'cluster_poses': self.cluster_poses,
            'clustering_method': self.clustering_method,
            'rmsd_cutoff': self.rmsd_cutoff,
            'analyze_interactions': self.analyze_interactions,
            'interaction_types': self.interaction_types,
            'classify_modes': self.classify_modes,
            'discover_modes': self.discover_modes,
            'n_modes': self.n_modes,
            'energy_decomposition': self.energy_decomposition,
            'per_residue_energy': self.per_residue_energy,
            'generate_analysis_report': self.generate_analysis_report,
            'analysis_report_format': self.analysis_report_format,
            'analysis_report_sections': self.analysis_report_sections
        }
    
    def from_dict(self, data: Dict[str, Any]):
        for key, value in data.items():
            if hasattr(self, key):
                setattr(self, key, value)


def create_default_config() -> DockingConfig:
    """Create a default configuration instance."""
    return DockingConfig()


def load_config_from_file(filepath: str) -> DockingConfig:
    """Load configuration from file."""
    config = DockingConfig()
    if os.path.exists(filepath):
        config.load(filepath)
    return config


def save_config_to_file(config: DockingConfig, filepath: str):
    """Save configuration to file."""
    config.save(filepath)
