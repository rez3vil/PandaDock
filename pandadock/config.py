# -*- coding: utf-8 -*-
"""
Configuration management for PandaDock
"""

import os
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class DockingConfig:
    """Configuration for docking parameters"""
    mode: str = "balanced"  # precise, balanced, fast
    num_poses: int = 10
    exhaustiveness: int = 8
    energy_range: float = 3.0
    seed: int = 42
    
    # Flexible docking
    flexible_residues: List[str] = field(default_factory=list)
    side_chain_flexibility: bool = False
    
    # Physics engine specific
    use_clash_detection: bool = True
    minimization_steps: int = 100
    vdw_scale: float = 1.0
    electrostatic_scale: float = 1.0
    
    # ML engine specific
    model_path: Optional[str] = None
    confidence_threshold: float = 0.5
    
    # GA engine specific
    population_size: int = 150
    generations: int = 27000
    mutation_rate: float = 0.02
    crossover_rate: float = 0.8
    elitism_rate: float = 0.1


@dataclass
class IOConfig:
    """Configuration for input/output handling"""
    protein_file: str = ""
    ligand_file: str = ""
    output_dir: str = "output"
    
    # Grid box parameters
    center_x: float = 0.0
    center_y: float = 0.0
    center_z: float = 0.0
    size_x: float = 22.5
    size_y: float = 22.5
    size_z: float = 22.5
    
    # Output formats
    save_poses: bool = True
    save_complex: bool = True
    save_report: bool = True
    report_format: str = "html"  # html, json, csv


@dataclass
class ScoringConfig:
    """Configuration for scoring functions"""
    scoring_function: str = "pandacore"  # pandacore, pandaml, pandaphysics
    include_solvation: bool = True
    include_entropy: bool = True
    use_ml_rescoring: bool = False
    
    # Weights for scoring terms
    vdw_weight: float = 1.0
    electrostatic_weight: float = 1.0
    hbond_weight: float = 1.0
    hydrophobic_weight: float = 1.0
    solvation_weight: float = 1.0
    entropy_weight: float = 1.0


@dataclass
class PandaDockConfig:
    """Main configuration class for PandaDock"""
    docking: DockingConfig = field(default_factory=DockingConfig)
    io: IOConfig = field(default_factory=IOConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    
    # Global settings
    verbose: bool = False
    debug: bool = False
    n_jobs: int = 1
    gpu_enabled: bool = False
    
    @classmethod
    def from_file(cls, config_file: str) -> 'PandaDockConfig':
        """Load configuration from JSON file"""
        import json
        
        with open(config_file, 'r', encoding='utf-8') as f:
            config_dict = json.load(f)
        
        return cls.from_dict(config_dict)
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PandaDockConfig':
        """Create configuration from dictionary"""
        config = cls()
        
        if 'docking' in config_dict:
            for key, value in config_dict['docking'].items():
                if hasattr(config.docking, key):
                    setattr(config.docking, key, value)
        
        if 'io' in config_dict:
            for key, value in config_dict['io'].items():
                if hasattr(config.io, key):
                    setattr(config.io, key, value)
        
        if 'scoring' in config_dict:
            for key, value in config_dict['scoring'].items():
                if hasattr(config.scoring, key):
                    setattr(config.scoring, key, value)
        
        # Global settings
        for key in ['verbose', 'debug', 'n_jobs', 'gpu_enabled']:
            if key in config_dict:
                setattr(config, key, config_dict[key])
        
        return config
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary"""
        return {
            'docking': self.docking.__dict__,
            'io': self.io.__dict__,
            'scoring': self.scoring.__dict__,
            'verbose': self.verbose,
            'debug': self.debug,
            'n_jobs': self.n_jobs,
            'gpu_enabled': self.gpu_enabled
        }
    
    def save(self, config_file: str):
        """Save configuration to JSON file"""
        import json
        
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    def validate(self) -> bool:
        """Validate configuration parameters"""
        errors = []
        
        # Validate docking mode
        if self.docking.mode not in ['precise', 'balanced', 'fast']:
            errors.append(f"Invalid docking mode: {self.docking.mode}")
        
        # Validate file paths
        if self.io.protein_file and not os.path.exists(self.io.protein_file):
            errors.append(f"Protein file not found: {self.io.protein_file}")
        
        if self.io.ligand_file and not os.path.exists(self.io.ligand_file):
            errors.append(f"Ligand file not found: {self.io.ligand_file}")
        
        # Validate scoring function
        if self.scoring.scoring_function not in ['pandacore', 'pandaml', 'pandaphysics']:
            errors.append(f"Invalid scoring function: {self.scoring.scoring_function}")
        
        # Validate numeric parameters
        if self.docking.num_poses <= 0:
            errors.append("Number of poses must be positive")
        
        if self.docking.exhaustiveness <= 0:
            errors.append("Exhaustiveness must be positive")
        
        if errors:
            for error in errors:
                print(f"Configuration error: {error}")
            return False
        
        return True


# Default configuration
DEFAULT_CONFIG = PandaDockConfig()