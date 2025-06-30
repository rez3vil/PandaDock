"""
Command-line argument parsing for PandaDock.

This module handles all CLI argument parsing with clear organization
and helpful error messages.
"""

import argparse
from pathlib import Path
from typing import Optional, List, Dict, Any


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser for PandaDock.
    
    Returns:
        Configured ArgumentParser instance
    """
    parser = argparse.ArgumentParser(
        description="🐼 PandaDock: Physics-based Molecular Docking 🚀",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic docking
  pandadock -p protein.pdb -l ligand.sdf -o results
  
  # Enhanced docking with GPU
  pandadock -p protein.pdb -l ligand.sdf -o results --enhanced --use-gpu
  
  # Physics-based docking with flexible residues
  pandadock -p protein.pdb -l ligand.sdf -o results --physics-based --auto-flex
  
  # Ensemble docking with reference validation
  pandadock -p protein.pdb -l ligand.sdf -o results --exhaustiveness 5 --reference ref.sdf

For more information, visit: https://github.com/pritampanda15/PandaDock
        """
    )
    
    # Add command subparsers
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Main docking command
    dock_parser = subparsers.add_parser('dock', help='Run molecular docking')
    _add_docking_arguments(dock_parser)
    
    # Analysis command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze docking results')
    _add_analysis_arguments(analyze_parser)
    
    # Preparation command
    prep_parser = subparsers.add_parser('prepare', help='Prepare molecules for docking')
    _add_preparation_arguments(prep_parser)
    
    # If no subcommand is provided, default to docking
    _add_docking_arguments(parser)
    
    return parser


def _add_docking_arguments(parser: argparse.ArgumentParser) -> None:
    """Add docking-specific arguments."""
    
    # Required arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument(
        '-p', '--protein', 
        required=True,
        type=str,
        help='Path to protein PDB file'
    )
    required.add_argument(
        '-l', '--ligand',
        required=True, 
        type=str,
        help='Path to ligand MOL/SDF file'
    )
    required.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        help='Output directory for results'
    )
    
    # Basic options
    basic = parser.add_argument_group('Basic Options')
    basic.add_argument(
        '-a', '--algorithm',
        choices=['genetic', 'random', 'monte-carlo', 'pandadock'],
        default='genetic',
        help='Docking algorithm to use (default: genetic)'
    )
    basic.add_argument(
        '-i', '--iterations',
        type=int,
        default=100,
        help='Number of iterations/generations (default: 100)'
    )
    basic.add_argument(
        '--exhaustiveness',
        type=int,
        default=1,
        help='Number of independent docking runs (default: 1)'
    )
    
    # Binding site options
    binding_site = parser.add_argument_group('Binding Site Options')
    binding_site.add_argument(
        '-s', '--site',
        nargs=3,
        type=float,
        metavar=('X', 'Y', 'Z'),
        help='Active site center coordinates'
    )
    binding_site.add_argument(
        '-r', '--radius',
        type=float,
        default=10.0,
        help='Active site radius in Angstroms (default: 10.0)'
    )
    binding_site.add_argument(
        '--detect-pockets',
        action='store_true',
        help='Automatically detect binding pockets'
    )
    
    # Hardware options
    hardware = parser.add_argument_group('Hardware Options')
    hardware.add_argument(
        '--use-gpu',
        action='store_true',
        help='Use GPU acceleration if available'
    )
    hardware.add_argument(
        '--gpu-id',
        type=int,
        default=0,
        help='GPU device ID to use (default: 0)'
    )
    hardware.add_argument(
        '--cpu-workers',
        type=int,
        default=None,
        help='Number of CPU workers for parallel processing'
    )
    
    # Scoring options
    scoring = parser.add_argument_group('Scoring Options')
    scoring.add_argument(
        '--enhanced-scoring',
        action='store_true',
        help='Use enhanced scoring function with electrostatics'
    )
    scoring.add_argument(
        '--physics-based',
        action='store_true',
        help='Use full physics-based scoring (very slow but most accurate)'
    )
    
    # Molecule preparation
    prep = parser.add_argument_group('Molecule Preparation')
    prep.add_argument(
        '--prepare-molecules',
        action='store_true',
        help='Prepare protein and ligand before docking'
    )
    prep.add_argument(
        '--ph',
        type=float,
        default=7.4,
        help='pH for protein preparation (default: 7.4)'
    )
    
    # Flexibility options
    flexibility = parser.add_argument_group('Protein Flexibility')
    flexibility.add_argument(
        '--flex-residues',
        nargs='+',
        help='Specify flexible residue IDs (e.g., A_42 B_57)'
    )
    flexibility.add_argument(
        '--auto-flex',
        action='store_true',
        help='Automatically detect flexible residues'
    )
    flexibility.add_argument(
        '--max-flex-residues',
        type=int,
        default=5,
        help='Maximum number of flexible residues (default: 5)'
    )
    
    # Advanced options
    advanced = parser.add_argument_group('Advanced Options')
    advanced.add_argument(
        '--local-opt',
        action='store_true',
        help='Enable local optimization on top poses'
    )
    advanced.add_argument(
        '--population-size',
        type=int,
        default=100,
        help='Population size for genetic algorithm (default: 100)'
    )
    advanced.add_argument(
        '--mc-steps',
        type=int,
        default=1000,
        help='Number of Monte Carlo steps (default: 1000)'
    )
    advanced.add_argument(
        '--temperature',
        type=float,
        default=300.0,
        help='Temperature for Monte Carlo simulation in Kelvin (default: 300K)'
    )
    
    # Reference docking
    reference = parser.add_argument_group('Reference Docking')
    reference.add_argument(
        '--reference',
        type=str,
        help='Reference ligand structure for validation'
    )
    reference.add_argument(
        '--tethered-docking',
        action='store_true',
        help='Use tethered scoring with reference structure'
    )
    reference.add_argument(
        '--exact-alignment',
        action='store_true',
        help='Align docked pose exactly to reference structure'
    )
    
    # Output options
    output = parser.add_argument_group('Output Options')
    output.add_argument(
        '--report-format',
        choices=['text', 'csv', 'json', 'html', 'all'],
        default='all',
        help='Report format (default: all)'
    )
    output.add_argument(
        '--detailed-energy',
        action='store_true',
        help='Include detailed energy analysis'
    )
    output.add_argument(
        '--skip-plots',
        action='store_true',
        help='Skip generating plots'
    )
    
    # Quick mode options
    modes = parser.add_argument_group('Quick Modes')
    modes.add_argument(
        '--fast-mode',
        action='store_true',
        help='Run with minimal enhancements for quick results'
    )
    modes.add_argument(
        '--enhanced',
        action='store_true',
        help='Use enhanced algorithms for better accuracy'
    )
    modes.add_argument(
        '--auto-algorithm',
        action='store_true',
        help='Automatically select best algorithm for your system'
    )


def _add_analysis_arguments(parser: argparse.ArgumentParser) -> None:
    """Add analysis-specific arguments."""
    
    parser.add_argument(
        'results_dir',
        type=str,
        help='Directory containing docking results to analyze'
    )
    
    # Analysis options
    analysis = parser.add_argument_group('Analysis Options')
    analysis.add_argument(
        '--cluster-poses',
        action='store_true',
        help='Perform clustering of docking poses'
    )
    analysis.add_argument(
        '--rmsd-cutoff',
        type=float,
        default=2.0,
        help='RMSD cutoff for pose clustering (default: 2.0)'
    )
    analysis.add_argument(
        '--analyze-interactions',
        action='store_true',
        help='Analyze protein-ligand interactions'
    )
    analysis.add_argument(
        '--energy-decomposition',
        action='store_true',
        help='Perform energy decomposition analysis'
    )
    
    # Output options
    output = parser.add_argument_group('Output Options')
    output.add_argument(
        '--output-format',
        choices=['html', 'pdf', 'txt'],
        default='html',
        help='Analysis report format (default: html)'
    )


def _add_preparation_arguments(parser: argparse.ArgumentParser) -> None:
    """Add preparation-specific arguments."""
    
    parser.add_argument(
        'input_file',
        type=str,
        help='Input molecule file to prepare'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help='Output file for prepared molecule'
    )
    
    # Preparation options
    prep = parser.add_argument_group('Preparation Options')
    prep.add_argument(
        '--molecule-type',
        choices=['protein', 'ligand'],
        required=True,
        help='Type of molecule to prepare'
    )
    prep.add_argument(
        '--add-hydrogens',
        action='store_true',
        help='Add hydrogen atoms'
    )
    prep.add_argument(
        '--remove-water',
        action='store_true',
        help='Remove water molecules'
    )
    prep.add_argument(
        '--ph',
        type=float,
        default=7.4,
        help='pH for protonation states (default: 7.4)'
    )


def validate_arguments(args: argparse.Namespace) -> None:
    """
    Validate command-line arguments.
    
    Args:
        args: Parsed arguments
        
    Raises:
        argparse.ArgumentError: If arguments are invalid
    """
    # Check file existence
    if hasattr(args, 'protein') and args.protein:
        if not Path(args.protein).exists():
            raise argparse.ArgumentError(None, f"Protein file not found: {args.protein}")
    
    if hasattr(args, 'ligand') and args.ligand:
        if not Path(args.ligand).exists():
            raise argparse.ArgumentError(None, f"Ligand file not found: {args.ligand}")
    
    if hasattr(args, 'reference') and args.reference:
        if not Path(args.reference).exists():
            raise argparse.ArgumentError(None, f"Reference file not found: {args.reference}")
    
    # Validate numerical ranges
    if hasattr(args, 'iterations') and args.iterations <= 0:
        raise argparse.ArgumentError(None, "Iterations must be positive")
    
    if hasattr(args, 'exhaustiveness') and args.exhaustiveness <= 0:
        raise argparse.ArgumentError(None, "Exhaustiveness must be positive")
    
    if hasattr(args, 'radius') and args.radius <= 0:
        raise argparse.ArgumentError(None, "Radius must be positive")
    
    if hasattr(args, 'temperature') and args.temperature <= 0:
        raise argparse.ArgumentError(None, "Temperature must be positive")
    
    # Validate combinations
    if hasattr(args, 'exact_alignment') and args.exact_alignment and not args.reference:
        raise argparse.ArgumentError(None, "--exact-alignment requires --reference")
    
    if hasattr(args, 'tethered_docking') and args.tethered_docking and not args.reference:
        raise argparse.ArgumentError(None, "--tethered-docking requires --reference")


def get_config_from_args(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Convert command-line arguments to configuration dictionary.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Configuration dictionary for the docking engine
    """
    config = {}
    
    # Basic parameters
    if hasattr(args, 'algorithm'):
        config['algorithm'] = args.algorithm
    if hasattr(args, 'iterations'):
        config['iterations'] = args.iterations
    if hasattr(args, 'exhaustiveness'):
        config['exhaustiveness'] = args.exhaustiveness
    
    # Hardware settings
    config['use_gpu'] = getattr(args, 'use_gpu', False)
    config['gpu_id'] = getattr(args, 'gpu_id', 0)
    config['cpu_workers'] = getattr(args, 'cpu_workers', None)
    
    # Binding site
    if hasattr(args, 'site') and args.site:
        config['binding_site_center'] = args.site
    if hasattr(args, 'radius'):
        config['binding_site_radius'] = args.radius
    config['detect_pockets'] = getattr(args, 'detect_pockets', False)
    
    # Scoring
    config['enhanced_scoring'] = getattr(args, 'enhanced_scoring', False)
    config['physics_based'] = getattr(args, 'physics_based', False)
    
    # Preparation
    config['prepare_molecules'] = getattr(args, 'prepare_molecules', False)
    config['ph'] = getattr(args, 'ph', 7.4)
    
    # Flexibility
    config['flexible_residues'] = getattr(args, 'flex_residues', None)
    config['auto_flex'] = getattr(args, 'auto_flex', False)
    config['max_flex_residues'] = getattr(args, 'max_flex_residues', 5)
    
    # Advanced options
    config['local_opt'] = getattr(args, 'local_opt', False)
    config['population_size'] = getattr(args, 'population_size', 100)
    config['mc_steps'] = getattr(args, 'mc_steps', 1000)
    config['temperature'] = getattr(args, 'temperature', 300.0)
    
    # Reference docking
    config['reference_ligand'] = getattr(args, 'reference', None)
    config['tethered_docking'] = getattr(args, 'tethered_docking', False)
    config['exact_alignment'] = getattr(args, 'exact_alignment', False)
    
    # Output
    config['report_format'] = getattr(args, 'report_format', 'all')
    config['detailed_energy'] = getattr(args, 'detailed_energy', False)
    config['skip_plots'] = getattr(args, 'skip_plots', False)
    
    # Mode flags
    if getattr(args, 'fast_mode', False):
        config.update({
            'enhanced_scoring': False,
            'physics_based': False,
            'local_opt': False,
            'population_size': 50,
            'prepare_molecules': False
        })
    
    if getattr(args, 'enhanced', False):
        config.update({
            'enhanced_scoring': True,
            'local_opt': True,
            'population_size': max(config.get('population_size', 100), 100)
        })
    
    return config