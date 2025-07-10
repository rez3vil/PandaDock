#!/usr/bin/env python3
"""
PandaDock: Main entry point for molecular docking
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional, List
import json
import logging

# Import from parent directory modules
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import PandaDockConfig
from docking.physics_engine import PhysicsEngine
from docking.ml_engine import MLEngine
from docking.ga_engine import GAEngine
from reports.html_report import HTMLReportGenerator


def setup_logging(verbose: bool, debug: bool):
    """Setup logging configuration"""
    level = logging.DEBUG if debug else (logging.INFO if verbose else logging.WARNING)
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('pandadock.log')
        ]
    )


def parse_flexible_residues(residues_str: str) -> List[str]:
    """Parse flexible residues string into list"""
    if not residues_str:
        return []
    return [res.strip() for res in residues_str.split(',')]


def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='PandaDock: Modular, Multi-Strategy, High-Performance Docking Software',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic docking with balanced mode
  pandadock --protein receptor.pdb --ligand ligand.sdf --mode balanced

  # Physics-based docking with flexible residues
  pandadock --protein receptor.pdb --ligand ligand.sdf --mode precise \\
           --flexible-residues "HIS57,SER195" --out report.html

  # Fast virtual screening
  pandadock --protein receptor.pdb --screen ligands.smi --mode fast \\
           --num-poses 5 --exhaustiveness 16

  # Use custom configuration
  pandadock --config config.json
        """
    )
    
    # Input files
    parser.add_argument('--protein', '-p', type=str, required=True,
                       help='Protein receptor file (PDB format)')
    parser.add_argument('--ligand', '-l', type=str,
                       help='Ligand file (SDF, MOL2, or SMILES format)')
    parser.add_argument('--screen', '-s', type=str,
                       help='File containing multiple ligands for virtual screening')
    
    # Docking mode
    parser.add_argument('--mode', '-m', type=str, default='balanced',
                       choices=['precise', 'balanced', 'fast'],
                       help='Docking mode: precise (physics-based), balanced (ML-based), fast (GA-based)')
    
    # Grid box parameters
    parser.add_argument('--center', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       help='Grid box center coordinates')
    parser.add_argument('--size', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                       default=[22.5, 22.5, 22.5],
                       help='Grid box size (default: 22.5 22.5 22.5)')
    
    # Docking parameters
    parser.add_argument('--num-poses', '-n', type=int, default=10,
                       help='Number of poses to generate (default: 10)')
    parser.add_argument('--exhaustiveness', '-e', type=int, default=8,
                       help='Exhaustiveness of search (default: 8)')
    parser.add_argument('--energy-range', type=float, default=3.0,
                       help='Energy range for pose selection (default: 3.0)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    
    # Flexible docking
    parser.add_argument('--flexible-residues', type=str,
                       help='Comma-separated list of flexible residues (e.g., HIS57,SER195)')
    parser.add_argument('--side-chain-flexibility', action='store_true',
                       help='Enable side chain flexibility')
    
    # Scoring
    parser.add_argument('--scoring', type=str, default='vina',
                       choices=['vina', 'glide', 'ml'],
                       help='Scoring function (default: vina)')
    parser.add_argument('--ml-rescoring', action='store_true',
                       help='Use ML-based rescoring')
    
    # Output
    parser.add_argument('--out', '-o', type=str, default='output',
                       help='Output directory or file (default: output)')
    parser.add_argument('--save-poses', action='store_true', default=True,
                       help='Save docked poses')
    parser.add_argument('--save-complex', action='store_true', default=True,
                       help='Save protein-ligand complexes')
    parser.add_argument('--report-format', type=str, default='html',
                       choices=['html', 'json', 'csv'],
                       help='Report format (default: html)')
    
    # Configuration
    parser.add_argument('--config', '-c', type=str,
                       help='Configuration file (JSON format)')
    parser.add_argument('--save-config', type=str,
                       help='Save current configuration to file')
    
    # Performance
    parser.add_argument('--n-jobs', '-j', type=int, default=1,
                       help='Number of parallel jobs (default: 1)')
    parser.add_argument('--gpu', action='store_true',
                       help='Use GPU acceleration if available')
    
    # Logging
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    parser.add_argument('--debug', action='store_true',
                       help='Debug output')
    
    return parser


def create_config_from_args(args: argparse.Namespace) -> PandaDockConfig:
    """Create configuration object from command line arguments"""
    if args.config:
        config = PandaDockConfig.from_file(args.config)
    else:
        config = PandaDockConfig()
    
    # Override with command line arguments
    config.docking.mode = args.mode
    config.docking.num_poses = args.num_poses
    config.docking.exhaustiveness = args.exhaustiveness
    config.docking.energy_range = args.energy_range
    config.docking.seed = args.seed
    
    if args.flexible_residues:
        config.docking.flexible_residues = parse_flexible_residues(args.flexible_residues)
    config.docking.side_chain_flexibility = args.side_chain_flexibility
    
    config.io.protein_file = args.protein
    config.io.ligand_file = args.ligand or ""
    config.io.output_dir = args.out
    config.io.save_poses = args.save_poses
    config.io.save_complex = args.save_complex
    config.io.report_format = args.report_format
    
    if args.center:
        config.io.center_x, config.io.center_y, config.io.center_z = args.center
    if args.size:
        config.io.size_x, config.io.size_y, config.io.size_z = args.size
    
    config.scoring.scoring_function = args.scoring
    config.scoring.use_ml_rescoring = args.ml_rescoring
    
    config.verbose = args.verbose
    config.debug = args.debug
    config.n_jobs = args.n_jobs
    config.gpu_enabled = args.gpu
    
    return config


def get_docking_engine(config: PandaDockConfig):
    """Get appropriate docking engine based on configuration"""
    if config.docking.mode == 'precise':
        return PhysicsEngine(config)
    elif config.docking.mode == 'balanced':
        return MLEngine(config)
    elif config.docking.mode == 'fast':
        return GAEngine(config)
    else:
        raise ValueError(f"Unknown docking mode: {config.docking.mode}")


def run_docking(config: PandaDockConfig, ligand_file: str = None):
    """Run docking with given configuration"""
    logger = logging.getLogger(__name__)
    
    # Create output directory
    os.makedirs(config.io.output_dir, exist_ok=True)
    
    # Initialize docking engine
    logger.info(f"Initializing {config.docking.mode} docking engine")
    engine = get_docking_engine(config)
    
    # Use specified ligand file or default from config
    ligand_path = ligand_file or config.io.ligand_file
    
    try:
        # Run docking
        logger.info(f"Starting docking: {config.io.protein_file} + {ligand_path}")
        results = engine.dock(config.io.protein_file, ligand_path)
        
        # Generate report
        logger.info("Generating report")
        report_generator = HTMLReportGenerator(config)
        report_generator.generate_report(results)
        
        logger.info(f"Docking completed successfully. Results saved to {config.io.output_dir}")
        return results
        
    except Exception as e:
        logger.error(f"Docking failed: {str(e)}")
        raise


def run_virtual_screening(config: PandaDockConfig, ligand_library: str):
    """Run virtual screening with multiple ligands"""
    logger = logging.getLogger(__name__)
    
    # Read ligand library
    with open(ligand_library, 'r') as f:
        ligands = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Starting virtual screening with {len(ligands)} ligands")
    
    all_results = []
    for i, ligand in enumerate(ligands):
        logger.info(f"Processing ligand {i+1}/{len(ligands)}: {ligand}")
        
        # Create temporary ligand file or use SMILES directly
        try:
            results = run_docking(config, ligand)
            all_results.extend(results)
        except Exception as e:
            logger.warning(f"Failed to dock ligand {ligand}: {str(e)}")
            continue
    
    logger.info(f"Virtual screening completed. Processed {len(all_results)} poses")
    return all_results


def main():
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose, args.debug)
    logger = logging.getLogger(__name__)
    
    try:
        # Create configuration
        config = create_config_from_args(args)
        
        # Save configuration if requested
        if args.save_config:
            config.save(args.save_config)
            logger.info(f"Configuration saved to {args.save_config}")
        
        # Validate configuration
        if not config.validate():
            logger.error("Configuration validation failed")
            sys.exit(1)
        
        # Run docking or virtual screening
        if args.screen:
            results = run_virtual_screening(config, args.screen)
        else:
            if not args.ligand:
                logger.error("Either --ligand or --screen must be specified")
                sys.exit(1)
            results = run_docking(config)
        
        logger.info("PandaDock completed successfully")
        
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()