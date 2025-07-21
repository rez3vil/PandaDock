#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import pandas as pd

# Add PandaDock to path if running as script
if __name__ == '__main__':
    project_root = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, project_root)

from pandadock.config import PandaDockConfig
from pandadock.docking.physics_engine import PhysicsEngine
from pandadock.docking.ml_engine import MLEngine
from pandadock.docking.ga_engine import GAEngine
from pandadock.reports.html_report import HTMLReportGenerator
from pandadock.version import __version__, check_and_notify_updates


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
  # Basic docking with PandaCore scoring
  pandadock --protein receptor.pdb --ligand ligand.sdf --mode balanced --scoring pandacore

  # Physics-based docking with comprehensive plots
  pandadock --protein receptor.pdb --ligand ligand.sdf --mode precise \\
           --flexible-residues "HIS57,SER195" --scoring pandaphysics \\
           --plots --interaction-maps --out results

  # Docking with publication-ready outputs
  pandadock --protein beta-2_alpha-1.pdb --ligand propofol.pdb --mode balanced \\
           --scoring pandacore --flexible-residues "ASN265" \\
           --center -15.7 -17.7 8.18 --size 40 40 40 \\
           --all-outputs --protein-name "GABA_A Receptor" --ligand-name "Propofol"

  # Professional interaction analysis with PandaMap
  pandadock --protein receptor.pdb --ligand compound.sdf --mode precise \\
           --pandamap --pandamap-3d --plots --out results

  # Fast virtual screening with PandaML scoring
  pandadock --protein receptor.pdb --screen ligands.smi --mode fast \\
           --num-poses 5 --scoring pandaml --exhaustiveness 16 --master-plot

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
                       help='Docking mode: precise (PandaPhysics), balanced (PandaML), fast (PandaCore)')
    
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
    parser.add_argument('--scoring', type=str, default='pandacore',
                       choices=['pandacore', 'pandaml', 'pandaphysics'],
                       help='PandaDock scoring function (default: pandacore)')
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
    
    # Plotting and visualization
    parser.add_argument('--plots', action='store_true',
                       help='Generate comprehensive plots (binding metrics, IC50/EC50, master publication plot)')
    parser.add_argument('--interaction-maps', action='store_true',
                       help='Generate 2D interaction maps for top poses')
    parser.add_argument('--master-plot', action='store_true',
                       help='Generate publication-ready master plot')
    parser.add_argument('--txt-report', action='store_true',
                       help='Generate detailed TXT report with algorithm and command info')
    parser.add_argument('--pandamap', action='store_true',
                       help='Use PandaMap for professional protein-ligand interaction visualization (2D + 3D)')
    parser.add_argument('--pandamap-3d', action='store_true',
                       help='Generate interactive 3D visualizations using PandaMap')
    parser.add_argument('--all-outputs', action='store_true',
                       help='Generate all output formats (HTML, CSV, JSON, plots, interaction maps, TXT)')
    parser.add_argument('--protein-name', type=str,
                       help='Protein name for plot titles and reports')
    parser.add_argument('--ligand-name', type=str,
                       help='Ligand name for plot titles and reports')
    
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
    
    # Version
    parser.add_argument('--version', '-v', action='version', 
                       version=f'PandaDock {__version__}',
                       help='Show version and exit')
    
    # Logging
    parser.add_argument('--verbose', action='store_true',
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
        config.center = args.center  # Also set for base engine compatibility
    if args.size:
        config.io.size_x, config.io.size_y, config.io.size_z = args.size
    
    config.scoring.scoring_function = args.scoring
    config.scoring.use_ml_rescoring = args.ml_rescoring
    
    # Plotting and visualization settings
    config.enable_plots = args.plots or args.all_outputs
    config.enable_interaction_maps = args.interaction_maps or args.all_outputs
    config.enable_master_plot = args.master_plot or args.all_outputs
    config.enable_txt_report = args.txt_report or args.all_outputs
    config.enable_pandamap = args.pandamap or args.all_outputs
    config.enable_pandamap_3d = args.pandamap_3d or args.pandamap or args.all_outputs
    config.comprehensive_output = args.all_outputs
    config.protein_name = args.protein_name
    config.ligand_name = args.ligand_name
    
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
        print(f"\n🐼 PandaDock {config.docking.mode.upper()} Mode")
        print(f"📁 Protein: {config.io.protein_file}")
        print(f"💊 Ligand: {ligand_path}")
        print(f"🎯 Target poses: {config.docking.num_poses}")
        print("─" * 60)
        results = engine.dock(config.io.protein_file, ligand_path)
        print("─" * 60)
        
        # Save poses if requested
        if config.io.save_poses:
            poses_dir = os.path.join(config.io.output_dir, "poses")
            print(f"💾 Saving {len(results)} poses to {poses_dir}...")
            logger.info(f"Saving poses to {poses_dir}")
            engine.save_poses(results, poses_dir)
        
        # Generate comprehensive report with plots if requested
        report_generator = HTMLReportGenerator(config)
        
        # Check if comprehensive plotting is enabled
        if (hasattr(config, 'enable_plots') and config.enable_plots) or \
           (hasattr(config, 'comprehensive_output') and config.comprehensive_output):
            
            print("🎨 Generating comprehensive plots and reports...")
            logger.info("Generating comprehensive plots and reports")
            
            # Prepare names for plots
            protein_name = getattr(config, 'protein_name', None) or Path(config.io.protein_file).stem
            ligand_name = getattr(config, 'ligand_name', None) or Path(ligand_path).stem
            
            # Prepare algorithm and command information
            algorithm_info = {
                'algorithm': 'PandaDock',
                'version': '1.0.0',
                'scoring_function': config.scoring.scoring_function.upper(),
                'engine': config.docking.mode.title(),
                'mode': config.docking.mode
            }
            
            # Build command from config
            command_parts = ['python', '-m', 'pandadock']
            command_parts.extend(['--protein', config.io.protein_file])
            command_parts.extend(['--ligand', ligand_path])
            command_parts.extend(['--mode', config.docking.mode])
            command_parts.extend(['--scoring', config.scoring.scoring_function])
            
            if config.docking.flexible_residues:
                command_parts.extend(['--flexible-residues', ','.join(config.docking.flexible_residues)])
            if hasattr(config.io, 'center_x'):
                command_parts.extend(['--center', str(config.io.center_x), str(config.io.center_y), str(config.io.center_z)])
            if hasattr(config.io, 'size_x'):
                command_parts.extend(['--size', str(config.io.size_x), str(config.io.size_y), str(config.io.size_z)])
            command_parts.extend(['--out', config.io.output_dir])
            
            command_info = {
                'command': ' '.join(command_parts),
                'protein': config.io.protein_file,
                'ligand': ligand_path,
                'center': f"{getattr(config.io, 'center_x', 'auto')} {getattr(config.io, 'center_y', 'auto')} {getattr(config.io, 'center_z', 'auto')}",
                'size': f"{getattr(config.io, 'size_x', 22.5)} {getattr(config.io, 'size_y', 22.5)} {getattr(config.io, 'size_z', 22.5)}",
                'exhaustiveness': config.docking.exhaustiveness
            }
            
            # Generate comprehensive report
            generated_files = report_generator.generate_comprehensive_report(
                results=results,
                output_dir=config.io.output_dir,
                protein_name=protein_name,
                ligand_name=ligand_name,
                algorithm_info=algorithm_info,
                command_info=command_info
            )
            
            # Generate PandaMap visualizations if requested
            if (hasattr(config, 'enable_pandamap') and config.enable_pandamap):
                print("🗺️  Generating PandaMap professional interaction visualizations...")
                logger.info("Generating PandaMap visualizations")
                
                try:
                    from pandadock.reports.pandamap_integration import create_pandamap_visualizations
                    
                    # Create poses DataFrame
                    poses_data = []
                    for i, pose in enumerate(results):
                        binding_affinity = pose.get_binding_affinity()
                        ic50_um = pose.get_ic50(units='uM')
                        
                        poses_data.append({
                            'Rank': i + 1,
                            'Pose_ID': pose.pose_id,
                            'Binding_Affinity': binding_affinity,
                            'IC50_uM': f"{ic50_um:.2e}" if ic50_um != float('inf') else ''
                        })
                    
                    poses_df = pd.DataFrame(poses_data)
                    
                    # Generate PandaMap visualizations
                    pandamap_files = create_pandamap_visualizations(
                        poses_df=poses_df,
                        poses_dir=os.path.join(config.io.output_dir, "poses"),
                        output_dir=config.io.output_dir,
                        top_n=3,
                        generate_3d=getattr(config, 'enable_pandamap_3d', False)
                    )
                    
                    # Add to generated files
                    if pandamap_files:
                        generated_files.update(pandamap_files)
                        
                        map_count = len(pandamap_files.get('2d_maps', []))
                        viz_count = len(pandamap_files.get('3d_maps', []))
                        print(f"   🗺️  Generated {map_count} PandaMap 2D interaction maps")
                        if viz_count > 0:
                            print(f"   🌐 Generated {viz_count} interactive 3D visualizations")
                    
                except ImportError:
                    print("   ⚠️  PandaMap not available - install with: pip install pandamap")
                    logger.warning("PandaMap not available for visualization")
                except Exception as e:
                    print(f"   ⚠️  PandaMap visualization failed: {e}")
                    logger.warning(f"PandaMap visualization failed: {e}")
            
            # Report generated files
            plot_count = sum(1 for k, v in generated_files.items() if 'plot' in k.lower() or 'map' in k.lower())
            print(f"✨ Generated {len(generated_files)} output files:")
            print(f"   📊 Plots: {plot_count}")
            print(f"   📄 Reports: {len(generated_files) - plot_count}")
            
            if 'master_publication' in generated_files:
                print(f"   🏆 Master publication plot: {Path(generated_files['master_publication']).name}")
            
        else:
            # Standard report generation
            if config.io.report_format == 'json':
                json_output_path = os.path.join(config.io.output_dir, "pandadock_report.json")
                print("📊 Generating JSON report...")
                logger.info("Generating JSON report")
                report_generator.export_data(results, format='json', output_path=json_output_path)
            elif config.io.report_format == 'csv':
                csv_output_path = os.path.join(config.io.output_dir, "pandadock_report.csv")
                print("📊 Generating CSV report...")
                logger.info("Generating CSV report")
                report_generator.export_data(results, format='csv', output_path=csv_output_path)
            else: # Default to HTML
                print("📊 Generating interactive HTML report...")
                logger.info("Generating HTML report")
                report_generator.generate_report(results)
        
        print(f"✅ Docking completed successfully!")
        print(f"📂 Results saved to: {config.io.output_dir}")
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
    
    # Check for updates (non-blocking)
    try:
        check_and_notify_updates(silent=False)
    except Exception:
        pass  # Don't let update check fail the main program
    
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