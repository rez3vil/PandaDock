"""
New main entry point for PandaDock with clean, modular architecture.

This replaces the complex main.py with a simple, clear interface that
delegates to the appropriate modules.
"""

import sys
import logging
from pathlib import Path

from .cli import create_argument_parser, validate_arguments, get_config_from_args
from .core import DockingEngine
from .utils_new import setup_logging, print_welcome_message


def main():
    """Main entry point for PandaDock."""
    
    try:
        # Parse arguments
        parser = create_argument_parser()
        args = parser.parse_args()
        
        # Show welcome message
        print_welcome_message()
        
        # Validate arguments
        validate_arguments(args)
        
        # Setup logging
        logger = setup_logging(args.output if hasattr(args, 'output') else '.')
        logger.info("Starting PandaDock with new modular architecture")
        
        # Handle different commands
        if args.command == 'analyze':
            return handle_analyze_command(args, logger)
        elif args.command == 'prepare':
            return handle_prepare_command(args, logger)
        else:
            # Default to docking
            return handle_dock_command(args, logger)
        
    except KeyboardInterrupt:
        print("\n🐼 PandaDock interrupted by user")
        return 1
    except Exception as e:
        print(f"🐼❌ PandaDock Error: {str(e)}")
        return 1


def handle_dock_command(args, logger):
    """Handle the docking command."""
    
    logger.info(f"Running docking: {Path(args.protein).name} + {Path(args.ligand).name}")
    
    # Convert arguments to configuration
    config = get_config_from_args(args)
    
    # Create docking engine
    engine = DockingEngine(config)
    
    try:
        # Run docking
        result = engine.run_docking(
            protein_path=args.protein,
            ligand_path=args.ligand,
            output_dir=args.output
        )
        
        if result['success']:
            logger.info("✅ Docking completed successfully!")
            print_success_message(result)
            return 0
        else:
            logger.error(f"❌ Docking failed: {result['error']}")
            print_failure_message(result)
            return 1
            
    finally:
        engine.cleanup()


def handle_analyze_command(args, logger):
    """Handle the analysis command."""
    
    logger.info(f"Analyzing results in: {args.results_dir}")
    
    from .analysis import AnalysisEngine
    
    # Create analysis engine
    analyzer = AnalysisEngine()
    
    # Run analysis
    result = analyzer.analyze_results(
        results_dir=args.results_dir,
        cluster_poses=args.cluster_poses,
        analyze_interactions=args.analyze_interactions,
        energy_decomposition=args.energy_decomposition,
        output_format=args.output_format
    )
    
    if result['success']:
        logger.info("✅ Analysis completed successfully!")
        return 0
    else:
        logger.error(f"❌ Analysis failed: {result['error']}")
        return 1


def handle_prepare_command(args, logger):
    """Handle the preparation command."""
    
    logger.info(f"Preparing molecule: {args.input_file}")
    
    from .molecules import StructurePreparation
    
    # Create preparation engine
    prep = StructurePreparation()
    
    # Prepare molecule
    if args.molecule_type == 'protein':
        result = prep.prepare_protein_file(
            input_file=args.input_file,
            output_file=args.output,
            add_hydrogens=args.add_hydrogens,
            remove_water=args.remove_water,
            ph=args.ph
        )
    else:
        result = prep.prepare_ligand_file(
            input_file=args.input_file,
            output_file=args.output,
            add_hydrogens=args.add_hydrogens
        )
    
    if result['success']:
        logger.info("✅ Molecule preparation completed successfully!")
        return 0
    else:
        logger.error(f"❌ Preparation failed: {result['error']}")
        return 1


def print_success_message(result):
    """Print success message with results summary."""
    
    report = result.get('report', {})
    
    print("\n🎉 Docking Completed Successfully! 🎉")
    print("=" * 50)
    
    if report.get('total_poses', 0) > 0:
        print(f"📊 Generated {report['total_poses']} poses")
        if report.get('best_score') is not None:
            print(f"🏆 Best score: {report['best_score']:.2f}")
    
    print(f"⏱️  Total time: {result['elapsed_time']:.1f} seconds")
    print(f"💻 Device used: {report.get('device_used', 'Unknown')}")
    print(f"🧬 Algorithm: {report.get('algorithm_used', 'Unknown')}")
    
    print("\n🐼 Thank you for using PandaDock! 🐼")


def print_failure_message(result):
    """Print failure message with helpful information."""
    
    print("\n❌ Docking Failed")
    print("=" * 30)
    print(f"Error: {result.get('error', 'Unknown error')}")
    print(f"Time elapsed: {result.get('elapsed_time', 0):.1f} seconds")
    
    # Print performance report if available
    if 'performance_report' in result:
        print("\nPerformance Information:")
        print(result['performance_report'])
    
    print("\n💡 Troubleshooting tips:")
    print("• Check that input files are valid")
    print("• Try --fast-mode for quicker results")
    print("• Use --use-gpu false if GPU issues occur")
    print("• Check the log files for detailed error information")


if __name__ == "__main__":
    sys.exit(main())