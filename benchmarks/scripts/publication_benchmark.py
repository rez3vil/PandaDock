#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PDBbind Benchmark for PandaDock

This script benchmarks PandaDock's docking capabilities against the PDBbind
database. It supports various docking engines, including physics-based,
ML-based, and specialized metal docking. The script generates
publication-ready plots and a comprehensive performance analysis.

Usage:
    python publication_benchmark.py --pdbbind_dir /path/to/pdbbind --output_dir results/
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from scipy.stats import pearsonr, spearmanr
import warnings
from biopandas.pdb import PandasPdb

# Add PandaDock to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(project_root)

from pandadock.docking.metal_docking_engine import MetalDockingEngine, MetalDockingConfig
from pandadock.docking.physics_engine import PhysicsEngine
from pandadock.docking.ml_engine import MLEngine
from pandadock.docking.base_engine import Pose
from pandadock.scoring.metal_scoring import MetalScoringFunction
from pandadock.utils.ic50_calculator import IC50Calculator
from pandadock.config import PandaDockConfig

warnings.filterwarnings('ignore')

# Set up plotting style for publication
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 11,
    'figure.titlesize': 16,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.format': 'png'
})


@dataclass
class BenchmarkEntry:
    """Represents a single benchmark entry"""
    pdb_code: str
    binding_affinity: float
    affinity_type: str
    resolution: float
    year: int
    has_metal: bool
    metal_types: List[str]
    ligand_file: Path
    receptor_file: Path
    crystal_pose_file: Path
    binding_site: Dict[str, float]
    molecular_weight: float
    heavy_atoms: int


@dataclass
class DockingResult:
    """Represents docking results for a single complex"""
    pdb_code: str
    predicted_score: float
    predicted_affinity: float
    experimental_affinity: float
    rmsd_best_pose: float
    success_rate: float
    docking_time: float
    num_poses: int
    metal_coordination: bool
    coordination_score: float
    engine_type: str
    confidence: float


class PDBbindProcessor:
    """Processes PDBbind database for docking benchmark"""

    def __init__(self, pdbbind_dir: str):
        self.pdbbind_dir = Path(pdbbind_dir)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.index_file = self.pdbbind_dir / "index" / "INDEX_demo_PL_data.2021"
        self.metal_elements = {'ZN', 'FE', 'MG', 'CA', 'MN', 'CU', 'NI', 'CO', 'MO', 'W', 'V', 'CR'}

    def load_pdbbind_index(self) -> pd.DataFrame:
        """Load PDBbind index file with binding affinity data"""
        self.logger.info("Loading PDBbind index file...")
        if not self.index_file.exists():
            self.logger.error(f"PDBbind index file not found at: {self.index_file}")
            raise FileNotFoundError(f"PDBbind index file not found at: {self.index_file}")

        entries = []
        with open(self.index_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 5:
                    pdb_code, resolution, year, affinity_data = parts[0], parts[1], parts[2], parts[3]
                    if 'Kd' in affinity_data or 'Ki' in affinity_data or 'IC50' in affinity_data:
                        try:
                            value_str = affinity_data.split('=')[1]
                            affinity_value = float(value_str[:-2])
                            affinity_type = value_str[-2:]
                            entries.append({
                                'pdb_code': pdb_code,
                                'resolution': float(resolution) if resolution != 'NMR' else 2.5,
                                'year': int(year),
                                'binding_affinity': affinity_value,
                                'affinity_type': affinity_type
                            })
                        except (ValueError, IndexError):
                            continue
        df = pd.DataFrame(entries)
        self.logger.info(f"Loaded {len(df)} entries from PDBbind index")
        return df

    def prepare_benchmark_entries(self, df: pd.DataFrame) -> List[BenchmarkEntry]:
        """Prepare benchmark entries with structure files and binding sites"""
        self.logger.info("Preparing benchmark entries...")
        entries = []
        for _, row in df.iterrows():
            pdb_code = row['pdb_code']
            pdb_dir = self.pdbbind_dir / pdb_code
            if not pdb_dir.exists():
                continue

            receptor_file = pdb_dir / f"{pdb_code}_protein.pdb"
            ligand_file = pdb_dir / f"{pdb_code}_ligand.sdf"
            crystal_pose_file = pdb_dir / f"{pdb_code}_ligand.mol2"

            if not receptor_file.exists() or not ligand_file.exists() or not crystal_pose_file.exists():
                continue

            has_metal, metal_types = self._detect_metal(receptor_file)
            
            # Placeholder for binding site, molecular weight and heavy atoms
            binding_site = {'center': [0, 0, 0], 'size': [20, 20, 20]}
            molecular_weight = 300.0
            heavy_atoms = 20

            entry = BenchmarkEntry(
                pdb_code=pdb_code,
                binding_affinity=row['binding_affinity'],
                affinity_type=row['affinity_type'],
                resolution=row['resolution'],
                year=row['year'],
                has_metal=has_metal,
                metal_types=metal_types,
                ligand_file=ligand_file,
                receptor_file=receptor_file,
                crystal_pose_file=crystal_pose_file,
                binding_site=binding_site,
                molecular_weight=molecular_weight,
                heavy_atoms=heavy_atoms
            )
            entries.append(entry)
        self.logger.info(f"Prepared {len(entries)} benchmark entries")
        return entries

    def _detect_metal(self, pdb_file: Path) -> Tuple[bool, List[str]]:
        """Detects metal ions in a PDB file"""
        try:
            pdb = PandasPdb().read_pdb(str(pdb_file))
            metals = pdb.df['HETATM'][pdb.df['HETATM']['residue_name'].isin(self.metal_elements)]
            if not metals.empty:
                metal_types = metals['residue_name'].unique().tolist()
                return True, metal_types
        except Exception as e:
            self.logger.warning(f"Could not parse {pdb_file} for metal detection: {e}")
        return False, []


class DockingBenchmark:
    """Main benchmark class for evaluating docking performance"""

    def __init__(self, output_dir: str, config_file: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.config = PandaDockConfig.from_file(config_file) if config_file else self._create_default_config()
        self.engines = self._initialize_engines()
        self.results: List[DockingResult] = []
        self.benchmark_data: List[BenchmarkEntry] = []

    def _create_default_config(self):
        """Creates a default configuration for the benchmark"""
        from easydict import EasyDict
        default_pandadock_config = PandaDockConfig()
        return EasyDict(default_pandadock_config.to_dict())

    def _initialize_engines(self):
        """Initialize different docking engines for comparison"""
        self.logger.info("Initializing docking engines...")
        engines = {}
        
        # Physics Engine
        try:
            engines['physics'] = PhysicsEngine(self.config)
        except Exception as e:
            self.logger.warning(f"Could not initialize PhysicsEngine: {e}")

        # ML Engine
        try:
            from pandadock.docking.ml_engine import MLEngine
            engines['ml'] = MLEngine(self.config)
        except ImportError:
            self.logger.warning("Could not import MLEngine. ML engine will be disabled.")
        except Exception as e:
            self.logger.warning(f"Could not initialize MLEngine: {e}")

        # Metal Docking Engines
        metal_configs = {
            "metal_strict": MetalDockingConfig(require_metal_coordination=True),
            "metal_flexible": MetalDockingConfig(require_metal_coordination=False, geometric_constraint_weight=1.5)
        }
        for name, metal_config in metal_configs.items():
            try:
                engines[name] = MetalDockingEngine(self.config, metal_config)
            except Exception as e:
                self.logger.warning(f"Could not initialize {name}: {e}")
        
        self.logger.info(f"Initialized engines: {list(engines.keys())}")
        return engines

    def run_benchmark(self, benchmark_entries: List[BenchmarkEntry], max_entries: Optional[int] = None):
        """Run benchmark on all entries"""
        self.benchmark_data = benchmark_entries
        if max_entries:
            benchmark_entries = benchmark_entries[:max_entries]
        
        self.logger.info(f"Starting benchmark on {len(benchmark_entries)} entries...")
        for entry in benchmark_entries:
            for engine_name, engine in self.engines.items():
                try:
                    result = self._dock_single_complex(entry, engine, engine_name)
                    if result:
                        self.results.append(result)
                except Exception as e:
                    self.logger.error(f"Failed to dock {entry.pdb_code} with {engine_name}: {e}")
        self.logger.info(f"Benchmark completed. Generated {len(self.results)} results.")

    def _dock_single_complex(self, entry: BenchmarkEntry, engine, engine_name: str) -> Optional[DockingResult]:
        """Dock a single protein-ligand complex"""
        import time
        start_time = time.time()
        
        poses = engine.dock(str(entry.receptor_file), str(entry.ligand_file))
        docking_time = time.time() - start_time

        if not poses:
            return None

        best_pose = poses[0]
        rmsd = self.calculate_rmsd(best_pose, entry.crystal_pose_file)
        
        return DockingResult(
            pdb_code=entry.pdb_code,
            predicted_score=best_pose.score,
            predicted_affinity=-best_pose.score,  # Assuming score correlates with affinity
            experimental_affinity=entry.binding_affinity,
            rmsd_best_pose=rmsd,
            success_rate=float(rmsd < 2.0),
            docking_time=docking_time,
            num_poses=len(poses),
            metal_coordination=hasattr(best_pose, 'coordinating_atoms') and best_pose.coordinating_atoms,
            coordination_score=getattr(best_pose, 'coordination_quality', {}).get('coordination_score', 0.0) if hasattr(best_pose, 'coordination_quality') else 0.0,
            engine_type=engine_name,
            confidence=best_pose.confidence if hasattr(best_pose, 'confidence') else 0.0
        )

    def calculate_rmsd(self, pose: Pose, crystal_pose_file: Path) -> float:
        """Calculate RMSD between docked pose and crystal pose"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            docked_mol = self._create_rdkit_mol_from_pose(pose)
            if docked_mol is None: return 10.0

            if str(crystal_pose_file).endswith('.mol2'):
                ref_mol = Chem.MolFromMol2File(str(crystal_pose_file), removeHs=False)
            else:
                ref_mol = Chem.MolFromPDBFile(str(crystal_pose_file), removeHs=False)

            if ref_mol is None: return 10.0
            
            return AllChem.GetBestRMS(docked_mol, ref_mol)
        except Exception as e:
            self.logger.warning(f"RMSD calculation failed: {e}")
            return 10.0

    def _create_rdkit_mol_from_pose(self, pose: Pose) -> Optional[Any]:
        """Creates an RDKit molecule from a pose object"""
        try:
            from rdkit import Chem
            
            # This is a simplified conversion. A real implementation would need atom types.
            mol = Chem.RWMol()
            for i, coord in enumerate(pose.coordinates):
                atom = Chem.Atom(6) # Assume Carbon
                mol.AddAtom(atom)
            
            conf = Chem.Conformer(mol.GetNumAtoms())
            for i, coord in enumerate(pose.coordinates):
                conf.SetAtomPosition(i, tuple(coord))
            mol.AddConformer(conf)
            return mol.GetMol()
        except Exception as e:
            self.logger.warning(f"Could not create RDKit molecule from pose: {e}")
            return None

    def analyze_and_plot(self):
        """Analyze results and generate all plots and reports"""
        if not self.results:
            self.logger.warning("No results to analyze.")
            return

        df = self.get_results_dataframe()
        analysis = self.analyze_results(df)
        self.generate_publication_plots(df, analysis)
        self.generate_benchmark_report(analysis)

    def get_results_dataframe(self) -> pd.DataFrame:
        """Converts results to a pandas DataFrame"""
        results_list = [r.__dict__ for r in self.results]
        df = pd.DataFrame(results_list)
        
        metal_map = {entry.pdb_code: entry.has_metal for entry in self.benchmark_data}
        df['has_metal'] = df['pdb_code'].map(metal_map)
        return df

    def analyze_results(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Comprehensive analysis of benchmark results"""
        self.logger.info("Analyzing benchmark results...")
        analysis = {'by_engine': {}}
        for engine in df['engine_type'].unique():
            engine_data = df[df['engine_type'] == engine]
            if len(engine_data) < 2: continue
            
            pearson_r, _ = pearsonr(engine_data['experimental_affinity'], engine_data['predicted_affinity'])
            spearman_r, _ = spearmanr(engine_data['experimental_affinity'], engine_data['predicted_affinity'])
            rmse = np.sqrt(np.mean((engine_data['experimental_affinity'] - engine_data['predicted_affinity'])**2))
            
            analysis['by_engine'][engine] = {
                'pearson_r': pearson_r,
                'spearman_r': spearman_r,
                'rmse': rmse,
                'rmsd_success_rate': (engine_data['rmsd_best_pose'] < 2.0).mean(),
                'mean_rmsd': engine_data['rmsd_best_pose'].mean(),
                'mean_docking_time': engine_data['docking_time'].mean()
            }
        return analysis

    def generate_publication_plots(self, df: pd.DataFrame, analysis: Dict[str, Any]):
        """Generate all publication-ready plots"""
        self._plot_correlation_analysis(df)
        self._plot_rmsd_analysis(df)
        self._plot_engine_performance(df, analysis)
        self._plot_metal_vs_nonmetal(df)
        self._create_summary_figure(df, analysis)

    def _plot_correlation_analysis(self, df: pd.DataFrame):
        """Plot experimental vs predicted binding affinity correlations"""
        g = sns.FacetGrid(df, col="engine_type", col_wrap=3, hue="has_metal", height=4)
        g.map(sns.scatterplot, "experimental_affinity", "predicted_affinity", alpha=0.7)
        g.add_legend()
        g.set_axis_labels("Experimental Affinity (pKa)", "Predicted Affinity (pKa)")
        g.set_titles("{col_name}")
        for ax in g.axes.flat:
            lims = [
                np.min([ax.get_xlim(), ax.get_ylim()]),
                np.max([ax.get_xlim(), ax.get_ylim()]),
            ]
            ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
            ax.set_aspect('equal', adjustable='box')
        plt.savefig(self.output_dir / "correlation_analysis.png")
        plt.close()

    def _plot_rmsd_analysis(self, df: pd.DataFrame):
        """Plot RMSD distribution and success rates"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        sns.kdeplot(data=df, x="rmsd_best_pose", hue="engine_type", ax=ax1, fill=True)
        ax1.axvline(2.0, color='r', linestyle='--', label='Success Threshold (2Å)')
        ax1.set_title("RMSD Distribution by Engine")
        ax1.set_xlabel("RMSD (Å)")
        
        success_rates = df.groupby('engine_type')['rmsd_best_pose'].apply(lambda x: (x < 2.0).mean()).reset_index()
        sns.barplot(data=success_rates, x='engine_type', y='rmsd_best_pose', ax=ax2)
        ax2.set_ylabel("Success Rate (RMSD < 2Å)")
        ax2.set_title("Docking Success Rate")
        ax2.set_ylim(0, 1)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(self.output_dir / "rmsd_analysis.png")
        plt.close()

    def _plot_engine_performance(self, df: pd.DataFrame, analysis: Dict[str, Any]):
        """Plot comprehensive engine performance comparison"""
        perf_data = pd.DataFrame(analysis['by_engine']).T.reset_index().rename(columns={'index': 'engine_type'})
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        sns.barplot(data=perf_data, x='engine_type', y='pearson_r', ax=axes[0, 0])
        axes[0, 0].set_title("Pearson Correlation")
        sns.barplot(data=perf_data, x='engine_type', y='rmse', ax=axes[0, 1])
        axes[0, 1].set_title("RMSE")
        sns.barplot(data=perf_data, x='engine_type', y='mean_rmsd', ax=axes[1, 0])
        axes[1, 0].set_title("Mean RMSD")
        sns.barplot(data=perf_data, x='engine_type', y='mean_docking_time', ax=axes[1, 1])
        axes[1, 1].set_title("Mean Docking Time")
        for ax in axes.flat:
            ax.tick_params(axis='x', rotation=45)
        plt.tight_layout()
        plt.savefig(self.output_dir / "engine_performance.png")
        plt.close()

    def _plot_metal_vs_nonmetal(self, df: pd.DataFrame):
        """Plot performance for metal vs non-metal complexes"""
        g = sns.catplot(data=df, x="engine_type", y="rmsd_best_pose", hue="has_metal", kind="bar", height=6, aspect=1.5)
        g.set_axis_labels("Engine", "Mean RMSD (Å)")
        g.fig.suptitle("Metal vs Non-metal Performance")
        plt.tight_layout()
        plt.savefig(self.output_dir / "metal_vs_nonmetal.png")
        plt.close()

    def _create_summary_figure(self, df: pd.DataFrame, analysis: Dict[str, Any]):
        """Create a comprehensive summary figure"""
        fig = plt.figure(figsize=(20, 15))
        gs = fig.add_gridspec(2, 2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])
        
        # Correlation plot
        sns.scatterplot(data=df, x="experimental_affinity", y="predicted_affinity", hue="engine_type", style="has_metal", ax=ax1, alpha=0.7)
        ax1.set_title("Affinity Prediction")
        
        # RMSD plot
        sns.boxplot(data=df, x="engine_type", y="rmsd_best_pose", ax=ax2)
        ax2.axhline(2.0, color='r', linestyle='--')
        ax2.set_title("Pose Accuracy (RMSD)")
        
        # Summary table
        summary_df = pd.DataFrame(analysis['by_engine']).T
        summary_df = summary_df[['pearson_r', 'rmse', 'rmsd_success_rate', 'mean_rmsd', 'mean_docking_time']].round(3)
        ax3.axis('off')
        table = ax3.table(cellText=summary_df.values, colLabels=summary_df.columns, rowLabels=summary_df.index, cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        ax3.set_title("Performance Summary")
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "summary_figure.png")
        plt.close()

    def generate_benchmark_report(self, analysis: Dict[str, Any]):
        """Generate comprehensive benchmark report in Markdown"""
        report_path = self.output_dir / "benchmark_report.md"
        report_content = [
            "# PandaDock Benchmark Report\n",
            "## Overall Performance\n",
            pd.DataFrame(analysis['by_engine']).T.to_markdown(),
            "\n\n## Plots\n",
            "![Correlation Analysis](correlation_analysis.png)\n",
            "![RMSD Analysis](rmsd_analysis.png)\n",
            "![Engine Performance](engine_performance.png)\n",
            "![Metal vs Non-metal](metal_vs_nonmetal.png)\n",
            "![Summary Figure](summary_figure.png)\n"
        ]
        with open(report_path, 'w') as f:
            f.write("\n".join(report_content))
        self.logger.info(f"Benchmark report saved to {report_path}")


def main():
    """Main function to run the benchmark"""
    parser = argparse.ArgumentParser(description='PandaDock Benchmark')
    parser.add_argument('--pdbbind_dir', type=str, required=True, help='Path to PDBbind database directory')
    parser.add_argument('--output_dir', type=str, default='benchmark_results', help='Output directory for results')
    parser.add_argument('--config', type=str, help='Path to configuration file')
    parser.add_argument('--max_entries', type=int, default=None, help='Maximum number of entries to process')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    processor = PDBbindProcessor(args.pdbbind_dir)
    benchmark_df = processor.load_pdbbind_index()
    benchmark_entries = processor.prepare_benchmark_entries(benchmark_df)

    benchmark = DockingBenchmark(args.output_dir, args.config)
    benchmark.run_benchmark(benchmark_entries, args.max_entries)
    benchmark.analyze_and_plot()

    print(f"Benchmark finished. Results are in {args.output_dir}")


if __name__ == "__main__":
    main()
