#!/usr/bin/env python3
"""
PDBbind Benchmark for PandaDock

This script performs comprehensive benchmarking of PandaDock using the PDBbind dataset.
It compares:
1. CPU vs GPU performance
2. Predicted binding affinities vs experimental IC50/Kd values
3. Different algorithms and scoring functions

Usage:
    python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --index_file INDEX_refined_set.2020 --output benchmark_results

Requirements:
    - PDBbind dataset with refined set
    - PandaDock installed
    - GPU support (optional)
"""

import os
import sys
import argparse
import json
import time
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
import traceback

# Add PandaDock to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: matplotlib/seaborn not available. Plots will be skipped.")

try:
    from scipy import stats
    from sklearn.metrics import mean_squared_error, r2_score
    HAS_STATS = True
except ImportError:
    HAS_STATS = False
    print("Warning: scipy/sklearn not available. Advanced statistics will be limited.")


@dataclass
class PDBBindEntry:
    """Represents a single entry from PDBbind dataset."""
    pdb_code: str
    resolution: float
    release_year: int
    binding_data_raw: str
    binding_value: float  # Converted to standard units (M)
    binding_type: str  # 'Kd', 'Ki', 'IC50'
    reference: str
    ligand_name: str
    protein_file: Optional[str] = None
    ligand_file: Optional[str] = None
    pocket_file: Optional[str] = None


@dataclass
class BenchmarkResult:
    """Results from a single docking run."""
    pdb_code: str
    algorithm: str
    device: str  # 'CPU' or 'GPU'
    scoring_function: str
    
    # Performance metrics
    total_time: float
    setup_time: float
    docking_time: float
    
    # Docking results
    best_score: float
    poses_generated: int
    success: bool
    error_message: Optional[str] = None
    
    # Binding affinity predictions
    predicted_delta_g: Optional[float] = None
    predicted_kd: Optional[float] = None
    predicted_ic50: Optional[float] = None
    predicted_ki: Optional[float] = None
    
    # Experimental values
    experimental_value: Optional[float] = None
    experimental_type: Optional[str] = None
    
    # Correlation metrics (filled during analysis)
    correlation_r: Optional[float] = None
    rmse: Optional[float] = None


class PDBBindParser:
    """Parser for PDBbind dataset index files and directories."""
    
    def __init__(self, pdbbind_dir: str, index_file: str):
        """
        Initialize parser.
        
        Args:
            pdbbind_dir: Path to PDBbind dataset directory
            index_file: Name of index file (e.g., INDEX_refined_set.2020)
        """
        self.pdbbind_dir = Path(pdbbind_dir)
        self.index_file = self.pdbbind_dir / index_file
        self.logger = logging.getLogger(__name__)
        
        if not self.pdbbind_dir.exists():
            raise FileNotFoundError(f"PDBbind directory not found: {pdbbind_dir}")
        if not self.index_file.exists():
            raise FileNotFoundError(f"Index file not found: {self.index_file}")
    
    def parse_binding_data(self, binding_str: str) -> Tuple[float, str]:
        """
        Parse binding data string and convert to standard units (M).
        
        Args:
            binding_str: String like "Kd=49uM", "Ki=190uM", "IC50=1.2nM"
            
        Returns:
            Tuple of (value_in_M, binding_type)
        """
        # Remove spaces and make uppercase for easier parsing
        binding_str = binding_str.strip().upper()
        
        # Extract binding type and value
        if '=' not in binding_str:
            raise ValueError(f"Invalid binding data format: {binding_str}")
        
        binding_type, value_str = binding_str.split('=', 1)
        binding_type = binding_type.strip()
        
        # Parse value and unit
        value_pattern = r'([0-9.]+)\s*([A-Z]+)'
        match = re.match(value_pattern, value_str.strip())
        
        if not match:
            raise ValueError(f"Could not parse value from: {value_str}")
        
        value = float(match.group(1))
        unit = match.group(2)
        
        # Convert to Molar
        unit_conversions = {
            'M': 1.0,
            'MM': 1e-3,  # millimolar
            'UM': 1e-6,  # micromolar
            'NM': 1e-9,  # nanomolar
            'PM': 1e-12, # picomolar
            'FM': 1e-15  # femtomolar
        }
        
        if unit not in unit_conversions:
            raise ValueError(f"Unknown unit: {unit}")
        
        value_in_m = value * unit_conversions[unit]
        
        return value_in_m, binding_type
    
    def parse_index_file(self) -> List[PDBBindEntry]:
        """Parse the PDBbind index file and return list of entries."""
        entries = []
        
        with open(self.index_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                try:
                    # Parse line format: PDB resolution year binding_data reference ligand_name
                    parts = line.split(None, 5)  # Split on whitespace, max 6 parts
                    
                    if len(parts) < 5:
                        self.logger.warning(f"Skipping malformed line {line_num}: {line}")
                        continue
                    
                    pdb_code = parts[0].lower()
                    resolution = float(parts[1])
                    release_year = int(parts[2])
                    binding_data_raw = parts[3]
                    reference = parts[4] if len(parts) > 4 else ""
                    ligand_name = parts[5] if len(parts) > 5 else ""
                    
                    # Parse binding data
                    try:
                        binding_value, binding_type = self.parse_binding_data(binding_data_raw)
                    except ValueError as e:
                        self.logger.warning(f"Could not parse binding data for {pdb_code}: {e}")
                        continue
                    
                    # Check if files exist
                    entry_dir = self.pdbbind_dir / pdb_code
                    if not entry_dir.exists():
                        self.logger.warning(f"Directory not found for {pdb_code}")
                        continue
                    
                    protein_file = self._find_file(entry_dir, f"{pdb_code}_protein.pdb")
                    ligand_file = self._find_file(entry_dir, f"{pdb_code}_ligand.sdf") or \
                                 self._find_file(entry_dir, f"{pdb_code}_ligand.mol2")
                    pocket_file = self._find_file(entry_dir, f"{pdb_code}_pocket.pdb")
                    
                    if not protein_file or not ligand_file:
                        self.logger.warning(f"Missing files for {pdb_code}")
                        continue
                    
                    entry = PDBBindEntry(
                        pdb_code=pdb_code,
                        resolution=resolution,
                        release_year=release_year,
                        binding_data_raw=binding_data_raw,
                        binding_value=binding_value,
                        binding_type=binding_type,
                        reference=reference,
                        ligand_name=ligand_name,
                        protein_file=str(protein_file),
                        ligand_file=str(ligand_file),
                        pocket_file=str(pocket_file) if pocket_file else None
                    )
                    
                    entries.append(entry)
                    
                except (ValueError, IndexError) as e:
                    self.logger.warning(f"Error parsing line {line_num}: {e}")
                    continue
        
        self.logger.info(f"Parsed {len(entries)} valid entries from {self.index_file}")
        return entries
    
    def _find_file(self, directory: Path, filename: str) -> Optional[Path]:
        """Find file in directory, case-insensitive."""
        file_path = directory / filename
        if file_path.exists():
            return file_path
        
        # Try case-insensitive search
        for file in directory.iterdir():
            if file.name.lower() == filename.lower():
                return file
        
        return None


class PandaDockBenchmark:
    """Main benchmark class for running PandaDock benchmarks."""
    
    def __init__(self, output_dir: str, grid_centers_file: Optional[str] = None):
        """Initialize benchmark runner."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.logger = self._setup_logging()
        
        # Results storage
        self.results: List[BenchmarkResult] = []
        
        # Load grid centers if provided
        self.grid_centers = {}
        if grid_centers_file:
            self.grid_centers = self._load_grid_centers(grid_centers_file)
    
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration."""
        logger = logging.getLogger('pandadock_benchmark')
        logger.setLevel(logging.INFO)
        
        # File handler
        log_file = self.output_dir / f"benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # Formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        
        return logger
    
    def _load_grid_centers(self, grid_centers_file: str) -> Dict[str, Tuple[float, float, float]]:
        """
        Load grid centers from CSV file.
        
        Args:
            grid_centers_file: Path to CSV file with columns: ProteinID, X, Y, Z
            
        Returns:
            Dictionary mapping protein IDs to (x, y, z) coordinates
        """
        grid_centers = {}
        
        try:
            df = pd.read_csv(grid_centers_file)
            
            # Validate required columns
            required_cols = ['ProteinID', 'X', 'Y', 'Z']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing required columns in grid centers file: {missing_cols}")
            
            # Load grid centers
            for _, row in df.iterrows():
                protein_id = str(row['ProteinID']).lower()  # Convert to lowercase for matching
                x, y, z = float(row['X']), float(row['Y']), float(row['Z'])
                grid_centers[protein_id] = (x, y, z)
            
            self.logger.info(f"Loaded {len(grid_centers)} grid centers from {grid_centers_file}")
            
            # Log first few entries for verification
            for i, (protein_id, coords) in enumerate(list(grid_centers.items())[:5]):
                self.logger.info(f"  {protein_id}: ({coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f})")
            
            if len(grid_centers) > 5:
                self.logger.info(f"  ... and {len(grid_centers) - 5} more")
            
        except Exception as e:
            self.logger.error(f"Failed to load grid centers from {grid_centers_file}: {e}")
            raise
        
        return grid_centers
    
    def extract_binding_site_from_pocket(self, pocket_file: str) -> Tuple[List[float], float]:
        """
        Extract binding site center and radius from pocket PDB file.
        
        Args:
            pocket_file: Path to pocket PDB file
            
        Returns:
            Tuple of (center_coords, radius)
        """
        try:
            coords = []
            
            with open(pocket_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Extract coordinates from PDB format
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append([x, y, z])
            
            if not coords:
                raise ValueError("No coordinates found in pocket file")
            
            coords = np.array(coords)
            
            # Calculate center as mean of coordinates
            center = np.mean(coords, axis=0).tolist()
            
            # Calculate radius as max distance from center + some padding
            distances = np.linalg.norm(coords - center, axis=1)
            radius = float(np.max(distances)) + 2.0  # Add 2Å padding
            
            return center, radius
            
        except Exception as e:
            self.logger.warning(f"Could not extract binding site from {pocket_file}: {e}")
            # Return default values
            return [0.0, 0.0, 0.0], 10.0
    
    def run_single_docking(self, entry: PDBBindEntry, algorithm: str = "genetic",
                          device: str = "CPU", scoring: str = "standard",
                          iterations: int = 50) -> BenchmarkResult:
        """
        Run a single docking job and return results.
        
        Args:
            entry: PDBbind entry to dock
            algorithm: Docking algorithm ('genetic', 'random', 'pandadock')
            device: Computing device ('CPU' or 'GPU')
            scoring: Scoring function ('standard', 'enhanced', 'physics-based')
            iterations: Number of iterations
            
        Returns:
            BenchmarkResult object
        """
        start_time = time.time()
        
        # Prepare output directory
        run_output = self.output_dir / f"{entry.pdb_code}_{algorithm}_{device}_{scoring}"
        run_output.mkdir(exist_ok=True)
        
        try:
            # Use grid centers if available, otherwise extract from pocket file
            if entry.pdb_code in self.grid_centers:
                center = list(self.grid_centers[entry.pdb_code])
                radius = 10.0  # Default radius for grid centers
                self.logger.debug(f"Using grid center for {entry.pdb_code}: {center}")
            elif entry.pocket_file:
                center, radius = self.extract_binding_site_from_pocket(entry.pocket_file)
                self.logger.debug(f"Extracted binding site from pocket file for {entry.pdb_code}")
            else:
                # Use default binding site
                center, radius = [0.0, 0.0, 0.0], 15.0
                self.logger.warning(f"No grid center or pocket file for {entry.pdb_code}, using default binding site")
            
            setup_time = time.time() - start_time
            
            # Build PandaDock command
            cmd = [
                "python", "-m", "pandadock",
                "-p", entry.protein_file,
                "-l", entry.ligand_file,
                "-o", str(run_output),
                "-a", algorithm,
                "-i", str(iterations),
                "-s"] + [str(c) for c in center] + [
                "-r", str(radius),
                "--fast-mode"
            ]
            
            # Add device-specific options
            if device == "GPU":
                cmd.extend(["--use-gpu"])
            
            # Add scoring function options
            if scoring == "enhanced":
                cmd.extend(["--enhanced-scoring"])
            elif scoring == "physics-based":
                cmd.extend(["--physics-based"])
            
            # Run docking
            self.logger.info(f"Running {entry.pdb_code} with {algorithm}/{device}/{scoring}")
            docking_start = time.time()
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )
            
            docking_time = time.time() - docking_start
            total_time = time.time() - start_time
            
            if result.returncode == 0:
                # Parse results
                return self._parse_docking_results(
                    entry, run_output, algorithm, device, scoring,
                    total_time, setup_time, docking_time
                )
            else:
                # Docking failed
                error_msg = result.stderr or "Unknown error"
                self.logger.error(f"Docking failed for {entry.pdb_code}: {error_msg}")
                
                return BenchmarkResult(
                    pdb_code=entry.pdb_code,
                    algorithm=algorithm,
                    device=device,
                    scoring_function=scoring,
                    total_time=total_time,
                    setup_time=setup_time,
                    docking_time=docking_time,
                    best_score=float('inf'),
                    poses_generated=0,
                    success=False,
                    error_message=error_msg,
                    experimental_value=entry.binding_value,
                    experimental_type=entry.binding_type
                )
        
        except subprocess.TimeoutExpired:
            error_msg = f"Docking timeout after 10 minutes"
            self.logger.error(f"Timeout for {entry.pdb_code}")
            
            return BenchmarkResult(
                pdb_code=entry.pdb_code,
                algorithm=algorithm,
                device=device,
                scoring_function=scoring,
                total_time=600.0,
                setup_time=setup_time if 'setup_time' in locals() else 0.0,
                docking_time=600.0,
                best_score=float('inf'),
                poses_generated=0,
                success=False,
                error_message=error_msg,
                experimental_value=entry.binding_value,
                experimental_type=entry.binding_type
            )
            
        except Exception as e:
            error_msg = f"Unexpected error: {str(e)}"
            self.logger.error(f"Error for {entry.pdb_code}: {error_msg}")
            
            return BenchmarkResult(
                pdb_code=entry.pdb_code,
                algorithm=algorithm,
                device=device,
                scoring_function=scoring,
                total_time=time.time() - start_time,
                setup_time=0.0,
                docking_time=0.0,
                best_score=float('inf'),
                poses_generated=0,
                success=False,
                error_message=error_msg,
                experimental_value=entry.binding_value,
                experimental_type=entry.binding_type
            )
    
    def _parse_docking_results(self, entry: PDBBindEntry, output_dir: Path,
                              algorithm: str, device: str, scoring: str,
                              total_time: float, setup_time: float, 
                              docking_time: float) -> BenchmarkResult:
        """Parse results from a completed docking run."""
        try:
            # Read JSON results
            json_file = output_dir / "results_summary.json"
            if json_file.exists():
                with open(json_file, 'r') as f:
                    results_data = json.load(f)
                
                poses = results_data.get('poses', [])
                stats = results_data.get('statistics', {})
                
                best_score = stats.get('best_score', float('inf'))
                poses_generated = len(poses)
                
                # Try to read binding affinity report
                affinity_data = self._read_binding_affinity_results(output_dir)
                
                return BenchmarkResult(
                    pdb_code=entry.pdb_code,
                    algorithm=algorithm,
                    device=device,
                    scoring_function=scoring,
                    total_time=total_time,
                    setup_time=setup_time,
                    docking_time=docking_time,
                    best_score=best_score,
                    poses_generated=poses_generated,
                    success=True,
                    predicted_delta_g=affinity_data.get('delta_g'),
                    predicted_kd=affinity_data.get('kd'),
                    predicted_ic50=affinity_data.get('ic50'),
                    predicted_ki=affinity_data.get('ki'),
                    experimental_value=entry.binding_value,
                    experimental_type=entry.binding_type
                )
            else:
                raise FileNotFoundError("Results summary not found")
                
        except Exception as e:
            self.logger.warning(f"Could not parse results for {entry.pdb_code}: {e}")
            
            return BenchmarkResult(
                pdb_code=entry.pdb_code,
                algorithm=algorithm,
                device=device,
                scoring_function=scoring,
                total_time=total_time,
                setup_time=setup_time,
                docking_time=docking_time,
                best_score=float('inf'),
                poses_generated=0,
                success=False,
                error_message=f"Could not parse results: {str(e)}",
                experimental_value=entry.binding_value,
                experimental_type=entry.binding_type
            )
    
    def _read_binding_affinity_results(self, output_dir: Path) -> Dict[str, float]:
        """Read binding affinity predictions from output directory."""
        try:
            csv_file = output_dir / "binding_affinity_report.csv"
            if csv_file.exists():
                df = pd.read_csv(csv_file)
                if len(df) > 0:
                    # Get best pose (first row)
                    row = df.iloc[0]
                    return {
                        'delta_g': float(row['Delta_G_kcal_mol']),
                        'kd': float(row['Kd_M'].replace('e', 'E')) if 'e' in str(row['Kd_M']) else float(row['Kd_M']),
                        'ic50': float(row['IC50_M'].replace('e', 'E')) if 'e' in str(row['IC50_M']) else float(row['IC50_M']),
                        'ki': float(row['Ki_M'].replace('e', 'E')) if 'e' in str(row['Ki_M']) else float(row['Ki_M'])
                    }
        except Exception as e:
            self.logger.debug(f"Could not read binding affinity data: {e}")
        
        return {}
    
    def run_benchmark(self, entries: List[PDBBindEntry], 
                     algorithms: List[str] = ["genetic"],
                     devices: List[str] = ["CPU"], 
                     scoring_functions: List[str] = ["standard"],
                     iterations: int = 50,
                     max_entries: Optional[int] = None,
                     parallel_jobs: int = 1) -> List[BenchmarkResult]:
        """
        Run comprehensive benchmark on multiple entries.
        
        Args:
            entries: List of PDBbind entries
            algorithms: List of algorithms to test
            devices: List of devices to test  
            scoring_functions: List of scoring functions to test
            iterations: Number of iterations per run
            max_entries: Maximum entries to process (for testing)
            parallel_jobs: Number of parallel jobs
            
        Returns:
            List of benchmark results
        """
        if max_entries:
            entries = entries[:max_entries]
        
        self.logger.info(f"Starting benchmark with {len(entries)} entries")
        self.logger.info(f"Testing: algorithms={algorithms}, devices={devices}, scoring={scoring_functions}")
        
        # Generate all combinations
        tasks = []
        for entry in entries:
            for algorithm in algorithms:
                for device in devices:
                    for scoring in scoring_functions:
                        tasks.append((entry, algorithm, device, scoring, iterations))
        
        self.logger.info(f"Total benchmark tasks: {len(tasks)}")
        
        # Run tasks
        if parallel_jobs > 1:
            results = self._run_parallel_benchmark(tasks, parallel_jobs)
        else:
            results = self._run_sequential_benchmark(tasks)
        
        self.results.extend(results)
        
        # Save results
        self._save_results(results)
        
        return results
    
    def _run_sequential_benchmark(self, tasks: List[Tuple]) -> List[BenchmarkResult]:
        """Run benchmark tasks sequentially."""
        results = []
        
        for i, (entry, algorithm, device, scoring, iterations) in enumerate(tasks, 1):
            self.logger.info(f"Task {i}/{len(tasks)}: {entry.pdb_code} {algorithm}/{device}/{scoring}")
            
            result = self.run_single_docking(entry, algorithm, device, scoring, iterations)
            results.append(result)
            
            # Log progress
            if i % 10 == 0:
                success_rate = sum(1 for r in results if r.success) / len(results) * 100
                self.logger.info(f"Progress: {i}/{len(tasks)} ({success_rate:.1f}% success rate)")
        
        return results
    
    def _run_parallel_benchmark(self, tasks: List[Tuple], parallel_jobs: int) -> List[BenchmarkResult]:
        """Run benchmark tasks in parallel."""
        results = []
        
        with ProcessPoolExecutor(max_workers=parallel_jobs) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(self.run_single_docking, entry, algorithm, device, scoring, iterations): 
                (entry, algorithm, device, scoring) 
                for entry, algorithm, device, scoring, iterations in tasks
            }
            
            # Collect results
            for i, future in enumerate(as_completed(future_to_task), 1):
                try:
                    result = future.result(timeout=60)
                    results.append(result)
                    
                    if i % 10 == 0:
                        success_rate = sum(1 for r in results if r.success) / len(results) * 100
                        self.logger.info(f"Progress: {i}/{len(tasks)} ({success_rate:.1f}% success rate)")
                        
                except Exception as e:
                    entry, algorithm, device, scoring = future_to_task[future]
                    self.logger.error(f"Task failed: {entry.pdb_code} {algorithm}/{device}/{scoring} - {e}")
        
        return results
    
    def _save_results(self, results: List[BenchmarkResult]) -> None:
        """Save benchmark results to files."""
        # Save as JSON
        json_file = self.output_dir / "benchmark_results.json"
        with open(json_file, 'w') as f:
            json.dump([asdict(r) for r in results], f, indent=2)
        
        # Save as CSV
        csv_file = self.output_dir / "benchmark_results.csv"
        df = pd.DataFrame([asdict(r) for r in results])
        df.to_csv(csv_file, index=False)
        
        self.logger.info(f"Results saved to {json_file} and {csv_file}")
    
    def analyze_results(self, results: Optional[List[BenchmarkResult]] = None) -> Dict[str, Any]:
        """
        Analyze benchmark results and generate statistics.
        
        Args:
            results: Results to analyze (default: self.results)
            
        Returns:
            Dictionary with analysis results
        """
        if results is None:
            results = self.results
        
        if not results:
            self.logger.warning("No results to analyze")
            return {}
        
        analysis = {}
        
        # Basic statistics
        total_runs = len(results)
        successful_runs = sum(1 for r in results if r.success)
        success_rate = successful_runs / total_runs * 100 if total_runs > 0 else 0
        
        analysis['basic_stats'] = {
            'total_runs': total_runs,
            'successful_runs': successful_runs,
            'success_rate': success_rate,
            'average_runtime': np.mean([r.total_time for r in results if r.success]),
            'median_runtime': np.median([r.total_time for r in results if r.success])
        }
        
        # Performance comparison (CPU vs GPU)
        cpu_results = [r for r in results if r.device == 'CPU' and r.success]
        gpu_results = [r for r in results if r.device == 'GPU' and r.success]
        
        if cpu_results and gpu_results:
            cpu_times = [r.total_time for r in cpu_results]
            gpu_times = [r.total_time for r in gpu_results]
            
            analysis['performance_comparison'] = {
                'cpu_avg_time': np.mean(cpu_times),
                'gpu_avg_time': np.mean(gpu_times),
                'speedup_factor': np.mean(cpu_times) / np.mean(gpu_times) if gpu_times else 1.0,
                'cpu_std': np.std(cpu_times),
                'gpu_std': np.std(gpu_times)
            }
        
        # Binding affinity correlation
        if HAS_STATS:
            analysis['affinity_correlation'] = self._analyze_affinity_correlation(results)
        
        # Algorithm comparison
        analysis['algorithm_comparison'] = self._analyze_algorithm_performance(results)
        
        # Save analysis
        analysis_file = self.output_dir / "benchmark_analysis.json"
        with open(analysis_file, 'w') as f:
            json.dump(analysis, f, indent=2)
        
        self.logger.info(f"Analysis saved to {analysis_file}")
        
        return analysis
    
    def _analyze_affinity_correlation(self, results: List[BenchmarkResult]) -> Dict[str, Any]:
        """Analyze correlation between predicted and experimental binding affinities."""
        # Filter results with both predicted and experimental values
        valid_results = [
            r for r in results 
            if r.success and r.predicted_kd is not None and r.experimental_value is not None
        ]
        
        if len(valid_results) < 10:
            return {'error': 'Insufficient data for correlation analysis'}
        
        # Extract values for correlation
        experimental_values = [r.experimental_value for r in valid_results]
        predicted_kd = [r.predicted_kd for r in valid_results]
        predicted_ic50 = [r.predicted_ic50 for r in valid_results if r.predicted_ic50 is not None]
        
        # Convert to log scale for better correlation
        log_exp = np.log10(experimental_values)
        log_pred_kd = np.log10(predicted_kd)
        log_pred_ic50 = np.log10(predicted_ic50) if predicted_ic50 else []
        
        correlation_results = {}
        
        # Kd correlation
        if len(log_pred_kd) > 0:
            r_kd, p_kd = stats.pearsonr(log_exp, log_pred_kd)
            rmse_kd = np.sqrt(mean_squared_error(log_exp, log_pred_kd))
            r2_kd = r2_score(log_exp, log_pred_kd)
            
            correlation_results['kd_correlation'] = {
                'pearson_r': float(r_kd),
                'p_value': float(p_kd),
                'rmse': float(rmse_kd),
                'r_squared': float(r2_kd),
                'n_samples': len(log_pred_kd)
            }
        
        # IC50 correlation (if available)
        if len(log_pred_ic50) > 0:
            r_ic50, p_ic50 = stats.pearsonr(log_exp[:len(log_pred_ic50)], log_pred_ic50)
            rmse_ic50 = np.sqrt(mean_squared_error(log_exp[:len(log_pred_ic50)], log_pred_ic50))
            r2_ic50 = r2_score(log_exp[:len(log_pred_ic50)], log_pred_ic50)
            
            correlation_results['ic50_correlation'] = {
                'pearson_r': float(r_ic50),
                'p_value': float(p_ic50),
                'rmse': float(rmse_ic50),
                'r_squared': float(r2_ic50),
                'n_samples': len(log_pred_ic50)
            }
        
        return correlation_results
    
    def _analyze_algorithm_performance(self, results: List[BenchmarkResult]) -> Dict[str, Any]:
        """Analyze performance differences between algorithms."""
        algorithm_stats = {}
        
        algorithms = set(r.algorithm for r in results)
        
        for algorithm in algorithms:
            alg_results = [r for r in results if r.algorithm == algorithm]
            successful = [r for r in alg_results if r.success]
            
            if successful:
                algorithm_stats[algorithm] = {
                    'total_runs': len(alg_results),
                    'successful_runs': len(successful),
                    'success_rate': len(successful) / len(alg_results) * 100,
                    'avg_runtime': np.mean([r.total_time for r in successful]),
                    'avg_score': np.mean([r.best_score for r in successful if r.best_score != float('inf')]),
                    'avg_poses': np.mean([r.poses_generated for r in successful])
                }
        
        return algorithm_stats
    
    def generate_plots(self, results: Optional[List[BenchmarkResult]] = None) -> None:
        """Generate visualization plots for benchmark results."""
        if not HAS_PLOTTING:
            self.logger.warning("Plotting libraries not available, skipping plots")
            return
        
        if results is None:
            results = self.results
        
        if not results:
            self.logger.warning("No results to plot")
            return
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create plots directory
        plots_dir = self.output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)
        
        # Plot 1: Performance comparison (CPU vs GPU)
        self._plot_performance_comparison(results, plots_dir)
        
        # Plot 2: Algorithm comparison
        self._plot_algorithm_comparison(results, plots_dir)
        
        # Plot 3: Binding affinity correlation
        self._plot_binding_affinity_correlation(results, plots_dir)
        
        # Plot 4: Success rate by category
        self._plot_success_rates(results, plots_dir)
        
        # Plot 5: Runtime distribution
        self._plot_runtime_distribution(results, plots_dir)
        
        self.logger.info(f"Plots saved to {plots_dir}")
    
    def _plot_performance_comparison(self, results: List[BenchmarkResult], plots_dir: Path) -> None:
        """Plot CPU vs GPU performance comparison."""
        cpu_results = [r for r in results if r.device == 'CPU' and r.success]
        gpu_results = [r for r in results if r.device == 'GPU' and r.success]
        
        if not cpu_results or not gpu_results:
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Runtime comparison
        cpu_times = [r.total_time for r in cpu_results]
        gpu_times = [r.total_time for r in gpu_results]
        
        ax1.boxplot([cpu_times, gpu_times], labels=['CPU', 'GPU'])
        ax1.set_ylabel('Runtime (seconds)')
        ax1.set_title('Runtime Comparison: CPU vs GPU')
        ax1.grid(True, alpha=0.3)
        
        # Speedup histogram
        if len(cpu_times) == len(gpu_times):
            speedups = [c/g for c, g in zip(cpu_times, gpu_times)]
            ax2.hist(speedups, bins=20, alpha=0.7, edgecolor='black')
            ax2.axvline(np.mean(speedups), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(speedups):.1f}x')
            ax2.set_xlabel('Speedup Factor (CPU time / GPU time)')
            ax2.set_ylabel('Frequency')
            ax2.set_title('GPU Speedup Distribution')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'performance_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_algorithm_comparison(self, results: List[BenchmarkResult], plots_dir: Path) -> None:
        """Plot algorithm performance comparison."""
        successful_results = [r for r in results if r.success]
        
        if not successful_results:
            return
        
        algorithms = list(set(r.algorithm for r in successful_results))
        
        if len(algorithms) < 2:
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Runtime by algorithm
        runtime_data = []
        labels = []
        for alg in algorithms:
            alg_results = [r for r in successful_results if r.algorithm == alg]
            runtime_data.append([r.total_time for r in alg_results])
            labels.append(f"{alg}\n(n={len(alg_results)})")
        
        ax1.boxplot(runtime_data, labels=labels)
        ax1.set_ylabel('Runtime (seconds)')
        ax1.set_title('Runtime by Algorithm')
        ax1.grid(True, alpha=0.3)
        
        # Success rate by algorithm
        success_rates = []
        for alg in algorithms:
            alg_total = len([r for r in results if r.algorithm == alg])
            alg_success = len([r for r in results if r.algorithm == alg and r.success])
            success_rates.append(alg_success / alg_total * 100 if alg_total > 0 else 0)
        
        ax2.bar(algorithms, success_rates, alpha=0.7)
        ax2.set_ylabel('Success Rate (%)')
        ax2.set_title('Success Rate by Algorithm')
        ax2.set_ylim(0, 100)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'algorithm_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_binding_affinity_correlation(self, results: List[BenchmarkResult], plots_dir: Path) -> None:
        """Plot correlation between predicted and experimental binding affinities."""
        valid_results = [
            r for r in results 
            if r.success and r.predicted_kd is not None and r.experimental_value is not None
        ]
        
        if len(valid_results) < 10:
            return
        
        experimental = [r.experimental_value for r in valid_results]
        predicted = [r.predicted_kd for r in valid_results]
        
        # Convert to log scale
        log_exp = np.log10(experimental)
        log_pred = np.log10(predicted)
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Scatter plot
        ax.scatter(log_exp, log_pred, alpha=0.6, s=50)
        
        # Add diagonal line (perfect correlation)
        min_val = min(min(log_exp), min(log_pred))
        max_val = max(max(log_exp), max(log_pred))
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.8, 
                label='Perfect correlation')
        
        # Calculate and display correlation
        if HAS_STATS:
            r, p = stats.pearsonr(log_exp, log_pred)
            ax.text(0.05, 0.95, f'R = {r:.3f}\np = {p:.3e}', 
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        ax.set_xlabel('log10(Experimental Kd/Ki/IC50 [M])')
        ax.set_ylabel('log10(Predicted Kd [M])')
        ax.set_title('Predicted vs Experimental Binding Affinity')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'binding_affinity_correlation.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_success_rates(self, results: List[BenchmarkResult], plots_dir: Path) -> None:
        """Plot success rates by different categories."""
        categories = ['device', 'algorithm', 'scoring_function']
        
        fig, axes = plt.subplots(1, len(categories), figsize=(15, 5))
        if len(categories) == 1:
            axes = [axes]
        
        for i, category in enumerate(categories):
            category_values = list(set(getattr(r, category) for r in results))
            success_rates = []
            
            for value in category_values:
                category_results = [r for r in results if getattr(r, category) == value]
                success_count = sum(1 for r in category_results if r.success)
                total_count = len(category_results)
                success_rate = success_count / total_count * 100 if total_count > 0 else 0
                success_rates.append(success_rate)
            
            axes[i].bar(category_values, success_rates, alpha=0.7)
            axes[i].set_ylabel('Success Rate (%)')
            axes[i].set_title(f'Success Rate by {category.replace("_", " ").title()}')
            axes[i].set_ylim(0, 100)
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'success_rates.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_runtime_distribution(self, results: List[BenchmarkResult], plots_dir: Path) -> None:
        """Plot runtime distribution."""
        successful_results = [r for r in results if r.success]
        
        if not successful_results:
            return
        
        runtimes = [r.total_time for r in successful_results]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Histogram
        ax1.hist(runtimes, bins=30, alpha=0.7, edgecolor='black')
        ax1.axvline(np.mean(runtimes), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(runtimes):.1f}s')
        ax1.axvline(np.median(runtimes), color='orange', linestyle='--', 
                   label=f'Median: {np.median(runtimes):.1f}s')
        ax1.set_xlabel('Runtime (seconds)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Runtime Distribution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Box plot by device (if available)
        devices = list(set(r.device for r in successful_results))
        if len(devices) > 1:
            runtime_by_device = []
            for device in devices:
                device_runtimes = [r.total_time for r in successful_results if r.device == device]
                runtime_by_device.append(device_runtimes)
            
            ax2.boxplot(runtime_by_device, labels=devices)
            ax2.set_ylabel('Runtime (seconds)')
            ax2.set_title('Runtime by Device')
            ax2.grid(True, alpha=0.3)
        else:
            ax2.text(0.5, 0.5, 'Insufficient device\nvariability for\ncomparison', 
                    ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Runtime by Device')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'runtime_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()


def main():
    """Main entry point for the benchmark script."""
    parser = argparse.ArgumentParser(
        description="PDBbind Benchmark for PandaDock",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic benchmark (CPU only, 10 entries)
  python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --max_entries 10

  # Using predefined grid centers
  python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --grid_centers grid_centers.csv --max_entries 50

  # CPU vs GPU comparison with grid centers
  python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --grid_centers grid_centers.csv \\
                               --devices CPU GPU --max_entries 50

  # Algorithm comparison with grid centers
  python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --grid_centers grid_centers.csv \\
                               --algorithms genetic random --max_entries 100

  # Full benchmark with grid centers
  python pdbbind_benchmark.py --pdbbind_dir /path/to/pdbbind --grid_centers grid_centers.csv \\
                               --devices CPU GPU --algorithms genetic random pandadock \\
                               --scoring standard enhanced physics-based \\
                               --max_entries 500 --parallel_jobs 4
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--pdbbind_dir', 
        required=True,
        help='Path to PDBbind dataset directory'
    )
    
    # Optional arguments
    parser.add_argument(
        '--index_file', 
        default='INDEX_refined_set.2020',
        help='PDBbind index file name (default: INDEX_refined_set.2020)'
    )
    
    parser.add_argument(
        '--output', 
        default='benchmark_results',
        help='Output directory for results (default: benchmark_results)'
    )
    
    parser.add_argument(
        '--grid_centers', 
        help='CSV file with grid centers (columns: ProteinID, X, Y, Z)'
    )
    
    parser.add_argument(
        '--algorithms', 
        nargs='+',
        default=['genetic'],
        choices=['genetic', 'random', 'pandadock'],
        help='Algorithms to benchmark (default: genetic)'
    )
    
    parser.add_argument(
        '--devices', 
        nargs='+',
        default=['CPU'],
        choices=['CPU', 'GPU'],
        help='Devices to benchmark (default: CPU)'
    )
    
    parser.add_argument(
        '--scoring', 
        nargs='+',
        default=['standard'],
        choices=['standard', 'enhanced', 'physics-based'],
        help='Scoring functions to benchmark (default: standard)'
    )
    
    parser.add_argument(
        '--iterations', 
        type=int,
        default=50,
        help='Number of iterations per docking run (default: 50)'
    )
    
    parser.add_argument(
        '--max_entries', 
        type=int,
        help='Maximum number of entries to process (for testing)'
    )
    
    parser.add_argument(
        '--parallel_jobs', 
        type=int,
        default=1,
        help='Number of parallel jobs (default: 1)'
    )
    
    parser.add_argument(
        '--skip_plots', 
        action='store_true',
        help='Skip generating plots'
    )
    
    parser.add_argument(
        '--verbose', 
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    try:
        # Parse PDBbind dataset
        print("Parsing PDBbind dataset...")
        parser = PDBBindParser(args.pdbbind_dir, args.index_file)
        entries = parser.parse_index_file()
        
        if not entries:
            print("No valid entries found in PDBbind dataset")
            return 1
        
        print(f"Found {len(entries)} valid entries")
        
        # Initialize benchmark
        benchmark = PandaDockBenchmark(args.output, args.grid_centers)
        
        # Run benchmark
        print(f"Starting benchmark with {len(args.algorithms)} algorithms, "
              f"{len(args.devices)} devices, {len(args.scoring)} scoring functions")
        
        results = benchmark.run_benchmark(
            entries=entries,
            algorithms=args.algorithms,
            devices=args.devices,
            scoring_functions=args.scoring,
            iterations=args.iterations,
            max_entries=args.max_entries,
            parallel_jobs=args.parallel_jobs
        )
        
        # Analyze results
        print("Analyzing results...")
        analysis = benchmark.analyze_results(results)
        
        # Generate plots
        if not args.skip_plots:
            print("Generating plots...")
            benchmark.generate_plots(results)
        
        # Print summary
        print("\n" + "="*60)
        print("BENCHMARK SUMMARY")
        print("="*60)
        
        if 'basic_stats' in analysis:
            stats = analysis['basic_stats']
            print(f"Total runs: {stats['total_runs']}")
            print(f"Successful runs: {stats['successful_runs']}")
            print(f"Success rate: {stats['success_rate']:.1f}%")
            print(f"Average runtime: {stats['average_runtime']:.1f} seconds")
        
        if 'performance_comparison' in analysis:
            perf = analysis['performance_comparison']
            print(f"\nCPU average time: {perf['cpu_avg_time']:.1f} seconds")
            print(f"GPU average time: {perf['gpu_avg_time']:.1f} seconds")
            print(f"GPU speedup: {perf['speedup_factor']:.1f}x")
        
        if 'affinity_correlation' in analysis:
            corr = analysis['affinity_correlation']
            if 'kd_correlation' in corr:
                kd_corr = corr['kd_correlation']
                print(f"\nBinding affinity correlation (Kd):")
                print(f"  Pearson R: {kd_corr['pearson_r']:.3f}")
                print(f"  R²: {kd_corr['r_squared']:.3f}")
                print(f"  RMSE: {kd_corr['rmse']:.3f}")
        
        print(f"\nResults saved to: {args.output}")
        print("="*60)
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())