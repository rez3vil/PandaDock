"""
Binding affinity calculations for molecular docking.

This module provides comprehensive binding affinity analysis including:
- Free energy calculations (ΔG)
- Dissociation constant (Kd) estimation
- IC50 value estimation
- Ki value calculation
- Thermodynamic analysis
"""

import numpy as np
import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
import math

# Physical constants
R = 1.98720425864083e-3  # Gas constant in kcal/(mol·K)
T_STANDARD = 298.15  # Standard temperature in Kelvin


@dataclass
class BindingAffinityResult:
    """Results from binding affinity analysis."""
    
    pose_rank: int
    docking_score: float
    delta_g: float  # Free energy change (kcal/mol)
    kd: float  # Dissociation constant (M)
    ic50: float  # IC50 value (M)
    ki: float  # Inhibition constant (M)
    efficiency_indices: Dict[str, float]
    confidence: str  # High, Medium, Low
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for export."""
        return {
            'pose_rank': self.pose_rank,
            'docking_score': self.docking_score,
            'delta_g_kcal_mol': self.delta_g,
            'kd_M': self.kd,
            'ic50_M': self.ic50,
            'ki_M': self.ki,
            'ligand_efficiency': self.efficiency_indices.get('LE', 0.0),
            'size_independent_LE': self.efficiency_indices.get('SILE', 0.0),
            'binding_efficiency_index': self.efficiency_indices.get('BEI', 0.0),
            'confidence': self.confidence
        }


class BindingAffinityCalculator:
    """
    Calculator for comprehensive binding affinity analysis.
    
    Provides methods to estimate binding affinity from docking scores
    using established thermodynamic relationships and empirical correlations.
    """
    
    def __init__(self, temperature: float = T_STANDARD, 
                 protein_concentration: float = 1e-9,
                 ligand_concentration: float = 1e-6):
        """
        Initialize binding affinity calculator.
        
        Args:
            temperature: Temperature in Kelvin (default: 298.15K)
            protein_concentration: Protein concentration in M (default: 1 nM)
            ligand_concentration: Ligand concentration in M (default: 1 μM)
        """
        self.temperature = temperature
        self.protein_concentration = protein_concentration
        self.ligand_concentration = ligand_concentration
        self.logger = logging.getLogger(__name__)
    
    def calculate_binding_affinity(self, docking_score: float, pose_rank: int,
                                 molecular_weight: Optional[float] = None,
                                 heavy_atom_count: Optional[int] = None) -> BindingAffinityResult:
        """
        Calculate comprehensive binding affinity metrics.
        
        Args:
            docking_score: Docking score in kcal/mol
            pose_rank: Rank of the pose
            molecular_weight: Molecular weight of ligand (Da)
            heavy_atom_count: Number of heavy atoms in ligand
            
        Returns:
            BindingAffinityResult with all calculated metrics
        """
        # Calculate free energy change (ΔG)
        delta_g = self._score_to_delta_g(docking_score)
        
        # Calculate dissociation constant (Kd)
        kd = self._delta_g_to_kd(delta_g)
        
        # Calculate IC50 and Ki
        ic50 = self._calculate_ic50(kd)
        ki = self._calculate_ki(ic50)
        
        # Calculate efficiency indices
        efficiency_indices = self._calculate_efficiency_indices(
            delta_g, molecular_weight, heavy_atom_count
        )
        
        # Assess confidence
        confidence = self._assess_confidence(docking_score, pose_rank)
        
        return BindingAffinityResult(
            pose_rank=pose_rank,
            docking_score=docking_score,
            delta_g=delta_g,
            kd=kd,
            ic50=ic50,
            ki=ki,
            efficiency_indices=efficiency_indices,
            confidence=confidence
        )
    
    def _score_to_delta_g(self, docking_score: float) -> float:
        """
        Convert docking score to free energy change.
        
        For most scoring functions, the docking score approximates ΔG.
        This can be refined with empirical corrections.
        """
        # Basic conversion (can be refined with calibration)
        # Most scoring functions are designed to approximate ΔG
        delta_g = docking_score
        
        # Apply empirical correction if needed
        # This could be calibrated against experimental data
        # delta_g = self._apply_empirical_correction(docking_score)
        
        return delta_g
    
    def _delta_g_to_kd(self, delta_g: float) -> float:
        """
        Calculate dissociation constant from free energy change.
        
        Uses the fundamental thermodynamic relationship:
        ΔG = RT ln(Kd)
        Therefore: Kd = exp(ΔG / RT)
        """
        try:
            kd = math.exp(delta_g / (R * self.temperature))
            return kd
        except OverflowError:
            # Handle very large positive ΔG values
            return float('inf')
    
    def _calculate_ic50(self, kd: float) -> float:
        """
        Calculate IC50 from Kd.
        
        Uses the Cheng-Prusoff equation:
        IC50 = Kd * (1 + [L]/Km)
        
        For competitive inhibition with [L] << Km:
        IC50 ≈ Kd * (1 + [L]/Km) ≈ 2 * Kd
        """
        # Simplified assumption: IC50 ≈ 2 * Kd for competitive inhibition
        # This can be refined with more detailed protein-ligand information
        ic50 = 2.0 * kd
        
        # Alternative more detailed calculation if protein concentration known
        if self.protein_concentration > 0:
            # More accurate relationship considering protein concentration
            ic50 = kd * (1 + self.ligand_concentration / kd)
        
        return ic50
    
    def _calculate_ki(self, ic50: float) -> float:
        """
        Calculate inhibition constant (Ki) from IC50.
        
        For competitive inhibition:
        Ki = IC50 / (1 + [S]/Km)
        
        Simplified: Ki ≈ IC50 / 2 for many cases
        """
        # Simplified relationship
        ki = ic50 / 2.0
        return ki
    
    def _calculate_efficiency_indices(self, delta_g: float, 
                                    molecular_weight: Optional[float] = None,
                                    heavy_atom_count: Optional[int] = None) -> Dict[str, float]:
        """
        Calculate ligand efficiency indices.
        
        Returns:
            Dictionary with efficiency metrics:
            - LE: Ligand Efficiency = ΔG / heavy_atom_count
            - SILE: Size-Independent Ligand Efficiency
            - BEI: Binding Efficiency Index
        """
        efficiency_indices = {}
        
        # Ligand Efficiency (LE)
        if heavy_atom_count and heavy_atom_count > 0:
            le = -delta_g / heavy_atom_count  # Negative delta_g for favorable binding
            efficiency_indices['LE'] = le
        else:
            efficiency_indices['LE'] = 0.0
        
        # Size-Independent Ligand Efficiency (SILE)
        if molecular_weight and molecular_weight > 0:
            sile = -delta_g / (molecular_weight ** 0.3)
            efficiency_indices['SILE'] = sile
        else:
            efficiency_indices['SILE'] = 0.0
        
        # Binding Efficiency Index (BEI)
        if molecular_weight and molecular_weight > 0:
            bei = -delta_g / (molecular_weight / 1000.0)  # Normalized by kDa
            efficiency_indices['BEI'] = bei
        else:
            efficiency_indices['BEI'] = 0.0
        
        return efficiency_indices
    
    def _assess_confidence(self, docking_score: float, pose_rank: int) -> str:
        """
        Assess confidence in binding affinity prediction.
        
        Based on docking score quality and pose rank.
        """
        # High confidence: Strong binding score and top-ranked pose
        if docking_score <= -10.0 and pose_rank <= 3:
            return "High"
        
        # Medium confidence: Moderate binding score or good rank
        elif docking_score <= -7.0 and pose_rank <= 5:
            return "Medium"
        
        # Low confidence: Weak binding or low-ranked pose
        else:
            return "Low"
    
    def calculate_batch_affinities(self, results: List[Dict[str, Any]]) -> List[BindingAffinityResult]:
        """
        Calculate binding affinities for multiple poses.
        
        Args:
            results: List of docking results with scores and ranks
            
        Returns:
            List of BindingAffinityResult objects
        """
        affinity_results = []
        
        for result in results:
            score = result.get('score', 0.0)
            rank = result.get('rank', 1)
            mw = result.get('molecular_weight')
            heavy_atoms = result.get('heavy_atom_count')
            
            affinity = self.calculate_binding_affinity(
                docking_score=score,
                pose_rank=rank,
                molecular_weight=mw,
                heavy_atom_count=heavy_atoms
            )
            
            affinity_results.append(affinity)
        
        return affinity_results
    
    def generate_affinity_report(self, affinity_results: List[BindingAffinityResult]) -> str:
        """
        Generate a formatted binding affinity report.
        
        Args:
            affinity_results: List of binding affinity results
            
        Returns:
            Formatted text report
        """
        report_lines = []
        
        # Header
        report_lines.append("=" * 80)
        report_lines.append("    Binding Affinity Report")
        report_lines.append("=" * 80)
        report_lines.append("")
        
        # Summary statistics
        if affinity_results:
            scores = [r.docking_score for r in affinity_results]
            kd_values = [r.kd for r in affinity_results if r.kd != float('inf')]
            ic50_values = [r.ic50 for r in affinity_results if r.ic50 != float('inf')]
            
            report_lines.append("SUMMARY STATISTICS")
            report_lines.append("-" * 18)
            report_lines.append(f"Number of poses: {len(affinity_results)}")
            report_lines.append(f"Best docking score: {min(scores):.2f} kcal/mol")
            report_lines.append(f"Average docking score: {np.mean(scores):.2f} kcal/mol")
            
            if kd_values:
                report_lines.append(f"Best Kd: {min(kd_values):.2e} M")
                report_lines.append(f"Average Kd: {np.mean(kd_values):.2e} M")
            
            if ic50_values:
                report_lines.append(f"Best IC50: {min(ic50_values):.2e} M")
                report_lines.append(f"Average IC50: {np.mean(ic50_values):.2e} M")
            
            report_lines.append("")
        
        # Detailed results table
        report_lines.append("DETAILED BINDING AFFINITY RESULTS")
        report_lines.append("-" * 35)
        report_lines.append("")
        report_lines.append("Pose    ΔG (kcal/mol)    Kd (M)        IC50 (M)      Ki (M)        LE     Confidence")
        report_lines.append("-" * 85)
        
        for result in affinity_results:
            kd_str = f"{result.kd:.2e}" if result.kd != float('inf') else "N/A"
            ic50_str = f"{result.ic50:.2e}" if result.ic50 != float('inf') else "N/A"
            ki_str = f"{result.ki:.2e}" if result.ki != float('inf') else "N/A"
            le_str = f"{result.efficiency_indices.get('LE', 0.0):.3f}"
            
            report_lines.append(
                f"{result.pose_rank:4d}    {result.delta_g:8.2f}      {kd_str:>9s}     "
                f"{ic50_str:>9s}     {ki_str:>9s}     {le_str:>5s}    {result.confidence}"
            )
        
        report_lines.append("")
        report_lines.append("=" * 80)
        report_lines.append("")
        report_lines.append("NOTES:")
        report_lines.append("- ΔG: Free energy change (more negative = stronger binding)")
        report_lines.append("- Kd: Dissociation constant (lower = stronger binding)")
        report_lines.append("- IC50: Half-maximal inhibitory concentration")
        report_lines.append("- Ki: Inhibition constant")
        report_lines.append("- LE: Ligand efficiency (binding energy per heavy atom)")
        report_lines.append("- Confidence: Based on score quality and pose rank")
        report_lines.append("")
        
        return "\n".join(report_lines)


def calculate_binding_affinities(docking_results: List[Dict[str, Any]], 
                               temperature: float = T_STANDARD) -> List[BindingAffinityResult]:
    """
    Convenience function to calculate binding affinities from docking results.
    
    Args:
        docking_results: List of docking result dictionaries
        temperature: Temperature in Kelvin
        
    Returns:
        List of BindingAffinityResult objects
    """
    calculator = BindingAffinityCalculator(temperature=temperature)
    return calculator.calculate_batch_affinities(docking_results)


def export_affinity_results_csv(affinity_results: List[BindingAffinityResult], 
                               filename: str) -> None:
    """
    Export binding affinity results to CSV file.
    
    Args:
        affinity_results: List of binding affinity results
        filename: Output CSV filename
    """
    import csv
    
    if not affinity_results:
        return
    
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = ['Pose', 'Docking_Score_kcal_mol', 'Delta_G_kcal_mol', 
                     'Kd_M', 'IC50_M', 'Ki_M', 'Ligand_Efficiency', 
                     'SILE', 'BEI', 'Confidence']
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in affinity_results:
            row = {
                'Pose': result.pose_rank,
                'Docking_Score_kcal_mol': f"{result.docking_score:.3f}",
                'Delta_G_kcal_mol': f"{result.delta_g:.3f}",
                'Kd_M': f"{result.kd:.2e}" if result.kd != float('inf') else "N/A",
                'IC50_M': f"{result.ic50:.2e}" if result.ic50 != float('inf') else "N/A",
                'Ki_M': f"{result.ki:.2e}" if result.ki != float('inf') else "N/A",
                'Ligand_Efficiency': f"{result.efficiency_indices.get('LE', 0.0):.3f}",
                'SILE': f"{result.efficiency_indices.get('SILE', 0.0):.3f}",
                'BEI': f"{result.efficiency_indices.get('BEI', 0.0):.3f}",
                'Confidence': result.confidence
            }
            writer.writerow(row)