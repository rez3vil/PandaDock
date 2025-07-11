"""
IC50 and binding affinity calculator
"""

import numpy as np
from typing import Dict, Any, Optional, List, Tuple
import logging


class IC50Calculator:
    """
    Calculator for IC50, binding affinity, and related thermodynamic properties
    
    Features:
    - ΔG to IC50 conversion
    - Ligand efficiency calculation
    - Binding kinetics estimation
    - Thermodynamic property calculation
    - Drug-like property assessment
    """
    
    def __init__(self, temperature: float = 298.15):
        self.logger = logging.getLogger(__name__)
        
        # Thermodynamic constants
        self.temperature = temperature  # Kelvin
        self.R = 8.314  # J/(mol·K) - Gas constant
        self.R_kcal = 1.987e-3  # kcal/(mol·K) - Gas constant
        self.kT = self.R * self.temperature  # J/mol
        self.kT_kcal = self.R_kcal * self.temperature  # kcal/mol
        
        # Standard reference concentrations
        self.standard_conc = 1.0  # M (for standard state)
        self.protein_conc = 1e-9  # M (typical protein concentration)
        
        # Conversion factors
        self.kcal_to_joule = 4184  # J/kcal
        self.joule_to_kcal = 1 / 4184  # kcal/J
        
        self.logger.info(f"Initialized IC50Calculator at T={temperature} K")
    
    def delta_g_to_ic50(self, delta_g: float, protein_conc: float = None) -> float:
        """
        Convert binding free energy to IC50
        
        Args:
            delta_g: Binding free energy in kcal/mol
            protein_conc: Protein concentration in M (optional)
            
        Returns:
            IC50 in M
        """
        if protein_conc is None:
            protein_conc = self.protein_conc
        
        # ΔG = -RT ln(Ka)
        # Ka = 1/Kd
        # For competitive inhibition: IC50 ≈ Kd (assuming Ki ≈ Kd)
        
        if delta_g >= 0:
            return float('inf')  # No binding
        
        # Calculate dissociation constant
        kd = np.exp(delta_g / self.kT_kcal)
        
        # IC50 approximation for competitive inhibition
        ic50 = kd * (1 + protein_conc / kd)
        
        return ic50
    
    def ic50_to_delta_g(self, ic50: float, protein_conc: float = None) -> float:
        """
        Convert IC50 to binding free energy
        
        Args:
            ic50: IC50 in M
            protein_conc: Protein concentration in M (optional)
            
        Returns:
            Binding free energy in kcal/mol
        """
        if protein_conc is None:
            protein_conc = self.protein_conc
        
        if ic50 <= 0:
            return float('-inf')  # Infinite binding
        
        # For competitive inhibition: IC50 ≈ Kd
        kd = ic50 / (1 + protein_conc / ic50)
        
        # ΔG = -RT ln(1/Kd) = RT ln(Kd)
        delta_g = -self.kT_kcal * np.log(1 / kd)
        
        return delta_g
    
    def calculate_ligand_efficiency(self, delta_g: float, num_heavy_atoms: int) -> float:
        """
        Calculate ligand efficiency (LE)
        
        Args:
            delta_g: Binding free energy in kcal/mol
            num_heavy_atoms: Number of heavy atoms in ligand
            
        Returns:
            Ligand efficiency in kcal/mol per heavy atom
        """
        if num_heavy_atoms <= 0:
            return 0.0
        
        le = delta_g / num_heavy_atoms
        return le
    
    def calculate_lipophilic_efficiency(self, delta_g: float, logp: float) -> float:
        """
        Calculate lipophilic efficiency (LipE)
        
        Args:
            delta_g: Binding free energy in kcal/mol
            logp: LogP value
            
        Returns:
            Lipophilic efficiency
        """
        # Convert ΔG to pIC50 (assuming IC50 ≈ Kd)
        ic50 = self.delta_g_to_ic50(delta_g)
        
        if ic50 <= 0:
            return float('inf')
        
        pic50 = -np.log10(ic50)
        lipe = pic50 - logp
        
        return lipe
    
    def calculate_size_independent_ligand_efficiency(self, delta_g: float, num_heavy_atoms: int) -> float:
        """
        Calculate size-independent ligand efficiency (SILE)
        
        Args:
            delta_g: Binding free energy in kcal/mol
            num_heavy_atoms: Number of heavy atoms in ligand
            
        Returns:
            Size-independent ligand efficiency
        """
        if num_heavy_atoms <= 0:
            return 0.0
        
        # SILE = ΔG / (N^0.3) where N is number of heavy atoms
        sile = delta_g / (num_heavy_atoms ** 0.3)
        
        return sile
    
    def calculate_fit_quality(self, delta_g: float, num_heavy_atoms: int) -> float:
        """
        Calculate fit quality (FQ)
        
        Args:
            delta_g: Binding free energy in kcal/mol
            num_heavy_atoms: Number of heavy atoms in ligand
            
        Returns:
            Fit quality score
        """
        if num_heavy_atoms <= 0:
            return 0.0
        
        # FQ = LE / LE_scale where LE_scale = 1.37 - 0.116 * log(N)
        le = self.calculate_ligand_efficiency(delta_g, num_heavy_atoms)
        le_scale = 1.37 - 0.116 * np.log(num_heavy_atoms)
        
        fq = le / le_scale if le_scale > 0 else 0.0
        
        return fq
    
    def calculate_binding_kinetics(self, delta_g: float, delta_h: float = None, 
                                  delta_s: float = None) -> Dict[str, float]:
        """
        Calculate binding kinetics parameters
        
        Args:
            delta_g: Binding free energy in kcal/mol
            delta_h: Enthalpy change in kcal/mol (optional)
            delta_s: Entropy change in kcal/(mol·K) (optional)
            
        Returns:
            Dictionary of kinetic parameters
        """
        kinetics = {}
        
        # Dissociation constant
        if delta_g < 0:
            kd = np.exp(-delta_g / self.kT_kcal)
            kinetics['kd'] = kd
            kinetics['ka'] = 1 / kd
        else:
            kinetics['kd'] = float('inf')
            kinetics['ka'] = 0.0
        
        # Estimate association and dissociation rate constants
        # These are rough estimates based on typical protein-ligand interactions
        
        if delta_g < 0:
            # Typical association rate constant (M⁻¹s⁻¹)
            k_on_typical = 1e6  # M⁻¹s⁻¹
            
            # Estimate based on binding affinity
            # Stronger binding often correlates with slower dissociation
            k_off = k_on_typical * kd
            
            kinetics['k_on'] = k_on_typical
            kinetics['k_off'] = k_off
            kinetics['residence_time'] = 1 / k_off if k_off > 0 else float('inf')
        else:
            kinetics['k_on'] = 0.0
            kinetics['k_off'] = float('inf')
            kinetics['residence_time'] = 0.0
        
        # Thermodynamic parameters
        if delta_h is not None:
            kinetics['delta_h'] = delta_h
            
            if delta_s is not None:
                kinetics['delta_s'] = delta_s
            else:
                # Calculate entropy from G and H
                kinetics['delta_s'] = (delta_h - delta_g) / self.temperature
        
        elif delta_s is not None:
            kinetics['delta_s'] = delta_s
            kinetics['delta_h'] = delta_g + self.temperature * delta_s
        
        return kinetics
    
    def calculate_cooperativity(self, ic50_values: List[float], hill_coefficient: float = 1.0) -> Dict[str, float]:
        """
        Calculate cooperativity parameters
        
        Args:
            ic50_values: List of IC50 values for multiple binding sites
            hill_coefficient: Hill coefficient (default: 1.0)
            
        Returns:
            Dictionary of cooperativity parameters
        """
        cooperativity = {}
        
        if len(ic50_values) < 2:
            return cooperativity
        
        # Calculate individual binding constants
        kd_values = [ic50 for ic50 in ic50_values if ic50 > 0]
        
        if len(kd_values) < 2:
            return cooperativity
        
        # Geometric mean of binding constants
        geometric_mean_kd = np.power(np.prod(kd_values), 1/len(kd_values))
        
        # Cooperativity index
        cooperativity['geometric_mean_kd'] = geometric_mean_kd
        cooperativity['hill_coefficient'] = hill_coefficient
        
        # Positive cooperativity: hill_coefficient > 1
        # Negative cooperativity: hill_coefficient < 1
        if hill_coefficient > 1:
            cooperativity['cooperativity_type'] = 'positive'
        elif hill_coefficient < 1:
            cooperativity['cooperativity_type'] = 'negative'
        else:
            cooperativity['cooperativity_type'] = 'independent'
        
        return cooperativity
    
    def calculate_selectivity_index(self, ic50_target: float, ic50_offtarget: float) -> float:
        """
        Calculate selectivity index
        
        Args:
            ic50_target: IC50 for target protein
            ic50_offtarget: IC50 for off-target protein
            
        Returns:
            Selectivity index (ratio of IC50 values)
        """
        if ic50_target <= 0:
            return 0.0
        
        if ic50_offtarget <= 0:
            return float('inf')
        
        selectivity = ic50_offtarget / ic50_target
        return selectivity
    
    def assess_drug_likeness(self, delta_g: float, molecular_weight: float, 
                           logp: float, hbd: int, hba: int, num_heavy_atoms: int) -> Dict[str, Any]:
        """
        Assess drug-likeness based on binding affinity and molecular properties
        
        Args:
            delta_g: Binding free energy in kcal/mol
            molecular_weight: Molecular weight in Da
            logp: LogP value
            hbd: Number of hydrogen bond donors
            hba: Number of hydrogen bond acceptors
            num_heavy_atoms: Number of heavy atoms
            
        Returns:
            Dictionary of drug-likeness assessment
        """
        assessment = {}
        
        # Calculate derived properties
        ic50 = self.delta_g_to_ic50(delta_g)
        le = self.calculate_ligand_efficiency(delta_g, num_heavy_atoms)
        lipe = self.calculate_lipophilic_efficiency(delta_g, logp)
        
        assessment['ic50'] = ic50 * 1e9  # Convert from M to nM
        assessment['delta_g'] = delta_g
        assessment['ligand_efficiency'] = le
        assessment['lipophilic_efficiency'] = lipe
        
        # Lipinski's Rule of Five
        lipinski_violations = 0
        
        if molecular_weight > 500:
            lipinski_violations += 1
        if logp > 5:
            lipinski_violations += 1
        if hbd > 5:
            lipinski_violations += 1
        if hba > 10:
            lipinski_violations += 1
        
        assessment['lipinski_violations'] = lipinski_violations
        assessment['lipinski_compliant'] = lipinski_violations <= 1
        
        # Veber's rules
        veber_compliant = True
        if num_heavy_atoms > 10:  # Simplified rotatable bond count
            veber_compliant = False
        
        # Estimate TPSA (topological polar surface area)
        tpsa_estimate = (hbd + hba) * 20  # Simplified estimate
        if tpsa_estimate > 140:
            veber_compliant = False
        
        assessment['veber_compliant'] = veber_compliant
        assessment['tpsa_estimate'] = tpsa_estimate
        
        # Ligand efficiency assessment
        if le < -0.3:
            assessment['le_assessment'] = 'excellent'
        elif le < -0.24:
            assessment['le_assessment'] = 'good'
        elif le < -0.15:
            assessment['le_assessment'] = 'moderate'
        else:
            assessment['le_assessment'] = 'poor'
        
        # Lipophilic efficiency assessment
        if lipe > 5:
            assessment['lipe_assessment'] = 'excellent'
        elif lipe > 3:
            assessment['lipe_assessment'] = 'good'
        elif lipe > 1:
            assessment['lipe_assessment'] = 'moderate'
        else:
            assessment['lipe_assessment'] = 'poor'
        
        # Overall drug-likeness score
        score = 0
        
        if assessment['lipinski_compliant']:
            score += 2
        if assessment['veber_compliant']:
            score += 2
        if assessment['le_assessment'] in ['excellent', 'good']:
            score += 2
        if assessment['lipe_assessment'] in ['excellent', 'good']:
            score += 2
        if ic50 < 1e-6:  # nM range
            score += 2
        
        assessment['drug_likeness_score'] = score
        assessment['max_score'] = 10
        
        # Drug-likeness category
        if score >= 8:
            assessment['drug_likeness_category'] = 'excellent'
        elif score >= 6:
            assessment['drug_likeness_category'] = 'good'
        elif score >= 4:
            assessment['drug_likeness_category'] = 'moderate'
        else:
            assessment['drug_likeness_category'] = 'poor'
        
        return assessment
    
    def calculate_dose_response_curve(self, ic50: float, hill_coefficient: float = 1.0, 
                                     max_response: float = 100.0, 
                                     min_response: float = 0.0) -> Dict[str, Any]:
        """
        Calculate dose-response curve parameters
        
        Args:
            ic50: IC50 value in M
            hill_coefficient: Hill coefficient
            max_response: Maximum response (%)
            min_response: Minimum response (%)
            
        Returns:
            Dictionary of dose-response parameters
        """
        curve_params = {
            'ic50': ic50,
            'hill_coefficient': hill_coefficient,
            'max_response': max_response,
            'min_response': min_response
        }
        
        # Calculate concentrations for dose-response curve
        concentrations = np.logspace(-12, -3, 100)  # 1 pM to 1 mM
        
        # Hill equation: Response = min + (max - min) / (1 + (IC50/[L])^n)
        responses = []
        for conc in concentrations:
            if conc > 0:
                response = min_response + (max_response - min_response) / (
                    1 + (ic50 / conc) ** hill_coefficient
                )
            else:
                response = min_response
            responses.append(response)
        
        curve_params['concentrations'] = concentrations
        curve_params['responses'] = np.array(responses)
        
        # Calculate EC50 (concentration for 50% response)
        ec50 = ic50 * ((max_response - min_response) / (50 - min_response) - 1) ** (1 / hill_coefficient)
        curve_params['ec50'] = ec50
        
        return curve_params
    
    def convert_units(self, value: float, from_unit: str, to_unit: str) -> float:
        """
        Convert between different concentration and energy units
        
        Args:
            value: Value to convert
            from_unit: Original unit
            to_unit: Target unit
            
        Returns:
            Converted value
        """
        # Define conversion factors
        concentration_conversions = {
            'M': 1.0,
            'mM': 1e-3,
            'µM': 1e-6,
            'nM': 1e-9,
            'pM': 1e-12
        }
        
        energy_conversions = {
            'kcal/mol': 1.0,
            'kJ/mol': self.kcal_to_joule,
            'J/mol': self.kcal_to_joule,
            'eV': 0.0433634  # kcal/mol to eV
        }
        
        # Convert concentrations
        if from_unit in concentration_conversions and to_unit in concentration_conversions:
            return value * concentration_conversions[from_unit] / concentration_conversions[to_unit]
        
        # Convert energies
        elif from_unit in energy_conversions and to_unit in energy_conversions:
            return value * energy_conversions[from_unit] / energy_conversions[to_unit]
        
        else:
            raise ValueError(f"Cannot convert from {from_unit} to {to_unit}")
    
    def get_thermodynamic_summary(self, delta_g: float, num_heavy_atoms: int, 
                                 molecular_weight: float, logp: float) -> Dict[str, Any]:
        """
        Get comprehensive thermodynamic summary
        
        Args:
            delta_g: Binding free energy in kcal/mol
            num_heavy_atoms: Number of heavy atoms
            molecular_weight: Molecular weight in Da
            logp: LogP value
            
        Returns:
            Comprehensive thermodynamic summary
        """
        summary = {
            'thermodynamics': {
                'delta_g': delta_g,
                'temperature': self.temperature,
                'ic50': self.delta_g_to_ic50(delta_g),
                'kd': np.exp(-delta_g / self.kT_kcal) if delta_g < 0 else float('inf'),
                'ka': np.exp(delta_g / self.kT_kcal) if delta_g < 0 else 0.0
            },
            'efficiency_metrics': {
                'ligand_efficiency': self.calculate_ligand_efficiency(delta_g, num_heavy_atoms),
                'lipophilic_efficiency': self.calculate_lipophilic_efficiency(delta_g, logp),
                'size_independent_le': self.calculate_size_independent_ligand_efficiency(delta_g, num_heavy_atoms),
                'fit_quality': self.calculate_fit_quality(delta_g, num_heavy_atoms)
            },
            'kinetics': self.calculate_binding_kinetics(delta_g),
            'molecular_properties': {
                'num_heavy_atoms': num_heavy_atoms,
                'molecular_weight': molecular_weight,
                'logp': logp
            }
        }
        
        return summary