#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metal Docking Examples

This module provides comprehensive examples of metal docking using PandaDock's
advanced metalloprotein support. Includes examples for zinc fingers, 
iron-containing enzymes, calcium-binding proteins, and more.
"""

import numpy as np
import os
import sys
from typing import List, Dict, Any, Optional
import logging

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docking.metal_docking_engine import MetalDockingEngine, MetalDockingConfig, MetalPose
from docking.metal_coordination import MetalCenter, MetalType, CoordinationGeometry
from scoring.metal_scoring import MetalScoringFunction, MetalScoringParameters
from utils.metal_constraints import MetalConstraintManager, ConstraintSetPresets
from utils.ic50_calculator import IC50Calculator

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MetalDockingExamples:
    """Collection of metal docking examples"""
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.ic50_calc = IC50Calculator()
        
    def example_1_zinc_finger_docking(self):
        """Example 1: Docking to a zinc finger protein"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 1: Zinc Finger Metal Docking")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates docking to a zinc finger protein.
        Zinc fingers typically have tetrahedral coordination with 
        2 cysteine and 2 histidine residues.
        """)
        
        # Create mock configuration (in real use, load from config file)
        config = self._create_mock_config()
        
        # Configure metal docking parameters
        metal_config = MetalDockingConfig(
            detect_metals_automatically=True,
            enforce_geometric_constraints=True,
            coordination_focused_sampling=True,
            use_metal_scoring=True,
            geometric_constraint_weight=2.0,
            metal_focused_exhaustiveness=25,
            require_metal_coordination=True,
            min_coordinating_atoms=1,
            filter_non_coordinating_poses=True
        )
        
        # Initialize metal docking engine
        engine = MetalDockingEngine(config, metal_config)
        
        # Create example zinc finger metal center
        zinc_center = self._create_zinc_finger_center()
        engine.metal_centers = [zinc_center]
        
        # Create example ligand (small molecule that could interact with zinc)
        ligand_data = self._create_zinc_binding_ligand()
        engine.ligand = ligand_data
        
        # Initialize metal scoring
        engine.metal_scorer = MetalScoringFunction([zinc_center])
        
        print(f"Zinc center at: {zinc_center.coordinates}")
        print(f"Coordination geometry: {zinc_center.geometry.value}")
        print(f"Coordination number: {zinc_center.coordination_number}")
        
        # Perform metal docking
        poses = self._mock_metal_docking(engine, "zinc_finger")
        
        # Analyze results
        self._analyze_metal_poses(poses, "Zinc Finger Docking")
        
        # Generate detailed report
        report = engine.get_metal_docking_report(poses)
        self._print_docking_report(report)
        
        return poses, report
    
    def example_2_iron_enzyme_docking(self):
        """Example 2: Docking to an iron-containing enzyme"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 2: Iron Enzyme Metal Docking")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates docking to an iron-containing enzyme.
        Iron centers often have octahedral coordination and can exist
        in multiple oxidation states (Fe2+/Fe3+).
        """)
        
        config = self._create_mock_config()
        
        # Configure for iron enzyme
        metal_config = MetalDockingConfig(
            detect_metals_automatically=True,
            enforce_geometric_constraints=True,
            coordination_focused_sampling=True,
            use_metal_scoring=True,
            geometric_constraint_weight=1.5,  # Slightly more flexible for Fe
            metal_focused_exhaustiveness=30,
            require_metal_coordination=True,
            min_coordinating_atoms=1,
            max_coordinating_atoms=6,  # Fe can be 6-coordinate
            filter_non_coordinating_poses=True
        )
        
        engine = MetalDockingEngine(config, metal_config)
        
        # Create iron center (octahedral, Fe3+)
        iron_center = self._create_iron_enzyme_center()
        engine.metal_centers = [iron_center]
        
        # Create substrate/inhibitor ligand
        ligand_data = self._create_iron_binding_ligand()
        engine.ligand = ligand_data
        
        # Initialize scoring with iron-specific parameters
        metal_params = MetalScoringParameters()
        metal_params.coordination_strength['O'] = -4.5  # Strong Fe-O interaction
        metal_params.coordination_strength['N'] = -4.0  # Strong Fe-N interaction
        engine.metal_scorer = MetalScoringFunction([iron_center], metal_params)
        
        print(f"Iron center at: {iron_center.coordinates}")
        print(f"Coordination geometry: {iron_center.geometry.value}")
        print(f"Oxidation state: +{iron_center.oxidation_state}")
        
        poses = self._mock_metal_docking(engine, "iron_enzyme")
        self._analyze_metal_poses(poses, "Iron Enzyme Docking")
        
        return poses
    
    def example_3_calcium_binding_protein(self):
        """Example 3: Docking to a calcium-binding protein"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 3: Calcium-Binding Protein Docking")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates docking to a calcium-binding protein.
        Calcium typically has 6-8 coordination with oxygen-rich ligands
        and more flexible geometry than transition metals.
        """)
        
        config = self._create_mock_config()
        
        # Configure for calcium binding (more flexible)
        metal_config = MetalDockingConfig(
            detect_metals_automatically=True,
            enforce_geometric_constraints=False,  # Ca is more flexible
            coordination_focused_sampling=True,
            use_metal_scoring=True,
            geometric_constraint_weight=0.5,  # Low weight for Ca
            metal_focused_exhaustiveness=20,
            require_metal_coordination=True,
            min_coordinating_atoms=2,
            max_coordinating_atoms=8,  # Ca can have high coordination
            filter_non_coordinating_poses=False,  # More permissive
            distance_tolerance=0.5,  # Larger tolerance for Ca
            angle_tolerance=20.0
        )
        
        engine = MetalDockingEngine(config, metal_config)
        
        # Create calcium center
        calcium_center = self._create_calcium_center()
        engine.metal_centers = [calcium_center]
        
        # Create calcium-chelating ligand
        ligand_data = self._create_calcium_chelator()
        engine.ligand = ligand_data
        
        # Calcium-specific scoring
        metal_params = MetalScoringParameters()
        metal_params.coordination_strength['O'] = -2.5  # Weaker than transition metals
        metal_params.geometric_penalty_weight = 0.5    # Less strict geometry
        metal_params.distance_tolerance = 0.4          # More flexible distances
        engine.metal_scorer = MetalScoringFunction([calcium_center], metal_params)
        
        print(f"Calcium center at: {calcium_center.coordinates}")
        print(f"Coordination geometry: {calcium_center.geometry.value}")
        print(f"Expected coordination number: {calcium_center.coordination_number}")
        
        poses = self._mock_metal_docking(engine, "calcium_protein")
        self._analyze_metal_poses(poses, "Calcium-Binding Protein Docking")
        
        return poses
    
    def example_4_multi_metal_cluster(self):
        """Example 4: Docking to a multi-metal cluster"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 4: Multi-Metal Cluster Docking")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates docking to a protein with multiple
        metal centers, such as a [2Fe-2S] cluster or binuclear enzyme.
        """)
        
        config = self._create_mock_config()
        
        # Configure for multi-metal system
        metal_config = MetalDockingConfig(
            detect_metals_automatically=True,
            enforce_geometric_constraints=True,
            coordination_focused_sampling=True,
            use_metal_scoring=True,
            geometric_constraint_weight=2.0,
            metal_focused_exhaustiveness=35,  # More sampling for complex system
            require_metal_coordination=True,
            min_coordinating_atoms=1,
            filter_non_coordinating_poses=True
        )
        
        engine = MetalDockingEngine(config, metal_config)
        
        # Create binuclear zinc cluster
        metal_centers = self._create_binuclear_zinc_cluster()
        engine.metal_centers = metal_centers
        
        # Create bridging ligand
        ligand_data = self._create_bridging_ligand()
        engine.ligand = ligand_data
        
        # Multi-metal scoring
        engine.metal_scorer = MetalScoringFunction(metal_centers)
        
        print(f"Multi-metal cluster with {len(metal_centers)} centers:")
        for i, center in enumerate(metal_centers):
            print(f"  Metal {i+1}: {center.metal_type.value} at {center.coordinates}")
        
        poses = self._mock_metal_docking(engine, "multi_metal")
        self._analyze_metal_poses(poses, "Multi-Metal Cluster Docking")
        
        # Analyze inter-metal distances in poses
        self._analyze_inter_metal_interactions(poses, metal_centers)
        
        return poses
    
    def example_5_constraint_optimization(self):
        """Example 5: Constraint-driven metal docking optimization"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 5: Constraint-Driven Optimization")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates how to use explicit geometric
        constraints to optimize metal coordination during docking.
        """)
        
        # Create zinc finger example
        config = self._create_mock_config()
        metal_config = MetalDockingConfig()
        engine = MetalDockingEngine(config, metal_config)
        
        zinc_center = self._create_zinc_finger_center()
        ligand_data = self._create_zinc_binding_ligand()
        
        print("Testing different constraint presets...")
        
        # Test different constraint presets
        presets = ["standard", "strict", "flexible", "distance_only"]
        results = {}
        
        for preset in presets:
            print(f"\n--- {preset.upper()} Constraints ---")
            
            # Apply constraints
            from utils.metal_constraints import apply_metal_constraints_to_pose
            
            initial_coords = ligand_data['coordinates']
            atom_types = ligand_data['atom_types']
            
            optimized_coords, constraint_results = apply_metal_constraints_to_pose(
                initial_coords, atom_types, [zinc_center], preset
            )
            
            print(f"Constraint satisfaction: {constraint_results['constraint_satisfaction']:.3f}")
            print(f"Total violations: {constraint_results['total_violations']}")
            print(f"Final penalty: {constraint_results['final_penalty']:.3f}")
            
            # Show coordinate changes
            coord_change = np.linalg.norm(optimized_coords - initial_coords)
            print(f"Coordinate RMSD: {coord_change:.3f} Å")
            
            results[preset] = constraint_results
        
        return results
    
    def example_6_drug_design_workflow(self):
        """Example 6: Complete drug design workflow for metalloprotein"""
        self.logger.info("="*60)
        self.logger.info("EXAMPLE 6: Metalloprotein Drug Design Workflow")
        self.logger.info("="*60)
        
        print("""
        This example demonstrates a complete drug design workflow
        for targeting a metalloprotein, including lead optimization
        and affinity prediction.
        """)
        
        # Target: Matrix metalloproteinase (MMP) with zinc
        config = self._create_mock_config()
        metal_config = MetalDockingConfig(
            coordination_focused_sampling=True,
            use_metal_scoring=True,
            require_metal_coordination=True,
            min_coordinating_atoms=1
        )
        
        engine = MetalDockingEngine(config, metal_config)
        
        # Create MMP-like zinc center
        mmp_zinc = self._create_mmp_zinc_center()
        engine.metal_centers = [mmp_zinc]
        engine.metal_scorer = MetalScoringFunction([mmp_zinc])
        
        print("Target: Matrix Metalloproteinase (MMP)")
        print(f"Zinc center: {mmp_zinc.coordinates}")
        print(f"Geometry: {mmp_zinc.geometry.value}")
        
        # Test multiple drug candidates
        candidates = self._create_drug_candidates()
        
        results = []
        for i, candidate in enumerate(candidates):
            print(f"\n--- Testing Candidate {i+1}: {candidate['name']} ---")
            
            engine.ligand = candidate
            poses = self._mock_metal_docking(engine, f"candidate_{i+1}")
            
            if poses:
                best_pose = poses[0]
                
                # Calculate drug-like properties
                analysis = self._analyze_drug_candidate(best_pose, candidate, mmp_zinc)
                results.append(analysis)
                
                print(f"Docking score: {best_pose.score:.3f}")
                print(f"Metal coordination: {len(best_pose.coordinating_atoms)} atoms")
                print(f"Predicted IC50: {analysis['predicted_ic50']:.1f} nM")
                print(f"Ligand efficiency: {analysis['ligand_efficiency']:.3f}")
                print(f"Drug-likeness: {analysis['drug_likeness']}")
        
        # Rank candidates
        self._rank_drug_candidates(results)
        
        return results
    
    def _create_mock_config(self):
        """Create mock configuration object"""
        class MockConfig:
            def __init__(self):
                self.docking = self._create_docking_config()
                self.io = self._create_io_config()
            
            def _create_docking_config(self):
                class DockingConfig:
                    num_poses = 10
                    exhaustiveness = 8
                    energy_range = 3.0
                return DockingConfig()
            
            def _create_io_config(self):
                class IOConfig:
                    center_x = 25.0
                    center_y = 30.0
                    center_z = 15.0
                    size_x = 20.0
                    size_y = 20.0
                    size_z = 20.0
                return IOConfig()
            
            def to_dict(self):
                return {'mock': 'config'}
        
        return MockConfig()
    
    def _create_zinc_finger_center(self) -> MetalCenter:
        """Create example zinc finger metal center"""
        # Typical zinc finger coordination
        coordinating_atoms = [
            {'coordinates': np.array([25.2, 29.8, 15.1]), 'element': 'S', 'residue': 'CYS'},
            {'coordinates': np.array([24.8, 30.2, 14.9]), 'element': 'S', 'residue': 'CYS'},
            {'coordinates': np.array([25.3, 30.1, 15.3]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([24.7, 29.9, 14.7]), 'element': 'N', 'residue': 'HIS'}
        ]
        
        return MetalCenter(
            metal_type=MetalType.ZN,
            coordinates=np.array([25.0, 30.0, 15.0]),
            coordination_number=4,
            geometry=CoordinationGeometry.TETRAHEDRAL,
            coordinating_atoms=coordinating_atoms,
            oxidation_state=2,
            charge=2.0
        )
    
    def _create_iron_enzyme_center(self) -> MetalCenter:
        """Create example iron enzyme center"""
        coordinating_atoms = [
            {'coordinates': np.array([25.1, 29.9, 15.2]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([24.9, 30.1, 14.8]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([25.2, 30.2, 15.1]), 'element': 'O', 'residue': 'GLU'},
            {'coordinates': np.array([24.8, 29.8, 14.9]), 'element': 'O', 'residue': 'ASP'},
            {'coordinates': np.array([25.0, 30.0, 16.0]), 'element': 'O', 'residue': 'H2O'},
            {'coordinates': np.array([25.0, 30.0, 14.0]), 'element': 'O', 'residue': 'H2O'}
        ]
        
        return MetalCenter(
            metal_type=MetalType.FE,
            coordinates=np.array([25.0, 30.0, 15.0]),
            coordination_number=6,
            geometry=CoordinationGeometry.OCTAHEDRAL,
            coordinating_atoms=coordinating_atoms,
            oxidation_state=3,
            charge=3.0
        )
    
    def _create_calcium_center(self) -> MetalCenter:
        """Create example calcium center"""
        coordinating_atoms = [
            {'coordinates': np.array([25.1, 29.8, 15.3]), 'element': 'O', 'residue': 'ASP'},
            {'coordinates': np.array([24.7, 30.3, 14.8]), 'element': 'O', 'residue': 'ASP'},
            {'coordinates': np.array([25.4, 30.1, 15.2]), 'element': 'O', 'residue': 'GLU'},
            {'coordinates': np.array([24.6, 29.7, 14.9]), 'element': 'O', 'residue': 'GLU'},
            {'coordinates': np.array([25.2, 30.4, 15.5]), 'element': 'O', 'residue': 'H2O'},
            {'coordinates': np.array([24.8, 29.6, 14.5]), 'element': 'O', 'residue': 'H2O'},
            {'coordinates': np.array([25.5, 30.2, 14.7]), 'element': 'O', 'residue': 'H2O'}
        ]
        
        return MetalCenter(
            metal_type=MetalType.CA,
            coordinates=np.array([25.0, 30.0, 15.0]),
            coordination_number=7,
            geometry=CoordinationGeometry.OCTAHEDRAL,  # Simplified
            coordinating_atoms=coordinating_atoms,
            oxidation_state=2,
            charge=2.0
        )
    
    def _create_binuclear_zinc_cluster(self) -> List[MetalCenter]:
        """Create binuclear zinc cluster"""
        # Zinc 1
        zn1_coords = [
            {'coordinates': np.array([24.5, 29.8, 15.1]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([24.3, 30.2, 14.9]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([24.7, 30.0, 15.3]), 'element': 'S', 'residue': 'CYS'},
            {'coordinates': np.array([24.2, 29.9, 14.7]), 'element': 'O', 'residue': 'bridging'}
        ]
        
        zn1 = MetalCenter(
            metal_type=MetalType.ZN,
            coordinates=np.array([24.5, 30.0, 15.0]),
            coordination_number=4,
            geometry=CoordinationGeometry.TETRAHEDRAL,
            coordinating_atoms=zn1_coords,
            oxidation_state=2,
            charge=2.0
        )
        
        # Zinc 2
        zn2_coords = [
            {'coordinates': np.array([25.5, 30.2, 14.9]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([25.7, 29.8, 15.1]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([25.3, 30.0, 14.7]), 'element': 'S', 'residue': 'CYS'},
            {'coordinates': np.array([25.8, 30.1, 15.3]), 'element': 'O', 'residue': 'bridging'}
        ]
        
        zn2 = MetalCenter(
            metal_type=MetalType.ZN,
            coordinates=np.array([25.5, 30.0, 15.0]),
            coordination_number=4,
            geometry=CoordinationGeometry.TETRAHEDRAL,
            coordinating_atoms=zn2_coords,
            oxidation_state=2,
            charge=2.0
        )
        
        return [zn1, zn2]
    
    def _create_mmp_zinc_center(self) -> MetalCenter:
        """Create MMP-like zinc center"""
        coordinating_atoms = [
            {'coordinates': np.array([25.1, 29.9, 15.2]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([24.9, 30.1, 14.8]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([25.2, 30.0, 15.1]), 'element': 'N', 'residue': 'HIS'},
            {'coordinates': np.array([25.0, 30.0, 14.0]), 'element': 'O', 'residue': 'H2O'}
        ]
        
        return MetalCenter(
            metal_type=MetalType.ZN,
            coordinates=np.array([25.0, 30.0, 15.0]),
            coordination_number=4,
            geometry=CoordinationGeometry.TETRAHEDRAL,
            coordinating_atoms=coordinating_atoms,
            oxidation_state=2,
            charge=2.0
        )
    
    def _create_zinc_binding_ligand(self) -> Dict[str, Any]:
        """Create example zinc-binding ligand"""
        # Small molecule with nitrogen and sulfur coordination
        coordinates = np.array([
            [26.0, 30.0, 15.0],  # N - potential coordinator
            [26.8, 30.2, 15.3],  # C
            [27.5, 29.8, 15.1],  # C
            [28.0, 30.0, 14.8],  # S - potential coordinator
            [26.5, 29.5, 14.7],  # C
            [25.8, 29.7, 14.9],  # C
        ])
        
        atom_types = ['N', 'C', 'C', 'S', 'C', 'C']
        bonds = [(0, 1, 'single'), (1, 2, 'single'), (2, 3, 'single'), 
                (1, 4, 'single'), (4, 5, 'single'), (5, 0, 'single')]
        
        return {
            'coordinates': coordinates,
            'atom_types': atom_types,
            'bonds': bonds,
            'name': 'zinc_chelator'
        }
    
    def _create_iron_binding_ligand(self) -> Dict[str, Any]:
        """Create example iron-binding ligand"""
        # Molecule with multiple oxygen coordinators
        coordinates = np.array([
            [26.0, 30.0, 15.0],  # O - potential coordinator
            [26.5, 30.3, 15.2],  # C
            [27.0, 29.8, 15.1],  # O - potential coordinator  
            [27.5, 30.1, 14.9],  # C
            [26.2, 29.6, 14.8],  # O - potential coordinator
            [25.7, 29.9, 15.1],  # C
        ])
        
        atom_types = ['O', 'C', 'O', 'C', 'O', 'C']
        bonds = [(0, 1, 'single'), (1, 2, 'single'), (2, 3, 'single'),
                (1, 4, 'single'), (4, 5, 'single')]
        
        return {
            'coordinates': coordinates,
            'atom_types': atom_types,
            'bonds': bonds,
            'name': 'iron_chelator'
        }
    
    def _create_calcium_chelator(self) -> Dict[str, Any]:
        """Create calcium-chelating ligand"""
        # EDTA-like structure with multiple oxygens
        coordinates = np.array([
            [26.0, 30.0, 15.0],  # O
            [26.3, 30.2, 15.3],  # C
            [26.6, 29.8, 15.1],  # O
            [27.0, 30.1, 14.8],  # C
            [26.8, 29.6, 15.2],  # O
            [25.5, 29.8, 14.9],  # C
            [25.2, 30.2, 15.1],  # O
            [25.8, 30.4, 15.4],  # C
        ])
        
        atom_types = ['O', 'C', 'O', 'C', 'O', 'C', 'O', 'C']
        bonds = [(0, 1, 'single'), (1, 2, 'single'), (2, 3, 'single'),
                (3, 4, 'single'), (4, 5, 'single'), (5, 6, 'single'),
                (6, 7, 'single'), (7, 0, 'single')]
        
        return {
            'coordinates': coordinates,
            'atom_types': atom_types,
            'bonds': bonds,
            'name': 'calcium_chelator'
        }
    
    def _create_bridging_ligand(self) -> Dict[str, Any]:
        """Create ligand that can bridge multiple metals"""
        coordinates = np.array([
            [24.8, 30.0, 15.0],  # O - coordinates to Zn1
            [25.0, 30.1, 15.1],  # C
            [25.2, 30.0, 15.0],  # O - coordinates to Zn2
            [25.1, 29.8, 14.8],  # C
            [24.9, 29.9, 14.9],  # N
        ])
        
        atom_types = ['O', 'C', 'O', 'C', 'N']
        bonds = [(0, 1, 'single'), (1, 2, 'single'), (2, 3, 'single'),
                (3, 4, 'single'), (4, 0, 'single')]
        
        return {
            'coordinates': coordinates,
            'atom_types': atom_types,
            'bonds': bonds,
            'name': 'bridging_ligand'
        }
    
    def _create_drug_candidates(self) -> List[Dict[str, Any]]:
        """Create multiple drug candidate molecules"""
        candidates = []
        
        # Candidate 1: Hydroxamic acid (common MMP inhibitor)
        candidate1 = {
            'coordinates': np.array([
                [26.0, 30.0, 15.0],  # O (hydroxamic)
                [26.2, 30.2, 15.3],  # N
                [26.5, 29.8, 15.1],  # O (hydroxamic)
                [27.0, 30.1, 14.8],  # C
                [26.8, 29.6, 15.2],  # C (aromatic)
                [25.5, 29.8, 14.9],  # C
            ]),
            'atom_types': ['O', 'N', 'O', 'C', 'C', 'C'],
            'bonds': [(0, 1, 'single'), (1, 2, 'single'), (1, 3, 'single'),
                     (3, 4, 'aromatic'), (4, 5, 'aromatic')],
            'name': 'hydroxamic_acid_inhibitor',
            'molecular_weight': 180.0,
            'heavy_atoms': 6
        }
        
        # Candidate 2: Carboxylate inhibitor
        candidate2 = {
            'coordinates': np.array([
                [26.0, 30.0, 15.0],  # O (carboxylate)
                [26.3, 30.2, 15.3],  # C
                [26.6, 29.8, 15.1],  # O (carboxylate)
                [27.0, 30.1, 14.8],  # C
                [26.8, 29.6, 15.2],  # N
                [25.5, 29.8, 14.9],  # C
            ]),
            'atom_types': ['O', 'C', 'O', 'C', 'N', 'C'],
            'bonds': [(0, 1, 'single'), (1, 2, 'double'), (1, 3, 'single'),
                     (3, 4, 'single'), (4, 5, 'single')],
            'name': 'carboxylate_inhibitor',
            'molecular_weight': 165.0,
            'heavy_atoms': 6
        }
        
        # Candidate 3: Phosphonate inhibitor
        candidate3 = {
            'coordinates': np.array([
                [26.0, 30.0, 15.0],  # P
                [26.3, 30.2, 15.3],  # O
                [26.6, 29.8, 15.1],  # O
                [27.0, 30.1, 14.8],  # O
                [26.8, 29.6, 15.2],  # C
                [25.5, 29.8, 14.9],  # C
            ]),
            'atom_types': ['P', 'O', 'O', 'O', 'C', 'C'],
            'bonds': [(0, 1, 'single'), (0, 2, 'double'), (0, 3, 'single'),
                     (0, 4, 'single'), (4, 5, 'single')],
            'name': 'phosphonate_inhibitor',
            'molecular_weight': 195.0,
            'heavy_atoms': 6
        }
        
        candidates.extend([candidate1, candidate2, candidate3])
        return candidates
    
    def _mock_metal_docking(self, engine: MetalDockingEngine, prefix: str) -> List[MetalPose]:
        """Mock metal docking with realistic results"""
        poses = []
        
        for i in range(5):  # Generate 5 example poses
            # Create realistic coordinates near metal center
            if engine.metal_centers:
                metal_pos = engine.metal_centers[0].coordinates
                # Place ligand near metal with some randomness
                offset = np.random.normal(0, 1.0, 3)
                base_coords = engine.ligand['coordinates'] + offset
                
                # Ensure at least one atom is close to metal
                coordinating_idx = 0  # First atom
                desired_distance = 2.1  # Typical coordination distance
                direction = (base_coords[coordinating_idx] - metal_pos)
                direction = direction / np.linalg.norm(direction)
                base_coords[coordinating_idx] = metal_pos + direction * desired_distance
            else:
                base_coords = engine.ligand['coordinates']
            
            # Create pose with realistic metal properties
            pose = MetalPose(
                coordinates=base_coords,
                score=-8.0 + i * 0.5 + np.random.normal(0, 0.3),  # Realistic scores
                energy=-8.5 + i * 0.6 + np.random.normal(0, 0.4),
                pose_id=f"{prefix}_pose_{i+1}",
                ligand_name=engine.ligand.get('name', 'unknown')
            )
            
            # Add metal-specific properties
            if engine.metal_centers:
                metal_center = engine.metal_centers[0]
                
                # Mock coordination analysis
                coordinating_atoms = [0, 2] if len(engine.ligand['atom_types']) > 2 else [0]
                pose.coordinating_atoms = coordinating_atoms
                
                # Mock distances
                pose.coordination_distances = [
                    np.linalg.norm(base_coords[idx] - metal_center.coordinates)
                    for idx in coordinating_atoms
                ]
                
                # Mock interactions
                pose.metal_interactions = [
                    {
                        'type': 'metal_coordination',
                        'subtype': 'direct_coordination',
                        'metal_type': metal_center.metal_type.value,
                        'ligand_atom': idx,
                        'ligand_element': engine.ligand['atom_types'][idx],
                        'distance': pose.coordination_distances[j],
                        'energy': -3.0 + np.random.normal(0, 0.5)
                    }
                    for j, idx in enumerate(coordinating_atoms)
                ]
                
                # Mock quality scores
                pose.coordination_quality = {
                    'coordination_score': 0.8 - i * 0.1 + np.random.normal(0, 0.1),
                    'geometric_score': 0.7 - i * 0.08 + np.random.normal(0, 0.1),
                    'overall_validity': True
                }
                
                # Mock energy breakdown
                pose.metal_energy_breakdown = {
                    'coordination_energy': -4.0 + np.random.normal(0, 0.5),
                    'geometric_penalty': i * 0.2 + np.random.normal(0, 0.1),
                    'electrostatic_energy': -2.0 + np.random.normal(0, 0.3),
                    'solvation_energy': 1.0 + np.random.normal(0, 0.2),
                    'total_metal_energy': pose.energy
                }
            
            # Calculate confidence
            pose.confidence = max(0.1, 0.9 - i * 0.15 + np.random.normal(0, 0.05))
            
            poses.append(pose)
        
        # Sort by score (lower is better)
        poses.sort(key=lambda x: x.score)
        return poses
    
    def _analyze_metal_poses(self, poses: List[MetalPose], title: str):
        """Analyze and display metal pose results"""
        print(f"\n{title} Results:")
        print("="*50)
        
        if not poses:
            print("No poses generated!")
            return
        
        print(f"Generated {len(poses)} poses")
        print(f"Best pose score: {poses[0].score:.3f}")
        print(f"Score range: {poses[0].score:.3f} to {poses[-1].score:.3f}")
        
        # Coordination analysis
        coordination_counts = [len(pose.coordinating_atoms) for pose in poses]
        avg_coordination = np.mean(coordination_counts) if coordination_counts else 0
        print(f"Average coordination: {avg_coordination:.1f} atoms")
        
        # Quality analysis
        confidence_scores = [pose.confidence for pose in poses]
        avg_confidence = np.mean(confidence_scores)
        print(f"Average confidence: {avg_confidence:.3f}")
        
        # Best pose detailed analysis
        best_pose = poses[0]
        print(f"\nBest Pose Analysis:")
        print(f"  Pose ID: {best_pose.pose_id}")
        print(f"  Score: {best_pose.score:.3f}")
        print(f"  Confidence: {best_pose.confidence:.3f}")
        print(f"  Coordinating atoms: {len(best_pose.coordinating_atoms)}")
        
        if best_pose.coordination_distances:
            print(f"  Coordination distances: {', '.join(f'{d:.2f}' for d in best_pose.coordination_distances)} Å")
        
        if best_pose.metal_interactions:
            print(f"  Metal interactions: {len(best_pose.metal_interactions)}")
            for interaction in best_pose.metal_interactions[:3]:  # Show first 3
                print(f"    {interaction.get('ligand_element', 'X')}-{interaction.get('metal_type', 'M')}: "
                      f"{interaction.get('distance', 0):.2f} Å, "
                      f"E = {interaction.get('energy', 0):.2f} kcal/mol")
    
    def _analyze_drug_candidate(self, pose: MetalPose, candidate: Dict, metal_center: MetalCenter) -> Dict[str, Any]:
        """Analyze drug candidate properties"""
        # Calculate binding affinity
        binding_affinity = pose.get_metal_binding_affinity()
        
        # Calculate IC50
        ic50 = self.ic50_calc.delta_g_to_ic50(binding_affinity)
        
        # Calculate ligand efficiency
        num_heavy_atoms = candidate.get('heavy_atoms', len(candidate['coordinates']))
        le = self.ic50_calc.calculate_ligand_efficiency(binding_affinity, num_heavy_atoms)
        
        # Drug-likeness assessment
        mol_weight = candidate.get('molecular_weight', 200.0)
        drug_likeness = "Good" if mol_weight < 500 and le > 0.3 else "Poor"
        
        return {
            'name': candidate['name'],
            'pose_score': pose.score,
            'binding_affinity': binding_affinity,
            'predicted_ic50': ic50,
            'ligand_efficiency': le,
            'molecular_weight': mol_weight,
            'coordination_atoms': len(pose.coordinating_atoms),
            'confidence': pose.confidence,
            'drug_likeness': drug_likeness
        }
    
    def _rank_drug_candidates(self, results: List[Dict[str, Any]]):
        """Rank and display drug candidates"""
        print(f"\n{'='*60}")
        print("DRUG CANDIDATE RANKING")
        print(f"{'='*60}")
        
        # Sort by composite score (lower IC50 + higher LE)
        def composite_score(result):
            ic50_score = 1.0 / (result['predicted_ic50'] + 1)  # Lower IC50 is better
            le_score = result['ligand_efficiency']  # Higher LE is better
            return ic50_score + le_score
        
        results.sort(key=composite_score, reverse=True)
        
        print(f"{'Rank':<4} {'Name':<25} {'IC50 (nM)':<10} {'LE':<6} {'MW':<6} {'Score':<8}")
        print("-" * 70)
        
        for i, result in enumerate(results, 1):
            print(f"{i:<4} {result['name']:<25} {result['predicted_ic50']:<10.1f} "
                  f"{result['ligand_efficiency']:<6.3f} {result['molecular_weight']:<6.1f} "
                  f"{result['pose_score']:<8.3f}")
        
        print(f"\nRecommendation: {results[0]['name']} shows the best overall profile")
    
    def _analyze_inter_metal_interactions(self, poses: List[MetalPose], metal_centers: List[MetalCenter]):
        """Analyze interactions between multiple metal centers"""
        print(f"\nMulti-Metal Analysis:")
        print("-" * 30)
        
        # Calculate inter-metal distances
        if len(metal_centers) >= 2:
            metal_distance = np.linalg.norm(
                metal_centers[0].coordinates - metal_centers[1].coordinates
            )
            print(f"Inter-metal distance: {metal_distance:.2f} Å")
        
        # Analyze bridging interactions
        bridging_poses = []
        for pose in poses:
            # Check if ligand coordinates to multiple metals
            metal_coords = [center.coordinates for center in metal_centers]
            ligand_coords = pose.coordinates
            
            coordinating_metals = 0
            for metal_pos in metal_coords:
                min_distance = min(
                    np.linalg.norm(ligand_pos - metal_pos) 
                    for ligand_pos in ligand_coords
                )
                if min_distance < 3.5:
                    coordinating_metals += 1
            
            if coordinating_metals > 1:
                bridging_poses.append(pose)
        
        print(f"Poses with bridging interactions: {len(bridging_poses)}")
        
        if bridging_poses:
            best_bridging = bridging_poses[0]
            print(f"Best bridging pose score: {best_bridging.score:.3f}")
    
    def _print_docking_report(self, report: Dict[str, Any]):
        """Print detailed docking report"""
        print(f"\n{'='*60}")
        print("DETAILED DOCKING REPORT")
        print(f"{'='*60}")
        
        summary = report.get('docking_summary', {})
        print(f"Total poses: {summary.get('total_poses', 0)}")
        print(f"Metal centers detected: {summary.get('metal_centers_detected', 0)}")
        print(f"Poses with coordination: {summary.get('poses_with_coordination', 0)}")
        print(f"Average coordination count: {summary.get('average_coordination_count', 0):.1f}")
        print(f"Best score: {summary.get('best_score', 0):.3f}")
        
        # Metal center details
        metal_centers = report.get('metal_centers', [])
        for metal in metal_centers:
            print(f"\nMetal Center:")
            print(f"  Type: {metal.get('metal_type', 'Unknown')}")
            print(f"  Geometry: {metal.get('geometry', 'Unknown')}")
            print(f"  Coordination number: {metal.get('coordination_number', 0)}")
            print(f"  Poses coordinating: {metal.get('poses_coordinating', 0)}")


def main():
    """Main function to run all metal docking examples"""
    print("PandaDock Metal Docking Examples")
    print("="*50)
    
    examples = MetalDockingExamples()
    
    # Run all examples
    try:
        print("Running Example 1: Zinc Finger Docking")
        poses1, report1 = examples.example_1_zinc_finger_docking()
        
        print("\nRunning Example 2: Iron Enzyme Docking")
        poses2 = examples.example_2_iron_enzyme_docking()
        
        print("\nRunning Example 3: Calcium-Binding Protein")
        poses3 = examples.example_3_calcium_binding_protein()
        
        print("\nRunning Example 4: Multi-Metal Cluster")
        poses4 = examples.example_4_multi_metal_cluster()
        
        print("\nRunning Example 5: Constraint Optimization")
        results5 = examples.example_5_constraint_optimization()
        
        print("\nRunning Example 6: Drug Design Workflow")
        drug_results = examples.example_6_drug_design_workflow()
        
        print("\n" + "="*60)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("="*60)
        
    except Exception as e:
        logger.error(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()