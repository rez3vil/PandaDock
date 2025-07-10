Metal Docking Tutorial
=====================

This tutorial will guide you through performing metal docking with PandaDock, specifically targeting metalloproteins that contain metal ions in their active sites.

Prerequisites
-------------

- Completed the :doc:`basic_docking` tutorial
- Understanding of metalloprotein structure and coordination chemistry
- PDB file containing metal ions
- Ligand that can potentially coordinate to metals

Tutorial Overview
-----------------

You will learn to:

1. Set up metal docking configurations
2. Detect and analyze metal centers
3. Perform metal-aware docking
4. Analyze coordination interactions
5. Optimize poses with geometric constraints

Step 1: Setting Up Metal Docking
---------------------------------

Start by importing the metal docking modules:

.. code-block:: python

   from pandadock.docking.metal_docking_engine import MetalDockingEngine, MetalDockingConfig
   from pandadock.docking.metal_coordination import MetalType, CoordinationGeometry
   from pandadock.scoring.metal_scoring import MetalScoringFunction
   
   # Configure metal docking parameters
   metal_config = MetalDockingConfig(
       detect_metals_automatically=True,      # Auto-detect metals in PDB
       enforce_geometric_constraints=True,    # Apply coordination geometry
       coordination_focused_sampling=True,    # Focus sampling near metals
       use_metal_scoring=True,               # Use metal-specific scoring
       require_metal_coordination=True,       # Require coordination
       min_coordinating_atoms=1,             # Minimum coordination
       geometric_constraint_weight=2.0       # Constraint penalty weight
   )

Step 2: Initialize Metal Docking Engine
----------------------------------------

Create the metal docking engine with your configuration:

.. code-block:: python

   # Create basic config (adapt to your needs)
   class BasicConfig:
       def __init__(self):
           self.docking = type('obj', (object,), {'num_poses': 10, 'exhaustiveness': 8})
           self.io = type('obj', (object,), {
               'center_x': 25.0, 'center_y': 30.0, 'center_z': 15.0,
               'size_x': 20.0, 'size_y': 20.0, 'size_z': 20.0
           })
       def to_dict(self): return {}
   
   config = BasicConfig()
   
   # Initialize metal docking engine
   metal_docker = MetalDockingEngine(config, metal_config)
   
   print("Metal docking engine initialized")
   print(f"Auto-detect metals: {metal_config.detect_metals_automatically}")
   print(f"Geometric constraints: {metal_config.enforce_geometric_constraints}")

Step 3: Prepare Metalloprotein Structure
-----------------------------------------

Load and analyze the metalloprotein structure:

.. code-block:: python

   # Load protein structure (example with zinc finger)
   protein_file = "zinc_finger.pdb"  # Replace with your PDB file
   
   # Prepare receptor and detect metal centers
   metal_docker.prepare_receptor(protein_file)
   
   # Display detected metal centers
   print(f"Detected {len(metal_docker.metal_centers)} metal centers:")
   
   for i, metal_center in enumerate(metal_docker.metal_centers):
       print(f"\nMetal Center {i+1}:")
       print(f"  Type: {metal_center.metal_type.value}")
       print(f"  Coordinates: {metal_center.coordinates}")
       print(f"  Geometry: {metal_center.geometry.value}")
       print(f"  Coordination number: {metal_center.coordination_number}")
       print(f"  Oxidation state: +{metal_center.oxidation_state}")
       
       # Show coordinating residues
       print(f"  Coordinating residues:")
       for coord_atom in metal_center.coordinating_atoms:
           print(f"    {coord_atom['residue']} {coord_atom['element']}: "
                 f"{coord_atom['coordinates']} "
                 f"(d = {np.linalg.norm(coord_atom['coordinates'] - metal_center.coordinates):.2f} Å)")

Step 4: Prepare Metal-Binding Ligand
-------------------------------------

Prepare a ligand that can potentially bind to the metal:

.. code-block:: python

   # Example: Create a simple zinc chelator
   # In practice, load from SDF file: metal_docker.prepare_ligand("chelator.sdf")
   
   import numpy as np
   
   # Simple bidentate ligand with N and O coordination
   chelator_coords = np.array([
       [26.0, 30.0, 15.0],  # N - potential zinc coordinator
       [26.5, 30.3, 15.2],  # C
       [27.0, 29.8, 15.1],  # O - potential zinc coordinator
       [27.5, 30.1, 14.9],  # C
       [26.2, 29.6, 14.8],  # C
       [25.7, 29.9, 15.1],  # C
   ])
   
   chelator_atoms = ['N', 'C', 'O', 'C', 'C', 'C']
   chelator_bonds = [
       (0, 1, 'single'), (1, 2, 'single'), (2, 3, 'single'),
       (1, 4, 'single'), (4, 5, 'single'), (5, 0, 'single')
   ]
   
   # Set ligand data
   metal_docker.ligand = {
       'coordinates': chelator_coords,
       'atom_types': chelator_atoms,
       'bonds': chelator_bonds,
       'name': 'zinc_chelator'
   }
   
   print(f"Prepared ligand: {metal_docker.ligand['name']}")
   print(f"Atoms: {len(chelator_coords)}")
   print(f"Potential coordinators: {[i for i, atom in enumerate(chelator_atoms) if atom in ['N', 'O', 'S', 'P']]}")

Step 5: Perform Metal Docking
------------------------------

Run the metal-aware docking:

.. code-block:: python

   # Perform metal docking
   print("\nStarting metal docking...")
   metal_poses = metal_docker.dock(protein_file, "ligand_data")  # Mock call
   
   # In practice, the dock method would be called like this:
   # metal_poses = metal_docker.dock(protein_file, ligand_file)
   
   # For this tutorial, we'll simulate the results
   # Create mock poses for demonstration
   mock_poses = []
   for i in range(5):
       from pandadock.docking.metal_docking_engine import MetalPose
       
       # Create example pose
       pose_coords = chelator_coords + np.random.normal(0, 0.3, chelator_coords.shape)
       
       pose = MetalPose(
           coordinates=pose_coords,
           score=-8.0 + i * 0.8,
           energy=-8.5 + i * 0.9,
           pose_id=f"metal_pose_{i+1}",
           ligand_name="zinc_chelator",
           coordinating_atoms=[0, 2],  # N and O atoms
           coordination_distances=[2.1, 2.0],  # Typical Zn-N and Zn-O distances
           confidence=0.9 - i * 0.1
       )
       
       # Add metal interaction data
       pose.metal_interactions = [
           {
               'type': 'metal_coordination',
               'metal_type': 'Zn',
               'ligand_atom': 0,
               'ligand_element': 'N',
               'distance': 2.1,
               'energy': -3.5,
               'subtype': 'direct_coordination'
           },
           {
               'type': 'metal_coordination', 
               'metal_type': 'Zn',
               'ligand_atom': 2,
               'ligand_element': 'O', 
               'distance': 2.0,
               'energy': -3.2,
               'subtype': 'direct_coordination'
           }
       ]
       
       # Add coordination quality scores
       pose.coordination_quality = {
           'coordination_score': 0.85 - i * 0.1,
           'geometric_score': 0.8 - i * 0.08,
           'overall_validity': True
       }
       
       mock_poses.append(pose)
   
   metal_poses = mock_poses
   
   print(f"Generated {len(metal_poses)} metal-coordinated poses")
   print(f"Best pose score: {metal_poses[0].score:.3f}")

Step 6: Analyze Metal Coordination
-----------------------------------

Analyze the coordination interactions in detail:

.. code-block:: python

   print("\nMetal Coordination Analysis:")
   print("=" * 50)
   
   # Analyze the best pose
   best_pose = metal_poses[0]
   
   print(f"Best Pose: {best_pose.pose_id}")
   print(f"  Overall Score: {best_pose.score:.3f}")
   print(f"  Confidence: {best_pose.confidence:.3f}")
   print(f"  Coordinating Atoms: {len(best_pose.coordinating_atoms)}")
   
   # Display coordination details
   print(f"\nCoordination Details:")
   for i, atom_idx in enumerate(best_pose.coordinating_atoms):
       atom_type = chelator_atoms[atom_idx]
       distance = best_pose.coordination_distances[i]
       print(f"  Atom {atom_idx} ({atom_type}): {distance:.2f} Å from metal")
   
   # Display metal interactions
   print(f"\nMetal-Ligand Interactions:")
   for interaction in best_pose.metal_interactions:
       print(f"  {interaction['ligand_element']}-{interaction['metal_type']}:")
       print(f"    Distance: {interaction['distance']:.2f} Å")
       print(f"    Energy: {interaction['energy']:.2f} kcal/mol")
       print(f"    Type: {interaction['subtype']}")
   
   # Coordination quality assessment
   print(f"\nCoordination Quality:")
   quality = best_pose.coordination_quality
   print(f"  Coordination Score: {quality['coordination_score']:.3f}")
   print(f"  Geometric Score: {quality['geometric_score']:.3f}")
   print(f"  Overall Validity: {quality['overall_validity']}")

Step 7: Generate Metal Docking Report
--------------------------------------

Generate a comprehensive report of the metal docking results:

.. code-block:: python

   # Generate metal docking report
   report = metal_docker.get_metal_docking_report(metal_poses)
   
   print("\nMetal Docking Report:")
   print("=" * 50)
   
   # Overall summary
   summary = report['docking_summary']
   print(f"Total Poses: {summary['total_poses']}")
   print(f"Metal Centers Detected: {summary['metal_centers_detected']}")
   print(f"Poses with Coordination: {summary['poses_with_coordination']}")
   print(f"Average Coordination Count: {summary['average_coordination_count']:.1f}")
   print(f"Best Score: {summary['best_score']:.3f}")
   
   # Metal center analysis
   print(f"\nMetal Center Analysis:")
   for metal in report['metal_centers']:
       print(f"  Metal {metal['metal_id']}: {metal['metal_type']}")
       print(f"    Geometry: {metal['geometry']}")
       print(f"    Poses Coordinating: {metal['poses_coordinating']}")
   
   # Best pose analysis
   if report['best_pose']:
       best = report['best_pose']
       print(f"\nBest Pose Detailed Analysis:")
       print(f"  Pose ID: {best['pose_id']}")
       print(f"  Score: {best['score']:.3f}")
       print(f"  Coordinating Atoms: {best['coordinating_atoms']}")
       print(f"  Coordination Distances: {[f'{d:.2f}' for d in best['coordination_distances']]}")
       print(f"  Number of Interactions: {len(best['interactions'])}")

Step 8: Constraint Optimization
--------------------------------

Optimize poses using geometric constraints:

.. code-block:: python

   from pandadock.utils.metal_constraints import apply_metal_constraints_to_pose
   
   print("\nConstraint Optimization:")
   print("=" * 40)
   
   # Test different constraint presets
   constraint_presets = ["standard", "strict", "flexible"]
   
   for preset in constraint_presets:
       print(f"\nTesting {preset.upper()} constraints:")
       
       # Apply constraints to best pose
       initial_coords = best_pose.coordinates
       atom_types = metal_docker.ligand['atom_types']
       
       optimized_coords, constraint_results = apply_metal_constraints_to_pose(
           initial_coords,
           atom_types, 
           metal_docker.metal_centers,
           constraint_preset=preset
       )
       
       # Analyze results
       satisfaction = constraint_results['constraint_satisfaction']
       violations = constraint_results['total_violations']
       penalty = constraint_results['final_penalty']
       
       print(f"  Constraint Satisfaction: {satisfaction:.3f}")
       print(f"  Total Violations: {violations}")
       print(f"  Final Penalty: {penalty:.3f}")
       
       # Calculate coordinate change
       coord_rmsd = np.sqrt(np.mean((optimized_coords - initial_coords)**2))
       print(f"  Coordinate RMSD: {coord_rmsd:.3f} Å")
       
       # Show individual violations
       if constraint_results['violations']:
           print(f"  Violations:")
           for violation in constraint_results['violations'][:3]:  # Show first 3
               print(f"    - {violation['type']}: "
                     f"expected {violation['expected_value']:.2f}, "
                     f"got {violation['actual_value']:.2f}")

Step 9: Calculate Binding Properties
------------------------------------

Calculate metal-corrected binding properties:

.. code-block:: python

   from pandadock.utils.ic50_calculator import IC50Calculator
   
   ic50_calc = IC50Calculator()
   
   print("\nBinding Properties Analysis:")
   print("=" * 40)
   
   for i, pose in enumerate(metal_poses[:3]):  # Analyze top 3 poses
       print(f"\nPose {i+1}:")
       
       # Standard binding affinity
       std_affinity = pose.get_binding_affinity()
       
       # Metal-corrected binding affinity
       metal_affinity = pose.get_metal_binding_affinity()
       
       # Calculate IC50 values
       std_ic50 = ic50_calc.delta_g_to_ic50(std_affinity)
       metal_ic50 = ic50_calc.delta_g_to_ic50(metal_affinity)
       
       print(f"  Standard ΔG: {std_affinity:.2f} kcal/mol")
       print(f"  Metal-corrected ΔG: {metal_affinity:.2f} kcal/mol")
       print(f"  Standard IC50: {std_ic50:.1f} nM")
       print(f"  Metal-corrected IC50: {metal_ic50:.1f} nM")
       
       # Ligand efficiency
       num_heavy_atoms = len([a for a in chelator_atoms if a != 'H'])
       le = ic50_calc.calculate_ligand_efficiency(metal_affinity, num_heavy_atoms)
       print(f"  Ligand Efficiency: {le:.3f} kcal/mol per heavy atom")
       
       # Coordination contribution
       coord_bonus = metal_affinity - std_affinity
       print(f"  Coordination Bonus: {coord_bonus:.2f} kcal/mol")

Step 10: Save Results
---------------------

Save the metal docking results:

.. code-block:: python

   import os
   
   # Create output directory
   output_dir = "metal_docking_results"
   os.makedirs(output_dir, exist_ok=True)
   
   # Save metal poses with additional metal information
   metal_docker.save_metal_poses(metal_poses, output_dir)
   
   print(f"\nResults saved to: {output_dir}/")
   print("Generated files:")
   print("  - metal_docking_report.json: Comprehensive metal analysis")
   print("  - metal_interactions.csv: Detailed interaction data")
   print("  - poses_summary.csv: Pose scores and properties")
   print("  - all.ligands.sdf: All poses in SDF format")
   print("  - Individual pose files: pose_*.pdb and pose_*.sdf")

Advanced Tips
-------------

**1. Customizing Metal Parameters:**

.. code-block:: python

   from pandadock.scoring.metal_scoring import MetalScoringParameters
   
   # Custom parameters for specific metals
   custom_params = MetalScoringParameters(
       coordination_strength={
           'N': -4.5,  # Stronger Zn-N interaction
           'O': -4.0,  # Stronger Zn-O interaction
           'S': -3.5   # Moderate Zn-S interaction
       },
       geometric_penalty_weight=3.0,  # Stricter geometry
       distance_tolerance=0.15,       # Tighter distance tolerance
       angle_tolerance=8.0            # Tighter angle tolerance
   )
   
   # Use in scoring function
   custom_scorer = MetalScoringFunction(metal_docker.metal_centers, custom_params)

**2. Multi-Metal Systems:**

.. code-block:: python

   # For proteins with multiple metal centers
   if len(metal_docker.metal_centers) > 1:
       print(f"Multi-metal system detected: {len(metal_docker.metal_centers)} centers")
       
       # Calculate inter-metal distances
       for i in range(len(metal_docker.metal_centers)):
           for j in range(i+1, len(metal_docker.metal_centers)):
               metal1 = metal_docker.metal_centers[i]
               metal2 = metal_docker.metal_centers[j]
               distance = np.linalg.norm(metal1.coordinates - metal2.coordinates)
               print(f"  {metal1.metal_type.value}-{metal2.metal_type.value} distance: {distance:.2f} Å")

**3. Validation Against Experimental Data:**

.. code-block:: python

   # If you have experimental IC50 data
   experimental_ic50 = 50.0  # nM
   predicted_ic50 = ic50_calc.delta_g_to_ic50(best_pose.get_metal_binding_affinity())
   
   fold_error = max(predicted_ic50/experimental_ic50, experimental_ic50/predicted_ic50)
   print(f"Prediction accuracy: {fold_error:.1f}-fold error")
   
   if fold_error < 3.0:
       print("✓ Good prediction accuracy")
   elif fold_error < 10.0:
       print("⚠ Moderate prediction accuracy")
   else:
       print("✗ Poor prediction accuracy - consider parameter tuning")

Troubleshooting
---------------

**Common Issues:**

1. **No metal centers detected:**
   - Check that PDB file contains HETATM records for metals
   - Verify metal ion names are standard (ZN, FE, etc.)
   - Check metal_detection_distance parameter

2. **No coordination detected:**
   - Ensure ligand has coordinating atoms (N, O, S, P)
   - Check if ligand is positioned near metal center
   - Adjust coordination distance thresholds

3. **High constraint violations:**
   - Use more flexible constraint presets
   - Increase distance and angle tolerances
   - Check if metal geometry assignment is correct

4. **Poor poses:**
   - Increase metal_focused_exhaustiveness
   - Enable coordination_focused_sampling
   - Check ligand conformer generation

This tutorial provides a comprehensive introduction to metal docking with PandaDock. The metal docking capabilities enable accurate modeling of metalloprotein-ligand interactions, which is crucial for drug discovery targeting this important class of proteins.