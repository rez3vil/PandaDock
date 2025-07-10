Basic Docking Tutorial
======================

This tutorial will guide you through performing basic molecular docking using PandaDock.

Prerequisites
-------------

Before starting, ensure you have:

- PandaDock installed
- A protein structure (PDB format)
- A ligand structure (SDF, MOL2, or SMILES format)
- Basic understanding of molecular docking concepts

Tutorial Overview
-----------------

In this tutorial, you will learn to:

1. Prepare protein and ligand structures
2. Configure docking parameters
3. Run docking calculations
4. Analyze results
5. Visualize poses

Step 1: Preparing Your Files
----------------------------

Create a new directory for this tutorial:

.. code-block:: bash

   mkdir docking_tutorial
   cd docking_tutorial

Download example files or use your own:

.. code-block:: bash

   # Download example receptor
   wget https://files.rcsb.org/download/1HSG.pdb
   
   # Or use your own protein structure
   cp /path/to/your/protein.pdb receptor.pdb

For the ligand, we'll use a SMILES string in this example:

.. code-block:: python

   # Example: Benzene ring with hydroxyl group
   ligand_smiles = "C1=CC=C(C=C1)O"

Step 2: Basic Docking Script
----------------------------

Create a Python script called ``basic_docking.py``:

.. code-block:: python

   from pandadock import PandaDock
   import os
   
   # Initialize PandaDock with PandaCore algorithm
   docker = PandaDock(
       engine='pandacore',
       exhaustiveness=8,
       num_poses=10,
       energy_range=3.0
   )
   
   # Define docking parameters
   receptor_file = "receptor.pdb"
   ligand_smiles = "C1=CC=C(C=C1)O"  # Phenol
   
   # Define binding site (center coordinates and size)
   binding_site = {
       'center': [25.0, 30.0, 15.0],  # Adjust based on your protein
       'size': [20.0, 20.0, 20.0]    # Search box size in Angstroms
   }
   
   # Run docking
   print("Starting docking calculation...")
   results = docker.dock(
       receptor=receptor_file,
       ligand=ligand_smiles,
       center=binding_site['center'],
       size=binding_site['size']
   )
   
   # Print results
   print(f"Docking completed in {results.runtime:.2f} seconds")
   print(f"Number of poses found: {len(results.poses)}")
   print(f"Best pose score: {results.best_pose.score:.3f}")
   print(f"Best pose energy: {results.best_pose.energy:.2f} kcal/mol")

Step 3: Finding the Binding Site
---------------------------------

If you don't know the binding site coordinates, you can use PandaDock's cavity detection:

.. code-block:: python

   from pandadock.utils import cavity_detection
   
   # Automatically detect binding cavities
   cavities = cavity_detection.find_cavities(receptor_file)
   
   print(f"Found {len(cavities)} potential binding sites:")
   for i, cavity in enumerate(cavities):
       print(f"Cavity {i+1}:")
       print(f"  Center: {cavity.center}")
       print(f"  Volume: {cavity.volume:.1f} Ų")
       print(f"  Druggability score: {cavity.druggability:.2f}")
   
   # Use the most druggable cavity
   best_cavity = max(cavities, key=lambda x: x.druggability)
   binding_site = {
       'center': best_cavity.center,
       'size': [20.0, 20.0, 20.0]
   }

Step 4: Advanced Configuration
------------------------------

For more control over the docking process:

.. code-block:: python

   # Advanced configuration with PandaPhysics algorithm
   docker = PandaDock(
       engine='pandaphysics',
       config={
           'exhaustiveness': 16,          # Higher for better accuracy
           'num_poses': 20,               # More poses
           'energy_range': 4.0,           # Wider energy range
           'force_field': 'amber',        # Force field choice
           'minimization_steps': 2000,    # Energy minimization
           'rigid_receptor': True,        # Keep receptor rigid
           'add_hydrogens': True,         # Add missing hydrogens
           'assign_charges': True,        # Assign partial charges
           'pH': 7.4                      # Physiological pH
       }
   )

Step 5: Analyzing Results
-------------------------

Examine the docking results in detail:

.. code-block:: python

   # Analyze all poses
   print("\nDetailed pose analysis:")
   print("-" * 50)
   
   for i, pose in enumerate(results.poses):
       print(f"Pose {i+1}:")
       print(f"  Score: {pose.score:.3f}")
       print(f"  Energy: {pose.energy:.2f} kcal/mol")
       print(f"  RMSD from best: {pose.rmsd:.2f} Å")
       
       # Energy decomposition
       if hasattr(pose, 'energy_terms'):
           print(f"  Energy breakdown:")
           for term, value in pose.energy_terms.items():
               print(f"    {term}: {value:.2f} kcal/mol")
       
       # Interactions
       if hasattr(pose, 'interactions'):
           print(f"  Interactions:")
           print(f"    H-bonds: {len(pose.interactions.hbonds)}")
           print(f"    Hydrophobic: {len(pose.interactions.hydrophobic)}")
           print(f"    Salt bridges: {len(pose.interactions.salt_bridges)}")
       
       print()

Step 6: Saving Results
----------------------

Save your results for further analysis:

.. code-block:: python

   # Save poses in SDF format
   results.save_poses("docking_poses.sdf")
   
   # Save detailed report
   results.save_report("docking_report.html")
   
   # Save CSV summary
   results.save_csv("pose_summary.csv")
   
   # Save individual pose files
   for i, pose in enumerate(results.poses):
       pose.save(f"pose_{i+1}.pdb")

Step 7: Visualization
---------------------

Visualize the results (requires PyMOL):

.. code-block:: python

   from pandadock.visualization import PyMOLVisualizer
   
   # Create visualizer
   viz = PyMOLVisualizer()
   
   # Load protein and best pose
   viz.load_receptor(receptor_file)
   viz.load_poses(results.poses[:5])  # Show top 5 poses
   
   # Customize visualization
   viz.show_binding_site(binding_site)
   viz.show_interactions(results.best_pose)
   viz.color_by_score()
   
   # Save visualization
   viz.save_image("docking_result.png")
   viz.save_session("docking_session.pse")

Step 8: Validating Results
--------------------------

If you have experimental data, validate your results:

.. code-block:: python

   # Compare with known active site
   if os.path.exists("reference_ligand.sdf"):
       reference_pose = docker.load_pose("reference_ligand.sdf")
       
       # Calculate RMSD to reference
       rmsd_to_ref = results.best_pose.rmsd(reference_pose)
       print(f"RMSD to reference structure: {rmsd_to_ref:.2f} Å")
       
       # Success criteria (typically < 2.0 Å)
       if rmsd_to_ref < 2.0:
           print("✓ Docking successful (RMSD < 2.0 Å)")
       else:
           print("⚠ Docking may need refinement")

Complete Example Script
-----------------------

Here's the complete working script:

.. code-block:: python

   #!/usr/bin/env python3
   """
   PandaDock Basic Docking Tutorial
   
   This script demonstrates basic molecular docking using PandaDock.
   """
   
   from pandadock import PandaDock
   from pandadock.utils import cavity_detection
   import os
   import sys
   
   def main():
       # Check if receptor file exists
       receptor_file = "receptor.pdb"
       if not os.path.exists(receptor_file):
           print(f"Error: {receptor_file} not found!")
           print("Please download a PDB file or adjust the filename.")
           sys.exit(1)
       
       # Initialize PandaDock with PandaCore algorithm
       print("Initializing PandaDock...")
       docker = PandaDock(
           engine='pandacore',
           exhaustiveness=8,
           num_poses=10,
           energy_range=3.0
       )
       
       # Define ligand (using SMILES)
       ligand_smiles = "C1=CC=C(C=C1)O"  # Phenol
       
       # Detect binding sites
       print("Detecting binding sites...")
       try:
           cavities = cavity_detection.find_cavities(receptor_file)
           if cavities:
               best_cavity = max(cavities, key=lambda x: x.druggability)
               binding_site = {
                   'center': best_cavity.center,
                   'size': [20.0, 20.0, 20.0]
               }
               print(f"Using cavity at {best_cavity.center}")
           else:
               # Fallback to manual coordinates
               binding_site = {
                   'center': [0.0, 0.0, 0.0],
                   'size': [20.0, 20.0, 20.0]
               }
               print("No cavities detected, using default coordinates")
       except Exception as e:
           print(f"Cavity detection failed: {e}")
           # Use default binding site
           binding_site = {
               'center': [0.0, 0.0, 0.0],
               'size': [20.0, 20.0, 20.0]
           }
       
       # Run docking
       print("Starting docking calculation...")
       try:
           results = docker.dock(
               receptor=receptor_file,
               ligand=ligand_smiles,
               center=binding_site['center'],
               size=binding_site['size']
           )
           
           # Print results
           print(f"\nDocking completed successfully!")
           print(f"Runtime: {results.runtime:.2f} seconds")
           print(f"Number of poses: {len(results.poses)}")
           print(f"Best pose score: {results.best_pose.score:.3f}")
           print(f"Best pose energy: {results.best_pose.energy:.2f} kcal/mol")
           
           # Save results
           print("\nSaving results...")
           results.save_poses("docking_poses.sdf")
           results.save_report("docking_report.html")
           results.save_csv("pose_summary.csv")
           
           print("Results saved:")
           print("- docking_poses.sdf: All poses")
           print("- docking_report.html: Interactive report")
           print("- pose_summary.csv: Summary table")
           
       except Exception as e:
           print(f"Docking failed: {e}")
           sys.exit(1)
   
   if __name__ == "__main__":
       main()

Running the Tutorial
--------------------

1. Save the script as ``basic_docking.py``
2. Ensure you have a protein structure file
3. Run the script:

.. code-block:: bash

   python basic_docking.py

Expected Output
---------------

You should see output similar to:

.. code-block:: text

   Initializing PandaDock...
   Detecting binding sites...
   Using cavity at [12.3, 15.7, 22.1]
   Starting docking calculation...
   
   Docking completed successfully!
   Runtime: 45.23 seconds
   Number of poses: 10
   Best pose score: -7.8
   Best pose energy: -8.2 kcal/mol
   
   Saving results...
   Results saved:
   - docking_poses.sdf: All poses
   - docking_report.html: Interactive report
   - pose_summary.csv: Summary table

Troubleshooting
---------------

**Common Issues:**

1. **"Receptor file not found"**
   - Ensure the PDB file exists in the current directory
   - Check file permissions

2. **"No cavities detected"**
   - Manually specify binding site coordinates
   - Use a different cavity detection method

3. **"Docking failed"**
   - Check ligand SMILES format
   - Verify binding site coordinates are reasonable
   - Increase exhaustiveness for difficult cases

4. **Poor docking results**
   - Increase exhaustiveness (16-32)
   - Try PandaML algorithm for better affinity prediction
   - Use PandaPhysics for metal complexes
   - Check protein preparation

**Getting Help:**

- Check the :doc:`../user_guide/troubleshooting` guide
- Visit the GitHub issues page
- Join the PandaDock community forum

Next Steps
----------

After completing this tutorial, you can:

- Try the :doc:`ml_enhanced_docking` tutorial
- Learn about :doc:`virtual_screening` 
- Explore :doc:`../examples/flexible_docking`
- Read about :doc:`../user_guide/scoring_functions`

Congratulations on completing your first docking calculation with PandaDock!