Quick Start Guide
================

This guide will get you up and running with PandaDock in just a few minutes.

Basic Docking Workflow
-----------------------

Here's a complete example of performing molecular docking with PandaDock:

.. code-block:: python

   from pandadock import PandaDock
   import os
   
   # Initialize the docking engine
   docker = PandaDock(
       engine='ml',           # Use ML-enhanced docking
       scoring='vina',        # AutoDock Vina scoring function
       exhaustiveness=8,      # Search thoroughness
       num_poses=10          # Number of poses to generate
   )
   
   # Perform docking
   results = docker.dock(
       receptor='examples/protein.pdb',     # Protein structure
       ligand='examples/ligand.sdf',        # Ligand structure
       center=[25.0, 30.0, 15.0],          # Binding site center
       size=[20.0, 20.0, 20.0],            # Search box size
       output_dir='docking_results/'       # Output directory
   )
   
   # Print results
   print(f"Generated {len(results.poses)} poses")
   print(f"Best score: {results.best_pose.score:.3f}")
   print(f"Best binding affinity: {results.best_pose.binding_affinity:.2f} kcal/mol")

Input File Formats
------------------

PandaDock supports multiple input formats:

**Protein (Receptor) Formats:**
- PDB (.pdb) - Preferred format
- MOL2 (.mol2)
- SDF (.sdf) - For small proteins

**Ligand Formats:**
- SDF (.sdf) - Preferred format
- MOL2 (.mol2)
- PDB (.pdb)
- SMILES (.smi) - Converted to 3D automatically

Example with different formats:

.. code-block:: python

   # Using SMILES input
   results = docker.dock(
       receptor='protein.pdb',
       ligand='CCO',  # Ethanol SMILES
       center=[25.0, 30.0, 15.0],
       size=[15.0, 15.0, 15.0]
   )
   
   # Using MOL2 format
   results = docker.dock(
       receptor='protein.mol2',
       ligand='ligand.mol2',
       center=[25.0, 30.0, 15.0],
       size=[15.0, 15.0, 15.0]
   )

Configuration Options
---------------------

Customize docking behavior with various options:

.. code-block:: python

   # Advanced configuration
   docker = PandaDock(
       engine='physics',              # Available: 'physics', 'ml', 'ga'
       scoring='vina',               # Available: 'vina', 'chemplp', 'custom'
       exhaustiveness=16,            # Higher = more thorough (1-32)
       num_poses=20,                # Number of output poses (1-100)
       energy_range=3.0,            # kcal/mol range for poses
       seed=42,                     # Random seed for reproducibility
       cpu_threads=8,               # Number of CPU threads
       gpu_acceleration=True        # Use GPU if available
   )

Docking Engines
---------------

Choose the appropriate docking engine for your needs:

**Physics Engine (Default)**
- Fast and reliable
- Based on AutoDock Vina algorithm
- Good for general-purpose docking

.. code-block:: python

   docker = PandaDock(engine='physics')

**Machine Learning Engine**
- State-of-the-art accuracy
- Slower but more precise
- Best for drug discovery projects

.. code-block:: python

   docker = PandaDock(engine='ml')

**Genetic Algorithm Engine**
- Excellent for flexible receptors
- Good conformational sampling
- Best for induced-fit docking

.. code-block:: python

   docker = PandaDock(engine='ga')

Working with Results
--------------------

Access and analyze docking results:

.. code-block:: python

   # Get the best pose
   best_pose = results.best_pose
   print(f"Score: {best_pose.score}")
   print(f"Energy: {best_pose.energy} kcal/mol")
   print(f"IC50: {best_pose.ic50} nM")
   
   # Iterate through all poses
   for i, pose in enumerate(results.poses):
       print(f"Pose {i+1}: Score={pose.score:.3f}, Energy={pose.energy:.2f}")
   
   # Get binding site interactions
   interactions = best_pose.interactions
   print(f"H-bonds: {len(interactions.hbonds)}")
   print(f"Hydrophobic contacts: {len(interactions.hydrophobic)}")
   print(f"Salt bridges: {len(interactions.salt_bridges)}")

Saving Results
--------------

Save results in various formats:

.. code-block:: python

   # Save all poses as SDF
   results.save_poses('poses.sdf', format='sdf')
   
   # Save best pose as PDB
   results.best_pose.save('best_pose.pdb', format='pdb')
   
   # Generate HTML report
   results.generate_report('docking_report.html')
   
   # Export summary as CSV
   results.export_summary('summary.csv')

Batch Processing
----------------

Process multiple ligands efficiently:

.. code-block:: python

   import glob
   
   # Process all SDF files in a directory
   ligand_files = glob.glob('ligands/*.sdf')
   
   for ligand_file in ligand_files:
       print(f"Processing {ligand_file}...")
       
       results = docker.dock(
           receptor='protein.pdb',
           ligand=ligand_file,
           center=[25.0, 30.0, 15.0],
           size=[20.0, 20.0, 20.0],
           output_dir=f'results/{os.path.basename(ligand_file)[:-4]}/'
       )
       
       print(f"Best score: {results.best_pose.score:.3f}")

Virtual Screening
-----------------

Screen large compound libraries:

.. code-block:: python

   from pandadock import VirtualScreening
   
   # Initialize screening
   screening = VirtualScreening(
       receptor='protein.pdb',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0],
       engine='ml'
   )
   
   # Screen a compound library
   results = screening.screen(
       ligand_library='compounds.sdf',  # Multi-molecule SDF
       top_n=100,                       # Keep top 100 results
       output_dir='screening_results/'
   )
   
   # Analyze results
   print(f"Screened {results.total_compounds} compounds")
   print(f"Top compound: {results.top_hits[0].name}")
   print(f"Best score: {results.top_hits[0].score:.3f}")

Visualization
-------------

Visualize docking results:

.. code-block:: python

   # Generate 3D visualization
   results.visualize_3d('docking_viz.html')
   
   # Plot scoring distribution
   results.plot_scores('score_distribution.png')
   
   # Create interaction diagram
   best_pose.plot_interactions('interactions.png')

Command Line Interface
----------------------

Use PandaDock from the command line:

.. code-block:: bash

   # Basic docking
   pandadock dock \
     --receptor protein.pdb \
     --ligand ligand.sdf \
     --center 25,30,15 \
     --size 20,20,20 \
     --output results/
   
   # Virtual screening
   pandadock screen \
     --receptor protein.pdb \
     --library compounds.sdf \
     --center 25,30,15 \
     --size 20,20,20 \
     --top-n 50 \
     --output screening/
   
   # Generate report only
   pandadock report \
     --poses poses.sdf \
     --output report.html

Error Handling
--------------

Handle common errors gracefully:

.. code-block:: python

   from pandadock.exceptions import DockingError, InvalidInputError
   
   try:
       results = docker.dock(
           receptor='protein.pdb',
           ligand='ligand.sdf',
           center=[25.0, 30.0, 15.0],
           size=[20.0, 20.0, 20.0]
       )
   except InvalidInputError as e:
       print(f"Input file error: {e}")
   except DockingError as e:
       print(f"Docking failed: {e}")
   except Exception as e:
       print(f"Unexpected error: {e}")

Performance Tips
----------------

Optimize performance for your system:

.. code-block:: python

   # Use multiple CPU cores
   docker = PandaDock(cpu_threads=8)
   
   # Enable GPU acceleration
   docker = PandaDock(gpu_acceleration=True)
   
   # Adjust exhaustiveness based on needs
   docker = PandaDock(exhaustiveness=8)  # Fast
   docker = PandaDock(exhaustiveness=32) # Thorough
   
   # Limit memory usage
   docker = PandaDock(max_memory_gb=8)

Next Steps
----------

- Learn about :doc:`docking_modes` for specific use cases
- Explore :doc:`scoring_functions` to customize scoring
- Check out the :doc:`../tutorials/basic_docking` tutorial
- See :doc:`../examples/protein_ligand_docking` for detailed examples