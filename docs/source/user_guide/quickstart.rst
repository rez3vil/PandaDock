Quick Start Guide
================

This guide will get you up and running with PandaDock's comprehensive analysis and visualization capabilities in just a few minutes.

Command-Line Interface (Recommended)
------------------------------------

The easiest way to get started with PandaDock is using the command-line interface with comprehensive visualization:

**Basic Docking with All Outputs:**

.. code-block:: bash

   # Complete analysis with all visualizations
   pandadock --protein protein.pdb --ligand ligand.sdf --all-outputs
   
   # Generates:
   # - Master publication dashboard
   # - Binding metrics analysis
   # - Score distribution analysis
   # - IC50/EC50 drug discovery metrics
   # - Interactive HTML reports

**Professional Interaction Analysis:**

.. code-block:: bash

   # PandaMap integration for professional interaction visualization
   pandadock --protein protein.pdb --ligand ligand.sdf \
             --pandamap --pandamap-3d --all-outputs
   
   # Generates:
   # - 2D Discovery Studio-style interaction maps
   # - 3D interactive molecular models (HTML)
   # - Detailed interaction reports
   # - Complete analysis dashboard

**Algorithm-Specific Analysis:**

.. code-block:: bash

   # PandaML for superior affinity prediction
   pandadock --protein protein.pdb --ligand ligand.sdf \
             --scoring pandaml --master-plot --pandamap
   
   # PandaPhysics for metal complexes and detailed analysis
   pandadock --protein protein.pdb --ligand ligand.sdf \
             --scoring pandaphysics --flexible-residues "HIS57,SER195" \
             --pandamap --all-outputs
   
   # PandaCore for fast, reliable docking
   pandadock --protein protein.pdb --ligand ligand.sdf \
             --scoring pandacore --plots --master-plot

Python API Workflow
--------------------

For programmatic access and custom analysis:

.. code-block:: python

   from pandadock import PandaDock
   import os
   
   # Initialize the docking engine with comprehensive analysis
   docker = PandaDock(
       engine='pandaml',          # Use PandaML algorithm
       scoring='pandaml',         # PandaML scoring function
       exhaustiveness=8,          # Search thoroughness
       num_poses=10,             # Number of poses to generate
       generate_plots=True,       # Enable visualization
       pandamap_analysis=True     # Enable PandaMap integration
   )
   
   # Perform docking with comprehensive analysis
   results = docker.dock(
       receptor='examples/protein.pdb',     # Protein structure
       ligand='examples/ligand.sdf',        # Ligand structure
       center=[25.0, 30.0, 15.0],          # Binding site center
       size=[20.0, 20.0, 20.0],            # Search box size
       output_dir='docking_results/',       # Output directory
       analysis_level='comprehensive'       # Full analysis suite
   )
   
   # Print results with comprehensive metrics
   print(f"Generated {len(results.poses)} poses")
   print(f"Best score: {results.best_pose.score:.3f}")
   print(f"Best binding affinity: {results.best_pose.binding_affinity:.2f} kcal/mol")
   print(f"IC50 prediction: {results.best_pose.ic50:.2e} M")
   print(f"Confidence score: {results.best_pose.confidence:.3f}")
   
   # Access generated visualizations
   print("Generated visualizations:")
   for plot_file in results.visualization_files:
       print(f"  - {plot_file}")

Generated Outputs Overview
--------------------------

PandaDock generates comprehensive publication-ready outputs:

**Master Publication Dashboard:**

.. image:: /_static/master_publication.png
   :alt: Master Publication Dashboard
   :width: 100%
   :align: center

**Professional Interaction Maps:**

.. image:: /_static/pandamap_2d_ml_pose_1.png
   :alt: PandaMap Professional Interaction Analysis
   :width: 60%
   :align: center

**Complete Output Structure:**

.. code-block:: text

   docking_results/
   ├── master_publication.png           # Main analysis dashboard
   ├── binding_metrics_analysis.png     # Statistical validation
   ├── score_distribution_analysis.png  # Score validation
   ├── ic50_ec50_analysis.png           # Drug discovery metrics
   ├── pandamap_2d_pose_1.png          # 2D interaction maps
   ├── pandamap_3d_pose_1.html         # 3D interactive models
   ├── complex_interactions.png         # Interaction networks
   ├── pandadock_report.html           # Interactive HTML report
   ├── pandadock_report.json           # Structured data
   ├── detailed_analysis_report.txt     # Text summary
   └── poses/                          # Molecular structures
       ├── pose_1.pdb
       ├── pose_2.pdb
       └── complex_1.pdb

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
       engine='pandaphysics',         # Available: 'pandacore', 'pandaml', 'pandaphysics'
       scoring='pandaphysics',        # Available: 'pandacore', 'pandaml', 'pandaphysics'
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

**PandaCore Algorithm (Baseline)**
- Fast and reliable
- Robust baseline performance
- Good for general-purpose docking

.. code-block:: python

   docker = PandaDock(engine='pandacore')

**PandaML Algorithm (Advanced)**
- Superior affinity prediction (R² = 0.845)
- Machine learning-enhanced accuracy
- Best for drug discovery projects

.. code-block:: python

   docker = PandaDock(engine='pandaml')

**PandaPhysics Algorithm (Specialized)**
- Excellent for metal complexes
- Physics-based coordination modeling
- Best for metalloproteins and complex chemistry

.. code-block:: python

   docker = PandaDock(engine='pandaphysics')

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
       engine='pandaml'
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