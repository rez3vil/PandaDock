Docking Modes
=============

PandaDock offers multiple docking modes optimized for different use cases and computational requirements. This guide explains when and how to use each mode.

Overview of Docking Modes
--------------------------

.. list-table:: Docking Mode Comparison
   :header-rows: 1
   :widths: 20 15 15 15 15 20

   * - Mode
     - Speed
     - Accuracy
     - Use Case
     - Memory
     - Best For
   * - PandaPhysics
     - Fast
     - Good
     - General docking
     - Low
     - Most applications
   * - PandaML
     - Medium
     - Excellent
     - Drug discovery
     - Medium
     - High accuracy needs
   * - PandaCore
     - Slow
     - Very Good
     - Flexible receptors
     - High
     - Baseline algorithm
   * - Hybrid
     - Medium
     - Excellent
     - Best of both
     - Medium
     - Production use

PandaPhysics Algorithm
---------------------

The PandaPhysics algorithm uses classical force fields and energy minimization for pose prediction.

**Key Features:**
- Physics-based docking algorithm
- Fast and reliable
- Good balance of speed and accuracy
- Low memory requirements

**When to Use:**
- General-purpose molecular docking
- Virtual screening of large libraries
- Initial pose generation
- Educational purposes

**Configuration:**

.. code-block:: python

   from pandadock import PandaDock
   
   docker = PandaDock(
       engine='pandaphysics',
       exhaustiveness=8,        # Search thoroughness (1-32)
       num_poses=10,           # Number of output poses
       energy_range=3.0,       # Energy range for poses (kcal/mol)
       force_field='amber',    # Force field (amber, charmm, gaff)
       minimization_steps=1000 # Energy minimization steps
   )

**Example:**

.. code-block:: python

   # Standard docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   print(f"PandaPhysics docking completed in {results.runtime:.2f} seconds")
   print(f"Best score: {results.best_pose.score:.3f}")

**Performance Tips:**
- Use ``exhaustiveness=8`` for screening
- Use ``exhaustiveness=16`` for accurate docking
- Enable parallel processing for multiple ligands

PandaML Algorithm
----------------

PandaML algorithm uses deep learning models trained on experimental binding data.

**Key Features:**
- State-of-the-art accuracy
- Uncertainty quantification
- Ensemble predictions
- Chemical space aware

**When to Use:**
- Drug discovery projects
- Lead optimization
- Binding affinity prediction
- Cases requiring high accuracy

**Configuration:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',    # transformer, cnn, graph
       ensemble_size=5,            # Number of ensemble models
       uncertainty=True,           # Calculate prediction uncertainty
       confidence_threshold=0.8,   # Minimum confidence for poses
       augmentation=True          # Data augmentation during inference
   )

**Example:**

.. code-block:: python

   # ML-enhanced docking with uncertainty
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Access ML-specific results
   for pose in results.poses:
       print(f"Pose score: {pose.score:.3f}")
       print(f"Confidence: {pose.confidence:.3f}")
       print(f"Uncertainty: {pose.uncertainty:.3f}")
       print(f"Predicted IC50: {pose.predicted_ic50:.2e} nM")

**Advanced PandaML Configuration:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaml',
       ml_config={
           'model_path': 'models/protein_specific.pt',  # Custom model
           'feature_extraction': 'graph',              # graph, grid, sequence
           'attention_heads': 8,                       # For transformer models
           'hidden_size': 512,                         # Model hidden size
           'dropout': 0.1,                            # Dropout for uncertainty
           'temperature': 1.0,                        # Calibration temperature
           'batch_size': 16                           # Inference batch size
       }
   )

PandaCore Algorithm
------------------

PandaCore algorithm uses evolutionary algorithms for conformational sampling and optimization.

**Key Features:**
- Excellent conformational sampling
- Handles flexible receptors well
- Global optimization
- Customizable operators

**When to Use:**
- Induced fit docking
- Highly flexible systems
- Novel binding modes discovery
- Difficult docking problems

**Configuration:**

.. code-block:: python

   docker = PandaDock(
       engine='pandacore',
       population_size=150,        # GA population size
       generations=300,            # Number of generations
       mutation_rate=0.02,         # Mutation probability
       crossover_rate=0.8,         # Crossover probability
       selection='tournament',     # Selection method
       elitism=0.1,               # Elite fraction
       diversity_pressure=True     # Maintain population diversity
   )

**Example:**

.. code-block:: python

   # PandaCore docking with flexible receptor
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0],
       flexible_residues=['ARG123', 'TYR456']  # Flexible residues
   )
   
   # Access PandaCore-specific results
   print(f"Generations: {results.generations_completed}")
   print(f"Final population diversity: {results.final_diversity:.3f}")

**PandaCore Operators Customization:**

.. code-block:: python

   from pandadock.ga import custom_operators
   
   docker = PandaDock(
       engine='pandacore',
       ga_config={
           'mutation_operators': ['gaussian', 'uniform', 'adaptive'],
           'crossover_operators': ['uniform', 'single_point'],
           'selection_pressure': 2.0,
           'niching': True,           # Maintain diverse solutions
           'adaptive_parameters': True # Adapt PandaCore parameters during run
       }
   )

Hybrid Docking
--------------

Hybrid mode combines multiple approaches for optimal results.

**Key Features:**
- Best of all worlds
- Adaptive algorithm selection
- Multi-stage refinement
- Robust predictions

**Configuration:**

.. code-block:: python

   docker = PandaDock(
       engine='hybrid',
       hybrid_config={
           'stage1': 'pandaphysics',  # Initial sampling
           'stage2': 'pandaml',      # Refinement
           'stage3': 'pandaphysics', # Final optimization
           'consensus_scoring': True, # Use multiple scoring functions
           'pose_clustering': True   # Cluster similar poses
       }
   )

**Example:**

.. code-block:: python

   # Multi-stage hybrid docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Access stage-specific information
   for stage, stage_results in results.stages.items():
       print(f"Stage {stage}: {len(stage_results.poses)} poses")

Flexible Receptor Docking
--------------------------

Enable receptor flexibility in any docking mode.

**Configuration:**

.. code-block:: python

   # Specify flexible residues
   docker = PandaDock(
       engine='pandaml',  # Can be used with any engine
       flexible_residues=['ARG123', 'TYR456', 'ASP789'],
       flexibility_config={
           'backbone_flexibility': False,   # Keep backbone rigid
           'sidechain_flexibility': True,   # Allow sidechain movement
           'max_torsion_change': 30.0,     # Max angle change (degrees)
           'clash_tolerance': 0.5          # Flexibility clash tolerance
       }
   )

**Auto-Detection of Flexible Residues:**

.. code-block:: python

   docker = PandaDock(
       engine='pandacore',
       auto_flexible=True,
       flexibility_config={
           'detection_method': 'b_factor',  # b_factor, cavity, contact
           'flexibility_threshold': 30.0,   # B-factor threshold
           'max_flexible_residues': 5      # Limit number of flexible residues
       }
   )

Consensus Docking
-----------------

Use multiple docking modes and combine results.

.. code-block:: python

   from pandadock import ConsensusDocking
   
   # Define multiple docking protocols
   protocols = [
       {'engine': 'pandaphysics', 'exhaustiveness': 16},
       {'engine': 'pandaml', 'ensemble_size': 3},
       {'engine': 'pandacore', 'generations': 200}
   ]
   
   # Consensus docking
   consensus = ConsensusDocking(protocols=protocols)
   
   results = consensus.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Access consensus results
   print(f"Consensus score: {results.consensus_score:.3f}")
   print(f"Agreement between methods: {results.agreement:.3f}")

Virtual Screening Modes
-----------------------

Optimized configurations for different screening scenarios.

**High-Throughput Screening:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaphysics',
       screening_mode='hts',
       exhaustiveness=4,           # Fast screening
       num_poses=3,               # Few poses per ligand
       energy_range=2.0,          # Narrow energy window
       early_termination=True,    # Stop early for poor binders
       parallel_screening=True    # Process multiple ligands in parallel
   )

**Lead Optimization:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaml',
       screening_mode='lead_opt',
       exhaustiveness=16,          # Thorough sampling
       num_poses=10,              # More poses for analysis
       uncertainty=True,          # Include uncertainty estimates
       similarity_filtering=True   # Filter similar poses
   )

**Fragment Screening:**

.. code-block:: python

   docker = PandaDock(
       engine='pandacore',
       screening_mode='fragment',
       fragment_config={
           'min_heavy_atoms': 6,      # Minimum fragment size
           'max_heavy_atoms': 22,     # Maximum fragment size
           'growth_vectors': True,    # Identify growth vectors
           'hot_spots': True         # Find hot spot interactions
       }
   )

Performance Optimization
------------------------

**Memory-Efficient Mode:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaphysics',
       memory_efficient=True,
       memory_config={
           'max_memory_gb': 4,        # Memory limit
           'stream_processing': True, # Process poses one at a time
           'compress_intermediates': True,  # Compress intermediate data
           'minimal_storage': True    # Store only essential data
       }
   )

**GPU-Accelerated Mode:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaml',
       gpu_acceleration=True,
       gpu_config={
           'device': 'cuda:0',       # GPU device
           'mixed_precision': True,  # Use mixed precision
           'batch_size': 64,        # Large batch for GPU efficiency
           'memory_fraction': 0.8   # GPU memory fraction
       }
   )

Choosing the Right Mode
-----------------------

**Decision Tree:**

1. **Need highest accuracy?** → Use PandaML algorithm
2. **Have flexible receptor?** → Use PandaCore algorithm  
3. **Large-scale screening?** → Use PandaPhysics algorithm
4. **Production deployment?** → Use Hybrid mode
5. **Uncertain about target?** → Use Consensus docking

**Recommendations by Use Case:**

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Use Case
     - Recommended Mode
     - Configuration
   * - Virtual screening (>10K compounds)
     - PandaPhysics
     - ``exhaustiveness=4-8``
   * - Lead optimization (<100 compounds)
     - PandaML
     - ``ensemble_size=5, uncertainty=True``
   * - Novel target (no known binders)
     - PandaCore
     - ``generations=500, diversity=True``
   * - Allosteric sites
     - Consensus
     - ``[pandaphysics, pandaml, pandacore]``
   * - Fragment-based drug design
     - PandaCore
     - ``fragment_mode=True``
   * - Production pipeline
     - Hybrid
     - ``multi_stage=True``

Next Steps
----------

- Learn about :doc:`scoring_functions` for each mode
- See :doc:`../tutorials/ml_enhanced_docking` for ML mode details
- Check :doc:`../tutorials/virtual_screening` for screening workflows
- Explore :doc:`../examples/flexible_docking` for advanced flexibility