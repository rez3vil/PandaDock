Configuration
=============

PandaDock provides extensive configuration options to customize docking behavior, performance, and output. This guide covers all available configuration methods.

Configuration Methods
---------------------

PandaDock can be configured through:

1. **Python API** - Direct parameter passing
2. **Configuration Files** - YAML/JSON config files  
3. **Environment Variables** - System-level settings
4. **Command Line Arguments** - CLI parameter overrides

Python API Configuration
------------------------

Configure PandaDock directly in your Python code:

.. code-block:: python

   from pandadock import PandaDock
   
   # Basic configuration
   docker = PandaDock(
       engine='pandaml',
       scoring='pandaml',
       exhaustiveness=16,
       num_poses=10
   )
   
   # Advanced configuration
   docker = PandaDock(
       engine='pandaphysics',
       scoring='pandaphysics',
       exhaustiveness=8,
       num_poses=20,
       energy_range=3.0,
       seed=42,
       cpu_threads=8,
       gpu_acceleration=True,
       max_memory_gb=16,
       temp_dir='/tmp/pandadock',
       log_level='INFO'
   )

Configuration Files
-------------------

Use YAML configuration files for complex setups:

**config.yaml:**

.. code-block:: yaml

   # PandaDock Configuration File
   
   # Engine settings
   engine:
     type: "pandaml"               # pandacore, pandaml, pandaphysics
     exhaustiveness: 16            # Search thoroughness (1-32)
     num_poses: 10                # Number of output poses
     energy_range: 3.0            # Energy range for poses (kcal/mol)
     seed: 42                     # Random seed for reproducibility
   
   # Scoring function settings
   scoring:
     function: "pandaml"          # pandacore, pandaml, pandaphysics
     weights:
       vdw: 1.0                  # van der Waals weight
       electrostatic: 1.0        # Electrostatic weight
       hbond: 1.0               # Hydrogen bond weight
       hydrophobic: 1.0         # Hydrophobic weight
       solvation: 1.0           # Solvation weight
       entropy: 1.0             # Entropy penalty weight
   
   # Performance settings
   performance:
     cpu_threads: 8              # Number of CPU threads
     gpu_acceleration: true      # Enable GPU acceleration
     max_memory_gb: 16          # Maximum memory usage
     batch_size: 100            # Batch size for screening
   
   # Input/Output settings
   io:
     temp_dir: "/tmp/pandadock"  # Temporary directory
     keep_temp: false           # Keep temporary files
     output_formats: ["sdf", "pdb", "mol2"]  # Output formats
   
   # Logging settings
   logging:
     level: "INFO"              # DEBUG, INFO, WARNING, ERROR
     file: "pandadock.log"      # Log file path
     console: true              # Log to console
   
   # Validation settings
   validation:
     check_inputs: true         # Validate input files
     energy_threshold: 1000.0   # Maximum allowed energy
     rmsd_threshold: 2.0        # RMSD cutoff for clustering

Load configuration from file:

.. code-block:: python

   # Load from YAML
   docker = PandaDock.from_config('config.yaml')
   
   # Load from JSON
   docker = PandaDock.from_config('config.json')
   
   # Override specific parameters
   docker = PandaDock.from_config('config.yaml', exhaustiveness=32)

Environment Variables
---------------------

Set global defaults using environment variables:

.. code-block:: bash

   # Engine settings
   export PANDADOCK_ENGINE=pandaml
   export PANDADOCK_SCORING=pandaml
   export PANDADOCK_EXHAUSTIVENESS=16
   export PANDADOCK_NUM_POSES=10
   
   # Performance settings
   export PANDADOCK_CPU_THREADS=8
   export PANDADOCK_GPU_ACCELERATION=true
   export PANDADOCK_MAX_MEMORY_GB=16
   
   # Paths
   export PANDADOCK_TEMP_DIR=/tmp/pandadock
   export PANDADOCK_DATA_DIR=/data/pandadock
   
   # Logging
   export PANDADOCK_LOG_LEVEL=INFO
   export PANDADOCK_LOG_FILE=pandadock.log

These can be overridden in Python:

.. code-block:: python

   import os
   
   # Check environment settings
   print(f"Engine: {os.getenv('PANDADOCK_ENGINE', 'pandacore')}")
   
   # Override environment variable
   docker = PandaDock(engine='pandaphysics')  # Overrides PANDADOCK_ENGINE

Docking Engine Configuration
----------------------------

**Physics Engine:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaphysics',
       physics_config={
           'force_field': 'amber',      # amber, charmm, gaff
           'integrator': 'verlet',      # verlet, langevin
           'timestep': 0.002,           # Integration timestep (ps)
           'temperature': 300.0,        # Simulation temperature (K)
           'pressure': 1.0,             # Pressure (bar)
           'constraint_tolerance': 1e-6  # Constraint tolerance
       }
   )

**Machine Learning Engine:**

.. code-block:: python

   docker = PandaDock(
       engine='pandaml',
       ml_config={
           'model': 'transformer',      # transformer, cnn, graph
           'model_path': 'models/best.pt',  # Path to trained model
           'batch_size': 32,           # Inference batch size
           'uncertainty': True,         # Calculate uncertainty
           'ensemble_size': 5,         # Number of ensemble models
           'dropout_rate': 0.1         # Dropout for uncertainty
       }
   )

**Genetic Algorithm Engine:**

.. code-block:: python

   docker = PandaDock(
       engine='pandacore',
       ga_config={
           'population_size': 150,      # GA population size
           'generations': 300,          # Number of generations
           'mutation_rate': 0.02,       # Mutation probability
           'crossover_rate': 0.8,       # Crossover probability
           'selection': 'tournament',   # tournament, roulette
           'elitism': 0.1              # Fraction of elite individuals
       }
   )

Scoring Function Configuration
------------------------------

**AutoDock Vina Scoring:**

.. code-block:: python

   docker = PandaDock(
       scoring='pandacore',
       vina_config={
           'weights': {
               'gauss1': -0.035579,
               'gauss2': -0.005156,
               'repulsion': 0.840245,
               'hydrophobic': -0.035069,
               'hydrogen': -0.587439,
               'rot': 0.05846
           },
           'atom_types': ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'],
           'cutoff': 8.0               # Interaction cutoff (Å)
       }
   )

**ChemPLP Scoring:**

.. code-block:: python

   docker = PandaDock(
       scoring='pandaml',
       chemplp_config={
           'clash_threshold': 0.6,     # Clash detection threshold
           'metal_bonus': 1.0,         # Metal interaction bonus
           'hbond_penalty': 0.1,       # H-bond angle penalty
           'rotatable_penalty': 0.05   # Rotatable bond penalty
       }
   )

**Custom Scoring:**

.. code-block:: python

   def custom_scoring_function(pose, receptor, ligand):
       """Custom scoring function implementation"""
       # Your scoring logic here
       return score
   
   docker = PandaDock(
       scoring='custom',
       custom_scorer=custom_scoring_function,
       custom_config={
           'parameter1': 1.0,
           'parameter2': 'value'
       }
   )

Performance Tuning
------------------

**CPU Optimization:**

.. code-block:: python

   docker = PandaDock(
       cpu_threads=16,              # Use all available cores
       parallel_docking=True,       # Parallel ligand processing
       memory_efficient=True,       # Reduce memory usage
       optimization_level=3         # Compiler optimization level
   )

**GPU Acceleration:**

.. code-block:: python

   docker = PandaDock(
       gpu_acceleration=True,
       gpu_device='cuda:0',         # Specific GPU device
       gpu_memory_fraction=0.8,     # Fraction of GPU memory to use
       mixed_precision=True         # Use mixed precision for speed
   )

**Memory Management:**

.. code-block:: python

   docker = PandaDock(
       max_memory_gb=32,           # Maximum memory limit
       memory_pool_size=1024,      # Memory pool size (MB)
       garbage_collection=True,    # Enable aggressive GC
       cache_size=500             # Cache size (MB)
   )

Advanced Configuration
----------------------

**Flexible Receptor Docking:**

.. code-block:: python

   docker = PandaDock(
       flexible_residues=['ARG123', 'TYR456', 'ASP789'],
       flexibility_config={
           'backbone_flexibility': False,  # Allow backbone movement
           'sidechain_flexibility': True,  # Allow sidechain movement
           'torsion_amplitude': 30.0,      # Max torsion angle change (degrees)
           'clash_tolerance': 0.5          # Clash tolerance for flexibility
       }
   )

**Binding Site Configuration:**

.. code-block:: python

   docker = PandaDock(
       binding_site_config={
           'auto_detect': True,         # Auto-detect binding site
           'cavity_detection': 'fpocket', # fpocket, caver, sitemap
           'min_cavity_volume': 200.0,  # Minimum cavity volume (Å³)
           'probe_radius': 1.4,         # Probe radius for cavity detection
           'grid_spacing': 0.375        # Grid spacing (Å)
       }
   )

**Output Configuration:**

.. code-block:: python

   docker = PandaDock(
       output_config={
           'formats': ['sdf', 'pdb', 'mol2'],  # Output formats
           'include_receptor': True,            # Include receptor in output
           'compress_output': True,             # Compress output files
           'report_format': 'html',             # html, pdf, json
           'detailed_analysis': True,           # Include detailed analysis
           'interaction_plots': True           # Generate interaction plots
       }
   )

Configuration Validation
------------------------

Validate configuration before use:

.. code-block:: python

   from pandadock.config import validate_config
   
   config = {
       'engine': 'pandaml',
       'exhaustiveness': 16,
       'num_poses': 10
   }
   
   # Validate configuration
   is_valid, errors = validate_config(config)
   
   if not is_valid:
       for error in errors:
           print(f"Configuration error: {error}")
   else:
       docker = PandaDock(**config)

Default Configuration
---------------------

View current default configuration:

.. code-block:: python

   from pandadock.config import get_default_config
   
   # Get default configuration
   defaults = get_default_config()
   print(defaults)
   
   # Save default configuration to file
   with open('default_config.yaml', 'w') as f:
       yaml.dump(defaults, f)

Configuration Profiles
----------------------

Use predefined configuration profiles:

.. code-block:: python

   # Fast screening profile
   docker = PandaDock.from_profile('fast_screening')
   
   # High accuracy profile  
   docker = PandaDock.from_profile('high_accuracy')
   
   # GPU optimized profile
   docker = PandaDock.from_profile('gpu_optimized')
   
   # Custom profile with overrides
   docker = PandaDock.from_profile('fast_screening', num_poses=20)

Available profiles:
- ``fast_screening`` - Optimized for virtual screening
- ``high_accuracy`` - Maximum accuracy settings
- ``gpu_optimized`` - GPU acceleration settings
- ``memory_efficient`` - Low memory usage
- ``flexible_receptor`` - Induced fit docking

Configuration Examples
----------------------

**Virtual Screening Setup:**

.. code-block:: yaml

   engine:
     type: "pandaphysics"
     exhaustiveness: 8
     num_poses: 5
     energy_range: 2.0
   
   performance:
     cpu_threads: 16
     parallel_docking: true
     batch_size: 1000
   
   output:
     formats: ["sdf"]
     compress_output: true

**High-Accuracy Docking:**

.. code-block:: yaml

   engine:
     type: "pandaml"
     exhaustiveness: 32
     num_poses: 20
     energy_range: 5.0
   
   ml_config:
     ensemble_size: 10
     uncertainty: true
   
   output:
     detailed_analysis: true
     interaction_plots: true

Next Steps
----------

- Learn about specific :doc:`docking_modes`
- Explore :doc:`scoring_functions` customization
- See :doc:`../tutorials/custom_scoring` for advanced scoring
- Check :doc:`../examples/comparative_analysis` for benchmarking