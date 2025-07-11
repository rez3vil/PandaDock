API Reference
=============

This section provides detailed documentation for all PandaDock modules, classes, and functions.

.. toctree::
   :maxdepth: 2
   :caption: API Contents:

   docking
   scoring
   utils
   io

Module Overview
---------------

The PandaDock API is organized into several key modules:

**Core Modules**

* :doc:`docking` - Docking engines and pose management
* :doc:`scoring` - Scoring functions and energy calculations
* :doc:`utils` - Utility functions and helper classes
* :doc:`io` - Input/output and file format handling

Module Descriptions
-------------------

**Docking Module** (:doc:`docking`)
  Contains the core docking engines including PandaPhysics, PandaML, and PandaCore implementations. Also includes the base classes for poses and docking results.

**Scoring Module** (:doc:`scoring`)
  Provides various scoring functions for evaluating molecular poses, including classical force field terms, machine learning rescoring, and interaction detection methods.

**Utilities Module** (:doc:`utils`)
  Contains helper functions for IC50 calculations, mathematical operations, rotamer libraries, and other computational tools used throughout PandaDock.

**Input/Output Module** (:doc:`io`)
  Handles reading and writing of molecular structure files, ligand preparation, and format conversions between different molecular file formats.

API Usage Patterns
-------------------

**Basic Docking Workflow**

.. code-block:: python

   from pandadock import PandaDock
   from pandadock.docking.base_engine import Pose
   
   # Initialize docking engine
   docker = PandaDock(engine='pandaphysics')
   
   # Perform docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[0, 0, 0],
       size=[20, 20, 20]
   )
   
   # Access results
   best_pose = results.best_pose
   all_poses = results.poses

**Custom Scoring Function**

.. code-block:: python

   from pandadock.scoring.scoring_functions import ScoringFunctions
   
   class CustomScorer(ScoringFunctions):
       def calculate_custom_energy(self, coordinates):
           # Custom energy implementation
           return energy_value
   
   # Use custom scorer
   scorer = CustomScorer()
   energy = scorer.calculate_total_energy(coordinates)

**Ligand Preparation**

.. code-block:: python

   from pandadock.io.ligand_preparer import LigandPreparer
   
   preparer = LigandPreparer()
   ligand_data = preparer.prepare_ligand('compound.sdf')
   
   # Access prepared data
   coordinates = ligand_data['coordinates']
   atom_types = ligand_data['atom_types']
   bonds = ligand_data['bonds']

**IC50 Calculations**

.. code-block:: python

   from pandadock.utils.ic50_calculator import IC50Calculator
   
   calc = IC50Calculator()
   ic50 = calc.delta_g_to_ic50(binding_free_energy)
   ligand_efficiency = calc.calculate_ligand_efficiency(
       binding_free_energy, 
       num_heavy_atoms
   )

Class Hierarchies
-----------------

**Docking Engines**

.. code-block:: text

   DockingEngine (base class)
   ├── PandaPhysicsEngine
   ├── PandaMLEngine
   ├── PandaCoreEngine
   └── FlexibleDocking

**Scoring Functions**

.. code-block:: text

   ScoringFunctions (base class)
   ├── PandaCoreScoringFunction
   ├── PandaMLScoringFunction
   └── CustomScoringFunction

**Pose Classes**

.. code-block:: text

   Pose (base class)
   ├── FlexiblePose
   ├── EnsemblePose
   └── CovalentPose

Error Handling
--------------

PandaDock defines several custom exception types:

.. code-block:: python

   from pandadock.docking.base_engine import DockingError
   from pandadock.scoring.scoring_functions import ScoringError
   from pandadock.io.ligand_preparer import PreparationError
   
   try:
       results = docker.dock(receptor, ligand)
   except DockingError as e:
       print(f"Docking failed: {e}")
   except ScoringError as e:
       print(f"Scoring failed: {e}")
   except PreparationError as e:
       print(f"Ligand preparation failed: {e}")

Type Hints and Annotations
--------------------------

PandaDock uses comprehensive type hints for better IDE support:

.. code-block:: python

   from typing import List, Dict, Optional, Tuple
   from pandadock.docking.base_engine import Pose, DockingResults
   
   def analyze_poses(poses: List[Pose]) -> Dict[str, float]:
       """Analyze a list of poses and return statistics."""
       return {
           'mean_score': sum(p.score for p in poses) / len(poses),
           'best_score': min(p.score for p in poses),
           'score_std': calculate_std([p.score for p in poses])
       }

Configuration Objects
---------------------

Many PandaDock functions accept configuration dictionaries:

.. code-block:: python

   # Docking configuration
   docking_config = {
       'exhaustiveness': 16,
       'num_poses': 10,
       'energy_range': 3.0,
       'algorithm': 'pandaphysics'
   }
   
   # PandaML configuration
   pandaml_config = {
       'model_type': 'transformer',
       'ensemble_size': 5,
       'uncertainty': True,
       'batch_size': 32
   }
   
   # Flexibility configuration
   flexibility_config = {
       'backbone_flexibility': True,
       'max_torsion_change': 30.0,
       'clash_tolerance': 0.8
   }

Version Compatibility
---------------------

API compatibility information:

* **Python**: Requires Python 3.8+
* **Dependencies**: See requirements.txt for version constraints
* **Backward Compatibility**: API changes are documented in the changelog
* **Deprecation Policy**: Deprecated features are supported for 2 major versions

Performance Considerations
--------------------------

For optimal performance when using the API:

1. **Reuse Objects**: Create scoring functions and preparers once, use multiple times
2. **Batch Operations**: Use batch methods when processing multiple ligands
3. **Memory Management**: Clear large objects when no longer needed
4. **GPU Acceleration**: Enable GPU support for ML models when available

Contributing to the API
------------------------

When contributing new API features:

1. Follow existing naming conventions
2. Include comprehensive docstrings
3. Add type hints for all functions
4. Write unit tests for new functionality
5. Update documentation and examples

For detailed contribution guidelines, see the development documentation.