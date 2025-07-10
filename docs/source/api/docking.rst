Docking Module
==============

The docking module contains the core docking engines and pose management functionality.

Base Classes
------------

.. automodule:: docking.base_engine
   :members:
   :undoc-members:
   :show-inheritance:

Pose Class
^^^^^^^^^^

.. autoclass:: docking.base_engine.Pose
   :members:
   :special-members: __init__
   :show-inheritance:

DockingEngine
^^^^^^^^^^^^^

.. autoclass:: docking.base_engine.DockingEngine
   :members:
   :special-members: __init__
   :show-inheritance:

Physics Engine
--------------

Classical force field-based docking engine.

.. automodule:: docking.physics_engine
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking.physics_engine.PhysicsEngine
   :members:
   :special-members: __init__
   :show-inheritance:

Machine Learning Engine
-----------------------

Deep learning-enhanced docking engine.

.. automodule:: docking.ml_engine
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking.ml_engine.MLEngine
   :members:
   :special-members: __init__
   :show-inheritance:

Genetic Algorithm Engine
------------------------

Evolutionary algorithm-based docking engine.

.. automodule:: docking.ga_engine
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking.ga_engine.GAEngine
   :members:
   :special-members: __init__
   :show-inheritance:

Flexible Docking
-----------------

Enhanced docking with receptor flexibility.

.. automodule:: docking.flexible_docking
   :members:
   :undoc-members:
   :show-inheritance:

Pose Filtering
--------------

Post-docking pose filtering and clustering.

.. automodule:: docking.pose_filtering
   :members:
   :undoc-members:
   :show-inheritance:

Examples
--------

Basic Docking
^^^^^^^^^^^^^

.. code-block:: python

   from docking.physics_engine import PhysicsEngine
   from docking.base_engine import Pose
   import numpy as np
   
   # Initialize engine
   engine = PhysicsEngine(config)
   
   # Create a pose
   coordinates = np.random.rand(20, 3) * 10
   pose = Pose(coordinates=coordinates, score=0.5, energy=-5.2)
   
   # Score the pose
   score = engine.score(pose)
   print(f"Pose score: {score}")

Advanced Configuration
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from docking.ml_engine import MLEngine
   
   # Configure ML engine
   config = {
       'model_path': 'models/best_model.pt',
       'batch_size': 32,
       'uncertainty': True
   }
   
   engine = MLEngine(config)
   
   # Dock with uncertainty quantification
   results = engine.dock(ligand, receptor)
   
   for pose in results.poses:
       print(f"Score: {pose.score:.3f} ± {pose.uncertainty:.3f}")

Flexible Receptor Docking
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from docking.flexible_docking import FlexibleDocking
   
   # Configure flexible docking
   flexible_dock = FlexibleDocking(
       flexible_residues=['ARG123', 'TYR456'],
       backbone_flexibility=False,
       sidechain_flexibility=True
   )
   
   # Perform flexible docking
   results = flexible_dock.dock(ligand, receptor, binding_site)

Error Handling
--------------

.. code-block:: python

   from docking.base_engine import DockingError
   
   try:
       results = engine.dock(ligand, receptor)
   except DockingError as e:
       print(f"Docking failed: {e}")
       # Handle error appropriately

Performance Tips
----------------

1. **Use appropriate engine for your needs:**
   - PhysicsEngine: Fast, general-purpose
   - MLEngine: High accuracy, slower
   - GAEngine: Flexible receptors, thorough sampling

2. **Optimize parameters:**
   - Increase exhaustiveness for better accuracy
   - Reduce num_poses for faster screening
   - Use GPU acceleration when available

3. **Memory management:**
   - Process large libraries in batches
   - Clear poses from memory when not needed
   - Use memory-efficient modes for large-scale screening