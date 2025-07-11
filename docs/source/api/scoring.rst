Scoring Module
==============

The scoring module provides various scoring functions and energy calculation methods for evaluating molecular poses.

Scoring Functions
-----------------

Core scoring function implementations.

.. automodule:: scoring.scoring_functions
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: scoring.scoring_functions.ScoringFunctions
   :members:
   :special-members: __init__
   :show-inheritance:

Energy Terms
^^^^^^^^^^^^

Individual energy calculation methods:

.. automethod:: scoring.scoring_functions.ScoringFunctions.calculate_vdw_energy

.. automethod:: scoring.scoring_functions.ScoringFunctions.calculate_electrostatic_energy

.. automethod:: scoring.scoring_functions.ScoringFunctions.calculate_hbond_energy

.. automethod:: scoring.scoring_functions.ScoringFunctions.calculate_hydrophobic_energy

.. automethod:: scoring.scoring_functions.ScoringFunctions.calculate_solvation_energy

Interaction Detection
^^^^^^^^^^^^^^^^^^^^^

Methods for detecting molecular interactions:

.. automethod:: scoring.scoring_functions.ScoringFunctions.find_hbond_interactions

.. automethod:: scoring.scoring_functions.ScoringFunctions.find_hydrophobic_interactions

.. automethod:: scoring.scoring_functions.ScoringFunctions.find_salt_bridge_interactions

PandaML Rescorer
----------------

PandaML machine learning-based rescoring module.

.. automodule:: scoring.ml_rescorer
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: scoring.ml_rescorer.MLRescorer
   :members:
   :special-members: __init__
   :show-inheritance:

Examples
--------

Basic Scoring
^^^^^^^^^^^^^

.. code-block:: python

   from scoring.scoring_functions import ScoringFunctions
   import numpy as np
   
   # Initialize scoring function
   scorer = ScoringFunctions()
   
   # Calculate energy for coordinates
   coordinates = np.random.rand(20, 3) * 10
   total_energy = scorer.calculate_total_energy(coordinates)
   print(f"Total energy: {total_energy:.2f} kcal/mol")
   
   # Calculate individual energy terms
   vdw_energy = scorer.calculate_vdw_energy(coordinates)
   hbond_energy = scorer.calculate_hbond_energy(coordinates)
   
   print(f"vdW energy: {vdw_energy:.2f} kcal/mol")
   print(f"H-bond energy: {hbond_energy:.2f} kcal/mol")

Interaction Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Find molecular interactions
   hbonds = scorer.find_hbond_interactions(ligand_coords, protein_coords)
   hydrophobic = scorer.find_hydrophobic_interactions(ligand_coords, protein_coords)
   
   print(f"Found {len(hbonds)} hydrogen bonds:")
   for hbond in hbonds:
       print(f"  Donor: {hbond['donor_atom']}, Acceptor: {hbond['acceptor_atom']}")
       print(f"  Distance: {hbond['distance']:.2f} Å, Energy: {hbond['energy']:.2f}")
   
   print(f"Found {len(hydrophobic)} hydrophobic interactions")

Custom Scoring Function
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from scoring.scoring_functions import ScoringFunctions
   
   class CustomScorer(ScoringFunctions):
       def __init__(self, custom_weights=None):
           super().__init__()
           if custom_weights:
               self.weights.update(custom_weights)
       
       def calculate_custom_term(self, coordinates):
           """Custom energy term implementation"""
           # Your custom scoring logic here
           return energy_value
       
       def calculate_total_energy(self, coordinates):
           """Override total energy calculation"""
           total = super().calculate_total_energy(coordinates)
           custom = self.calculate_custom_term(coordinates)
           return total + custom * self.weights.get('custom', 1.0)
   
   # Use custom scorer
   custom_scorer = CustomScorer({'custom': 0.5})
   energy = custom_scorer.calculate_total_energy(coordinates)

PandaML Rescoring
^^^^^^^^^^^^^^^^^

.. code-block:: python

   from scoring.ml_rescorer import MLRescorer
   
   # Initialize PandaML rescorer
   rescorer = MLRescorer(model_path='models/pandaml_rescorer.pt')
   
   # Rescore poses
   poses = [pose1, pose2, pose3]  # List of poses
   rescored_poses = rescorer.rescore_poses(poses, receptor_features)
   
   for original, rescored in zip(poses, rescored_poses):
       print(f"Original score: {original.score:.3f}")
       print(f"PandaML score: {rescored.ml_score:.3f}")
       print(f"Confidence: {rescored.confidence:.3f}")

Feature Extraction
^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Extract features for PandaML scoring
   features = rescorer.extract_features(pose, receptor)
   
   print(f"Feature vector shape: {features.shape}")
   print(f"Feature names: {rescorer.feature_names}")
   
   # Get feature importance
   importance = rescorer.get_feature_importance()
   for name, imp in zip(rescorer.feature_names, importance):
       print(f"{name}: {imp:.4f}")

Scoring Configuration
---------------------

Configure scoring function parameters:

.. code-block:: python

   # Custom scoring weights
   custom_weights = {
       'vdw': 1.2,
       'electrostatic': 0.8,
       'hbond': 1.5,
       'hydrophobic': 1.0,
       'solvation': 0.9,
       'entropy': 1.1
   }
   
   scorer = ScoringFunctions(weights=custom_weights)
   
   # Custom interaction parameters
   scorer.hbond_cutoff = 3.0  # Å
   scorer.hydrophobic_cutoff = 4.0  # Å
   scorer.vdw_cutoff = 8.0  # Å

Energy Decomposition
--------------------

Analyze energy contributions:

.. code-block:: python

   def analyze_energy_components(scorer, coordinates):
       """Detailed energy analysis"""
       components = {
           'vdw': scorer.calculate_vdw_energy(coordinates),
           'electrostatic': scorer.calculate_electrostatic_energy(coordinates),
           'hbond': scorer.calculate_hbond_energy(coordinates),
           'hydrophobic': scorer.calculate_hydrophobic_energy(coordinates),
           'solvation': scorer.calculate_solvation_energy(coordinates),
           'entropy': scorer.calculate_entropy_penalty(coordinates)
       }
       
       total = sum(components.values())
       
       print("Energy Decomposition:")
       print("-" * 30)
       for term, energy in components.items():
           percentage = (energy / total) * 100 if total != 0 else 0
           print(f"{term:12s}: {energy:8.2f} kcal/mol ({percentage:5.1f}%)")
       print("-" * 30)
       print(f"{'Total':12s}: {total:8.2f} kcal/mol")
       
       return components

PandaCore Scoring
^^^^^^^^^^^^^^^^^

Use PandaCore scoring algorithm:

.. code-block:: python

   # Calculate PandaCore score
   pandacore_score = scorer.calculate_pandacore_score(coordinates)
   print(f"PandaCore score: {pandacore_score:.3f}")
   
   # Get detailed PandaCore terms
   pandacore_terms = scorer.get_pandacore_terms(coordinates)
   for term, value in pandacore_terms.items():
       print(f"{term}: {value:.4f}")

Performance Optimization
------------------------

Tips for efficient scoring:

.. code-block:: python

   # Pre-calculate distance matrices for multiple evaluations
   from scoring.scoring_functions import distance_matrix
   
   distances = distance_matrix(coordinates, coordinates)
   
   # Use distances for multiple energy calculations
   vdw_energy = scorer._calculate_vdw_from_distances(distances)
   hbond_energy = scorer._calculate_hbond_from_distances(distances)
   
   # Batch scoring for multiple poses
   energies = scorer.batch_score([pose1, pose2, pose3])

Error Handling
--------------

.. code-block:: python

   from scoring.scoring_functions import ScoringError
   
   try:
       energy = scorer.calculate_total_energy(coordinates)
   except ScoringError as e:
       print(f"Scoring failed: {e}")
   except Exception as e:
       print(f"Unexpected error: {e}")

Best Practices
--------------

1. **Choose appropriate cutoffs:**
   - Use larger cutoffs for more accurate but slower calculations
   - Use smaller cutoffs for faster screening

2. **Weight tuning:**
   - Adjust scoring weights based on your system
   - Use experimental data to optimize weights

3. **Memory efficiency:**
   - Clear large distance matrices when not needed
   - Use batch processing for multiple poses

4. **Validation:**
   - Always validate scoring results against known data
   - Check energy components for reasonableness