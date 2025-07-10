Metal Docking Guide
==================

PandaDock provides advanced support for molecular docking to metalloproteins, which contain metal ions as essential components of their active sites. This guide covers the specialized features and workflows for metal docking.

Introduction to Metal Docking
-----------------------------

Metalloproteins represent approximately 30% of all proteins and play crucial roles in:

- **Catalysis**: Enzymes like carbonic anhydrase (Zn), cytochrome c oxidase (Cu, Fe)
- **Electron Transfer**: Iron-sulfur clusters, copper centers
- **Structural Support**: Zinc fingers, calcium-binding proteins
- **Oxygen Transport**: Hemoglobin (Fe), hemocyanin (Cu)
- **DNA Binding**: Zinc finger transcription factors

Metal docking presents unique challenges:

1. **Coordination Geometry**: Metals have specific geometric preferences
2. **Electrostatic Effects**: Strong charge-charge interactions
3. **Constraint Satisfaction**: Rigid geometric requirements
4. **Multiple Binding Modes**: Direct coordination vs. allosteric effects

Key Features of PandaDock Metal Docking
---------------------------------------

**Automatic Metal Detection**
  - Identifies metal centers from PDB structures
  - Analyzes coordination environments
  - Determines coordination geometry and numbers

**Specialized Scoring Functions**
  - Metal-ligand coordination energies
  - Geometric constraint penalties
  - Electrostatic interactions with metal charges
  - Solvation effects for metal centers

**Constraint-Based Optimization**
  - Distance constraints for metal-ligand bonds
  - Angular constraints for coordination geometry
  - Coordination number requirements
  - Oxidation state considerations

**Multi-Metal Support**
  - Binuclear and polynuclear metal clusters
  - Metal-metal cooperative effects
  - Bridging ligand interactions

Quick Start
-----------

Basic metal docking example:

.. code-block:: python

   from pandadock import PandaDock
   from pandadock.docking import MetalDockingEngine, MetalDockingConfig
   
   # Configure metal docking
   metal_config = MetalDockingConfig(
       detect_metals_automatically=True,
       enforce_geometric_constraints=True,
       use_metal_scoring=True,
       coordination_focused_sampling=True
   )
   
   # Initialize metal docking engine
   docker = MetalDockingEngine(config, metal_config)
   
   # Perform metal docking
   poses = docker.dock(
       protein_file='metalloprotein.pdb',
       ligand_file='metal_chelator.sdf'
   )
   
   # Analyze results
   report = docker.get_metal_docking_report(poses)
   docker.save_metal_poses(poses, 'metal_results/')

Supported Metal Types
--------------------

PandaDock supports all biologically relevant metals:

**Transition Metals:**
  - Zinc (Zn²⁺) - Most common, tetrahedral coordination
  - Iron (Fe²⁺/Fe³⁺) - Octahedral, square pyramidal
  - Copper (Cu²⁺) - Square planar, tetrahedral
  - Manganese (Mn²⁺) - Octahedral
  - Nickel (Ni²⁺) - Square planar, octahedral
  - Cobalt (Co²⁺) - Tetrahedral, octahedral

**Alkaline Earth Metals:**
  - Magnesium (Mg²⁺) - Octahedral, flexible
  - Calcium (Ca²⁺) - 6-8 coordination, very flexible

**Other Metals:**
  - Molybdenum (Mo) - Variable coordination
  - Tungsten (W) - Similar to Mo
  - Vanadium (V) - Multiple oxidation states

Coordination Geometries
-----------------------

The system recognizes and handles various coordination geometries:

**Linear (CN=2)**
  - 180° angle between ligands
  - Common in Cu(I) complexes

**Trigonal Planar (CN=3)**
  - 120° angles between ligands
  - Some Cu(I) and Au(I) complexes

**Tetrahedral (CN=4)**
  - 109.5° angles, most common for Zn²⁺
  - Also common for Cu(I), Ni²⁺ (high spin)

**Square Planar (CN=4)**
  - 90° and 180° angles
  - Common for Cu²⁺, Ni²⁺ (low spin), Pt²⁺

**Octahedral (CN=6)**
  - 90° and 180° angles
  - Most common for Fe²⁺/Fe³⁺, Mg²⁺

**Trigonal Bipyramidal (CN=5)**
  - 90°, 120°, 180° angles
  - Some Fe complexes

Configuration Options
---------------------

MetalDockingConfig Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   metal_config = MetalDockingConfig(
       # Metal detection
       detect_metals_automatically=True,    # Auto-detect metals in structure
       metal_detection_distance=3.5,       # Max distance for coordination
       
       # Coordination requirements
       require_metal_coordination=True,     # Require ligand coordination
       min_coordinating_atoms=1,           # Minimum coordinating atoms
       max_coordinating_atoms=6,           # Maximum coordinating atoms
       
       # Geometric constraints
       enforce_geometric_constraints=True,  # Apply geometry constraints
       geometric_constraint_weight=2.0,     # Weight for violations
       distance_tolerance=0.3,             # Distance tolerance (Å)
       angle_tolerance=15.0,               # Angle tolerance (degrees)
       
       # Sampling parameters
       coordination_focused_sampling=True,  # Focus sampling near metals
       metal_focused_exhaustiveness=20,     # Enhanced sampling factor
       coordination_cone_angle=30.0,       # Sampling cone (degrees)
       
       # Scoring parameters
       use_metal_scoring=True,              # Use metal-specific scoring
       metal_scoring_weight=1.0,            # Weight for metal terms
       
       # Pose filtering
       filter_non_coordinating_poses=True,  # Remove non-coordinating poses
       min_coordination_score=0.3,         # Minimum coordination quality
       max_geometric_violations=2           # Max allowed violations
   )

Metal-Specific Scoring
----------------------

The metal scoring function includes several specialized terms:

**Coordination Energy**
  - Direct metal-ligand bond formation
  - Element-specific interaction strengths
  - Distance-dependent energy profiles

**Geometric Penalties**
  - Violations of ideal coordination angles
  - Deviations from expected distances
  - Coordination number mismatches

**Electrostatic Interactions**
  - Metal charge interactions with ligand charges
  - Distance-dependent dielectric effects
  - Polarization contributions

**Solvation Effects**
  - Metal desolvation penalties
  - Coordination sphere solvation
  - Entropy effects from coordination

**Example Scoring Configuration:**

.. code-block:: python

   from pandadock.scoring import MetalScoringParameters
   
   # Custom scoring parameters
   metal_params = MetalScoringParameters(
       coordination_strength={
           'N': -4.0,    # Strong coordination
           'O': -3.5,    # Moderate-strong
           'S': -3.0,    # Moderate
           'P': -2.5     # Weaker
       },
       geometric_penalty_weight=2.0,
       distance_tolerance=0.2,
       angle_tolerance=10.0
   )

Working with Constraints
------------------------

PandaDock provides flexible constraint handling for metal coordination:

**Automatic Constraints**
  - Generated based on metal type and geometry
  - Standard distance and angle requirements
  - Coordination number constraints

**Custom Constraints**
  - User-defined geometric requirements
  - Specific distance or angle constraints
  - Composite constraints for complex geometries

**Constraint Presets:**

.. code-block:: python

   from pandadock.utils.metal_constraints import ConstraintSetPresets
   
   # Strict constraints for precise docking
   strict_manager = ConstraintSetPresets.create_strict_coordination_constraints(
       metal_centers
   )
   
   # Flexible constraints for screening
   flexible_manager = ConstraintSetPresets.create_flexible_coordination_constraints(
       metal_centers
   )
   
   # Distance-only constraints
   distance_manager = ConstraintSetPresets.create_distance_only_constraints(
       metal_centers
   )

Practical Examples
------------------

Zinc Finger Docking
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Zinc finger proteins typically have Zn²⁺ coordinated by 2 Cys + 2 His
   metal_config = MetalDockingConfig(
       enforce_geometric_constraints=True,
       geometric_constraint_weight=3.0,    # Strict for Zn fingers
       require_metal_coordination=True,
       min_coordinating_atoms=1,
       distance_tolerance=0.2,             # Tight tolerance
       angle_tolerance=10.0
   )

Iron Enzyme Docking
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Iron enzymes often have flexible coordination
   metal_config = MetalDockingConfig(
       enforce_geometric_constraints=True,
       geometric_constraint_weight=1.5,    # More flexible for Fe
       max_coordinating_atoms=6,           # Fe can be 6-coordinate
       distance_tolerance=0.3,
       angle_tolerance=15.0,
       coordination_focused_sampling=True,
       metal_focused_exhaustiveness=30     # More sampling for Fe
   )

Calcium Binding Protein
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Calcium has very flexible coordination
   metal_config = MetalDockingConfig(
       enforce_geometric_constraints=False, # Ca is very flexible
       geometric_constraint_weight=0.5,     # Low penalty weight
       max_coordinating_atoms=8,            # Ca can coordinate many atoms
       distance_tolerance=0.5,              # Large tolerance
       angle_tolerance=25.0,
       filter_non_coordinating_poses=False  # Allow non-coordinating
   )

Multi-Metal Systems
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # For proteins with multiple metal centers
   metal_config = MetalDockingConfig(
       detect_metals_automatically=True,
       enforce_geometric_constraints=True,
       coordination_focused_sampling=True,
       metal_focused_exhaustiveness=40,     # More sampling needed
       require_metal_coordination=False,    # May coordinate to one metal
       min_coordinating_atoms=1
   )

Analysis and Interpretation
---------------------------

Metal Pose Analysis
^^^^^^^^^^^^^^^^^^^

Metal poses contain additional information:

.. code-block:: python

   # Analyze metal-specific properties
   for pose in poses:
       print(f"Pose ID: {pose.pose_id}")
       print(f"Coordinating atoms: {pose.coordinating_atoms}")
       print(f"Coordination distances: {pose.coordination_distances}")
       print(f"Metal interactions: {len(pose.metal_interactions)}")
       print(f"Coordination quality: {pose.coordination_quality}")
       print(f"Geometric violations: {len(pose.geometric_violations)}")

Detailed Reports
^^^^^^^^^^^^^^^^

.. code-block:: python

   # Generate comprehensive metal docking report
   report = engine.get_metal_docking_report(poses)
   
   # Access specific information
   summary = report['docking_summary']
   metal_centers = report['metal_centers']
   best_pose_analysis = report['best_pose']
   
   print(f"Poses with coordination: {summary['poses_with_coordination']}")
   print(f"Average coordination score: {summary['average_coordination_score']}")
   
   # Analyze individual metal centers
   for metal in metal_centers:
       print(f"Metal {metal['metal_id']}: {metal['metal_type']}")
       print(f"  Poses coordinating: {metal['poses_coordinating']}")

Interaction Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Analyze metal-ligand interactions
   best_pose = poses[0]
   
   for interaction in best_pose.metal_interactions:
       metal_type = interaction['metal_type']
       ligand_element = interaction['ligand_element']
       distance = interaction['distance']
       energy = interaction['energy']
       
       print(f"{ligand_element}-{metal_type} interaction:")
       print(f"  Distance: {distance:.2f} Å")
       print(f"  Energy: {energy:.2f} kcal/mol")
       print(f"  Type: {interaction['subtype']}")

Best Practices
--------------

**Protein Preparation**
  1. Ensure metal ions are properly placed in PDB file
  2. Include all coordinating residues
  3. Check for proper protonation states
  4. Verify metal oxidation states

**Ligand Design**
  - Include potential coordinating atoms (N, O, S, P)
  - Consider chelation effects
  - Account for coordination geometry preferences
  - Design for selectivity between metals

**Parameter Selection**
  - Use strict constraints for well-defined systems
  - Use flexible constraints for screening
  - Adjust tolerance based on metal type
  - Consider experimental data for validation

**Result Validation**
  - Check coordination distances against literature
  - Verify geometric parameters
  - Compare with crystal structures when available
  - Validate binding affinities experimentally

Troubleshooting
---------------

**Common Issues:**

1. **No coordination detected:**
   - Check metal detection distance
   - Verify ligand has coordinating atoms
   - Reduce distance/angle tolerances

2. **Too many constraint violations:**
   - Increase geometric tolerances
   - Use flexible constraint presets
   - Check metal center geometry assignment

3. **Poor pose quality:**
   - Increase sampling exhaustiveness
   - Use coordination-focused sampling
   - Check ligand conformer generation

4. **Unrealistic binding affinities:**
   - Verify metal charges and oxidation states
   - Check scoring function parameters
   - Validate against experimental data

**Performance Optimization:**

- Use distance-only constraints for screening
- Reduce exhaustiveness for initial filtering
- Focus sampling only when needed
- Use appropriate tolerance values

Advanced Features
-----------------

Custom Metal Centers
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from pandadock.docking.metal_coordination import MetalCenter, MetalType, CoordinationGeometry
   
   # Define custom metal center
   custom_metal = MetalCenter(
       metal_type=MetalType.ZN,
       coordinates=np.array([25.0, 30.0, 15.0]),
       coordination_number=4,
       geometry=CoordinationGeometry.TETRAHEDRAL,
       oxidation_state=2,
       charge=2.0
   )
   
   # Use in docking
   engine.metal_centers = [custom_metal]

Constraint Optimization
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from pandadock.utils.metal_constraints import apply_metal_constraints_to_pose
   
   # Optimize pose to satisfy constraints
   optimized_coords, results = apply_metal_constraints_to_pose(
       pose.coordinates,
       ligand_atom_types,
       metal_centers,
       constraint_preset="strict"
   )
   
   print(f"Constraint satisfaction: {results['constraint_satisfaction']:.3f}")
   print(f"Violations: {results['total_violations']}")

Integration with Standard Docking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Combine with standard PandaDock features
   from pandadock import PandaDock
   
   # Use MetalDockingEngine as engine
   docker = PandaDock(
       engine='metal',
       config=config,
       metal_config=metal_config
   )
   
   # Standard PandaDock interface
   results = docker.dock(
       receptor='metalloprotein.pdb',
       ligand='chelator.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )

This comprehensive metal docking system in PandaDock enables accurate modeling of metalloprotein-ligand interactions, supporting drug discovery efforts targeting this important class of proteins.