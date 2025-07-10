Utilities Module
================

The utils module provides utility functions and helper classes for various computational tasks.

IC50 Calculator
---------------

Calculate IC50 values and binding kinetics from binding affinities.

.. automodule:: utils.ic50_calculator
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: utils.ic50_calculator.IC50Calculator
   :members:
   :special-members: __init__
   :show-inheritance:

Mathematical Utilities
----------------------

Mathematical functions and computational tools.

.. automodule:: utils.math_utils
   :members:
   :undoc-members:
   :show-inheritance:

Rotamer Library
---------------

Rotamer library management and sampling.

.. automodule:: utils.rotamer_lib
   :members:
   :undoc-members:
   :show-inheritance:

Examples
--------

IC50 Calculations
^^^^^^^^^^^^^^^^^

.. code-block:: python

   from utils.ic50_calculator import IC50Calculator
   
   # Initialize calculator
   calc = IC50Calculator(temperature=298.15)  # Room temperature
   
   # Convert binding free energy to IC50
   delta_g = -8.5  # kcal/mol (favorable binding)
   ic50 = calc.delta_g_to_ic50(delta_g)
   print(f"IC50: {ic50:.2e} M ({ic50*1e9:.1f} nM)")
   
   # Convert IC50 to binding free energy
   ic50_measured = 1e-8  # M (10 nM)
   delta_g_calc = calc.ic50_to_delta_g(ic50_measured)
   print(f"ΔG: {delta_g_calc:.2f} kcal/mol")

Ligand Efficiency
^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Calculate ligand efficiency
   delta_g = -8.5  # kcal/mol
   num_heavy_atoms = 25
   
   le = calc.calculate_ligand_efficiency(delta_g, num_heavy_atoms)
   print(f"Ligand Efficiency: {le:.3f} kcal/mol per heavy atom")
   
   # Calculate lipophilic efficiency
   logp = 3.2
   lipe = calc.calculate_lipophilic_efficiency(delta_g, logp)
   print(f"Lipophilic Efficiency: {lipe:.2f}")

Drug-likeness Assessment
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Comprehensive drug-likeness assessment
   assessment = calc.assess_drug_likeness(
       delta_g=-8.5,
       molecular_weight=450.0,
       logp=3.2,
       hbd=2,  # Hydrogen bond donors
       hba=5,  # Hydrogen bond acceptors
       num_heavy_atoms=25
   )
   
   print(f"IC50: {assessment['ic50']:.2e} M")
   print(f"Ligand Efficiency: {assessment['ligand_efficiency']:.3f}")
   print(f"Lipinski compliant: {assessment['lipinski_compliant']}")
   print(f"Veber compliant: {assessment['veber_compliant']}")
   print(f"Drug-likeness category: {assessment['drug_likeness_category']}")

Binding Kinetics
^^^^^^^^^^^^^^^^

.. code-block:: python

   # Calculate binding kinetics
   kinetics = calc.calculate_binding_kinetics(
       delta_g=-8.5,
       delta_h=-12.0,  # Enthalpy change (kcal/mol)
       delta_s=-0.012  # Entropy change (kcal/mol/K)
   )
   
   print(f"Kd: {kinetics['kd']:.2e} M")
   print(f"Ka: {kinetics['ka']:.2e} M⁻¹")
   print(f"k_on: {kinetics['k_on']:.2e} M⁻¹s⁻¹")
   print(f"k_off: {kinetics['k_off']:.2e} s⁻¹")
   print(f"Residence time: {kinetics['residence_time']:.2f} s")

Mathematical Functions
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from utils.math_utils import distance_matrix, angle_between_vectors
   import numpy as np
   
   # Calculate distance matrix
   coords1 = np.random.rand(10, 3)
   coords2 = np.random.rand(15, 3)
   
   distances = distance_matrix(coords1, coords2)
   print(f"Distance matrix shape: {distances.shape}")
   
   # Calculate angles between vectors
   vec1 = np.array([1, 0, 0])
   vec2 = np.array([0, 1, 0])
   angle = angle_between_vectors(vec1, vec2)
   print(f"Angle: {np.degrees(angle):.1f} degrees")

Rotamer Library Usage
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from utils.rotamer_lib import RotamerLibrary
   
   # Load rotamer library
   rotlib = RotamerLibrary()
   
   # Get rotamers for a residue
   rotamers = rotlib.get_rotamers('ARG')
   print(f"Found {len(rotamers)} rotamers for ARG")
   
   for i, rotamer in enumerate(rotamers[:3]):
       print(f"Rotamer {i+1}: chi angles = {rotamer['chi_angles']}")
   
   # Sample random rotamer
   random_rotamer = rotlib.sample_rotamer('TYR')
   print(f"Random TYR rotamer: {random_rotamer}")

Advanced IC50 Calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Dose-response curve calculation
   ic50 = 1e-8  # M
   curve_params = calc.calculate_dose_response_curve(
       ic50=ic50,
       hill_coefficient=1.0,
       max_response=100.0,
       min_response=0.0
   )
   
   # Plot dose-response curve
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(8, 6))
   plt.semilogx(curve_params['concentrations'] * 1e9, curve_params['responses'])
   plt.xlabel('Concentration (nM)')
   plt.ylabel('Response (%)')
   plt.title('Dose-Response Curve')
   plt.grid(True)
   plt.show()

Cooperativity Analysis
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Analyze cooperativity for multiple binding sites
   ic50_values = [1e-8, 5e-9, 2e-8]  # IC50 values for different sites
   hill_coefficient = 1.2
   
   cooperativity = calc.calculate_cooperativity(ic50_values, hill_coefficient)
   
   print(f"Geometric mean Kd: {cooperativity['geometric_mean_kd']:.2e} M")
   print(f"Hill coefficient: {cooperativity['hill_coefficient']}")
   print(f"Cooperativity type: {cooperativity['cooperativity_type']}")

Selectivity Analysis
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Calculate selectivity index
   ic50_target = 1e-8     # Target protein IC50
   ic50_offtarget = 1e-6  # Off-target protein IC50
   
   selectivity = calc.calculate_selectivity_index(ic50_target, ic50_offtarget)
   print(f"Selectivity index: {selectivity:.1f}-fold")

Unit Conversions
^^^^^^^^^^^^^^^^

.. code-block:: python

   # Convert between concentration units
   ic50_m = 1e-8  # M
   ic50_nm = calc.convert_units(ic50_m, 'M', 'nM')
   print(f"IC50: {ic50_m:.2e} M = {ic50_nm:.1f} nM")
   
   # Convert between energy units
   delta_g_kcal = -8.5  # kcal/mol
   delta_g_kj = calc.convert_units(delta_g_kcal, 'kcal/mol', 'kJ/mol')
   print(f"ΔG: {delta_g_kcal:.1f} kcal/mol = {delta_g_kj:.1f} kJ/mol")

Thermodynamic Summary
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Get comprehensive thermodynamic summary
   summary = calc.get_thermodynamic_summary(
       delta_g=-8.5,
       num_heavy_atoms=25,
       molecular_weight=450.0,
       logp=3.2
   )
   
   print("Thermodynamic Summary:")
   print("-" * 40)
   print(f"ΔG: {summary['thermodynamics']['delta_g']:.2f} kcal/mol")
   print(f"IC50: {summary['thermodynamics']['ic50']:.2e} M")
   print(f"Kd: {summary['thermodynamics']['kd']:.2e} M")
   print(f"LE: {summary['efficiency_metrics']['ligand_efficiency']:.3f}")
   print(f"LipE: {summary['efficiency_metrics']['lipophilic_efficiency']:.2f}")

Mathematical Utilities Advanced
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from utils.math_utils import (
       rmsd, center_of_mass, radius_of_gyration,
       principal_components, align_structures
   )
   
   # Calculate RMSD between two structures
   coords1 = np.random.rand(20, 3)
   coords2 = coords1 + np.random.normal(0, 0.1, coords1.shape)
   
   rmsd_value = rmsd(coords1, coords2)
   print(f"RMSD: {rmsd_value:.3f} Å")
   
   # Calculate center of mass
   masses = np.ones(20)  # Assume unit masses
   com = center_of_mass(coords1, masses)
   print(f"Center of mass: {com}")
   
   # Calculate radius of gyration
   rog = radius_of_gyration(coords1, masses)
   print(f"Radius of gyration: {rog:.3f} Å")
   
   # Principal component analysis
   eigenvalues, eigenvectors = principal_components(coords1)
   print(f"Principal components: {eigenvalues}")

Error Handling
--------------

.. code-block:: python

   from utils.ic50_calculator import IC50Calculator
   
   calc = IC50Calculator()
   
   try:
       # Invalid input handling
       ic50 = calc.delta_g_to_ic50(5.0)  # Positive ΔG (unfavorable)
       print(f"IC50: {ic50}")  # Will be inf
   except ValueError as e:
       print(f"Error: {e}")
   
   try:
       # Unit conversion error
       value = calc.convert_units(1.0, 'invalid_unit', 'M')
   except ValueError as e:
       print(f"Unit conversion error: {e}")

Performance Tips
----------------

1. **Reuse calculator instances:**
   
   .. code-block:: python
   
      # Good: Reuse calculator
      calc = IC50Calculator()
      for delta_g in binding_energies:
          ic50 = calc.delta_g_to_ic50(delta_g)
   
2. **Batch calculations:**
   
   .. code-block:: python
   
      # Process multiple values efficiently
      delta_g_values = [-8.5, -7.2, -9.1, -6.8]
      ic50_values = [calc.delta_g_to_ic50(dg) for dg in delta_g_values]

3. **Avoid repeated expensive operations:**
   
   .. code-block:: python
   
      # Pre-calculate distance matrices for multiple operations
      distances = distance_matrix(coords1, coords2)
      # Use distances for multiple calculations

Constants and References
------------------------

The module provides standard thermodynamic constants:

.. code-block:: python

   calc = IC50Calculator(temperature=298.15)  # Room temperature
   print(f"Gas constant: {calc.R} J/(mol·K)")
   print(f"kT: {calc.kT_kcal:.4f} kcal/mol")
   print(f"Standard concentration: {calc.standard_conc} M")