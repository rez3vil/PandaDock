Input/Output Module
===================

The I/O module handles reading and writing of molecular structure files and data formats.

Ligand Preparer
---------------

Ligand preparation and file format conversion.

.. automodule:: io.ligand_preparer
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: io.ligand_preparer.LigandPreparer
   :members:
   :special-members: __init__
   :show-inheritance:

Examples
--------

Basic Ligand Preparation
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from io.ligand_preparer import LigandPreparer
   
   # Initialize preparer
   preparer = LigandPreparer()
   
   # Prepare ligand from SDF file
   ligand_data = preparer.prepare_ligand('ligand.sdf')
   
   print(f"Ligand name: {ligand_data['name']}")
   print(f"Number of atoms: {ligand_data['num_atoms']}")
   print(f"Number of bonds: {len(ligand_data['bonds'])}")
   print(f"Molecular weight: {ligand_data['molecular_weight']:.2f} Da")

SMILES to 3D Conversion
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Prepare ligand from SMILES string
   smiles = "CCO"  # Ethanol
   ligand_data = preparer.prepare_from_smiles(smiles)
   
   print(f"Generated {ligand_data['num_atoms']} atoms from SMILES")
   print(f"Coordinates shape: {ligand_data['coordinates'].shape}")
   print(f"Atom types: {ligand_data['atom_types']}")

File Format Support
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Supported input formats
   formats = ['.sdf', '.mol2', '.pdb', '.smi', '.smiles']
   
   for fmt in formats:
       if preparer.is_format_supported(fmt):
           print(f"Format {fmt} is supported")
   
   # Prepare from different formats
   sdf_data = preparer.prepare_ligand('compound.sdf')
   mol2_data = preparer.prepare_ligand('compound.mol2')
   pdb_data = preparer.prepare_ligand('compound.pdb')

Molecular Properties
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Access calculated molecular properties
   properties = ligand_data['properties']
   
   print(f"Molecular weight: {properties['molecular_weight']:.2f} Da")
   print(f"LogP: {properties['logp']:.2f}")
   print(f"H-bond donors: {properties['hbd']}")
   print(f"H-bond acceptors: {properties['hba']}")
   print(f"Rotatable bonds: {properties['num_rotatable_bonds']}")
   print(f"TPSA: {properties['tpsa']:.1f} Ų")
   print(f"Lipinski violations: {properties['lipinski_violations']}")

Structure Validation
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Validate prepared ligand structure
   is_valid = preparer.validate_ligand(ligand_data)
   
   if is_valid:
       print("Ligand structure is valid")
   else:
       print("Ligand structure validation failed")
   
   # Check specific validation criteria
   if ligand_data['molecular_weight'] > 500:
       print("Warning: High molecular weight (>500 Da)")
   
   if ligand_data['num_atoms'] > 100:
       print("Warning: Large number of atoms (>100)")

Save Prepared Ligand
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Save in different formats
   preparer.save_ligand(ligand_data, 'output.sdf', format='sdf')
   preparer.save_ligand(ligand_data, 'output.mol2', format='mol2')
   preparer.save_ligand(ligand_data, 'output.pdb', format='pdb')
   
   print("Ligand saved in multiple formats")

Advanced Features
^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Conformer generation configuration
   preparer = LigandPreparer()
   preparer.max_conformers = 1000
   preparer.energy_threshold = 10.0  # kcal/mol
   preparer.rmsd_threshold = 0.5     # Angstroms
   
   # Generate multiple conformers
   conformers = preparer.generate_conformers(smiles, num_conformers=10)
   
   print(f"Generated {len(conformers)} conformers")
   for i, conf in enumerate(conformers):
       print(f"Conformer {i+1}: Energy = {conf['energy']:.2f} kcal/mol")

Bond Inference
^^^^^^^^^^^^^^

.. code-block:: python

   # Infer bonds from coordinates (for PDB files)
   coordinates = ligand_data['coordinates']
   atom_types = ligand_data['atom_types']
   
   inferred_bonds = preparer.infer_bonds_from_coordinates(coordinates, atom_types)
   
   print(f"Inferred {len(inferred_bonds)} bonds")
   for bond in inferred_bonds[:5]:  # Show first 5 bonds
       atom1, atom2, bond_type = bond
       print(f"Bond: {atom1}-{atom2} ({bond_type})")

Partial Charges
^^^^^^^^^^^^^^^

.. code-block:: python

   # Calculate partial charges
   charges = preparer.calculate_partial_charges(smiles)
   
   print("Partial charges:")
   for i, (atom, charge) in enumerate(zip(ligand_data['atom_types'], charges)):
       print(f"  {i+1:2d} {atom:2s}: {charge:6.3f}")
   
   total_charge = sum(charges)
   print(f"Total charge: {total_charge:.3f}")

Rotatable Bonds
^^^^^^^^^^^^^^^

.. code-block:: python

   # Identify rotatable bonds
   bonds = ligand_data['bonds']
   rotatable_bonds = preparer.identify_rotatable_bonds(bonds)
   
   print(f"Found {len(rotatable_bonds)} rotatable bonds:")
   for bond_idx in rotatable_bonds:
       bond = bonds[bond_idx]
       atom1, atom2, bond_type = bond
       print(f"  Bond {bond_idx}: {atom1}-{atom2} ({bond_type})")

Geometric Constraints
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Apply geometric constraints to coordinates
   raw_coordinates = ligand_data['coordinates']
   constrained_coords = preparer.apply_geometric_constraints(raw_coordinates)
   
   # Check bond lengths
   from utils.math_utils import distance_matrix
   
   distances = distance_matrix(constrained_coords, constrained_coords)
   
   for bond in ligand_data['bonds']:
       atom1, atom2, _ = bond
       bond_length = distances[atom1, atom2]
       print(f"Bond {atom1}-{atom2}: {bond_length:.3f} Å")

Batch Processing
^^^^^^^^^^^^^^^^

.. code-block:: python

   import glob
   
   # Process multiple ligand files
   ligand_files = glob.glob('ligands/*.sdf')
   
   processed_ligands = []
   for file_path in ligand_files:
       try:
           ligand_data = preparer.prepare_ligand(file_path)
           processed_ligands.append(ligand_data)
           print(f"Processed: {file_path}")
       except Exception as e:
           print(f"Failed to process {file_path}: {e}")
   
   print(f"Successfully processed {len(processed_ligands)} ligands")

Custom Configuration
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Custom preparer configuration
   config = {
       'max_conformers': 500,
       'energy_threshold': 15.0,
       'rmsd_threshold': 0.3,
       'force_field': 'mmff94',
       'optimize_geometry': True,
       'add_hydrogens': True,
       'assign_stereochemistry': True
   }
   
   preparer = LigandPreparer(config=config)

Error Handling
^^^^^^^^^^^^^^

.. code-block:: python

   from io.ligand_preparer import LigandPreparer
   
   preparer = LigandPreparer()
   
   try:
       ligand_data = preparer.prepare_ligand('nonexistent.sdf')
   except FileNotFoundError:
       print("File not found")
   except ValueError as e:
       print(f"Invalid file format: {e}")
   except Exception as e:
       print(f"Preparation failed: {e}")

File Format Details
-------------------

**SDF Format Support:**

.. code-block:: python

   # Read SDF with multiple structures
   sdf_data = preparer.read_sdf_file('multi_structure.sdf')
   
   # Access SDF properties
   print(f"SDF properties: {sdf_data.get('properties', {})}")

**MOL2 Format Support:**

.. code-block:: python

   # Read MOL2 with atom types and charges
   mol2_data = preparer.read_mol2_file('ligand.mol2')
   
   # MOL2 includes detailed atom typing
   print(f"MOL2 atom types: {mol2_data['atom_types']}")

**PDB Format Support:**

.. code-block:: python

   # Read PDB (typically for small molecules/ligands)
   pdb_data = preparer.read_pdb_file('ligand.pdb')
   
   # PDB may require bond inference
   if not pdb_data['bonds']:
       print("No bonds found, inferring from coordinates...")
       pdb_data['bonds'] = preparer.infer_bonds_from_coordinates(
           pdb_data['coordinates'], 
           pdb_data['atom_types']
       )

Output Format Customization
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # Customize SDF output
   preparer.write_sdf(ligand_data, 'custom.sdf')
   
   # Customize MOL2 output with specific atom types
   preparer.write_mol2(ligand_data, 'custom.mol2')
   
   # Customize PDB output
   preparer.write_pdb(ligand_data, 'custom.pdb')

Performance Optimization
------------------------

.. code-block:: python

   # For large-scale processing
   preparer = LigandPreparer()
   
   # Disable expensive operations for screening
   preparer.optimize_geometry = False
   preparer.generate_conformers = False
   
   # Use simplified charge calculation
   preparer.charge_method = 'gasteiger'  # Faster than quantum methods
   
   # Process in parallel
   from multiprocessing import Pool
   
   def process_ligand(file_path):
       return preparer.prepare_ligand(file_path)
   
   with Pool(processes=8) as pool:
       results = pool.map(process_ligand, ligand_files)

Best Practices
--------------

1. **Always validate input files:**
   
   .. code-block:: python
   
      if not os.path.exists(ligand_file):
          raise FileNotFoundError(f"Ligand file not found: {ligand_file}")

2. **Handle different input formats appropriately:**
   
   .. code-block:: python
   
      file_ext = os.path.splitext(ligand_file)[1].lower()
      if file_ext not in preparer.supported_formats:
          raise ValueError(f"Unsupported format: {file_ext}")

3. **Check molecular properties:**
   
   .. code-block:: python
   
      if ligand_data['molecular_weight'] > 1000:
          print("Warning: Very large molecule")

4. **Use appropriate methods for your needs:**
   - Use SMILES input for virtual compound libraries
   - Use SDF for experimental compounds with known 3D structures
   - Use MOL2 for detailed atom typing information