Flexible Docking Examples
=========================

This section provides comprehensive examples of flexible docking using PandaDock, covering receptor flexibility, ligand conformational sampling, and induced fit docking.

Introduction
------------

Flexible docking accounts for the dynamic nature of protein-ligand interactions by allowing movement of receptor residues and exploring multiple ligand conformations. This is crucial for:

- Induced fit binding
- Allosteric site binding
- Novel binding mode discovery
- Accurate binding affinity prediction

Examples Overview
-----------------

1. **Basic Flexible Docking**: Simple receptor flexibility
2. **Induced Fit Docking**: Full induced fit protocol
3. **Ensemble Docking**: Multiple receptor conformations
4. **Flexible Loop Docking**: Specific loop flexibility
5. **Allosteric Site Docking**: Flexible allosteric binding

Example 1: Basic Flexible Docking
----------------------------------

Start with basic sidechain flexibility:

.. code-block:: python

   from pandadock import PandaDock
   
   # Configure basic flexible docking
   docker = PandaDock(
       engine='pandacore',  # PandaCore algorithm handles flexibility well
       flexible_residues=['ARG123', 'TYR456', 'ASP789'],
       flexibility_config={
           'backbone_flexibility': False,    # Keep backbone rigid
           'sidechain_flexibility': True,    # Allow sidechain movement
           'max_torsion_change': 30.0,      # Max angle change (degrees)
           'clash_tolerance': 0.8,          # Flexibility clash tolerance
           'minimization_steps': 1000       # Post-docking minimization
       }
   )
   
   # Run flexible docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[22.0, 22.0, 22.0]  # Slightly larger box for flexibility
   )
   
   # Analyze flexibility
   print(f"Flexible docking completed: {len(results.poses)} poses")
   
   for i, pose in enumerate(results.poses[:3]):
       print(f"\nPose {i+1}:")
       print(f"  Score: {pose.score:.3f}")
       print(f"  Receptor RMSD: {pose.receptor_rmsd:.2f} Å")
       print(f"  Ligand RMSD: {pose.ligand_rmsd:.2f} Å")
       
       # Check which residues moved
       if hasattr(pose, 'flexible_residues'):
           print(f"  Residue movements:")
           for res, rmsd in pose.flexible_residues.items():
               print(f"    {res}: {rmsd:.2f} Å")

Example 2: Auto-Detection of Flexible Residues
-----------------------------------------------

Let PandaDock automatically identify flexible residues:

.. code-block:: python

   # Auto-detect flexible residues
   docker = PandaDock(
       engine='pandacore',
       auto_flexible=True,
       flexibility_config={
           'detection_method': 'b_factor',     # Use B-factors
           'b_factor_threshold': 30.0,         # B-factor cutoff
           'max_flexible_residues': 5,         # Limit number
           'binding_site_radius': 8.0,         # Only near binding site
           'exclude_backbone': True,           # Exclude backbone atoms
           'include_waters': False,            # Exclude water molecules
           'clash_detection': True,            # Detect and resolve clashes
           'energy_minimization': True         # Minimize after movement
       }
   )
   
   # Run with auto-detection
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Check which residues were identified as flexible
   print("Auto-detected flexible residues:")
   for res_id, b_factor in results.flexible_residues_detected.items():
       print(f"  {res_id}: B-factor = {b_factor:.1f}")

Example 3: Induced Fit Docking Protocol
----------------------------------------

Implement full induced fit docking:

.. code-block:: python

   from pandadock import InducedFitDocking
   
   # Configure induced fit protocol
   ifd = InducedFitDocking(
       protocol='glide_like',  # Glide-like protocol
       stages={
           'initial_docking': {
               'engine': 'physics',
               'flexibility': 'none',
               'exhaustiveness': 8,
               'num_poses': 20
           },
           'receptor_refinement': {
               'method': 'prime_like',
               'flexible_residues': 'auto',
               'refinement_radius': 5.0,
               'minimization_steps': 5000,
               'constraint_weight': 10.0
           },
           'redocking': {
               'engine': 'pandaml',
               'flexibility': 'refined',
               'exhaustiveness': 16,
               'num_poses': 10
           }
       }
   )
   
   # Run induced fit docking
   results = ifd.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze induced fit results
   print("Induced Fit Docking Results:")
   print("=" * 35)
   
   for stage, stage_results in results.stages.items():
       print(f"\n{stage.title()}:")
       print(f"  Runtime: {stage_results.runtime:.2f} seconds")
       print(f"  Poses: {len(stage_results.poses)}")
       if hasattr(stage_results, 'receptor_rmsd'):
           print(f"  Receptor RMSD: {stage_results.receptor_rmsd:.2f} Å")
   
   # Best pose from final stage
   best_pose = results.best_pose
   print(f"\nBest Final Pose:")
   print(f"  Score: {best_pose.score:.3f}")
   print(f"  Induced fit score: {best_pose.if_score:.3f}")
   print(f"  Receptor conformational change: {best_pose.receptor_rmsd:.2f} Å")

Example 4: Ensemble Docking
----------------------------

Use multiple receptor conformations:

.. code-block:: python

   from pandadock import EnsembleDocking
   
   # Prepare ensemble of receptor conformations
   ensemble = EnsembleDocking(
       conformations=[
           'protein_conf1.pdb',
           'protein_conf2.pdb',  
           'protein_conf3.pdb'
       ],
       ensemble_config={
           'weighting': 'boltzmann',        # Boltzmann weighting
           'temperature': 298.15,           # Temperature for weighting
           'clustering': True,              # Cluster similar conformations
           'cluster_threshold': 1.0,        # RMSD threshold for clustering
           'max_conformations': 5,          # Max conformations to use
           'consensus_scoring': True         # Consensus across ensemble
       }
   )
   
   # Run ensemble docking
   results = ensemble.dock(
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze ensemble results
   print("Ensemble Docking Results:")
   print("=" * 30)
   
   for conf_id, conf_results in results.conformations.items():
       print(f"\nConformation {conf_id}:")
       print(f"  Weight: {conf_results.weight:.3f}")
       print(f"  Best score: {conf_results.best_score:.3f}")
       print(f"  Contribution: {conf_results.contribution:.3f}")
   
   # Consensus pose
   consensus_pose = results.consensus_pose
   print(f"\nConsensus Pose:")
   print(f"  Score: {consensus_pose.score:.3f}")
   print(f"  Ensemble agreement: {consensus_pose.agreement:.3f}")
   print(f"  Conformational diversity: {consensus_pose.diversity:.3f}")

Example 5: Flexible Loop Docking
---------------------------------

Target specific flexible loops:

.. code-block:: python

   # Configure flexible loop docking
   docker = PandaDock(
       engine='pandacore',
       flexible_loops=[
           {
               'residues': ['GLY45', 'ALA46', 'SER47', 'GLY48'],
               'flexibility': 'backbone',
               'max_deviation': 2.0,  # Max Å deviation from starting structure
               'anchor_residues': ['VAL44', 'PRO49']  # Fixed anchor points
           },
           {
               'residues': ['LYS156', 'GLU157', 'ASP158'],
               'flexibility': 'sidechain',
               'max_torsion_change': 45.0
           }
       ],
       loop_config={
           'loop_modeling': 'ab_initio',     # Ab initio loop modeling
           'num_loop_models': 10,           # Number of loop conformations
           'loop_refinement': True,         # Refine loop structures
           'clash_resolution': True,        # Resolve clashes
           'energy_minimization': True      # Minimize loop energy
       }
   )
   
   # Run flexible loop docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze loop flexibility
   print("Flexible Loop Docking Results:")
   print("=" * 35)
   
   for pose in results.poses[:3]:
       print(f"\nPose (Score: {pose.score:.3f}):")
       for loop_id, loop_data in pose.flexible_loops.items():
           print(f"  Loop {loop_id}:")
           print(f"    RMSD: {loop_data['rmsd']:.2f} Å")
           print(f"    Energy: {loop_data['energy']:.2f} kcal/mol")
           print(f"    Contacts: {loop_data['ligand_contacts']}")

Example 6: Allosteric Site Docking
-----------------------------------

Dock to allosteric sites with flexibility:

.. code-block:: python

   from pandadock.allosteric import AllostericDocking
   
   # Configure allosteric docking
   allo_docker = AllostericDocking(
       allosteric_sites=[
           {
               'center': [45.0, 30.0, 25.0],
               'size': [18.0, 18.0, 18.0],
               'type': 'allosteric',
               'communication_pathway': True  # Model pathway to active site
           }
       ],
       flexibility_config={
           'communication_residues': 'auto',  # Auto-detect pathway residues
           'pathway_flexibility': True,       # Allow pathway flexibility
           'cooperative_binding': True,       # Model cooperative effects
           'allosteric_networks': True        # Map allosteric networks
       }
   )
   
   # Run allosteric docking
   results = allo_docker.dock(
       receptor='protein.pdb',
       ligand='allosteric_ligand.sdf',
       orthosteric_ligand='active_site_ligand.sdf',  # If occupied
       cooperative_effects=True
   )
   
   # Analyze allosteric effects
   print("Allosteric Docking Results:")
   print("=" * 30)
   
   for pose in results.poses:
       print(f"\nPose (Score: {pose.score:.3f}):")
       print(f"  Allosteric effect: {pose.allosteric_effect:.3f}")
       print(f"  Pathway perturbation: {pose.pathway_perturbation:.3f}")
       print(f"  Cooperative binding: {pose.cooperative_score:.3f}")
       
       if hasattr(pose, 'communication_pathway'):
           print(f"  Communication pathway:")
           for residue in pose.communication_pathway:
               print(f"    {residue}")

Example 7: Multi-Domain Flexibility
------------------------------------

Handle multi-domain proteins with interdomain flexibility:

.. code-block:: python

   # Configure multi-domain flexibility
   docker = PandaDock(
       engine='pandacore',
       multi_domain=True,
       domain_config={
           'domains': [
               {
                   'name': 'domain1',
                   'residues': range(1, 150),
                   'flexibility': 'rigid_body',
                   'anchor_residues': [75, 76, 77]
               },
               {
                   'name': 'domain2', 
                   'residues': range(200, 350),
                   'flexibility': 'rigid_body',
                   'anchor_residues': [275, 276, 277]
               }
           ],
           'interdomain_flexibility': True,
           'hinge_regions': [
               {'residues': range(150, 200), 'flexibility': 'backbone'}
           ],
           'max_domain_movement': 5.0,  # Max Å movement between domains
           'interdomain_contacts': True  # Maintain key contacts
       }
   )
   
   # Run multi-domain docking
   results = docker.dock(
       receptor='multidomain_protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[25.0, 25.0, 25.0]
   )
   
   # Analyze domain movements
   print("Multi-Domain Docking Results:")
   print("=" * 35)
   
   for pose in results.poses:
       print(f"\nPose (Score: {pose.score:.3f}):")
       for domain_name, domain_data in pose.domain_movements.items():
           print(f"  {domain_name}:")
           print(f"    Translation: {domain_data['translation']:.2f} Å")
           print(f"    Rotation: {domain_data['rotation']:.2f}°")
           print(f"    Interdomain contacts: {domain_data['contacts']}")

Example 8: Membrane Protein Flexibility
----------------------------------------

Handle membrane protein flexibility:

.. code-block:: python

   from pandadock.membrane import MembraneFlexibleDocking
   
   # Configure membrane protein docking
   membrane_docker = MembraneFlexibleDocking(
       membrane_config={
           'membrane_type': 'lipid_bilayer',
           'lipid_composition': 'POPC',
           'membrane_thickness': 30.0,  # Å
           'membrane_center': [0, 0, 0],
           'membrane_normal': [0, 0, 1]
       },
       flexibility_config={
           'transmembrane_regions': 'auto',  # Auto-detect TM regions
           'tm_flexibility': 'limited',      # Limited TM flexibility
           'extracellular_loops': 'flexible', # Flexible EC loops
           'intracellular_loops': 'flexible', # Flexible IC loops
           'lipid_interactions': True,       # Model lipid interactions
           'membrane_insertion': True        # Allow membrane insertion
       }
   )
   
   # Run membrane protein docking
   results = membrane_docker.dock(
       receptor='membrane_protein.pdb',
       ligand='ligand.sdf',
       binding_site='extracellular',  # or 'intracellular', 'transmembrane'
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze membrane effects
   print("Membrane Protein Docking Results:")
   print("=" * 40)
   
   for pose in results.poses:
       print(f"\nPose (Score: {pose.score:.3f}):")
       print(f"  Membrane insertion depth: {pose.membrane_depth:.2f} Å")
       print(f"  Lipid interactions: {pose.lipid_interactions}")
       print(f"  Membrane orientation: {pose.membrane_orientation:.2f}°")

Example 9: Covalent Flexible Docking
-------------------------------------

Dock covalent inhibitors with flexibility:

.. code-block:: python

   from pandadock.covalent import CovalentFlexibleDocking
   
   # Configure covalent docking
   covalent_docker = CovalentFlexibleDocking(
       covalent_config={
           'target_residue': 'CYS145',
           'reaction_type': 'michael_addition',
           'warhead': 'acrylamide',
           'bond_formation': True,
           'covalent_geometry': 'tetrahedral'
       },
       flexibility_config={
           'target_residue_flexibility': True,
           'neighboring_residues': True,
           'flexible_radius': 5.0,  # Å around covalent bond
           'bond_constraint': True,  # Constrain covalent bond
           'reaction_coordinate': True  # Model reaction coordinate
       }
   )
   
   # Run covalent docking
   results = covalent_docker.dock(
       receptor='protein.pdb',
       ligand='covalent_ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze covalent binding
   print("Covalent Docking Results:")
   print("=" * 30)
   
   for pose in results.poses:
       print(f"\nPose (Score: {pose.score:.3f}):")
       print(f"  Covalent bond length: {pose.covalent_bond_length:.2f} Å")
       print(f"  Bond angle: {pose.covalent_bond_angle:.2f}°")
       print(f"  Reaction barrier: {pose.reaction_barrier:.2f} kcal/mol")
       print(f"  Covalent binding energy: {pose.covalent_energy:.2f} kcal/mol")

Example 10: Complete Flexible Docking Workflow
-----------------------------------------------

Comprehensive workflow combining multiple flexibility types:

.. code-block:: python

   #!/usr/bin/env python3
   """
   Complete Flexible Docking Workflow
   
   This script demonstrates a comprehensive flexible docking workflow
   that combines multiple types of flexibility.
   """
   
   from pandadock import PandaDock
   from pandadock.analysis import FlexibilityAnalysis
   from pandadock.visualization import FlexibilityVisualizer
   import numpy as np
   
   class FlexibleDockingWorkflow:
       def __init__(self, config):
           self.config = config
           self.results = None
           
       def analyze_receptor_flexibility(self, receptor_file):
           """Analyze receptor flexibility patterns"""
           print("Analyzing receptor flexibility...")
           
           analyzer = FlexibilityAnalysis()
           
           # B-factor analysis
           b_factors = analyzer.analyze_b_factors(receptor_file)
           
           # Cavity analysis
           cavities = analyzer.analyze_cavities(receptor_file)
           
           # Sequence analysis
           sequence_flexibility = analyzer.analyze_sequence_flexibility(receptor_file)
           
           # Combine analyses
           flexibility_profile = analyzer.combine_analyses(
               b_factors, cavities, sequence_flexibility
           )
           
           return flexibility_profile
           
       def configure_flexibility(self, flexibility_profile):
           """Configure flexibility based on analysis"""
           print("Configuring flexibility parameters...")
           
           # Select flexible residues
           flexible_residues = []
           for residue, score in flexibility_profile.items():
               if score > 0.7:  # High flexibility score
                   flexible_residues.append(residue)
           
           # Configure docker
           docker = PandaDock(
               engine='pandacore',
               flexible_residues=flexible_residues,
               flexibility_config={
                   'backbone_flexibility': True,
                   'sidechain_flexibility': True,
                   'max_torsion_change': 30.0,
                   'adaptive_flexibility': True,
                   'energy_minimization': True,
                   'clash_resolution': True
               },
               ga_config={
                   'population_size': 200,
                   'generations': 300,
                   'mutation_rate': 0.02,
                   'crossover_rate': 0.8,
                   'flexibility_mutation': True
               }
           )
           
           return docker
           
       def run_flexible_docking(self, docker):
           """Run flexible docking"""
           print("Running flexible docking...")
           
           results = docker.dock(
               receptor=self.config['receptor'],
               ligand=self.config['ligand'],
               center=self.config['center'],
               size=self.config['size']
           )
           
           return results
           
       def analyze_results(self, results):
           """Analyze flexibility results"""
           print("Analyzing flexibility results...")
           
           analysis = {
               'pose_diversity': [],
               'flexibility_utilization': [],
               'binding_mode_analysis': [],
               'energetic_analysis': []
           }
           
           for pose in results.poses:
               # Pose diversity
               diversity = self.calculate_pose_diversity(pose, results.poses)
               analysis['pose_diversity'].append(diversity)
               
               # Flexibility utilization
               utilization = self.calculate_flexibility_utilization(pose)
               analysis['flexibility_utilization'].append(utilization)
               
               # Binding mode analysis
               binding_mode = self.analyze_binding_mode(pose)
               analysis['binding_mode_analysis'].append(binding_mode)
               
               # Energetic analysis
               energetics = self.analyze_energetics(pose)
               analysis['energetic_analysis'].append(energetics)
           
           return analysis
           
       def calculate_pose_diversity(self, pose, all_poses):
           """Calculate pose diversity"""
           rmsds = []
           for other_pose in all_poses:
               if pose != other_pose:
                   rmsd = self.calculate_rmsd(pose, other_pose)
                   rmsds.append(rmsd)
           return np.mean(rmsds) if rmsds else 0.0
           
       def calculate_flexibility_utilization(self, pose):
           """Calculate how much flexibility was used"""
           if not hasattr(pose, 'flexible_residues'):
               return 0.0
           
           total_movement = 0.0
           count = 0
           
           for residue, rmsd in pose.flexible_residues.items():
               total_movement += rmsd
               count += 1
           
           return total_movement / count if count > 0 else 0.0
           
       def analyze_binding_mode(self, pose):
           """Analyze binding mode"""
           binding_mode = {
               'hbonds': len(pose.interactions.hbonds) if hasattr(pose, 'interactions') else 0,
               'hydrophobic': len(pose.interactions.hydrophobic) if hasattr(pose, 'interactions') else 0,
               'salt_bridges': len(pose.interactions.salt_bridges) if hasattr(pose, 'interactions') else 0,
               'binding_site_occupancy': getattr(pose, 'binding_site_occupancy', 0.0)
           }
           return binding_mode
           
       def analyze_energetics(self, pose):
           """Analyze energetic components"""
           energetics = {
               'total_energy': getattr(pose, 'energy', 0.0),
               'flexibility_penalty': getattr(pose, 'flexibility_penalty', 0.0),
               'binding_energy': getattr(pose, 'binding_energy', 0.0),
               'conformational_energy': getattr(pose, 'conformational_energy', 0.0)
           }
           return energetics
           
       def calculate_rmsd(self, pose1, pose2):
           """Calculate RMSD between poses"""
           # Implementation depends on pose structure
           return 0.0  # Placeholder
           
       def visualize_results(self, results, analysis):
           """Visualize flexible docking results"""
           print("Generating visualizations...")
           
           visualizer = FlexibilityVisualizer()
           
           # Flexibility heatmap
           visualizer.plot_flexibility_heatmap(
               results.poses,
               output_file='flexibility_heatmap.png'
           )
           
           # Binding mode diversity
           visualizer.plot_binding_mode_diversity(
               analysis['binding_mode_analysis'],
               output_file='binding_mode_diversity.png'
           )
           
           # Energy landscape
           visualizer.plot_energy_landscape(
               analysis['energetic_analysis'],
               output_file='energy_landscape.png'
           )
           
           # 3D visualization
           visualizer.create_3d_visualization(
               results.poses[:5],  # Top 5 poses
               output_file='flexible_docking_3d.pml'
           )
           
       def run_workflow(self):
           """Run complete flexible docking workflow"""
           try:
               # Step 1: Analyze receptor flexibility
               flexibility_profile = self.analyze_receptor_flexibility(
                   self.config['receptor']
               )
               
               # Step 2: Configure flexibility
               docker = self.configure_flexibility(flexibility_profile)
               
               # Step 3: Run flexible docking
               results = self.run_flexible_docking(docker)
               
               # Step 4: Analyze results
               analysis = self.analyze_results(results)
               
               # Step 5: Visualize results
               self.visualize_results(results, analysis)
               
               # Step 6: Generate report
               self.generate_report(results, analysis)
               
               print("Flexible docking workflow completed successfully!")
               return results, analysis
               
           except Exception as e:
               print(f"Workflow failed: {e}")
               return None, None
           
       def generate_report(self, results, analysis):
           """Generate comprehensive report"""
           print("Generating report...")
           
           report = f"""
           Flexible Docking Report
           =======================
           
           Configuration:
           - Receptor: {self.config['receptor']}
           - Ligand: {self.config['ligand']}
           - Flexible residues: {len(results.flexible_residues)}
           
           Results Summary:
           - Total poses: {len(results.poses)}
           - Best score: {results.best_pose.score:.3f}
           - Average flexibility utilization: {np.mean(analysis['flexibility_utilization']):.3f}
           - Pose diversity: {np.mean(analysis['pose_diversity']):.3f}
           
           Flexibility Analysis:
           - Backbone flexibility: {results.backbone_flexibility_used}
           - Sidechain flexibility: {results.sidechain_flexibility_used}
           - Induced fit score: {results.induced_fit_score:.3f}
           
           Binding Mode Analysis:
           - Average H-bonds: {np.mean([bm['hbonds'] for bm in analysis['binding_mode_analysis']]):.1f}
           - Average hydrophobic contacts: {np.mean([bm['hydrophobic'] for bm in analysis['binding_mode_analysis']]):.1f}
           - Average salt bridges: {np.mean([bm['salt_bridges'] for bm in analysis['binding_mode_analysis']]):.1f}
           """
           
           with open('flexible_docking_report.txt', 'w') as f:
               f.write(report)
   
   def main():
       # Example configuration
       config = {
           'receptor': 'protein.pdb',
           'ligand': 'ligand.sdf',
           'center': [25.0, 30.0, 15.0],
           'size': [22.0, 22.0, 22.0]
       }
       
       # Run workflow
       workflow = FlexibleDockingWorkflow(config)
       results, analysis = workflow.run_workflow()
       
       if results:
           print("Flexible docking completed successfully!")
           print(f"Best pose score: {results.best_pose.score:.3f}")
           print(f"Report saved to: flexible_docking_report.txt")
       else:
           print("Flexible docking failed!")
   
   if __name__ == "__main__":
       main()

Performance Considerations
--------------------------

**Computational Cost:**
- Flexible docking is 5-50x slower than rigid docking
- Cost scales with number of flexible residues
- Use GA engine for best flexibility handling

**Memory Usage:**
- Flexible docking requires more memory
- Consider memory-efficient modes for large systems
- Use checkpointing for long calculations

**Accuracy vs Speed:**
- More flexibility = higher accuracy but slower
- Start with limited flexibility and expand as needed
- Use ensemble methods for best results

Best Practices
--------------

1. **Start Simple**: Begin with sidechain flexibility only
2. **Gradual Expansion**: Add backbone flexibility carefully
3. **Validate Setup**: Use known flexible binders for validation
4. **Monitor Convergence**: Ensure sampling is adequate
5. **Post-Analysis**: Always analyze flexibility utilization

Troubleshooting
---------------

**Common Issues:**
- Excessive flexibility leads to poor convergence
- Insufficient sampling misses binding modes
- Clash resolution problems
- Memory limitations

**Solutions:**
- Limit flexibility to essential residues
- Increase sampling parameters
- Use appropriate clash tolerance
- Enable memory-efficient modes

This comprehensive guide provides the foundation for advanced flexible docking applications. Experiment with different flexibility configurations to find optimal settings for your specific systems!