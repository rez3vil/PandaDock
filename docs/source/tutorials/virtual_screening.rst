Virtual Screening Tutorial
=========================

This tutorial covers high-throughput virtual screening using PandaDock, from library preparation to hit analysis.

Introduction
------------

Virtual screening allows you to computationally evaluate thousands or millions of compounds against a target protein to identify potential drug candidates. PandaDock provides optimized workflows for efficient large-scale screening.

Prerequisites
-------------

- Completed :doc:`basic_docking` tutorial
- Understanding of drug discovery concepts
- Access to compound libraries (ZINC, ChEMBL, etc.)
- Sufficient computational resources

Tutorial Overview
-----------------

You will learn to:

1. Prepare compound libraries
2. Configure high-throughput screening
3. Run parallel screening campaigns
4. Analyze and filter results
5. Perform focused screening
6. Validate hits experimentally

Step 1: Library Preparation
---------------------------

Start by preparing your virtual compound library:

.. code-block:: python

   from pandadock import VirtualScreening
   from pandadock.utils import LibraryPreparation
   import pandas as pd
   
   # Initialize library preparation
   lib_prep = LibraryPreparation()
   
   # Load compound library (example: ZINC15)
   library_path = "zinc15_druglike.sdf"
   
   # Basic library statistics
   stats = lib_prep.analyze_library(library_path)
   print(f"Library contains {stats['num_compounds']} compounds")
   print(f"Average molecular weight: {stats['avg_mw']:.1f} Da")
   print(f"MW range: {stats['min_mw']:.1f} - {stats['max_mw']:.1f} Da")
   print(f"Lipinski violations: {stats['lipinski_violations']}")

Filter the library based on drug-likeness criteria:

.. code-block:: python

   # Define filtering criteria
   filters = {
       'molecular_weight': (150, 500),      # Da
       'logp': (-1, 5),                     # Partition coefficient
       'hbd': (0, 5),                       # H-bond donors
       'hba': (0, 10),                      # H-bond acceptors
       'rotatable_bonds': (0, 10),          # Flexibility
       'tpsa': (0, 140),                    # Polar surface area
       'formal_charge': (-2, 2),            # Charge range
       'heavy_atoms': (10, 50),             # Size range
       'aromatic_rings': (0, 5),            # Aromaticity
       'aliphatic_rings': (0, 5),           # Flexibility
       'heteroatoms': (0, 15),              # Heteroatom content
       'lipinski_violations': (0, 1),       # Allow 1 violation
       'veber_violations': (0, 0),          # Strict Veber compliance
       'reactive_groups': False,            # No reactive groups
       'pan_assay_interference': False      # No PAINS
   }
   
   # Apply filters
   filtered_library = lib_prep.filter_library(
       library_path, 
       filters=filters,
       output_file="filtered_library.sdf"
   )
   
   print(f"Filtered library: {len(filtered_library)} compounds")
   print(f"Filtering efficiency: {len(filtered_library)/stats['num_compounds']*100:.1f}%")

Step 2: High-Throughput Screening Setup
----------------------------------------

Configure PandaDock for high-throughput screening:

.. code-block:: python

   # Configure for high-throughput screening
   screener = VirtualScreening(
       engine='pandaphysics',               # Fast PandaPhysics engine
       screening_mode='hts',                # High-throughput mode
       config={
           'exhaustiveness': 4,             # Lower for speed
           'num_poses': 1,                  # Single pose per compound
           'energy_range': 2.0,             # Narrow energy window
           'early_termination': True,       # Stop poor binders early
           'score_threshold': -6.0,         # Pre-filter threshold
           'timeout_per_ligand': 60,        # Max time per ligand (seconds)
           'parallel_screening': True,      # Enable parallelization
           'num_workers': 8,                # Parallel workers
           'batch_size': 100,               # Compounds per batch
           'memory_limit': '16GB',          # Memory management
           'checkpoint_interval': 1000,     # Save progress every N compounds
           'restart_failed': True           # Retry failed compounds
       }
   )

Step 3: Defining the Screening Campaign
---------------------------------------

Set up the screening parameters:

.. code-block:: python

   # Define screening parameters
   screening_params = {
       'receptor': 'target_protein.pdb',
       'library': 'filtered_library.sdf',
       'binding_site': {
           'center': [25.0, 30.0, 15.0],
           'size': [20.0, 20.0, 20.0]
       },
       'output_dir': 'screening_results',
       'project_name': 'target_screening_2024'
   }
   
   # Advanced screening configuration
   advanced_config = {
       'diversity_selection': True,         # Maintain structural diversity
       'similarity_threshold': 0.7,         # Tanimoto similarity cutoff
       'pharmacophore_filtering': True,     # Use pharmacophore models
       'shape_screening': True,             # Include shape similarity
       'decoy_generation': True,            # Generate decoys for validation
       'multiple_conformers': False,        # Use single conformer for speed
       'binding_site_flexibility': False,   # Rigid receptor for speed
       'scoring_functions': ['pandacore', 'pandaml'], # Multiple scoring functions
       'consensus_scoring': True            # Combine multiple scores
   }

Step 4: Running the Screening Campaign
--------------------------------------

Execute the virtual screening:

.. code-block:: python

   # Start screening campaign
   print("Starting virtual screening campaign...")
   
   # Run screening with progress monitoring
   results = screener.run_campaign(
       **screening_params,
       advanced_config=advanced_config,
       verbose=True,
       log_file='screening.log'
   )
   
   # Monitor progress
   def progress_callback(completed, total, current_compound):
       progress = (completed / total) * 100
       print(f"Progress: {progress:.1f}% ({completed}/{total}) - {current_compound}")
   
   # Run with progress monitoring
   results = screener.run_campaign(
       **screening_params,
       progress_callback=progress_callback
   )
   
   print(f"Screening completed in {results.total_runtime:.2f} seconds")
   print(f"Screened {results.total_compounds} compounds")
   print(f"Success rate: {results.success_rate:.1f}%")

Step 5: Results Analysis
------------------------

Analyze the screening results:

.. code-block:: python

   # Load and analyze results
   results_df = pd.read_csv('screening_results/screening_results.csv')
   
   print("Screening Results Summary:")
   print("=" * 30)
   print(f"Total compounds screened: {len(results_df)}")
   print(f"Mean docking score: {results_df['score'].mean():.2f}")
   print(f"Standard deviation: {results_df['score'].std():.2f}")
   print(f"Best score: {results_df['score'].min():.2f}")
   print(f"Worst score: {results_df['score'].max():.2f}")
   
   # Score distribution analysis
   import matplotlib.pyplot as plt
   import numpy as np
   
   plt.figure(figsize=(15, 5))
   
   # Score histogram
   plt.subplot(1, 3, 1)
   plt.hist(results_df['score'], bins=50, alpha=0.7, edgecolor='black')
   plt.xlabel('Docking Score')
   plt.ylabel('Frequency')
   plt.title('Score Distribution')
   plt.axvline(results_df['score'].mean(), color='red', linestyle='--', label='Mean')
   plt.legend()
   
   # Score vs molecular weight
   plt.subplot(1, 3, 2)
   plt.scatter(results_df['molecular_weight'], results_df['score'], alpha=0.5)
   plt.xlabel('Molecular Weight (Da)')
   plt.ylabel('Docking Score')
   plt.title('Score vs Molecular Weight')
   
   # Score vs LogP
   plt.subplot(1, 3, 3)
   plt.scatter(results_df['logp'], results_df['score'], alpha=0.5)
   plt.xlabel('LogP')
   plt.ylabel('Docking Score')
   plt.title('Score vs LogP')
   
   plt.tight_layout()
   plt.savefig('screening_analysis.png', dpi=300)
   plt.show()

Step 6: Hit Identification and Filtering
-----------------------------------------

Identify and filter potential hits:

.. code-block:: python

   # Define hit criteria
   hit_criteria = {
       'score_threshold': -8.0,           # Docking score cutoff
       'efficiency_threshold': 0.3,       # Ligand efficiency
       'lipinski_compliant': True,        # Drug-likeness
       'similarity_filter': 0.7,          # Remove similar compounds
       'visual_inspection': True,         # Flag for visual inspection
       'interaction_requirements': {      # Required interactions
           'hbonds': 1,                   # Minimum H-bonds
           'hydrophobic': 2,              # Minimum hydrophobic contacts
           'aromatic': 1                  # Minimum aromatic interactions
       }
   }
   
   # Apply hit identification
   hits = screener.identify_hits(results_df, hit_criteria)
   
   print(f"Identified {len(hits)} potential hits")
   print(f"Hit rate: {len(hits)/len(results_df)*100:.2f}%")
   
   # Rank hits by multiple criteria
   hits_ranked = screener.rank_hits(
       hits,
       criteria=['score', 'efficiency', 'diversity', 'interactions'],
       weights=[0.4, 0.3, 0.2, 0.1]
   )
   
   # Display top hits
   print("\nTop 10 Hits:")
   print("-" * 50)
   for i, hit in enumerate(hits_ranked[:10]):
       print(f"{i+1:2d}. {hit['compound_id']:15s} Score: {hit['score']:6.2f} "
             f"LE: {hit['efficiency']:.3f} MW: {hit['molecular_weight']:6.1f}")

Step 7: Cluster Analysis
------------------------

Perform clustering to identify diverse scaffolds:

.. code-block:: python

   from pandadock.analysis import ClusterAnalysis
   
   # Perform clustering
   clusterer = ClusterAnalysis(
       method='butina',              # Butina clustering
       similarity_metric='tanimoto',  # Tanimoto coefficient
       threshold=0.7,                # Similarity threshold
       min_cluster_size=2            # Minimum cluster size
   )
   
   # Cluster hits
   clusters = clusterer.cluster_compounds(hits_ranked)
   
   print(f"Identified {len(clusters)} clusters")
   
   # Analyze clusters
   for i, cluster in enumerate(clusters):
       print(f"\nCluster {i+1}: {len(cluster)} compounds")
       print(f"  Representative: {cluster[0]['compound_id']}")
       print(f"  Score range: {min(c['score'] for c in cluster):.2f} - "
             f"{max(c['score'] for c in cluster):.2f}")
       print(f"  Scaffold diversity: {cluster[0]['scaffold_diversity']:.2f}")
   
   # Select diverse representatives
   diverse_hits = clusterer.select_diverse_representatives(
       clusters,
       selection_method='best_score',
       max_per_cluster=3
   )
   
   print(f"\nSelected {len(diverse_hits)} diverse hits for further analysis")

Step 8: Focused Screening
-------------------------

Perform focused screening around promising hits:

.. code-block:: python

   # Generate focused libraries around hits
   from pandadock.utils import FocusedLibraryGenerator
   
   generator = FocusedLibraryGenerator()
   
   # Generate analogs for top hits
   focused_libraries = []
   for hit in hits_ranked[:5]:  # Top 5 hits
       
       # Generate structural analogs
       analogs = generator.generate_analogs(
           hit['smiles'],
           methods=['fragment_replacement', 'functional_group_addition', 
                   'ring_expansion', 'bioisosteric_replacement'],
           max_analogs=100,
           similarity_range=(0.6, 0.9)
       )
       
       # Create focused library
       focused_lib = generator.create_focused_library(
           core_structure=hit['smiles'],
           analogs=analogs,
           diversity_filter=True,
           property_filter=filters
       )
       
       focused_libraries.append({
           'parent': hit['compound_id'],
           'library': focused_lib,
           'size': len(focused_lib)
       })
   
   # Screen focused libraries
   focused_results = []
   for lib_info in focused_libraries:
       print(f"Screening {lib_info['size']} analogs of {lib_info['parent']}...")
       
       # Configure for focused screening (higher accuracy)
       focused_screener = VirtualScreening(
           engine='pandaml',               # Use PandaML for higher accuracy
           screening_mode='focused',
           config={
               'exhaustiveness': 16,       # Higher exhaustiveness
               'num_poses': 5,             # More poses
               'uncertainty': True,        # Include uncertainty
               'affinity_prediction': True  # Predict binding affinity
           }
       )
       
       # Run focused screening
       focused_result = focused_screener.run_campaign(
           receptor=screening_params['receptor'],
           library=lib_info['library'],
           binding_site=screening_params['binding_site'],
           output_dir=f"focused_results/{lib_info['parent']}"
       )
       
       focused_results.append(focused_result)

Step 9: Experimental Validation Planning
-----------------------------------------

Plan experimental validation of top hits:

.. code-block:: python

   # Prioritize compounds for experimental testing
   validation_candidates = screener.prioritize_for_validation(
       hits_ranked,
       criteria={
           'commercial_availability': True,    # Must be purchasable
           'price_threshold': 1000,           # Max price per mg
           'purity_threshold': 95,            # Min purity %
           'stability_prediction': True,      # Predict stability
           'toxicity_prediction': True,       # Predict toxicity
           'synthetic_accessibility': True,   # Ease of synthesis
           'patent_check': True              # Check patent status
       },
       max_candidates=50
   )
   
   # Generate experimental design
   experimental_design = screener.design_validation_experiments(
       validation_candidates,
       assay_types=['binding', 'enzymatic', 'cellular'],
       controls=['positive', 'negative', 'vehicle'],
       replicates=3,
       concentration_range=(1e-8, 1e-4)  # M
   )
   
   # Export for experimental team
   experimental_design.to_csv('validation_experiments.csv')
   
   print(f"Designed validation experiments for {len(validation_candidates)} compounds")
   print("Experimental design saved to validation_experiments.csv")

Step 10: Results Reporting
--------------------------

Generate comprehensive screening reports:

.. code-block:: python

   # Generate comprehensive report
   report = screener.generate_screening_report(
       results_df,
       hits_ranked,
       clusters,
       focused_results,
       validation_candidates,
       include_plots=True,
       include_structures=True,
       format='html'
   )
   
   # Save report
   report.save('virtual_screening_report.html')
   
   # Generate summary for management
   summary = screener.generate_executive_summary(
       results_df,
       hits_ranked,
       validation_candidates,
       timeline=results.timeline,
       costs=results.computational_costs
   )
   
   print("Executive Summary:")
   print("=" * 20)
   print(summary)

Complete Screening Pipeline
---------------------------

Here's a complete pipeline script:

.. code-block:: python

   #!/usr/bin/env python3
   """
   Complete Virtual Screening Pipeline
   
   This script performs a complete virtual screening campaign
   from library preparation to hit identification.
   """
   
   import os
   import sys
   import pandas as pd
   from datetime import datetime
   
   from pandadock import VirtualScreening
   from pandadock.utils import LibraryPreparation, FocusedLibraryGenerator
   from pandadock.analysis import ClusterAnalysis
   
   class ScreeningPipeline:
       def __init__(self, config_file):
           self.config = self.load_config(config_file)
           self.screener = None
           self.results = None
           
       def load_config(self, config_file):
           """Load screening configuration"""
           # Implementation depends on your config format
           pass
           
       def prepare_library(self):
           """Prepare and filter compound library"""
           print("Step 1: Library Preparation")
           print("-" * 30)
           
           lib_prep = LibraryPreparation()
           
           # Analyze original library
           stats = lib_prep.analyze_library(self.config['library_path'])
           print(f"Original library: {stats['num_compounds']} compounds")
           
           # Apply filters
           filtered = lib_prep.filter_library(
               self.config['library_path'],
               filters=self.config['filters'],
               output_file='filtered_library.sdf'
           )
           
           print(f"Filtered library: {len(filtered)} compounds")
           return filtered
           
       def run_screening(self):
           """Run virtual screening campaign"""
           print("\nStep 2: Virtual Screening")
           print("-" * 30)
           
           # Configure screener
           self.screener = VirtualScreening(
               engine=self.config['engine'],
               screening_mode=self.config['screening_mode'],
               config=self.config['screening_config']
           )
           
           # Run campaign
           self.results = self.screener.run_campaign(
               receptor=self.config['receptor'],
               library='filtered_library.sdf',
               binding_site=self.config['binding_site'],
               output_dir='screening_results',
               verbose=True
           )
           
           print(f"Screening completed: {self.results.total_compounds} compounds")
           return self.results
           
       def analyze_results(self):
           """Analyze screening results"""
           print("\nStep 3: Results Analysis")
           print("-" * 30)
           
           # Load results
           results_df = pd.read_csv('screening_results/screening_results.csv')
           
           # Identify hits
           hits = self.screener.identify_hits(
               results_df, 
               self.config['hit_criteria']
           )
           
           # Rank hits
           hits_ranked = self.screener.rank_hits(
               hits,
               criteria=self.config['ranking_criteria'],
               weights=self.config['ranking_weights']
           )
           
           print(f"Identified {len(hits)} hits")
           return hits_ranked
           
       def cluster_hits(self, hits):
           """Perform clustering analysis"""
           print("\nStep 4: Clustering Analysis")
           print("-" * 30)
           
           clusterer = ClusterAnalysis(
               method='butina',
               threshold=0.7
           )
           
           clusters = clusterer.cluster_compounds(hits)
           diverse_hits = clusterer.select_diverse_representatives(clusters)
           
           print(f"Identified {len(clusters)} clusters")
           print(f"Selected {len(diverse_hits)} diverse hits")
           return diverse_hits
           
       def focused_screening(self, hits):
           """Perform focused screening"""
           print("\nStep 5: Focused Screening")
           print("-" * 30)
           
           generator = FocusedLibraryGenerator()
           
           # Generate focused libraries
           focused_results = []
           for hit in hits[:3]:  # Top 3 hits
               analogs = generator.generate_analogs(
                   hit['smiles'],
                   max_analogs=50
               )
               
               # Screen analogs
               focused_result = self.screener.run_campaign(
                   receptor=self.config['receptor'],
                   library=analogs,
                   binding_site=self.config['binding_site'],
                   output_dir=f"focused_results/{hit['compound_id']}"
               )
               
               focused_results.append(focused_result)
           
           return focused_results
           
       def generate_report(self, hits, clusters, focused_results):
           """Generate final report"""
           print("\nStep 6: Report Generation")
           print("-" * 30)
           
           report = self.screener.generate_screening_report(
               self.results,
               hits,
               clusters,
               focused_results,
               output_file='screening_report.html'
           )
           
           print("Report generated: screening_report.html")
           return report
           
       def run_pipeline(self):
           """Run complete screening pipeline"""
           start_time = datetime.now()
           
           try:
               # Run pipeline steps
               filtered_library = self.prepare_library()
               results = self.run_screening()
               hits = self.analyze_results()
               diverse_hits = self.cluster_hits(hits)
               focused_results = self.focused_screening(diverse_hits)
               report = self.generate_report(hits, diverse_hits, focused_results)
               
               # Calculate total time
               end_time = datetime.now()
               total_time = (end_time - start_time).total_seconds()
               
               print(f"\nPipeline completed in {total_time:.2f} seconds")
               print(f"Results saved to: screening_report.html")
               
               return True
               
           except Exception as e:
               print(f"Pipeline failed: {e}")
               return False
   
   def main():
       if len(sys.argv) != 2:
           print("Usage: python screening_pipeline.py config.yaml")
           sys.exit(1)
       
       config_file = sys.argv[1]
       
       # Run pipeline
       pipeline = ScreeningPipeline(config_file)
       success = pipeline.run_pipeline()
       
       if success:
           print("Virtual screening pipeline completed successfully!")
           sys.exit(0)
       else:
           print("Pipeline failed!")
           sys.exit(1)
   
   if __name__ == "__main__":
       main()

Performance Optimization
------------------------

**For Large Libraries (>100K compounds):**

1. **Use PandaPhysics Engine**: Fastest for initial screening
2. **Reduce Exhaustiveness**: Use 2-4 for initial screening
3. **Single Pose Output**: Reduce to 1 pose per compound
4. **Parallel Processing**: Use all available CPU cores
5. **Batch Processing**: Process in batches of 1000-10000
6. **Memory Management**: Set appropriate memory limits

**For High-Accuracy Screening:**

1. **Use PandaML Engine**: Best accuracy but slower
2. **Increase Exhaustiveness**: Use 16-32
3. **Multiple Poses**: Generate 5-10 poses per compound
4. **Ensemble Methods**: Use multiple models
5. **Uncertainty Quantification**: Include confidence estimates

Best Practices
--------------

1. **Library Preparation**: Always filter libraries for drug-likeness
2. **Validation**: Use known actives/inactives to validate setup
3. **Diversity**: Maintain structural diversity in hit sets
4. **Documentation**: Keep detailed records of parameters and results
5. **Experimental Validation**: Always validate computationally identified hits

Next Steps
----------

- Learn about :doc:`../examples/custom_libraries`
- Explore :doc:`../user_guide/result_analysis`
- Check out :doc:`structure_based_optimization`
- Try :doc:`../examples/consensus_scoring`

Virtual screening is a powerful tool for drug discovery when properly executed. Use this tutorial as a foundation for your own screening campaigns!