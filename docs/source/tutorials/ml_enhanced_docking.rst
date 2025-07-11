ML-Enhanced Docking Tutorial
============================

This tutorial demonstrates how to use PandaDock's PandaML algorithm for superior affinity prediction and high-accuracy molecular docking.

Introduction
------------

PandaDock's PandaML algorithm uses advanced machine learning models trained on experimental binding data to provide:

- Higher accuracy than traditional scoring functions
- Uncertainty quantification for predictions
- Binding affinity estimation
- Chemical space-aware predictions

Prerequisites
-------------

- Completed the :doc:`basic_docking` tutorial
- Understanding of machine learning concepts
- NVIDIA GPU recommended (but not required)

Tutorial Overview
-----------------

You will learn to:

1. Configure ML-enhanced docking
2. Use pre-trained models
3. Interpret uncertainty estimates
4. Perform binding affinity prediction
5. Customize ML models for your targets

Step 1: Basic ML Configuration
------------------------------

Start with a simple ML-enhanced docking setup:

.. code-block:: python

   from pandadock import PandaDock
   
   # Initialize with PandaML algorithm
   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',    # transformer, cnn, graph
       ensemble_size=3,            # Number of models in ensemble
       uncertainty=True,           # Enable uncertainty quantification
       gpu_acceleration=True       # Use GPU if available
   )

Step 2: Model Selection
-----------------------

Choose the appropriate ML model based on your needs:

.. code-block:: python

   # Graph neural network (best for novel scaffolds)
   docker_gnn = PandaDock(
       engine='pandaml',
       model_type='graph',
       ml_config={
           'architecture': 'graphsage',
           'num_layers': 4,
           'hidden_size': 256,
           'attention_heads': 8,
           'dropout': 0.1
       }
   )
   
   # Convolutional neural network (fast, good for screening)
   docker_cnn = PandaDock(
       engine='pandaml',
       model_type='cnn',
       ml_config={
           'architecture': '3dcnn',
           'voxel_size': 0.5,
           'grid_size': 48,
           'num_filters': [32, 64, 128],
           'kernel_size': 3
       }
   )
   
   # Transformer (best overall accuracy with PandaML)
   docker_transformer = PandaDock(
       engine='pandaml',
       model_type='transformer',
       ml_config={
           'architecture': 'bert',
           'num_layers': 12,
           'hidden_size': 512,
           'attention_heads': 16,
           'intermediate_size': 2048
       }
   )

Step 3: Advanced ML Docking
----------------------------

Run ML-enhanced docking with comprehensive analysis:

.. code-block:: python

   # Configure PandaML for high-accuracy docking
   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',
       ensemble_size=5,
       uncertainty=True,
       ml_config={
           'calibration_temperature': 1.2,  # For better uncertainty
           'monte_carlo_dropout': True,     # Additional uncertainty
           'num_mc_samples': 10,           # Monte Carlo samples
           'feature_attribution': True,    # Explainable AI
           'confidence_threshold': 0.7     # Minimum confidence
       }
   )
   
   # Run docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze ML-specific results
   print("PandaML Algorithm Results:")
   print("=" * 40)
   
   for i, pose in enumerate(results.poses):
       print(f"\nPose {i+1}:")
       print(f"  Score: {pose.score:.3f}")
       print(f"  ML Score: {pose.ml_score:.3f}")
       print(f"  Confidence: {pose.confidence:.3f}")
       print(f"  Uncertainty: {pose.uncertainty:.3f}")
       print(f"  Predicted pKd: {pose.predicted_pkd:.2f}")
       print(f"  Predicted IC50: {pose.predicted_ic50:.2e} nM")
       
       # Feature importance
       if hasattr(pose, 'feature_importance'):
           print(f"  Top contributing features:")
           for feature, importance in pose.feature_importance.items():
               print(f"    {feature}: {importance:.3f}")

Step 4: Uncertainty Analysis
-----------------------------

Understand and interpret uncertainty estimates:

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Extract uncertainty information
   scores = [pose.ml_score for pose in results.poses]
   uncertainties = [pose.uncertainty for pose in results.poses]
   confidences = [pose.confidence for pose in results.poses]
   
   # Plot score vs uncertainty
   plt.figure(figsize=(12, 4))
   
   plt.subplot(1, 3, 1)
   plt.scatter(scores, uncertainties, alpha=0.7)
   plt.xlabel('ML Score')
   plt.ylabel('Uncertainty')
   plt.title('Score vs Uncertainty')
   
   plt.subplot(1, 3, 2)
   plt.scatter(scores, confidences, alpha=0.7)
   plt.xlabel('ML Score')
   plt.ylabel('Confidence')
   plt.title('Score vs Confidence')
   
   plt.subplot(1, 3, 3)
   plt.hist(uncertainties, bins=20, alpha=0.7)
   plt.xlabel('Uncertainty')
   plt.ylabel('Count')
   plt.title('Uncertainty Distribution')
   
   plt.tight_layout()
   plt.savefig('uncertainty_analysis.png')
   plt.show()
   
   # Identify high-confidence predictions
   high_confidence_poses = [
       pose for pose in results.poses 
       if pose.confidence > 0.8
   ]
   
   print(f"\nHigh-confidence poses: {len(high_confidence_poses)}")
   for pose in high_confidence_poses[:3]:
       print(f"  Score: {pose.ml_score:.3f}, Confidence: {pose.confidence:.3f}")

Step 5: Binding Affinity Prediction
------------------------------------

Use ML models to predict binding affinities:

.. code-block:: python

   # Configure PandaML for affinity prediction
   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',
       affinity_prediction=True,
       ml_config={
           'affinity_model': 'trained_on_pdbbind',
           'units': 'log_molar',
           'temperature': 298.15,
           'ph': 7.4
       }
   )
   
   # Run docking with affinity prediction
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze affinity predictions
   print("Binding Affinity Predictions:")
   print("-" * 40)
   
   for i, pose in enumerate(results.poses):
       kd = pose.predicted_kd
       ic50 = pose.predicted_ic50
       ki = pose.predicted_ki
       
       print(f"Pose {i+1}:")
       print(f"  Kd: {kd:.2e} M ({kd*1e9:.1f} nM)")
       print(f"  IC50: {ic50:.2e} M ({ic50*1e9:.1f} nM)")
       print(f"  Ki: {ki:.2e} M ({ki*1e9:.1f} nM)")
       print(f"  pKd: {pose.predicted_pkd:.2f}")
       print(f"  pIC50: {pose.predicted_pic50:.2f}")
       print(f"  ΔG: {pose.predicted_delta_g:.2f} kcal/mol")
       
       # Affinity uncertainty
       if hasattr(pose, 'affinity_uncertainty'):
           print(f"  Affinity uncertainty: ±{pose.affinity_uncertainty:.2f} log units")

Step 6: Model Ensembling
------------------------

Use multiple models for robust predictions:

.. code-block:: python

   # Configure PandaML ensemble
   docker = PandaDock(
       engine='pandaml',
       ensemble_config={
           'models': [
               {'type': 'transformer', 'weight': 0.4},
               {'type': 'graph', 'weight': 0.3},
               {'type': 'cnn', 'weight': 0.3}
           ],
           'voting': 'weighted',           # weighted, majority, average
           'uncertainty_aggregation': 'variance'  # variance, entropy
       }
   )
   
   # Run ensemble docking
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze ensemble results
   for pose in results.poses:
       print(f"Ensemble score: {pose.ensemble_score:.3f}")
       print(f"Model agreement: {pose.model_agreement:.3f}")
       print(f"Individual scores: {pose.individual_scores}")

Step 7: Custom Model Training
------------------------------

Train models on your specific data:

.. code-block:: python

   from pandadock.ml import ModelTrainer
   
   # Prepare training data
   training_data = [
       {
           'receptor': 'receptor1.pdb',
           'ligand': 'ligand1.sdf',
           'affinity': 7.5,  # pKd or pIC50
           'pose': 'pose1.pdb'
       },
       # ... more training examples
   ]
   
   # Configure trainer
   trainer = ModelTrainer(
       model_type='transformer',
       training_config={
           'batch_size': 16,
           'learning_rate': 1e-4,
           'num_epochs': 100,
           'validation_split': 0.2,
           'early_stopping': True,
           'patience': 10
       }
   )
   
   # Train custom model
   model = trainer.train(
       training_data=training_data,
       validation_data=validation_data,
       output_dir='custom_model'
   )
   
   # Use custom model for docking
   docker = PandaDock(
       engine='pandaml',
       model_path='custom_model/best_model.pt'
   )

Step 8: Feature Attribution
----------------------------

Understand what the model is learning:

.. code-block:: python

   # Enable feature attribution
   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',
       explainable_ai=True,
       ml_config={
           'attribution_method': 'integrated_gradients',
           'attribution_baseline': 'zero',
           'num_attribution_steps': 50
       }
   )
   
   # Run docking with attribution
   results = docker.dock(
       receptor='protein.pdb',
       ligand='ligand.sdf',
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Analyze feature importance
   best_pose = results.best_pose
   attributions = best_pose.feature_attributions
   
   print("Feature Attributions:")
   print("-" * 30)
   
   # Sort by importance
   sorted_features = sorted(
       attributions.items(), 
       key=lambda x: abs(x[1]), 
       reverse=True
   )
   
   for feature, importance in sorted_features[:10]:
       print(f"{feature:20s}: {importance:+.4f}")
   
   # Visualize attributions
   from pandadock.visualization import plot_attributions
   
   plot_attributions(
       attributions, 
       output_file='feature_attributions.png'
   )

Step 9: Comparison with Traditional Scoring
--------------------------------------------

Compare ML and traditional scoring:

.. code-block:: python

   # Run both PandaML and PandaPhysics algorithms
   docker_ml = PandaDock(engine='pandaml', model_type='transformer')
   docker_physics = PandaDock(engine='pandaphysics')
   
   # Same docking parameters
   dock_params = {
       'receptor': 'protein.pdb',
       'ligand': 'ligand.sdf',
       'center': [25.0, 30.0, 15.0],
       'size': [20.0, 20.0, 20.0]
   }
   
   # Run both
   results_ml = docker_ml.dock(**dock_params)
   results_physics = docker_physics.dock(**dock_params)
   
   # Compare results
   print("Comparison of PandaML vs PandaPhysics Algorithms:")
   print("=" * 50)
   
   print(f"PandaML best score: {results_ml.best_pose.score:.3f}")
   print(f"PandaPhysics best score: {results_physics.best_pose.score:.3f}")
   
   print(f"PandaML runtime: {results_ml.runtime:.2f} seconds")
   print(f"PandaPhysics runtime: {results_physics.runtime:.2f} seconds")
   
   # Score correlation
   import numpy as np
   from scipy.stats import pearsonr
   
   ml_scores = [pose.score for pose in results_ml.poses]
   physics_scores = [pose.score for pose in results_physics.poses]
   
   # Align poses (assuming same number)
   if len(ml_scores) == len(physics_scores):
       correlation, p_value = pearsonr(ml_scores, physics_scores)
       print(f"Score correlation: {correlation:.3f} (p={p_value:.3f})")

Step 10: Production Deployment
-------------------------------

Deploy ML models in production:

.. code-block:: python

   # Production configuration with PandaML
   docker = PandaDock(
       engine='pandaml',
       model_type='transformer',
       production_config={
           'batch_size': 32,           # Optimize for throughput
           'mixed_precision': True,    # Faster inference
           'model_optimization': True, # Optimize model for inference
           'cache_features': True,     # Cache computed features
           'parallel_workers': 4,      # Parallel processing
           'memory_efficient': True    # Reduce memory usage
       }
   )
   
   # Batch processing multiple ligands
   ligands = ['ligand1.sdf', 'ligand2.sdf', 'ligand3.sdf']
   
   batch_results = docker.dock_batch(
       receptor='protein.pdb',
       ligands=ligands,
       center=[25.0, 30.0, 15.0],
       size=[20.0, 20.0, 20.0]
   )
   
   # Process results
   for ligand, results in batch_results.items():
       print(f"Ligand: {ligand}")
       print(f"  Best score: {results.best_pose.score:.3f}")
       print(f"  Confidence: {results.best_pose.confidence:.3f}")

Complete Example
----------------

Here's a comprehensive example combining all techniques:

.. code-block:: python

   #!/usr/bin/env python3
   """
   Complete ML-Enhanced Docking Example
   """
   
   from pandadock import PandaDock
   import numpy as np
   import matplotlib.pyplot as plt
   
   def ml_enhanced_docking(receptor_file, ligand_file, binding_site):
       """
       Perform ML-enhanced docking with comprehensive analysis
       """
       
       # Configure PandaML algorithm
       docker = PandaDock(
           engine='pandaml',
           model_type='transformer',
           ensemble_size=3,
           uncertainty=True,
           affinity_prediction=True,
           explainable_ai=True,
           ml_config={
               'calibration_temperature': 1.2,
               'monte_carlo_dropout': True,
               'num_mc_samples': 20,
               'confidence_threshold': 0.7
           }
       )
       
       # Run docking
       print("Running PandaML algorithm...")
       results = docker.dock(
           receptor=receptor_file,
           ligand=ligand_file,
           center=binding_site['center'],
           size=binding_site['size']
       )
       
       # Analyze results
       print(f"\nDocking completed in {results.runtime:.2f} seconds")
       print(f"Generated {len(results.poses)} poses")
       
       # Best pose analysis
       best_pose = results.best_pose
       print(f"\nBest Pose Analysis:")
       print(f"  PandaML Score: {best_pose.ml_score:.3f}")
       print(f"  Confidence: {best_pose.confidence:.3f}")
       print(f"  Uncertainty: {best_pose.uncertainty:.3f}")
       print(f"  Predicted pKd: {best_pose.predicted_pkd:.2f}")
       print(f"  Predicted IC50: {best_pose.predicted_ic50:.2e} nM")
       
       # High-confidence poses
       high_conf_poses = [
           pose for pose in results.poses 
           if pose.confidence > 0.8
       ]
       print(f"\nHigh-confidence poses: {len(high_conf_poses)}")
       
       # Feature importance
       if hasattr(best_pose, 'feature_attributions'):
           print("\nTop feature contributions:")
           sorted_features = sorted(
               best_pose.feature_attributions.items(),
               key=lambda x: abs(x[1]),
               reverse=True
           )
           for feature, importance in sorted_features[:5]:
               print(f"  {feature}: {importance:+.4f}")
       
       # Save results
       results.save_poses("pandaml_docking_poses.sdf")
       results.save_report("pandaml_docking_report.html")
       
       return results
   
   def main():
       # Example usage
       receptor_file = "protein.pdb"
       ligand_file = "ligand.sdf"
       binding_site = {
           'center': [25.0, 30.0, 15.0],
           'size': [20.0, 20.0, 20.0]
       }
       
       try:
           results = ml_enhanced_docking(receptor_file, ligand_file, binding_site)
           print("\nPandaML docking completed successfully!")
           
       except Exception as e:
           print(f"Error: {e}")
           return 1
       
       return 0
   
   if __name__ == "__main__":
       exit(main())

Performance Considerations
--------------------------

**GPU Acceleration:**
- Use NVIDIA GPU for 5-10x speedup
- Increase batch_size for GPU efficiency
- Enable mixed precision training

**Memory Optimization:**
- Use gradient checkpointing for large models
- Process ligands in batches
- Clear GPU cache between runs

**PandaML Model Selection:**
- Transformer: Best accuracy (R² = 0.845), slower
- Graph: Good for novel scaffolds
- CNN: Fastest, good for screening

Next Steps
----------

After completing this tutorial:

- Explore :doc:`virtual_screening` with PandaML algorithm
- Learn about :doc:`../examples/custom_ml_models`
- Check out :doc:`../user_guide/model_selection`
- Try :doc:`active_learning` for model improvement

The PandaML algorithm provides state-of-the-art accuracy for molecular docking applications with superior affinity prediction (R² = 0.878). Experiment with different models and configurations to find the best setup for your specific use case!