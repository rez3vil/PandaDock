Comprehensive Benchmark Tutorial
=================================

This tutorial demonstrates how to reproduce the comprehensive PandaDock benchmark results on the complete PDBbind database.

Overview
--------

The comprehensive benchmark evaluates all three PandaDock algorithms on 5,316 protein-ligand complexes from the complete PDBbind database, representing the most extensive molecular docking evaluation performed with PandaDock.

Prerequisites
-------------

- PandaDock installed with all algorithms
- Access to PDBbind database
- Sufficient computational resources (recommended: 8+ CPU cores, 16+ GB RAM)
- Python packages: pandas, numpy, matplotlib, seaborn

Running the Complete Benchmark
-------------------------------

Step 1: Setup
~~~~~~~~~~~~~

First, organize your PDBbind data:

.. code-block:: bash

   # Create benchmark directory
   mkdir comprehensive_benchmark_results
   cd comprehensive_benchmark_results

   # Ensure PDBbind structure
   ls /path/to/pdbbind/
   # Should contain: refined-set/, general-set/, index/

Step 2: Execute Comprehensive Benchmark
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Navigate to benchmark scripts
   cd /path/to/pandadock/benchmarks/scripts

   # Run comprehensive benchmark (this will take several hours)
   python comprehensive_benchmark.py \
       --pdbbind_dir /path/to/pdbbind \
       --output_dir comprehensive_benchmark_results \
       --algorithms pandacore,pandaml,pandaphysics \
       --num_poses 10 \
       --exhaustiveness 8

Step 3: Monitor Progress
~~~~~~~~~~~~~~~~~~~~~~~~

The benchmark will process all 5,316 complexes:

.. code-block:: bash

   # Monitor progress
   tail -f comprehensive_benchmark_results/benchmark.log

   # Check intermediate results
   ls comprehensive_benchmark_results/
   # Should show: raw_data/, plots/, reports/

Understanding the Results
-------------------------

Performance Metrics
~~~~~~~~~~~~~~~~~~~

The benchmark evaluates three key areas:

**1. Binding Affinity Prediction:**
   - R² (coefficient of determination)
   - Pearson correlation coefficient
   - RMSE (root mean square error)
   - MAE (mean absolute error)

**2. Pose Prediction Accuracy:**
   - RMSD to crystal structure
   - Success rate (RMSD < 2Å)
   - Median RMSD
   - RMSD distribution

**3. Computational Efficiency:**
   - Mean docking time per complex
   - Time per heavy atom
   - Scalability analysis

Generated Output Files
~~~~~~~~~~~~~~~~~~~~~~

The benchmark generates comprehensive output:

.. code-block:: text

   comprehensive_benchmark_results/
   ├── benchmark_raw_data.csv          # All raw results
   ├── benchmark_results.json          # Summary statistics
   ├── comprehensive_benchmark_report.md  # Detailed report
   ├── master_publication_figure.png   # Main results figure
   ├── correlation_analysis.png        # Affinity prediction plots
   ├── rmsd_analysis.png              # Pose accuracy analysis
   ├── engine_performance.png         # Algorithm comparison
   ├── performance_vs_properties.png  # Property dependence
   └── ligand_complexity_analysis.png # Size/complexity effects

Analyzing Specific Results
--------------------------

Algorithm Performance Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns

   # Load results
   data = pd.read_csv('comprehensive_benchmark_results/benchmark_raw_data.csv')

   # PandaML performance
   pandaml_data = data[data['engine_type'] == 'pandaml']
   print(f"PandaML R²: {pandaml_data['predicted_affinity'].corr(pandaml_data['experimental_affinity'])**2:.3f}")
   print(f"PandaML Success Rate: {(pandaml_data['rmsd_best_pose'] < 2.0).mean():.3f}")

   # Compare algorithms
   performance_summary = data.groupby('engine_type').agg({
       'predicted_affinity': lambda x: x.corr(data.loc[x.index, 'experimental_affinity'])**2,
       'rmsd_best_pose': ['mean', lambda x: (x < 2.0).mean()],
       'docking_time': 'mean'
   })
   print(performance_summary)

Visualizing Results
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Correlation plots
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))
   
   algorithms = ['pandacore', 'pandaml', 'pandaphysics']
   
   for i, alg in enumerate(algorithms):
       alg_data = data[data['engine_type'] == alg]
       axes[i].scatter(alg_data['experimental_affinity'], 
                      alg_data['predicted_affinity'], 
                      alpha=0.6)
       axes[i].plot([4, 11], [4, 11], 'k--')
       axes[i].set_title(f'{alg.upper()} Algorithm')
       axes[i].set_xlabel('Experimental Affinity (pKd/pKi)')
       axes[i].set_ylabel('Predicted Affinity')
   
   plt.tight_layout()
   plt.savefig('algorithm_comparison.png', dpi=300)

Custom Analysis Examples
------------------------

Affinity Range Analysis
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Analyze performance by affinity range
   data['affinity_range'] = pd.cut(data['experimental_affinity'], 
                                   bins=[0, 6, 8, 12], 
                                   labels=['Low', 'Medium', 'High'])

   range_performance = data.groupby(['engine_type', 'affinity_range']).agg({
       'predicted_affinity': lambda x: x.corr(data.loc[x.index, 'experimental_affinity'])**2,
       'rmsd_best_pose': lambda x: (x < 2.0).mean()
   })

   print("Performance by Affinity Range:")
   print(range_performance)

Ligand Size Dependence
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Performance vs ligand size
   data['ligand_size_range'] = pd.cut(data['ligand_atoms'], 
                                      bins=[0, 30, 50, 100], 
                                      labels=['Small', 'Medium', 'Large'])

   size_performance = data.groupby(['engine_type', 'ligand_size_range']).agg({
       'rmsd_best_pose': 'mean',
       'docking_time': 'mean'
   })

   print("Performance by Ligand Size:")
   print(size_performance)

Statistical Analysis
--------------------

Significance Testing
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from scipy import stats

   # Compare algorithm performance
   pandaml_rmsd = data[data['engine_type'] == 'pandaml']['rmsd_best_pose']
   pandacore_rmsd = data[data['engine_type'] == 'pandacore']['rmsd_best_pose']
   pandaphysics_rmsd = data[data['engine_type'] == 'pandaphysics']['rmsd_best_pose']

   # Wilcoxon rank-sum tests
   stat1, p1 = stats.ranksums(pandaml_rmsd, pandacore_rmsd)
   stat2, p2 = stats.ranksums(pandaml_rmsd, pandaphysics_rmsd)
   stat3, p3 = stats.ranksums(pandacore_rmsd, pandaphysics_rmsd)

   print(f"PandaML vs PandaCore: p = {p1:.4f}")
   print(f"PandaML vs PandaPhysics: p = {p2:.4f}")
   print(f"PandaCore vs PandaPhysics: p = {p3:.4f}")

Cross-Validation Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from sklearn.model_selection import KFold
   from sklearn.metrics import r2_score

   # 5-fold cross-validation for PandaML
   kf = KFold(n_splits=5, shuffle=True, random_state=42)
   pandaml_data = data[data['engine_type'] == 'pandaml']

   cv_scores = []
   for train_idx, test_idx in kf.split(pandaml_data):
       test_data = pandaml_data.iloc[test_idx]
       r2 = r2_score(test_data['experimental_affinity'], 
                     test_data['predicted_affinity'])
       cv_scores.append(r2)

   print(f"PandaML CV R²: {np.mean(cv_scores):.3f} ± {np.std(cv_scores):.3f}")

Reproducing Specific Figures
-----------------------------

Master Publication Figure
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Reproduce the main benchmark figure
   from benchmark_analysis import create_master_figure

   create_master_figure(
       data_file='comprehensive_benchmark_results/benchmark_raw_data.csv',
       output_file='master_figure_reproduction.png'
   )

Performance Comparison
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Create algorithm performance comparison
   algorithms = ['pandacore', 'pandaml', 'pandaphysics']
   metrics = []

   for alg in algorithms:
       alg_data = data[data['engine_type'] == alg]
       r2 = alg_data['predicted_affinity'].corr(alg_data['experimental_affinity'])**2
       success_rate = (alg_data['rmsd_best_pose'] < 2.0).mean()
       mean_rmsd = alg_data['rmsd_best_pose'].mean()
       mean_time = alg_data['docking_time'].mean()
       
       metrics.append({
           'Algorithm': alg.upper(),
           'R²': r2,
           'Success Rate': success_rate,
           'Mean RMSD': mean_rmsd,
           'Mean Time': mean_time
       })

   metrics_df = pd.DataFrame(metrics)
   print(metrics_df.round(3))

Best Practices
--------------

Resource Management
~~~~~~~~~~~~~~~~~~~

- **CPU cores**: Use all available cores for parallel processing
- **Memory**: Monitor RAM usage, especially with large datasets
- **Storage**: Ensure sufficient disk space for results (~10GB)
- **Time**: Allow 12-24 hours for complete benchmark on modern hardware

Quality Control
~~~~~~~~~~~~~~~

- **Data validation**: Verify PDBbind structure completeness
- **Progress monitoring**: Check intermediate results regularly
- **Error handling**: Review log files for any failed complexes
- **Result validation**: Compare with published benchmark values

Customization Options
---------------------

Algorithm Subset
~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Run only specific algorithms
   python comprehensive_benchmark.py \
       --algorithms pandaml,pandaphysics \
       --output_dir pandaml_pandaphysics_comparison

Custom Metrics
~~~~~~~~~~~~~~

.. code-block:: python

   # Add custom evaluation metrics
   def custom_metric(predicted, experimental):
       # Your custom metric implementation
       return metric_value

   # Integrate into analysis pipeline
   data['custom_metric'] = data.apply(
       lambda row: custom_metric(row['predicted_affinity'], 
                                row['experimental_affinity']), 
       axis=1
   )

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**1. Memory errors**
   - Reduce batch size in benchmark script
   - Process subsets and combine results

**2. Missing dependencies**
   - Install all required packages: ``pip install -r requirements.txt``

**3. PDBbind path issues**
   - Verify directory structure matches expected format
   - Check file permissions

**4. Long runtime**
   - Use subset for testing: ``--max_complexes 100``
   - Increase parallelization: ``--n_jobs 16``

Getting Help
~~~~~~~~~~~~

For benchmark-specific issues:
- Check the benchmark logs in ``results/benchmark.log``
- Review the GitHub issues for similar problems
- Contact the development team with specific error messages

This comprehensive benchmark provides definitive validation of PandaDock's performance across the complete chemical space represented in PDBbind.