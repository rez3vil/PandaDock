# PandaDock Comprehensive Benchmark Report

**Date:** 2025-07-10 14:21:52

**Total Complexes Evaluated:** 5316
**Total Docking Runs:** 15948
**Engines Evaluated:** pandacore, pandaml, pandaphysics

## Dataset Statistics

- **Experimental Affinity Range:** 4.00 - 10.50 pKd/pKi
- **Mean Experimental Affinity:** 7.22 ± 1.87
- **Ligand Size Range:** 15 - 79 heavy atoms
- **Mean Ligand Size:** 47.4 ± 18.9 heavy atoms

## Engine Performance Summary

### PANDACORE Engine

- **Number of complexes:** 5316
- **Affinity Prediction:**
  - Pearson correlation: 0.842
  - R²: 0.709
  - RMSE: 1.201
  - MAE: 0.954
- **Pose Prediction:**
  - Mean RMSD: 3.314 Å
  - Median RMSD: 2.200 Å
  - Success rate (RMSD < 2Å): 0.471
  - Success rate (RMSD < 3Å): 0.610
- **Computational Efficiency:**
  - Mean docking time: 33.7 seconds
  - Median docking time: 33.1 seconds
  - Time per heavy atom: 0.76 s/atom

### PANDAML Engine

- **Number of complexes:** 5316
- **Affinity Prediction:**
  - Pearson correlation: 0.919
  - R²: 0.845
  - RMSE: 0.803
  - MAE: 0.641
- **Pose Prediction:**
  - Mean RMSD: 3.111 Å
  - Median RMSD: 2.051 Å
  - Success rate (RMSD < 2Å): 0.490
  - Success rate (RMSD < 3Å): 0.634
- **Computational Efficiency:**
  - Mean docking time: 26.7 seconds
  - Median docking time: 25.2 seconds
  - Time per heavy atom: 0.62 s/atom

### PANDAPHYSICS Engine

- **Number of complexes:** 5316
- **Affinity Prediction:**
  - Pearson correlation: 0.877
  - R²: 0.769
  - RMSE: 1.018
  - MAE: 0.807
- **Pose Prediction:**
  - Mean RMSD: 3.203 Å
  - Median RMSD: 2.103 Å
  - Success rate (RMSD < 2Å): 0.483
  - Success rate (RMSD < 3Å): 0.624
- **Computational Efficiency:**
  - Mean docking time: 45.6 seconds
  - Median docking time: 45.5 seconds
  - Time per heavy atom: 1.01 s/atom

## Statistical Comparisons

### RMSD Comparisons (Wilcoxon Rank-Sum Test)

| Engine 1 | Engine 2 | p-value | Significant |
|----------|----------|---------|-------------|
| PANDACORE | PANDAML | 0.0031 | Yes |
| PANDACORE | PANDAPHYSICS | 0.0709 | No |
| PANDAML | PANDAPHYSICS | 0.2518 | No |

### Performance by Ligand Size

#### Very Small Ligands

**Size range:** 15-27 heavy atoms
**Number of complexes:** 1059

- **PANDACORE:** RMSD = 2.197 Å, Success = 0.624
- **PANDAML:** RMSD = 1.887 Å, Success = 0.665
- **PANDAPHYSICS:** RMSD = 2.013 Å, Success = 0.627

#### Small Ligands

**Size range:** 28-40 heavy atoms
**Number of complexes:** 1027

- **PANDACORE:** RMSD = 2.840 Å, Success = 0.497
- **PANDAML:** RMSD = 2.515 Å, Success = 0.553
- **PANDAPHYSICS:** RMSD = 2.812 Å, Success = 0.526

#### Medium Ligands

**Size range:** 41-53 heavy atoms
**Number of complexes:** 1073

- **PANDACORE:** RMSD = 3.169 Å, Success = 0.477
- **PANDAML:** RMSD = 3.121 Å, Success = 0.447
- **PANDAPHYSICS:** RMSD = 3.190 Å, Success = 0.463

#### Large Ligands

**Size range:** 54-66 heavy atoms
**Number of complexes:** 1039

- **PANDACORE:** RMSD = 3.817 Å, Success = 0.411
- **PANDAML:** RMSD = 3.679 Å, Success = 0.442
- **PANDAPHYSICS:** RMSD = 3.753 Å, Success = 0.428

#### Very Large Ligands

**Size range:** 67-79 heavy atoms
**Number of complexes:** 1118

- **PANDACORE:** RMSD = 4.480 Å, Success = 0.352
- **PANDAML:** RMSD = 4.278 Å, Success = 0.354
- **PANDAPHYSICS:** RMSD = 4.192 Å, Success = 0.375

## Generated Figures

- **Master Publication Figure:** ![Master](master_publication_figure.png)
- **Correlation Analysis:** ![Correlation](correlation_analysis.png)
- **RMSD Analysis:** ![RMSD](rmsd_analysis.png)
- **Engine Performance:** ![Performance](engine_performance.png)
- **Ligand Complexity Analysis:** ![Complexity](ligand_complexity_analysis.png)
- **Performance vs Properties:** ![Properties](performance_vs_properties.png)
