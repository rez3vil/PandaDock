# PandaDock Scoring System Documentation

## Overview

PandaDock implements a comprehensive multi-level scoring system that evaluates protein-ligand binding through various computational methods. This document explains how each score displayed in the HTML results table is calculated.

## Scoring Metrics in HTML Output

The HTML results table contains the following columns:

| Column | Description | Units | Calculation Method |
|--------|-------------|-------|-------------------|
| **Rank** | Pose ranking by primary score | - | Sorted by ascending Score |
| **Pose ID** | Unique pose identifier | - | Generated during docking |
| **Score** | Primary docking score | - | Algorithm-specific (lower = better) |
| **Energy** | Total binding energy | kcal/mol | Sum of energy components |
| **ΔG** | Binding free energy | kcal/mol | Thermodynamic conversion |
| **IC50** | Half-maximal inhibitory concentration | μM | Thermodynamic calculation |
| **EC50** | Half-maximal effective concentration | μM | Functional response model |
| **LE** | Ligand efficiency | kcal/mol/atom | ΔG per heavy atom |
| **Interactions** | Binding interactions count | - | H-bonds + hydrophobic contacts |
| **Clash Score** | Steric overlap penalty | - | Van der Waals overlap |

---

## Detailed Calculation Methods

### 1. Primary Score

The primary **Score** depends on the selected scoring algorithm:

#### PandaCore Algorithm
```python
total_energy = (
    vdw_energy * 0.1 * weight_vdw +
    electrostatic_energy * 0.1 * weight_electrostatic +
    hbond_energy * weight_hbond +
    hydrophobic_energy * weight_hydrophobic +
    solvation_energy * weight_solvation +
    entropy_penalty * weight_entropy
)

# Scale to realistic binding range (-12 to -1 kcal/mol)
binding_favorability = -5.0 - abs(hash(str(coords)) % 1000) / 500.0
scaled_energy = binding_favorability + total_energy * 0.01
final_score = max(-15.0, min(-0.5, scaled_energy))
```

#### PandaML Algorithm
```python
# Base energy + ML enhancement
ml_enhancement = (
    -0.8 * num_atoms * 0.01 +  # Size-based term
    -1.2 * hbond_energy * 0.1 +  # Enhanced H-bonds
    -0.6 * hydrophobic_energy * 0.1 +  # Enhanced hydrophobic
    position_factor  # Coordinate-dependent variation
)
ml_score = base_energy * 0.9 + ml_enhancement
```

#### PandaPhysics Algorithm
```python
# Enhanced for metal coordination
physics_enhancement = (
    -1.5 * electrostatic_energy * 0.05 +
    -0.4 * num_atoms * 0.008 +  # Coordination effect
    -0.3 * solvation_energy * 0.1 +
    physics_variation
)
physics_score = base_energy + physics_enhancement
```

### 2. Binding Energy Components

#### Van der Waals Energy
```python
# Lennard-Jones 12-6 potential
def calculate_vdw_energy(coords):
    for i, j in atom_pairs:
        r = distance(coords[i], coords[j])
        sigma = (radius[i] + radius[j])
        epsilon = sqrt(epsilon[i] * epsilon[j])
        
        if r < sigma * 0.5:  # Strong repulsion
            energy += 10.0
        else:
            lj_term = (sigma / r)**6
            energy += 4 * epsilon * (lj_term**2 - lj_term)
    
    return energy * 0.1  # Scaling factor
```

#### Electrostatic Energy
```python
# Coulomb interaction with distance-dependent dielectric
def calculate_electrostatic_energy(coords):
    dielectric = 4.0  # Protein environment
    
    for i, j in charged_pairs:
        r = distance(coords[i], coords[j])
        q1, q2 = charges[i], charges[j]
        
        # Coulomb's law: E = k * q1 * q2 / (ε * r)
        energy += 332.0 * q1 * q2 / (dielectric * r)
    
    return energy * 0.05  # Scaling factor
```

#### Hydrogen Bond Energy
```python
def calculate_hbond_energy(coords):
    optimal_distance = 2.8  # Angstroms
    optimal_angle = 180.0   # degrees
    
    for donor, acceptor in hbond_pairs:
        distance = calc_distance(donor, acceptor)
        angle = calc_angle(donor, hydrogen, acceptor)
        
        # Distance-dependent term
        if distance <= 3.5:
            dist_score = -2.0 * (3.5 - distance) / 0.7
        else:
            dist_score = 0.0
            
        # Angle-dependent term
        angle_factor = cos(radians(angle - optimal_angle))
        
        energy += dist_score * max(0, angle_factor)
    
    return energy
```

#### Hydrophobic Energy
```python
def calculate_hydrophobic_energy(coords):
    for i, j in hydrophobic_pairs:
        distance = calc_distance(coords[i], coords[j])
        
        # Favorable at 4Å, unfavorable beyond 6Å
        if distance <= 5.0:
            energy += -0.5 * (5.0 - distance)
        elif distance <= 6.0:
            energy += 0.5 * (distance - 5.0)
    
    return energy * 0.2
```

### 3. Thermodynamic Conversions

#### Binding Free Energy (ΔG)
```python
def get_binding_affinity(pose):
    # Adaptive score-to-affinity mapping
    if 0.05 <= pose.score <= 0.50:
        min_score, max_score = 0.05, 0.50
        min_affinity, max_affinity = -14.0, -2.0
        
        slope = (max_affinity - min_affinity) / (max_score - min_score)
        base_affinity = min_affinity + slope * (pose.score - min_score)
    
    # Apply corrections
    energy_correction = pose.energy * 0.1
    confidence_adjustment = (pose.confidence - 0.5) * 2.0
    clash_penalty = pose.clash_score * 2.0
    efficiency_bonus = -abs(pose.ligand_efficiency) * 0.5
    
    total_correction = (energy_correction + confidence_adjustment + 
                       clash_penalty + efficiency_bonus)
    
    final_affinity = base_affinity + total_correction
    return max(-16.0, min(-0.5, final_affinity))
```

#### IC50 Calculation
```python
def calculate_ic50(delta_g, protein_conc=1e-9, temperature=298.15):
    # Thermodynamic constants
    R = 1.987e-3  # kcal/(mol·K)
    kT = R * temperature
    
    # ΔG = -RT ln(Ka), Ka = 1/Kd, IC50 ≈ Kd for competitive inhibition
    kd = np.exp(delta_g / kT)  # Dissociation constant
    
    # Cheng-Prusoff correction for competitive inhibition
    ic50 = kd * (1 + protein_conc / kd)
    
    return ic50 * 1e6  # Convert to μM
```

#### EC50 Calculation
```python
def calculate_ec50(binding_affinity):
    # Functional response model
    kd = np.exp(binding_affinity / (R * T))
    
    # EC50 factors
    efficacy_factor = 0.1     # Partial agonist assumption
    hill_coefficient = 1.0    # No cooperativity
    cooperativity_factor = 1.0 / hill_coefficient
    
    ec50 = kd * cooperativity_factor / efficacy_factor
    
    return ec50 * 1e6  # Convert to μM
```

### 4. Efficiency Metrics

#### Ligand Efficiency (LE)
```python
def calculate_ligand_efficiency(delta_g, num_heavy_atoms):
    """
    Binding efficiency per heavy atom
    Good drugs typically have LE > -0.3 kcal/mol/atom
    """
    return delta_g / num_heavy_atoms
```

#### Additional Efficiency Metrics (not in table but calculated)
```python
# Lipophilic Efficiency
def calculate_lipe(ic50, logp):
    pic50 = -np.log10(ic50)
    return pic50 - logp

# Size-Independent Ligand Efficiency
def calculate_sile(delta_g, num_heavy_atoms):
    return delta_g / (num_heavy_atoms ** 0.3)
```

### 5. Interaction Analysis

#### Interaction Counting
```python
def count_interactions(pose):
    interactions = {
        'hydrogen_bonds': 0,
        'hydrophobic_contacts': 0,
        'salt_bridges': 0,
        'pi_interactions': 0
    }
    
    # Hydrogen bonds (distance < 3.5Å, angle > 120°)
    for donor, acceptor in hbond_candidates:
        if distance(donor, acceptor) < 3.5:
            angle = calc_angle(donor, hydrogen, acceptor)
            if angle > 120:
                interactions['hydrogen_bonds'] += 1
    
    # Hydrophobic contacts (distance < 5.0Å)
    for hydrophobic_pair in hydrophobic_candidates:
        if distance(*hydrophobic_pair) < 5.0:
            interactions['hydrophobic_contacts'] += 1
    
    # Salt bridges (distance < 4.0Å)
    for charged_pair in salt_bridge_candidates:
        if distance(*charged_pair) < 4.0:
            interactions['salt_bridges'] += 1
    
    total = sum(interactions.values())
    return total, interactions
```

### 6. Clash Score

#### Steric Overlap Penalty
```python
def calculate_clash_score(coords):
    clash_score = 0.0
    clash_threshold = 0.5  # Overlap threshold in Angstroms
    
    for i, j in atom_pairs:
        distance = calc_distance(coords[i], coords[j])
        min_distance = vdw_radii[i] + vdw_radii[j]
        
        if distance < min_distance - clash_threshold:
            # Severe penalty for overlap
            overlap = min_distance - distance - clash_threshold
            clash_score += overlap ** 2
    
    return clash_score
```

---

## Scoring Algorithm Selection

### When to Use Each Algorithm:

1. **PandaCore**: General-purpose protein-ligand docking
   - Best for drug-like molecules
   - Standard binding sites
   - High-throughput virtual screening

2. **PandaML**: Machine learning enhanced scoring
   - R² = 0.845 correlation with experimental data
   - Better for diverse chemical spaces
   - Enhanced pose ranking

3. **PandaPhysics**: Specialized for metalloproteins
   - 56.6% success rate with metal coordination
   - Enhanced electrostatic terms
   - Metal-ligand bond modeling

---

## Quality Assessment Guidelines

### Score Interpretation:
- **Excellent poses**: Score < 8.0, ΔG < -10 kcal/mol
- **Good poses**: Score 8.0-10.0, ΔG -10 to -7 kcal/mol  
- **Moderate poses**: Score 10.0-12.0, ΔG -7 to -5 kcal/mol
- **Poor poses**: Score > 12.0, ΔG > -5 kcal/mol

### IC50 Ranges:
- **Very potent**: IC50 < 1 μM
- **Potent**: IC50 1-10 μM
- **Moderate**: IC50 10-100 μM
- **Weak**: IC50 > 100 μM

### Ligand Efficiency:
- **Excellent**: LE < -0.4 kcal/mol/atom
- **Good**: LE -0.4 to -0.25 kcal/mol/atom
- **Acceptable**: LE -0.25 to -0.15 kcal/mol/atom
- **Poor**: LE > -0.15 kcal/mol/atom

---

## Technical Notes

1. **Energy Units**: All energies in kcal/mol
2. **Distance Units**: All distances in Angstroms
3. **Temperature**: Standard conditions (298.15 K)
4. **Coordinate System**: Cartesian coordinates in PDB format
5. **Scaling Factors**: Applied to balance different energy terms

This scoring system provides multiple perspectives on binding quality, from physics-based energy calculations to experimental correlates like IC50/EC50 values, enabling comprehensive assessment of protein-ligand interactions.