# MS-CASPT2 vs CASPT2 Comparison Report

## Executive Summary

This report presents a detailed comparison between Multi-State CASPT2 (MS-CASPT2) and standard CASPT2 calculations for two archetypal peptide geometries: Beta-Strand (φ=210°, ψ=135°) and Alpha-Helix (φ=240°, ψ=30°). The analysis focuses on transition energies, dipole moments, and absolute energies for four key electronic transitions.

## Key Findings

### 1. Overall Energy Differences

**Beta-Strand Geometry (φ=210°, ψ=135°):**
- Average MS-CASPT2 vs CASPT2 difference: **+0.144 ± 0.262 eV**
- Largest difference: **+0.572 eV** (pi_nb_to_pi_*_L transition)
- Smallest difference: **+0.015 eV** (n_to_pi_*_L transition)

**Alpha-Helix Geometry (φ=240°, ψ=30°):**
- Average MS-CASPT2 vs CASPT2 difference: **+0.065 ± 0.084 eV**
- Largest difference: **+0.169 eV** (n_to_pi_*_L transition)
- Smallest difference: **+0.032 eV** (pi_nb_to_pi_*_L transition)

### 2. Geometry-Dependent Effects

- **Beta-Strand geometry shows larger MS-CASPT2 vs CASPT2 differences** compared to Alpha-Helix
- **Geometry effect on MS-CASPT2 vs CASPT2 difference: -0.079 eV**
- **Most geometry-sensitive transition:** pi_nb_to_pi_*_L (0.700 eV difference between geometries)
- **Least geometry-sensitive transition:** n_to_pi_*_L (0.018 eV difference between geometries)

## Detailed Analysis by Geometry

### Beta-Strand Geometry (φ=210°, ψ=135°)

| Transition | CASPT2 (eV) | MS-CASPT2 (eV) | Difference (eV) | % Difference |
|------------|-------------|----------------|-----------------|--------------|
| pi_nb_to_pi_*_L | 6.989 | 7.561 | +0.572 | +8.18% |
| pi_nb_to_pi_*_R | 6.901 | 6.775 | -0.126 | -1.83% |
| n_to_pi_*_L | 5.807 | 5.822 | +0.015 | +0.26% |
| n_to_pi_*_R | 5.760 | 5.874 | +0.114 | +1.98% |

**Key Observations:**
- **pi_nb_to_pi_*_L** shows the largest difference (+0.572 eV), indicating significant multi-state effects
- **pi_nb_to_pi_*_R** shows a negative difference (-0.126 eV), suggesting state-specific effects
- **n_to_pi_* transitions** show relatively small differences, indicating less multi-state coupling

**Transition Dipole Moments:**
- Strongest transition: pi_nb_to_pi_*_L (4.047 D)
- Weakest transition: n_to_pi_*_L (0.178 D)
- Average dipole moment: 1.748 D

### Alpha-Helix Geometry (φ=240°, ψ=30°)

| Transition | CASPT2 (eV) | MS-CASPT2 (eV) | Difference (eV) | % Difference |
|------------|-------------|----------------|-----------------|--------------|
| pi_nb_to_pi_*_L | 6.829 | 6.861 | +0.032 | +0.47% |
| pi_nb_to_pi_*_R | 6.858 | 6.971 | +0.113 | +1.65% |
| n_to_pi_*_L | 5.635 | 5.804 | +0.169 | +3.00% |
| n_to_pi_*_R | 5.746 | 5.692 | -0.054 | -0.94% |

**Key Observations:**
- **n_to_pi_*_L** shows the largest difference (+0.169 eV), indicating significant multi-state effects
- **pi_nb_to_pi_*_L** shows the smallest difference (+0.032 eV), suggesting minimal multi-state coupling
- Overall smaller differences compared to Beta-Strand geometry

**Transition Dipole Moments:**
- Strongest transition: pi_nb_to_pi_*_L (4.119 D)
- Weakest transition: n_to_pi_*_R (0.147 D)
- Average dipole moment: 2.102 D

## Electronic Configuration Analysis

### Beta-Strand Geometry
- **pi_nb_to_pi_*_L**: Configuration `22222u0d`, Coefficient 0.811, Weight 0.658
- **pi_nb_to_pi_*_R**: Configuration `2222u2d0`, Coefficient -0.763, Weight 0.582
- **n_to_pi_*_L**: Configuration `222u220d`, Coefficient 0.957, Weight 0.915
- **n_to_pi_*_R**: Configuration `22u222d0`, Coefficient 0.921, Weight 0.849

### Alpha-Helix Geometry
- **pi_nb_to_pi_*_L**: Configuration `22u2220d`, Coefficient 0.766, Weight 0.587
- **pi_nb_to_pi_*_R**: Configuration `22222ud0`, Coefficient 0.907, Weight 0.822
- **n_to_pi_*_L**: Configuration `222u220d`, Coefficient 0.859, Weight 0.738
- **n_to_pi_*_R**: Configuration `2222u2d0`, Coefficient -0.954, Weight 0.910

## State-by-State Energy Analysis

### Beta-Strand Geometry
| State | CASPT2 (Hartree) | MS-CASPT2 (Hartree) | Difference (Hartree) |
|-------|------------------|---------------------|---------------------|
| 2 | -455.168000 | -455.170000 | -0.002000 |
| 3 | -455.169000 | -455.168000 | +0.001000 |
| 4 | -455.127000 | -455.135000 | -0.008000 |
| 6 | -455.124000 | -455.106000 | +0.018000 |

### Alpha-Helix Geometry
| State | CASPT2 (Hartree) | MS-CASPT2 (Hartree) | Difference (Hartree) |
|-------|------------------|---------------------|---------------------|
| 2 | -455.170000 | -455.174000 | -0.004000 |
| 3 | -455.174000 | -455.170000 | +0.004000 |
| 4 | -455.130000 | -455.131000 | -0.001000 |
| 5 | -455.129000 | -455.127000 | +0.002000 |

## Conclusions and Implications

### 1. Multi-State Effects
- **MS-CASPT2 shows systematic differences** from CASPT2, indicating the importance of multi-state coupling
- **Beta-Strand geometry exhibits larger multi-state effects** than Alpha-Helix geometry
- **pi_nb_to_pi_* transitions** show the most significant multi-state effects in Beta-Strand geometry

### 2. Geometry Dependence
- **Conformational changes significantly affect** the magnitude of MS-CASPT2 vs CASPT2 differences
- **Beta-Strand geometry** shows more pronounced multi-state coupling effects
- **Alpha-Helix geometry** shows more uniform and smaller differences

### 3. Transition-Specific Effects
- **pi_nb_to_pi_*_L transition** is most sensitive to geometry changes (0.700 eV difference)
- **n_to_pi_*_L transition** is least sensitive to geometry changes (0.018 eV difference)
- **Different transitions respond differently** to multi-state coupling effects

### 4. Practical Implications
- **MS-CASPT2 is essential** for accurate description of excited states in peptide systems
- **Geometry-dependent effects** should be considered when choosing between CASPT2 and MS-CASPT2
- **Beta-Strand conformations** require more careful treatment of multi-state effects
- **Transition-specific analysis** is crucial for understanding electronic structure

## Recommendations

1. **Use MS-CASPT2** for accurate excited state calculations in peptide systems
2. **Consider geometry effects** when interpreting MS-CASPT2 vs CASPT2 differences
3. **Pay special attention** to pi_nb_to_pi_* transitions in Beta-Strand geometries
4. **Perform transition-specific analysis** rather than relying on average differences
5. **Include multiple conformations** in computational studies to capture geometry-dependent effects

## Files Generated

1. **beta_strand_comparison.png** - Comprehensive comparison plot for Beta-Strand geometry
2. **alpha_helix_comparison.png** - Comprehensive comparison plot for Alpha-Helix geometry
3. **plot_comparison_for_single_geometries.py** - Script for generating comparison plots
4. **detailed_comparison_analysis.py** - Script for detailed numerical analysis

---

*Report generated from analysis of merged_results.csv data*
*Geometries analyzed: Beta-Strand (φ=210°, ψ=135°) and Alpha-Helix (φ=240°, ψ=30°)* 