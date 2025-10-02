# Lab Results Template - Localized Surface Plasmons

## Experiment Setup Summary

**Fixed Parameters:**
- Outer radius: a = 10 nm
- Core refractive index: nc = 1 (air)
- Background refractive index: nb = 1 (air)

**Date:** ___________

---

## Part 1: Initial Observation (Basic Example)

### Run: `example_mie_spectrum_3d_shell_eV.m`

**Parameters used:**
- Metal: ___________
- Volume fraction: f = ___________
- Outer radius: a = ___________ nm

### Observations:

#### Peak 1 (Lower Energy/Longer Wavelength):
- Energy: ___________ eV
- Wavelength: ___________ nm
- Absorption cross-section: ___________
- Scattering cross-section: ___________
- Ratio C_abs/C_sca: ___________

**Mode identification:**
- [ ] Dipole (l=1)
- [ ] Quadrupole (l=2)
- [ ] Higher order
- **Evidence:** (from angular scattering, field pattern)

**Physical character:**
- [ ] Bonding mode (symmetric)
- [ ] Antibonding mode (antisymmetric)
- **Field pattern description:**

#### Peak 2 (Higher Energy/Shorter Wavelength):
- Energy: ___________ eV
- Wavelength: ___________ nm
- Absorption cross-section: ___________
- Scattering cross-section: ___________
- Ratio C_abs/C_sca: ___________

**Mode identification:**
- [ ] Dipole (l=1)
- [ ] Quadrupole (l=2)
- [ ] Higher order
- **Evidence:**

**Physical character:**
- [ ] Bonding mode (symmetric)
- [ ] Antibonding mode (antisymmetric)
- **Field pattern description:**

#### Mode Splitting:
- Energy splitting: ΔE = ___________ eV
- Wavelength splitting: Δλ = ___________ nm
- Ratio ΔE/E_avg = ___________

### Near-Field Observations:

**Peak 1 fields:**
- Field hotspots located: [ ] Inner surface [ ] Outer surface [ ] Both
- Maximum enhancement |E|²/|E₀|²: ___________
- Power flow pattern:

**Peak 2 fields:**
- Field hotspots located: [ ] Inner surface [ ] Outer surface [ ] Both
- Maximum enhancement |E|²/|E₀|²: ___________
- Power flow pattern:

### Angular Scattering:

**Peak 1:**
- Forward scattering (0°): [ ] Strong [ ] Moderate [ ] Weak
- Backward scattering (180°): [ ] Strong [ ] Moderate [ ] Weak
- Pattern shape:

**Peak 2:**
- Forward scattering (0°): [ ] Strong [ ] Moderate [ ] Weak
- Backward scattering (180°): [ ] Strong [ ] Moderate [ ] Weak
- Pattern shape:

---

## Part 2: Shell Thickness Study

### Run: `example_shell_thickness_study.m`

**Fixed:** Metal = Au, a = 10 nm, nc = 1, nb = 1

### Results Table:

| f | Shell thickness (nm) | Peak 1 Energy (eV) | Peak 1 λ (nm) | Peak 2 Energy (eV) | Peak 2 λ (nm) | Splitting ΔE (eV) |
|---|---------------------|-------------------|---------------|-------------------|---------------|------------------|
| 0.1 | | | | | | |
| 0.3 | | | | | | |
| 0.5 | | | | | | |
| 0.7 | | | | | | |
| 0.9 | | | | | | |

**Note:** Shell thickness = a × [1 - (1-f)^(1/3)]

### Observations:

#### Trend 1: Peak Position vs Shell Thickness
- As f increases (thicker shell):
  - Peak 1 (bonding): [ ] Red-shifts [ ] Blue-shifts [ ] Stays same
  - Peak 2 (antibonding): [ ] Red-shifts [ ] Blue-shifts [ ] Stays same
  
#### Trend 2: Mode Splitting vs Shell Thickness
- As f increases: ΔE [ ] Increases [ ] Decreases [ ] Stays same
- Does this match ΔE ∝ exp(-d/R)? [ ] Yes [ ] No

#### Trend 3: Relative Peak Strengths
- For thin shells (f ≈ 0.1):
  - Which peak is stronger? ___________ 
  - Absorption vs scattering ratio: ___________

- For thick shells (f ≈ 0.9):
  - Which peak is stronger? ___________
  - Absorption vs scattering ratio: ___________

### Physical Interpretation:

**Plasmon hybridization model:**
- Thin shell behavior:

- Thick shell behavior:

- Convergence to solid sphere:

---

## Part 3: Metal Comparison Study

### Run: `example_metal_comparison.m`

**Fixed:** f = 0.5, a = 10 nm, nc = 1, nb = 1

### Results Table:

| Metal | Peak 1 Energy (eV) | Peak 1 λ (nm) | Peak 1 Width (eV) | Peak 2 Energy (eV) | Peak 2 λ (nm) | Peak 2 Width (eV) | Notes |
|-------|-------------------|---------------|-------------------|-------------------|---------------|-------------------|-------|
| **Au** | | | | | | | |
| **Ag** | | | | | | | |
| **Al** | | | | | | | |
| **Cu** | | | | | | | |

### Observations:

#### Peak Position Ranking (lowest to highest energy):
1. ___________ (most red-shifted)
2. ___________
3. ___________
4. ___________ (most blue-shifted)

**Physical explanation:**
(Relate to plasma frequency ωₚ)

#### Peak Width Ranking (narrowest to broadest):
1. ___________ (sharpest, best Q-factor)
2. ___________
3. ___________
4. ___________ (broadest, worst Q-factor)

**Physical explanation:**
(Relate to damping γ)

#### Peak Strength Observations:
- Which metal has strongest absorption peak? ___________
- Which metal has strongest scattering peak? ___________
- Which metal has best overall plasmonic response? ___________

### Metal Properties Correlation:

| Metal | Plasma freq ωₚ | Damping γ | Observable Effect |
|-------|---------------|-----------|-------------------|
| Au | ~9 eV | ~0.07 eV | |
| Ag | ~9 eV | ~0.02 eV | |
| Al | ~15 eV | ~0.5 eV | |
| Cu | ~9 eV | ~0.1 eV | |

**Conclusions:**
- Effect of ωₚ on peak position:

- Effect of γ on peak width:

- Trade-offs for applications:

---

## Part 4: Optional - Core Index Study

**If performed:**

### Parameters:
- Metal: ___________
- f: ___________
- Core index values tested: ___________

### Results:

| nc | Peak 1 (eV) | Peak 2 (eV) | Splitting ΔE | Notes |
|----|-------------|-------------|--------------|-------|
| 1.0 | | | | Air core |
| 1.5 | | | | Silica core |
| 2.0 | | | | High-index core |
| 2.5 | | | | Very high-index |

### Observations:
- Effect on bonding mode:
- Effect on antibonding mode:
- Physical explanation:

---

## Summary of Key Findings

### Finding 1: Mode Identification
- Peak 1 is: ___________
- Peak 2 is: ___________
- Evidence: ___________

### Finding 2: Physical Origin of Two Peaks
- Explanation using plasmon hybridization:

- Charge oscillation picture:

### Finding 3: Shell Thickness Effect
- Main trend:
- Physical explanation:
- Quantitative relationship:

### Finding 4: Metal Choice Effect
- Best metal for UV range: ___________
- Best metal for visible range: ___________
- Best metal for narrow linewidth: ___________
- Trade-offs: ___________

### Finding 5: Absorption vs Scattering
- Which mode absorbs more: ___________
- Which mode scatters more: ___________
- Size dependence explanation:

---

## Connections to Physical Models

### 1. Fröhlich Condition
- For solid sphere: ε_metal + 2ε_b = 0
- How does this explain your observations?

### 2. Plasmon Hybridization
- Cavity mode + Sphere mode → Bonding + Antibonding
- How does your data support this model?

### 3. Drude Model
- ε(ω) = 1 - ωₚ²/(ω² + iγω)
- How does this explain metal differences?

### 4. Quasi-static vs Retarded Regime
- For 10 nm particles, are we in quasi-static limit?
- Evidence from absorption/scattering ratio:

---

## Questions Answered

### Lab Question 1: "What's the difference between the two peaks?"
**Answer:**

### Lab Question 2: "Which modes are they (dipole or higher)?"
**Answer:**

### Lab Question 3: "How do parameters affect peak position, strength and character?"

**Geometry (shell thickness):**

**Metal choice:**

**Core index (if studied):**

### Lab Question 4: "How can results be related to physical models?"

**Charge distribution:**

**Plasmon hybridization:**

**Electromagnetic boundary conditions:**

---

## Figures to Include in Report

- [ ] Figure 1: Spectrum for basic example (2 peaks, absorption + scattering)
- [ ] Figure 2: Near-field plot for Peak 1 (bonding mode)
- [ ] Figure 3: Near-field plot for Peak 2 (antibonding mode)
- [ ] Figure 4: Angular scattering patterns
- [ ] Figure 5: Peak position vs shell thickness
- [ ] Figure 6: Mode splitting vs shell thickness
- [ ] Figure 7: Spectra for different metals (overlay)
- [ ] Figure 8: Peak position vs metal type
- [ ] Figure 9: Peak width vs metal type (Q-factor)
- [ ] Figure 10: (Optional) Core index effect

---

## Additional Notes & Observations

_(Use this space for any unexpected findings, interesting patterns, or ideas for discussion)_

---

## Checklist Before Submitting

- [ ] All figures properly labeled with axes, units, legends
- [ ] Physical interpretation provided for each trend
- [ ] Compared at least 2 metals (as required)
- [ ] Explored shell thickness effect (as required)
- [ ] Identified mode types (dipole/quadrupole)
- [ ] Explained difference between two peaks
- [ ] Related findings to physical models
- [ ] Checked calculations and units
- [ ] Proofread for clarity

---

**Lab completed by:** ___________________________
**Date:** _______________
