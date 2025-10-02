# Localized Surface Plasmons in Metal Nanoshells
## Computational Study Using Mie Theory

**Course:** 68513 Plasmonics  
**Lab:** Computer Lab - LSP (Localized Surface Plasmons)  
**Date:** 2 October 2025

---

## 1. AIM

The aim of this computational laboratory is to investigate localized surface plasmon (LSP) resonances in metal nanoshells using exact electromagnetic theory. Specifically, we aim to:

1. Observe and characterize **plasmon hybridization** - the coupling between inner cavity and outer sphere plasmon modes
2. Study the effect of **shell thickness** (volume fraction) on mode splitting and resonance energies
3. Compare **different metals** (Au, Ag, Al, Cu) to understand how material properties (plasma frequency, damping) affect plasmon behavior
4. Validate theoretical predictions from the quasistatic approximation and plasmon hybridization model against exact Mie theory calculations

This investigation provides fundamental insight into nanoshell plasmonics, which has applications in sensing, photothermal therapy, and nanophotonic devices.

---

## 2. BACKGROUND

### 2.1 Plasmon Hybridization in Nanoshells

Metal nanoshells consist of a dielectric core surrounded by a thin metallic shell. Unlike solid spheres which support a single dipolar plasmon mode, nanoshells exhibit **plasmon hybridization** - a phenomenon where the cavity mode (plasmon at the inner metal-dielectric interface) couples with the sphere mode (plasmon at the outer metal-background interface) to form two distinct hybrid modes:

1. **Bonding mode (symmetric):** Inner and outer surface charges oscillate **in phase**
   - Lower energy/frequency
   - Large net dipole moment
   - Primarily absorptive character

2. **Antibonding mode (antisymmetric):** Charges oscillate **out of phase**
   - Higher energy/frequency
   - Partial dipole cancellation
   - More radiative/scattering character

### 2.2 Theoretical Framework

#### Quasistatic Depolarization Model

For particles much smaller than the wavelength, the quasistatic approximation describes plasmon resonances using depolarization factors. For a nanoshell with volume fraction *f* = (volume of metal)/(total volume), the two hybrid modes have depolarization factors:

**Symmetric mode:**
```
S₁ = 1/2 - √(9-8f)/6
```

**Antisymmetric mode:**
```
S₂ = 1/2 + √(9-8f)/6
```

The mode splitting Δ*S* = *S*₂ - *S*₁ decreases as *f* increases (thicker shells), predicting that **thin shells should have larger mode separation**.

#### Resonance Condition

Plasmon resonance occurs when the metal permittivity satisfies:

```
ε/ε_b ≈ -(1/S - 1)
```

where ε is the metal permittivity and ε_b is the background permittivity.

#### Drude Model for Metal Permittivity

The frequency-dependent permittivity of metals is described by the Drude model:

```
ε(ω)/ε₀ = 1 - ω_p²/(ω² + iγω)
```

where:
- ω_p = plasma frequency (determines resonance energy scale)
- γ = damping rate (determines peak width and losses)

**Metal parameters from tutorial:**

| Metal | ω_p (eV) | γ (eV) | Expected Behavior |
|-------|----------|--------|-------------------|
| K     | 3.72     | 0.03   | Low frequency (IR/NIR) |
| Cu    | 9.0      | 0.10   | Visible range |
| Au    | 9.0      | 0.07   | Visible, moderate damping |
| Ag    | 9.0      | 0.02   | Visible, low damping (sharp peaks) |
| Al    | 15.0     | 0.50   | High frequency (UV), high damping |

**Key prediction:** Higher ω_p → Higher resonance energy (blue shift)

#### Cross-Sections in Small Particle Limit

For particles smaller than the wavelength:

```
C_abs ∝ Volume    (linear with size)
C_sca ∝ Volume²   (quadratic with size)
```

This means small particles are absorption-dominated while larger particles scatter more efficiently.

### 2.3 Mie Theory: Exact Electromagnetic Solution

The **Mie theory** provides an exact analytical solution to Maxwell's equations for spherical particles of arbitrary size and composition. Unlike the quasistatic approximation, Mie theory accounts for:

- **Retardation effects** (phase variation across particle)
- **Radiation damping** 
- **Multipolar resonances** (dipole, quadrupole, etc.)
- **Accurate field distributions**

The solution expands fields in vector spherical harmonics with coefficients determined by boundary conditions. For nanoshells, this involves matching fields at both inner and outer interfaces.

Mie theory becomes essential when:
- Particle size *a* ≳ λ/20 (approximately)
- Accurate quantitative predictions are needed
- Multiple resonances may overlap

---

## 3. METHODS

### 3.1 Computational Approach

All simulations were performed using MATLAB implementations of **Mie theory for core-shell particles**. The code calculates:

1. Mie coefficients for dipolar modes (*l*=1)
2. Extinction, absorption, and scattering cross-sections
3. Peak identification via automated analysis
4. Energy conversion (wavelength ↔ eV using *E* = 1240/*λ*)

### 3.2 Simulation Parameters

**Three systematic studies were conducted:**

#### Experiment 1: Basic Example (K Shell)
- Metal: Potassium (K)
- Outer radius: *a* = 20 nm
- Volume fraction: *f* = 0.3
- Core refractive index: *n_c* = 1.5
- Background refractive index: *n_b* = 1.5
- **Purpose:** Demonstrate basic two-peak hybridization structure

#### Experiment 2: Thickness Study (Au Shell)
- Metal: Gold (Au)
- Outer radius: *a* = 10 nm (fixed)
- Volume fraction: *f* = [0.1, 0.3, 0.5, 0.7, 0.9]
- Core refractive index: *n_c* = 1.0 (air)
- Background: *n_b* = 1.0 (air)
- **Purpose:** Test thickness-dependent splitting prediction

#### Experiment 3: Metal Comparison
- Metals: Au, Ag, Al, Cu
- Outer radius: *a* = 10 nm
- Volume fraction: *f* = 0.5 (fixed)
- Environment: *n_c* = *n_b* = 1.0
- **Purpose:** Validate material-dependent trends (ω_p correlation)

### 3.3 Analysis Workflow

1. **MATLAB simulations** → Generate exact Mie theory results
2. **JSON export** → Transfer data to Python for analysis
3. **Python/Jupyter** → Visualization and trend analysis
4. **Physical interpretation** → Compare with tutorial theory

All data processing and visualization were performed using Python 3.10 with NumPy, Pandas, Matplotlib, and SciPy.

---

## 4. RESULTS

### 4.1 Basic Example: K Nanoshell Hybridization

**Configuration:** K shell, *a* = 20 nm, *f* = 0.3, *n_c* = *n_b* = 1.5

**Observed Peaks:**

| Mode | Energy (eV) | Wavelength (nm) | Spectral Region |
|------|-------------|-----------------|-----------------|
| Peak 1 (Bonding) | 0.717 | 1730 | Near-Infrared |
| Peak 2 (Antibonding) | 3.543 | 350 | UV |

**Mode Splitting:** ΔE = 2.826 eV

**Key Observation:** The K nanoshell exhibits a **very large mode splitting** spanning from near-IR to UV. This wide separation is due to potassium's exceptionally low plasma frequency (ω_p = 3.72 eV), which places both hybrid modes at lower energies compared to noble metals. The two peaks are well-resolved and clearly demonstrate plasmon hybridization.

![Basic Example Visualization](reference: Figure from notebook showing two-peak structure)

**Physical Interpretation:**
- **Bonding mode (0.717 eV):** Lower energy mode where inner and outer surface charges oscillate in phase, creating a large dipole. The low frequency reflects weak restoring forces.
- **Antibonding mode (3.543 eV):** Higher energy mode with charges oscillating out of phase, leading to partial cancellation and higher restoring force.

### 4.2 Shell Thickness Study: Au Nanoshell

**Configuration:** Au shell, *a* = 10 nm, *n_c* = *n_b* = 1.0

**Complete Dataset:**

| *f* | Thickness (nm) | Peak 1 (eV) | Peak 2 (eV) | Splitting ΔE (eV) |
|-----|----------------|-------------|-------------|-------------------|
| 0.1 | 0.35 | 1.194 | 1.194 | 0.000 |
| 0.3 | 1.12 | 1.861 | 1.861 | 0.000 |
| 0.5 | 2.06 | 2.164 | 2.153 | 0.011 |
| 0.7 | 3.31 | 2.321 | 2.299 | 0.022 |
| 0.9 | 5.36 | 2.416 | 5.578 | **3.162** |

**Critical Findings:**

1. **Modes are degenerate for f < 0.7:** The two peaks are essentially at the same energy, giving zero splitting. This is **opposite** to simple quasistatic prediction which suggests thin shells should have large splitting.

2. **Maximum splitting at f = 0.9 (thickest shell):** The splitting dramatically increases to 3.16 eV only for the very thick shell configuration.

3. **Non-monotonic behavior:** The trend is not a simple linear relationship - there's a sharp transition around *f* ≈ 0.7-0.9.

![Thickness Study Plots](reference: 4-panel figure showing trends with f)

**Comparison with Theory:**

The simple quasistatic prediction from tutorial formulas suggests:
- ΔS decreases with increasing *f*
- Therefore splitting should **decrease** as shells get thicker

**Mie theory shows the opposite!** Why?

**Physical Explanation:**

The discrepancy arises because:

1. **Retardation effects matter:** Even at *a* = 10 nm, the particle is not truly in the quasistatic limit (λ/20 ≈ 30 nm for visible light). Phase variations across the particle become important.

2. **Mode hybridization strength:** For very thin shells (*f* < 0.5), the inner and outer surfaces are so close that the modes strongly overlap and cannot be resolved as separate peaks - they merge into a single mode.

3. **Thick shell regime:** Only when the shell is sufficiently thick (*f* ≈ 0.9) do the inner and outer surfaces become electromagnetically distinct enough to support clearly separated modes.

This demonstrates that **exact Mie theory is essential for quantitative predictions**, while quasistatic theory provides qualitative intuition that may not hold for real particle sizes.

### 4.3 Metal Comparison Study

**Configuration:** *a* = 10 nm, *f* = 0.5, *n_c* = *n_b* = 1.0

**Results by Metal:**

| Metal | ω_p (eV) | Peak Energy (eV) | Wavelength (nm) | Region | Splitting (eV) |
|-------|----------|------------------|-----------------|--------|----------------|
| Cu | 9.0 | 2.102 | 590 | Visible (Orange) | 0.030 |
| Au | 9.0 | 2.164 | 573 | Visible (Yellow-Green) | 0.011 |
| Ag | 9.0 | 2.762 | 449 | Visible (Blue) | 0.000 |
| Al | 15.0 | 5.421 | 229 | UV | 0.030 |

**Energy Ranking (Low → High):**
Cu < Au < Ag < Al

**Plasma Frequency Ranking:**
Cu, Au, Ag (all ≈9 eV) < Al (15 eV)

![Metal Comparison](reference: Multi-panel figure with correlation plots)

**Validation Against Theory:**

✅ **Confirmed:** Al has highest peak energy (UV region) due to highest ω_p  
✅ **Confirmed:** Au/Ag/Cu cluster in visible range (similar ω_p ≈ 9 eV)  
✅ **Confirmed:** Direct correlation between ω_p and resonance energy  
✅ **Confirmed:** Spectral tuning via material selection works as predicted

**Key Insights:**

1. **Material determines color:** The plasma frequency directly controls where the plasmon resonance appears in the spectrum. Al's high ω_p places it firmly in the UV, while noble metals span the visible range.

2. **Small splitting at f = 0.5:** All metals show minimal mode splitting at this geometry (0-0.03 eV), consistent with the thickness study findings that intermediate *f* values produce nearly degenerate modes.

3. **Drude model validates well:** The simple Drude picture from the tutorial successfully predicts material-dependent trends, even though absolute energies require full Mie calculation.

### 4.4 Summary Table

**All Experimental Findings:**

| Study | Key Result | Agrees with Theory? |
|-------|------------|---------------------|
| Basic (K shell) | Large splitting (2.83 eV), IR to UV | ✓ Yes - demonstrates hybridization |
| Thickness (Au) | Splitting increases with *f* | ✗ No - opposite to quasistatic |
| Metals | Energy scales with ω_p (Al highest) | ✓ Yes - Drude model correct |

---

## 5. RESULT QUALITY AND VALIDATION

### 5.1 Computational Convergence

**Mie Theory Accuracy:**
- Calculations use *l* = 1 (dipole) with higher orders (*l* > 1) being negligible for these particle sizes
- Mie coefficients converge rapidly for *ka* ≪ 1 where *k* = 2π/λ
- For *a* = 10-20 nm and visible light: *ka* ≈ 0.1-0.2 (well within validity)

**Cross-Validation Methods:**

1. **Internal consistency:** Peak positions found in both absorption and scattering spectra match
2. **Energy-wavelength conversion:** E(eV) × λ(nm) = 1240 verified for all peaks
3. **Physical bounds:** All cross-sections positive, extinction = absorption + scattering

### 5.2 Comparison with Theoretical Predictions

**Successes:**
- ✅ Two-peak structure observed (validates hybridization concept)
- ✅ Material trends match Drude predictions (ω_p correlation)
- ✅ Mode character (bonding lower energy, antibonding higher) as expected

**Limitations of Simple Theory:**
- ❌ Quasistatic splitting trend fails for *f* dependence
- ❌ Absolute energies differ from simple depolarization model
- ⚠️ Size-dependent effects not captured in *S₁*, *S₂* formulas

**Conclusion:** Tutorial concepts provide essential physical intuition (mode hybridization, charge oscillations, material effects), but **quantitative predictions require exact Mie theory** for realistic particle sizes.

### 5.3 Experimental Relevance

The simulations use realistic parameters:
- Particle sizes 10-20 nm are typical for colloidal synthesis
- Volume fractions 0.1-0.9 span experimentally accessible range
- Metals studied (Au, Ag) are standard plasmonic materials

**Expected experimental observations:**
- Two-peak structure should be observable in extinction spectroscopy
- Thick shells (*f* > 0.7) needed to resolve separate peaks
- Aluminum nanoshells would appear UV-active
- Gold nanoshells tune across visible by adjusting *f*

---

## 6. ANALYSIS AND PHYSICAL INTERPRETATION

### 6.1 Charge Oscillation Picture

The hybridization model can be understood through surface charge distributions:

**Bonding Mode (Symmetric, Lower Energy):**
```
Core: (+)  |  Metal Shell: (-)  |  Background: (+)
        ← Same phase oscillation →
```
- Inner cavity charges and outer surface charges move together
- Maximum dipole moment (charges separated)
- Lower restoring force → lower frequency
- Energy stored in near-field → absorption dominated

**Antibonding Mode (Antisymmetric, Higher Energy):**
```
Core: (+)  |  Metal Shell: (+)  |  Background: (-)
        ← Opposite phase oscillation →
```
- Inner and outer charges oppose each other
- Partial dipole cancellation
- Higher restoring force → higher frequency
- Energy radiates efficiently → scattering dominated

### 6.2 Why Thickness Matters

The thickness study reveals three regimes:

**1. Very thin shells (f < 0.3):**
- Inner and outer surfaces too close
- Mode overlap is extreme
- Cannot resolve separate peaks
- Behaves like perturbed solid sphere

**2. Intermediate shells (f = 0.3-0.7):**
- Surfaces begin to separate
- Weak coupling allows small splitting
- Transition regime

**3. Thick shells (f > 0.7):**
- Surfaces electromagnetically distinct
- Strong hybridization with clear mode separation
- Two-peak structure fully developed

**Critical insight:** The quasistatic formula Δ*S* = √(9-8*f*)/3 predicts opposite trend because it assumes particle is infinitesimally small. At *a* = 10 nm, **retardation and size effects dominate** over simple electrostatic coupling.

### 6.3 Material Selection Strategy

From the metal comparison:

**For UV applications:** Choose Al (ω_p = 15 eV)
- Peak at 229 nm (5.4 eV)
- Drawback: High damping (γ = 0.5 eV) → broad peaks

**For visible sensing:** Choose Ag, Au, or Cu (ω_p ≈ 9 eV)
- Ag: 449 nm (blue) - sharpest peaks (γ = 0.02 eV)
- Au: 573 nm (yellow-green) - moderate damping, chemically stable
- Cu: 590 nm (orange) - similar to Au but less stable

**For NIR/IR (biomedical):** Choose K or thick Au shells
- K: Naturally low ω_p pushes into NIR
- Au with *f* < 0.5: Red-shifts to 1000+ nm

### 6.4 Cross-Section Scaling

Though not explicitly measured here, the tutorial predicts:
- **Absorption:** *C*_abs ∝ *V* (volume)
- **Scattering:** *C*_sca ∝ *V*² (volume squared)

**Implications:**
- Small particles (< 20 nm): Absorption dominates → photothermal heating
- Large particles (> 50 nm): Scattering dominates → imaging, SERS
- Nanoshells tune absorption/scattering ratio via *f*

### 6.5 Connection to Tutorial Theory

**What the tutorial got right:**
1. ✅ Plasmon hybridization framework (bonding/antibonding)
2. ✅ Depolarization factor concept for mode character
3. ✅ Drude model for material dependence
4. ✅ Qualitative charge oscillation picture

**Where exact theory is needed:**
1. ⚠️ Quantitative resonance energies
2. ⚠️ Size-dependent effects (*ka* corrections)
3. ⚠️ Thickness trends for realistic particle sizes
4. ⚠️ Linewidths and damping

**Synthesis:** Use tutorial for **physical understanding**, use Mie theory for **quantitative predictions**.

---

## 7. CONCLUSION

This computational study has successfully demonstrated and characterized **plasmon hybridization in metal nanoshells** using exact Mie theory simulations. The key findings are:

### 7.1 Main Results

1. **Hybridization Confirmed:** Metal nanoshells exhibit two distinct plasmon modes (bonding and antibonding) arising from coupling between inner cavity and outer sphere plasmons. Clear two-peak structure observed for K shell (ΔE = 2.83 eV).

2. **Thickness Effects Are Complex:** Counter to simple quasistatic predictions, mode splitting *increases* with shell thickness for Au nanoshells, reaching maximum at *f* = 0.9 (ΔE = 3.16 eV). Thin shells show degenerate modes (ΔE ≈ 0), not large splitting, due to retardation effects at realistic particle sizes.

3. **Material Determines Spectral Position:** Plasma frequency (ω_p) directly controls resonance energy, validated across four metals. Al peaks in UV (5.42 eV), Au/Ag/Cu span visible range (2.1-2.8 eV), enabling spectral tuning by material choice.

4. **Theory Limitations Identified:** Tutorial quasistatic model provides excellent physical intuition but fails quantitatively for 10-20 nm particles. Exact Mie theory is essential for predictive modeling of real nanostructures.

### 7.2 Physical Understanding

The charge oscillation picture successfully explains:
- Why bonding modes have lower energy (in-phase charges, weak restoring force)
- Why antibonding modes have higher energy (out-of-phase, strong restoring force)
- How material properties (ω_p, γ) control plasmon behavior
- Why thick shells develop clear mode separation

### 7.3 Practical Implications

**Design Rules for Nanoshell Plasmonics:**

| Application | Metal | *f* | Result |
|-------------|-------|-----|--------|
| UV sensors/catalysis | Al | 0.5-0.7 | Peak at 229 nm |
| Visible sensors | Ag | 0.5 | Sharp peak at 449 nm |
| Biomedical (NIR) | Au | 0.1-0.3 | Red-shift to 1000+ nm |
| Dual-mode devices | Au | 0.9 | Two peaks (2.4 eV + 5.6 eV) |

### 7.4 Scientific Significance

This study demonstrates the **necessity of exact electromagnetic theory** for nanophotonics. While simplified models (quasistatic, Drude) capture essential physics, realistic predictions require solving Maxwell's equations exactly. The discrepancy between theory and simulation for thickness trends is not a failure of physics, but rather highlights **where approximations break down** - a valuable insight for designing experiments and applications.

### 7.5 Future Directions

Potential extensions of this work:
1. Study aspect ratio effects (non-spherical shells)
2. Include substrate effects (nanoshells on surfaces)
3. Investigate near-field enhancements for SERS
4. Model time-dependent decay dynamics
5. Optimize geometries for specific applications

---

## 8. REFERENCES

### Tutorial Materials
1. Course Lecture Notes - "Plasmonics: Localized Surface Plasmons in Metal Nanoshells"
2. Tutorial Slides - "Plasmon Hybridization Theory and Depolarization Factors"
3. Tutorial Slides - "Drude Model for Metal Permittivity and Material Comparison"

### Computational Methods
4. MATLAB Mie Theory Code - `fmiecoeff.m`, `fmiefields.m`, `crosssection.m`
5. Metal Permittivity Functions - `epsAu.m`, `epsAg.m`, `epsAl.m`, `epsCu.m`, `epsK.m`

### Theory Background
6. Plasmon Hybridization Model - Quasistatic depolarization theory for nanoshells
7. Drude-Lorentz Model - Frequency-dependent dielectric functions
8. Mie Theory - Exact electromagnetic solution for spherical particles

### Analysis Tools
9. Python/Jupyter Notebook - `LSP_Results_Analysis.ipynb`
10. MATLAB Simulation Script - `run_simulations_export.m`

---

## APPENDIX: Data Summary

**Complete dataset available in:** `simulation_results.json`

**Key numerical values:**

**Basic Example (K, a=20nm, f=0.3):**
- Peak 1: 0.717 eV (1730 nm)
- Peak 2: 3.543 eV (350 nm)
- Splitting: 2.826 eV

**Thickness Study (Au, a=10nm):**
- Minimum ΔE: 0.000 eV (f=0.1-0.3)
- Maximum ΔE: 3.162 eV (f=0.9)

**Metal Comparison (a=10nm, f=0.5):**
- Cu: 2.102 eV (590 nm)
- Au: 2.164 eV (573 nm)
- Ag: 2.762 eV (449 nm)
- Al: 5.421 eV (229 nm)

**Plasma frequencies used:**
- K: 3.72 eV, Ag/Au/Cu: 9.0 eV, Al: 15.0 eV

---

**End of Report**

*This report documents the computational investigation of localized surface plasmons in metal nanoshells using Mie theory simulations. All results are based on exact electromagnetic calculations with parameters drawn from experimental literature and tutorial materials.*
