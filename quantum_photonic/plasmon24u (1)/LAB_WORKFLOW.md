# Plasmonics Computer Lab - Step-by-Step Workflow

## Overview
This lab uses **existing Mie theory code** in MATLAB to simulate metal nanoshell plasmons.
All physics calculations are done by the provided functions - you just need to run them systematically and interpret results.

---

## Part 1: Understanding the Basic Example (Task 1)

### Step 1.1: Run the Original Example
```matlab
% In MATLAB, run:
example_mie_spectrum_3d_shell_eV
```

**What you should get:**
- **Figure 1**: Spectrum showing absorption (C_abs) and scattering (C_sca) vs energy
- **Figure 2+**: Near-field plots showing |E| field distribution with power flow (red lines)
- **Peaks**: You should see 2 distinct peaks in the spectrum

### Step 1.2: Understand What's Being Shown

**The Spectrum (Figure 1):**
- X-axis: Energy in eV (or wavelength)
- Y-axis: Cross-section normalized by volume and wavenumber
- **Two curves**: 
  - Absorption (solid line) - energy absorbed by metal
  - Scattering (dashed line) - energy scattered away

**Angular Scattering:**
- 0° = forward scattering (transmission direction)
- 180° = backward scattering (reflection direction)
- Shows which directions the light is scattered

**Near-field Plots (Figure 2+):**
- Color: |E| field magnitude
- Red streamlines: Poynting vector (power flow direction)
- Shows field distribution at each peak wavelength

### Step 1.3: Modify Parameters (Optional First Test)
Try changing in `example_mie_spectrum_3d_shell_eV.m`:
```matlab
a = 0.010;      % Try 10 nm instead of 20 nm
f = 0.5;        % Try different shell thickness
metal = 'Au';   % Try Au instead of K
nc = 1;         % Core = air
nb = 1;         % Background = air
```

---

## Part 2: Interpret the Results (Task 2)

### Questions to Answer:

#### Q1: What's the difference between the two peaks?

**Things to compare:**
1. **Absorption vs Scattering ratio**
   - Which peak has more absorption?
   - Which peak has more scattering?

2. **Field distributions**
   - Look at near-field plots for each peak
   - Where is field concentrated? (inside? outside? at surface?)
   - Are there hot spots?

3. **Power flow patterns**
   - Where does energy go?
   - Is energy flowing into particle or around it?

4. **Physical interpretation**
   - Lower energy peak: Bonding (symmetric) plasmon mode
   - Higher energy peak: Antibonding (antisymmetric) plasmon mode
   - This is **plasmon hybridization** - cavity mode couples with sphere mode

#### Q2: Which modes are they?

**How to determine:**
1. Look at angular scattering plot
   - Dipole mode: Mostly forward/backward (like p² pattern)
   - Quadrupole: More complex pattern with nulls
   - Higher order: Multiple lobes

2. Check field symmetry in near-field plots
   - Dipole: Two lobes (+ and -)
   - Quadrupole: Four lobes (+ - + -)

3. Look at expansion coefficients
   - Dominated by a₁, b₁? → Dipole (l=1)
   - Dominated by a₂, b₂? → Quadrupole (l=2)

---

## Part 3: Systematic Parameter Study (Task 3)

### Study 3.1: Shell Thickness Effect

**Use the script:** `example_shell_thickness_study.m`

**What it does:**
- Fixes: outer radius a=10nm, core index nc=1, background nb=1, metal=Au
- Varies: shell thickness through volume fraction f = 0.1, 0.3, 0.5, 0.7, 0.9
- Note: f = 1 means solid sphere, f → 0 means very thin shell

**To run:**
```matlab
example_shell_thickness_study
```

**What to observe:**
1. **Peak positions shift** with shell thickness
   - Thin shell (f=0.1): Large splitting between modes
   - Thick shell (f=0.9): Modes converge
   - Solid sphere (f→1): Single dipole mode

2. **Peak strengths change**
   - Which peak dominates for thin vs thick shells?

3. **Physical model (Plasmon Hybridization):**
   - Cavity mode (inner surface) + sphere mode (outer surface)
   - Thin shell → strong coupling → large splitting
   - Thick shell → weak coupling → small splitting

**Questions to answer:**
- How does peak position depend on f?
- How does mode splitting depend on f?
- Can you relate this to electromagnetic boundary conditions?

### Study 3.2: Metal Type Effect

**Use the script:** `example_metal_comparison.m`

**What it does:**
- Fixes: outer radius a=10nm, f=0.5, nc=1, nb=1
- Varies: Metal type (Au, Ag, Al, Cu)

**To run:**
```matlab
example_metal_comparison
```

**What to observe:**
1. **Peak positions** for different metals
   - Related to plasma frequency ωₚ of each metal
   - Related to screening parameter

2. **Peak widths** (Q-factor)
   - Related to damping (loss) in metal
   - Ag typically sharpest (lowest loss)
   - Al typically broadest (higher loss)

3. **Peak strengths**
   - Related to how negative ε can get
   - Related to quality of plasmonic response

**Physical background:**
- Drude model: ε(ω) = 1 - ωₚ²/(ω² + iγω)
- Different metals have different ωₚ and γ
- This shifts resonances and changes loss

**Questions to answer:**
- Which metal has peaks at highest energy? Why?
- Which metal has sharpest peaks? (best Q-factor)
- How does this relate to the Drude parameters?

### Study 3.3: Core Index Effect (Optional)

If you want to explore further:

**Modify either script to vary nc:**
```matlab
nc_values = [1.0, 1.5, 2.0, 2.5];  % Air, silica, high-index
```

**What to observe:**
- How does core dielectric affect mode coupling?
- Effect on bonding vs antibonding mode shifts

---

## Part 4: Relating to Physical Models

### Plasmon Hybridization Model

**Analogy to molecular orbitals:**
- Two plasmon modes: cavity (inner surface) + sphere (outer surface)
- Coupling creates: bonding (lower E) + antibonding (higher E)
- Splitting depends on overlap (shell thickness)

### Charge Distribution Picture

**For bonding mode (lower energy):**
- Charges on inner and outer surface oscillate in phase
- Creates large dipole moment
- Lower restoring force → lower frequency

**For antibonding mode (higher energy):**
- Charges on inner and outer surface oscillate out of phase
- Partially cancels dipole moment
- Higher restoring force → higher frequency

### Quasi-static Approximation

For small particles (a << λ):
- Fröhlich condition: ε_metal + 2ε_background = 0
- For shells: More complex, involves ε_core too
- Explains why peak position depends on all three dielectrics

---

## Key Functions Reference

**These functions do the actual physics:**

1. **`fmiecoeff(N, L, av, m, mu)`** - Calculates Mie coefficients
   - Uses spherical Bessel functions
   - Matches boundary conditions at each interface
   - Returns expansion coefficients aₙ, bₙ

2. **`crosssection(coeffs, L, av, m)`** - Calculates cross-sections
   - From Mie coefficients → total, absorption, scattering
   - Uses: C_sca = Σ(2n+1)(|aₙ|² + |bₙ|²)
   - C_abs = C_total - C_sca

3. **`epsAu(L)`, `epsAg(L)`, etc.** - Metal permittivity
   - Experimental data from literature
   - Returns complex ε(λ) for each metal
   - Real part: dispersion, Imaginary part: loss

4. **`fmiefields(coeffs, ...)`** - Calculates near fields
   - Evaluates spherical harmonics
   - Gives E and H fields in space

---

## Tips for Lab Report

### What to include:

1. **Introduction**
   - Brief background on localized surface plasmons
   - Mention Mie theory and plasmon hybridization

2. **Methods**
   - State that you used Mie theory code
   - Specify fixed parameters (a=10nm, nc=1, nb=1)
   - Describe what you varied

3. **Results**
   - Include spectrum plots
   - Include near-field plots (at least for one case)
   - Include parameter variation plots
   - Table of peak positions

4. **Discussion**
   - Identify modes (dipole, quadrupole, etc.)
   - Explain difference between bonding/antibonding
   - Relate peak shifts to physical models
   - Compare metals based on optical properties

5. **Conclusions**
   - Summary of main findings
   - Physical insights gained

### Common Questions to Address:

- **Why two peaks?** → Plasmon hybridization
- **Why peaks shift with thickness?** → Coupling strength changes
- **Why different metals give different peaks?** → Different ωₚ and γ
- **Absorption vs scattering?** → Depends on size, mode character, material loss
- **Field concentration where?** → At metal-dielectric interfaces (hot spots)

---

## Troubleshooting

**If you get errors:**
- Make sure you're in the correct directory
- Check that all .m files are in MATLAB path
- Verify metal functions exist (epsAu.m, epsAg.m, etc.)

**If plots look wrong:**
- Check units (eV vs wavelength)
- Verify normalization (C/V/k vs Q vs C)
- Make sure particle size is reasonable (1-100 nm typically)

**If simulations are slow:**
- Reduce number of wavelength points
- Set `do_fields=0` to skip field plots
- Reduce N (number of orders) if convergence is ok

---

## Next Steps

1. ✅ Read this guide
2. ✅ Run `example_mie_spectrum_3d_shell_eV.m` to see basic output
3. ✅ Understand what each figure shows
4. ✅ Run `example_shell_thickness_study.m` for systematic study
5. ✅ Run `example_metal_comparison.m` for metal comparison
6. ✅ Analyze results and answer questions
7. ✅ Write lab report with physical interpretation

Good luck with your lab! 🎯
