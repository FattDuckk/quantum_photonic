# Quick Start Guide - Plasmonics Lab

## Files Created for You

I've created several files to help you with the lab:

1. **LSP_Lab_Guide.md** - Comprehensive guide with theory and instructions
2. **example_shell_thickness_study.m** - Automated shell thickness parameter scan
3. **example_metal_comparison.m** - Automated metal comparison study

## Step-by-Step Instructions

### Step 1: Run the Original Example (5 minutes)

```matlab
% In MATLAB command window:
cd '/home/fatduck/git/quantum_photonic/quantum_photonic/plasmon24u (1)'
example_mie_spectrum_3d_shell_eV
```

**What you'll see:**
- Figure 1: Spectrum showing extinction (blue), absorption (orange), scattering (yellow)
- Figure 2: Angular scattering patterns (top) and near-field distributions (bottom) for two peaks

**What to do:**
- Take screenshots of both figures
- Note the peak positions (eV values)
- Look at the field patterns - count the lobes
- Observe the power flow (red streamlines)

---

### Step 2: Understand the Results (10 minutes)

**Look at the two peaks in Figure 2:**

**Peak 1 (usually lower energy):**
- Is it more absorbing or scattering?
- Field pattern: dipolar (2 lobes) or higher order (4+ lobes)?
- Power flow: circulating or radiating?
- Angular pattern: simple or complex?

**Peak 2 (usually higher energy):**
- Same questions as above
- How does it differ from Peak 1?

**Physical meaning:**
- These are the **bonding** and **antibonding** modes of the shell
- Like molecular orbitals, but for plasmons
- They arise from coupling between inner and outer surface plasmons

---

### Step 3: Shell Thickness Study (20-30 minutes)

Run the automated scan:

```matlab
example_shell_thickness_study
```

This will:
- Test 5 different shell thicknesses (f = 0.1, 0.3, 0.5, 0.7, 0.9)
- Plot all spectra together
- Create summary plots showing trends
- Print a table of peak positions

**What to observe:**
1. **Peak separation**: Do the two peaks move together or apart as shell gets thicker?
2. **Peak positions**: Do they red-shift or blue-shift?
3. **Peak intensities**: Which mode gets stronger/weaker?

**Expected behavior:**
- **Thin shell (f=0.1)**: Large energy splitting between modes
- **Thick shell (f=0.9)**: Small splitting, converging toward solid sphere
- This is **plasmon hybridization** - inner and outer surface modes couple

**For your report:**
- Include the spectrum plot (Figure 1)
- Include the summary plots (Figure 2)
- Paste the summary table
- Explain the trend using hybridization model

---

### Step 4: Metal Comparison (20-30 minutes)

Run the metal comparison:

```matlab
example_metal_comparison
```

This will:
- Compare Ag, Au, Al, Cu
- Show extinction, absorption, scattering for each
- Calculate peak widths and quality factors
- Create comparison plots

**What to observe:**
1. **Peak positions**: Which metal has highest/lowest energy resonances?
2. **Peak widths**: Which metal has sharpest peaks? (lower loss)
3. **Absorption vs scattering**: Which metal is more lossy?

**Expected ranking (quality):**
1. **Silver (Ag)**: Best - sharpest peaks, low loss, visible range
2. **Gold (Au)**: Good - works in red/NIR, more stable than Ag
3. **Copper (Cu)**: Moderate - similar to Au but less stable
4. **Aluminum (Al)**: High loss - but works in UV

**For your report:**
- Include comparison plots
- Paste the summary table
- Explain why Ag is best (filled d-bands, low damping)
- Explain why different metals have different resonance energies (plasma frequency)

---

### Step 5: Manual Exploration (Optional, 30+ minutes)

Open the original file and manually change parameters:

```matlab
edit example_mie_spectrum_3d_shell_eV
```

**Try these:**

#### Vary outer radius (Line 16):
```matlab
a = 0.005;  % 5 nm
a = 0.010;  % 10 nm (default)
a = 0.020;  % 20 nm
a = 0.050;  % 50 nm
```
Effect: Larger particles → red-shift, more scattering

#### Vary background index (Line 19):
```matlab
nb = 1.0;   % Air/vacuum
nb = 1.33;  % Water
nb = 1.5;   % Glass
nb = 2.0;   % High-index material
```
Effect: Higher index → red-shift (plasmons "see" effective longer wavelength)

#### Vary core index (Line 18):
```matlab
nc = 1.0;   % Air core (hollow)
nc = 1.5;   % Silica core
nc = 3.5;   % Semiconductor core
```
Effect: Similar to background, but affects mode splitting

---

## What to Include in Your Report

### Part 1: Baseline Results (Run original example)
- Screenshot of spectrum (Figure 1)
- Screenshots of angular scattering and fields (Figure 2)
- Description of the two peaks:
  - Peak positions (eV and nm)
  - Field patterns (dipole/quadrupole/etc.)
  - Absorption vs scattering character
  - Power flow patterns

### Part 2: Shell Thickness Study
- Spectrum plot showing all thicknesses
- Summary plots (peak positions vs thickness)
- Table of peak data
- **Interpretation**: Explain using plasmon hybridization:
  - Why does splitting change?
  - What's the physical picture? (Draw charge distribution if helpful)

### Part 3: Metal Comparison
- Comparison plots
- Table of peak positions and widths
- **Interpretation**: Explain:
  - Why different metals have different resonance energies
  - Why some have sharper peaks (relate to damping, d-band transitions)
  - Which is best for plasmonics and why

### Part 4: Physical Models
- **Plasmon hybridization**: Explain how shell modes arise
- **Quasi-static model**: For small particles, resonance when ε_metal = -2ε_background
- **Size effects**: Retardation, radiation damping for larger particles
- **Material effects**: Plasma frequency, interband transitions, damping

---

## Common Questions & Answers

**Q: Why are there two peaks?**
A: The shell has two surfaces (inner and outer). Each can support a surface plasmon. They couple to form bonding (lower energy, symmetric) and antibonding (higher energy, antisymmetric) modes.

**Q: Which peak is dipolar?**
A: Usually both are dipolar (l=1). Look at the field pattern - should see two main lobes. Higher-order modes (l=2, 3...) are weaker for small particles.

**Q: Why does one peak absorb more than scatter?**
A: Depends on mode character. Modes with field concentrated in the lossy metal → more absorption. Modes with field outside → more scattering.

**Q: How does this relate to the Fröhlich condition?**
A: For a sphere, resonance occurs when ε(ω) ≈ -2ε_outside. For a shell, it's modified by geometry, but same basic idea - surface charge oscillation is resonant at specific frequency.

**Q: Why does peak red-shift with size?**
A: For small increases: quasi-static still applies. For large particles: retardation effects (phase delay across particle) shift resonance. Also radiation damping increases.

**Q: Why is silver better than gold?**
A: Ag has filled d-bands further from visible range → less interband absorption. Lower damping → sharper resonances → stronger field enhancement.

---

## Tips for Success

1. ✅ **Run code first, understand later** - see what happens before diving deep
2. ✅ **Change one thing at a time** - easier to see cause and effect
3. ✅ **Take lots of screenshots** - you'll need them for your report
4. ✅ **Make tables** - organize your peak position/width data
5. ✅ **Plot trends** - peak position vs parameter is clearer than words
6. ✅ **Think physically** - every observation has a physical explanation
7. ✅ **Use the hybridization model** - it explains shell behavior beautifully

---

## Troubleshooting

**Problem: Code is slow**
```matlab
do_fields = 0;  % Skip field calculations (Line 27)
```

**Problem: No peaks visible**
- Try different metal or different size
- Check energy range covers resonance

**Problem: Peaks outside energy range**
```matlab
eV_min = 0.5;   % Lower energy bound
eV_max = 6.0;   % Higher energy bound
eV = linspace(eV_min, eV_max, 500)';
L = 1.24./eV;
```

**Problem: Error about material data**
- Wavelength might be outside tabulated range
- Use different metal or adjust wavelength range

---

## Additional Resources

### Recommended Reading:
1. **Plasmon Hybridization**: Prodan et al., Science 302, 419 (2003)
2. **Mie Theory**: Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
3. **Plasmonics Basics**: Maier, "Plasmonics: Fundamentals and Applications"

### Physical Picture:
```
Thin Shell (weak coupling):        Thick Shell (strong coupling):
    Inner    Outer                      Inner    Outer
     +++      ---                         +++      ---
  
  ω_-  ω_+  (large splitting)       ω_-  ≈  ω_+  (small splitting)
  
  Antibonding  Bonding              Both converge to dipole
```

---

## Time Management

- Step 1-2: 15 minutes (run and understand baseline)
- Step 3: 30 minutes (shell thickness study)
- Step 4: 30 minutes (metal comparison)
- Step 5: 30 minutes (optional explorations)
- **Total: ~2 hours of computation**
- **Write-up: 2-3 hours**

---

Good luck! The code is well-written and should work smoothly. Focus on understanding the physics, not debugging. 🚀
