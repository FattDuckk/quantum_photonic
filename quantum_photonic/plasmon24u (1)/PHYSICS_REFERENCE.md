# Physics Reference - Localized Surface Plasmons in Nanoshells

## Quick Theory Summary

### What is Mie Theory?
- **Exact** solution to Maxwell's equations for spherical particles
- Expands fields in vector spherical harmonics
- Matches boundary conditions at each interface
- Valid for **any** size (not just small particles)

### Expansion:
```
E, H = Σ [aₙ M_n + bₙ N_n]  (spherical harmonics)
      n=1,2,3...

n=1 → Dipole mode
n=2 → Quadrupole mode  
n=3 → Octupole mode
...
```

---

## Plasmon Resonance Condition

### For a solid sphere (simple case):
**Fröhlich condition:** ε_metal(ω) + 2ε_background = 0

For Drude metal: ε = 1 - ωₚ²/ω²

This gives resonance at: **ω_res = ωₚ/√(1 + 2ε_b)**

### For a nanoshell (more complex):
Two resonances arise from coupling of:
1. **Cavity mode** - plasmon on inner surface (core/metal interface)
2. **Sphere mode** - plasmon on outer surface (metal/background interface)

**Hybridization creates:**
- **Bonding mode** (ω₋): Lower energy, symmetric charge oscillation
- **Antibonding mode** (ω₊): Higher energy, antisymmetric charge oscillation

**Mode splitting:** Δω ∝ exp(-d/R) where d = shell thickness

---

## Physical Picture: Charge Oscillations

### Bonding Mode (Lower Energy Peak)
```
Core  |  Metal  |  Background
      |         |
  +   |    -    |    +         ← Charges in phase
      |         |
```
- Inner and outer charges oscillate **together**
- Large net dipole moment
- Lower restoring force → **lower frequency**
- Typically has more **absorption** (energy stored in oscillation)
- **Red-shifted** (lower energy, longer wavelength)

### Antibonding Mode (Higher Energy Peak)
```
Core  |  Metal  |  Background
      |         |
  +   |    +    |    -         ← Charges opposite phase
      |         |
```
- Inner and outer charges oscillate **opposite**
- Partial cancellation of dipole
- Higher restoring force → **higher frequency**
- Typically has more **scattering** (radiates efficiently)
- **Blue-shifted** (higher energy, shorter wavelength)

---

## Shell Thickness Effects

### Thin Shell (f → 0, small volume fraction)
- Strong coupling between inner/outer surfaces
- **Large mode splitting:** Δω is large
- Bonding mode: Very red-shifted
- Antibonding mode: Very blue-shifted
- Field penetrates through shell → strong interaction

### Thick Shell (f → 1, approaches solid sphere)
- Weak coupling between surfaces (inner surface too far away)
- **Small mode splitting:** Δω → 0
- Both modes converge to solid sphere dipole resonance
- Effectively just outer surface plasmon dominates

### Quantitative: Plasmon Hybridization Model
**Energy splitting:** 

Δω/ω₀ ≈ exp(-d/R)

where:
- d = shell thickness
- R = inner radius
- ω₀ = characteristic plasmon frequency

---

## Metal Type Effects

### Drude Model:
ε(ω) = ε∞ - ωₚ²/(ω² + iγω)

Where:
- **ωₚ**: Plasma frequency (sets energy scale)
- **γ**: Damping rate (sets linewidth/loss)
- **ε∞**: High-frequency permittivity (background polarization)

### Metal Comparison:

| Metal | ωₚ (eV) | γ (eV) | Notes |
|-------|---------|--------|-------|
| **Ag** | ~9 | ~0.02 | Lowest loss, sharpest peaks |
| **Au** | ~9 | ~0.07 | Moderate loss, red-visible peaks |
| **Al** | ~15 | ~0.5 | High loss, UV peaks |
| **Cu** | ~9 | ~0.1 | Similar to Au but oxidizes |
| **K** | ~4 | ~0.03 | Low ωₚ, IR peaks |

### How Metal Affects Plasmons:

1. **Peak Position** (wavelength/energy)
   - Higher ωₚ → higher frequency peaks (blue-shift)
   - Ag, Au, Cu: Similar ωₚ → visible range
   - Al: Higher ωₚ → UV range
   - K: Lower ωₚ → IR range

2. **Peak Width** (Q-factor)
   - Q ≈ ω_res/γ
   - Lower γ → sharper peaks
   - Ag has lowest γ → best Q-factor
   - Al has high γ → broad peaks

3. **Peak Strength**
   - Depends on how negative Re(ε) can get
   - Depends on Im(ε) - higher loss reduces peak
   - Trade-off between resonance strength and loss

---

## Absorption vs Scattering

### Size Dependence:
- **Small particles** (d << λ): Absorption dominates
  - Quasi-static limit
  - C_abs ∝ V (volume)
  - C_sca ∝ V² (volume squared)
  
- **Large particles** (d ~ λ): Scattering dominates
  - Retardation effects important
  - Radiation damping increases

### For 10 nm particle:
- Typically in quasi-static or slightly retarded regime
- Both absorption and scattering present
- Bonding mode: More absorption (trapped energy)
- Antibonding mode: More scattering (radiates)

---

## Field Enhancement

### Where are "hot spots"?
- **At metal-dielectric interfaces**
- Especially at points of high curvature
- For shells: Both inner and outer surfaces

### Why important?
- SERS (Surface Enhanced Raman Scattering)
- Enhanced absorption for photovoltaics
- Nonlinear optics
- Hot electron generation

### Typical enhancement factors:
- |E|²/|E₀|² can reach 10² - 10⁴
- Depends on:
  - Particle shape (sharper → stronger)
  - Material loss (lower loss → stronger)
  - Distance from surface (decays exponentially)

---

## Angular Scattering Patterns

### Dipole mode (l=1):
```
     ↑ k (incident)
     
  /--|--\    Mostly forward/backward
 |   •   |   
  \--|--/    Pattern: cos²θ like
```

### Quadrupole mode (l=2):
```
     ↑ k (incident)
     
  /     \    Multiple lobes
 | • - • |   Nulls at certain angles
  \     /    More complex pattern
```

### Higher order modes:
- More lobes
- More nulls
- Increasingly complex

---

## Normalization & Units

### Cross-sections:
- **C_ext** = C_abs + C_sca (total extinction)
- Units: area (e.g., nm²)

### Efficiency:
- **Q** = C / (πa²) (dimensionless)
- Normalized by geometric cross-section

### Normalized by volume:
- **C/V/k** = C / (Volume × wavenumber)
- Makes comparison easier between different sizes
- Removes trivial size dependence

---

## Key Equations

### Mie Coefficients (for homogeneous sphere):
```
       ψₙ(mx)ψₙ'(x) - mψₙ(x)ψₙ'(mx)
aₙ = ─────────────────────────────────
     ψₙ(mx)ξₙ'(x) - mξₙ(x)ψₙ'(mx)

       mψₙ(mx)ψₙ'(x) - ψₙ(x)ψₙ'(mx)
bₙ = ─────────────────────────────────
     mψₙ(mx)ξₙ'(x) - ξₙ(x)ψₙ'(mx)
```
where:
- ψₙ, ξₙ are Riccati-Bessel functions
- x = 2πa/λ (size parameter)
- m = n_particle/n_medium (relative refractive index)

### Cross-sections from coefficients:
```
         2π   ∞
C_sca = ───── Σ (2n+1)(|aₙ|² + |bₙ|²)
         k²  n=1

         2π   ∞
C_ext = ───── Σ (2n+1)Re(aₙ + bₙ)
         k²  n=1

C_abs = C_ext - C_sca
```

### For multilayer spheres:
- Recursively apply boundary conditions
- Start from innermost layer
- Each layer adds complexity to aₙ, bₙ formulas
- **The code handles this automatically!**

---

## Interpretation Checklist

When analyzing your results, ask:

### About the spectrum:
- [ ] How many peaks do I see?
- [ ] Which peak is absorption-dominated?
- [ ] Which peak is scattering-dominated?
- [ ] What are the peak wavelengths/energies?
- [ ] How wide are the peaks? (Q-factor)

### About the modes:
- [ ] Are they dipole, quadrupole, or higher?
- [ ] What's the mode splitting (energy difference)?
- [ ] How does splitting change with shell thickness?

### About the fields:
- [ ] Where is |E| strongest? (hot spots)
- [ ] Is field concentrated inside, outside, or at interface?
- [ ] What's the field pattern? (symmetric? antisymmetric?)
- [ ] Where does power flow? (Poynting vector)

### About physical models:
- [ ] Can I explain peak positions using Fröhlich condition?
- [ ] Can I explain mode splitting using hybridization model?
- [ ] Can I explain metal differences using Drude parameters?
- [ ] Do my observations match charge oscillation picture?

---

## Common Pitfalls

❌ **"The code calculates plasmons"** 
✅ The code solves Maxwell's equations exactly. Plasmons emerge as resonances.

❌ **"Thin shell means small f"**
✅ **Correct!** f is volume fraction. Thin shell → small f → large splitting.

❌ **"Red-shift means lower frequency"**
✅ **Correct!** Lower energy = lower frequency = longer wavelength = red-shift.

❌ **"Absorption and scattering are independent"**
✅ They're coupled through energy conservation: C_ext = C_abs + C_sca

❌ **"Bonding mode always absorbs more"**
✅ Usually, but depends on parameters. Check your actual results!

---

## Further Reading

- **Bohren & Huffman** - *Absorption and Scattering of Light by Small Particles* (Mie theory bible)
- **Prodan et al., Science 302, 419 (2003)** - Plasmon hybridization model
- **Maier** - *Plasmonics: Fundamentals and Applications*
- **Jackson** - *Classical Electrodynamics* (Multipole expansion, Ch. 9)

---

Good luck interpreting your results! 🔬✨
