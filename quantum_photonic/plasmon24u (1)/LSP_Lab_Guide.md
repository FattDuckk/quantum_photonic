# Plasmonics Computer Lab - LSP (Localized Surface Plasmons) Guide

## Overview of the Code

The main simulation file is `example_mie_spectrum_3d_shell_eV.m`, which simulates a **metal shell nanoparticle** illuminated by a plane wave using **Mie theory** - an exact analytical solution for electromagnetic scattering by spheres.

### Key Physics Concepts

**Mie Theory**: Solves Maxwell's equations exactly for spherical particles by expanding fields in spherical harmonics and matching boundary conditions at interfaces.

**Shell Structure**: 
- Inner core: dielectric (silica, air, etc.)
- Metal shell: plasmonic material (Au, Ag, Cu, Al, K)
- Outer medium: background (water, glass, air, etc.)

---

## Understanding the Current Setup

### Current Parameters (Line 16-24):
```matlab
a = 0.020;        % Outer sphere radius = 20 nm
f = 0.3;          % Volume fraction of metal = 30%
nc = 1.5;         % Core refractive index (silica-like)
nb = 1.5;         % Background refractive index  
metal = 'K';      % Metal type (Potassium)
```

**Geometry Calculation**:
- Outer radius: `a = 20 nm`
- Volume fraction metal: `f = 0.3` means 30% is metal
- Inner radius: `a*(1-f)^(1/3) = 20*(0.7)^(1/3) ≈ 17.7 nm`
- Shell thickness: `20 - 17.7 = 2.3 nm`

---

## What the Code Does

### 1. **Spectrum Calculation** (Figure 1)
- Calculates extinction, absorption, and scattering cross-sections vs photon energy (eV)
- **Extinction (ext)**: Total energy removed from incident beam (abs + sca)
- **Absorption (abs)**: Energy converted to heat in the metal
- **Scattering (sca)**: Energy re-radiated as light

### 2. **Peak Finding**
- Automatically identifies peaks in absorption and scattering
- These peaks correspond to **plasmon resonances**

### 3. **Angular Scattering** (Figure 2, top panels)
- Polar plots showing how light scatters in different directions
- 0° = forward, 180° = backward
- Three polarizations: unpolarized (u), perpendicular (s), parallel (p)

### 4. **Near-Field Plots** (Figure 2, bottom panels)
- Shows electric field magnitude |E| around the particle
- Red streamlines show **Poynting vector** (power flow)
- Reveals field enhancement and energy localization

---

## Task 1: Run the Example

### Expected Output:

**Figure 1**: Spectrum showing extinction (blue), absorption (orange), scattering (yellow) curves

**Figure 2**: For the two main peaks:
- Top: Angular scattering patterns
- Bottom: Near-field distributions with power flow

---

## Task 2: Interpret the Results

### Questions to Answer:

#### **What's the difference between the two peaks?**

Look for these characteristics:

1. **Absorption vs Scattering Dominance**:
   - Which peak has higher absorption relative to scattering?
   - Which is more "lossy" (more absorption)?

2. **E-field Patterns**:
   - **Dipolar mode**: Field pattern looks like classic dipole (+ on one side, - on other)
   - **Higher-order modes**: More complex patterns (quadrupole, octupole, etc.)
   - Count the number of field lobes/nodes

3. **Power Flow**:
   - Does power flow around the particle or radiate away?
   - Circulating flow suggests more localized (absorptive) mode
   - Radial flow suggests scattering mode

4. **Angular Distribution**:
   - Dipole: Forward-backward symmetric scattering
   - Higher orders: More complex angular patterns

### Physical Interpretation:

**Localized Surface Plasmon Resonances (LSPRs)** occur when:
- Conduction electrons in the metal oscillate collectively
- Resonance condition depends on geometry, metal properties, surrounding medium
- Different resonances correspond to different multipole moments (l=1: dipole, l=2: quadrupole, etc.)

**Two main resonances in a shell**:
1. **Bonding mode** (lower energy): Symmetric charge oscillation, typically more absorptive
2. **Antibonding mode** (higher energy): Antisymmetric charge oscillation, typically more scattering

---

## Task 3: Explore Parameters

### Recommended Fixed Parameters:
```matlab
nc = 1;           % Core permittivity (air)
nb = 1;           % Background permittivity (vacuum/air)  
a = 0.010;        % Outer radius = 10 nm (0.010 μm)
```

### A. Effect of Shell Thickness

Try these volume fractions:
```matlab
f = 0.1;    % Thin shell (9.7 nm inner radius, 0.3 nm shell)
f = 0.3;    % Medium shell (8.9 nm inner radius, 1.1 nm shell)
f = 0.5;    % Thick shell (7.9 nm inner radius, 2.1 nm shell)
f = 0.7;    % Very thick shell (6.7 nm inner radius, 3.3 nm shell)
f = 0.9;    % Almost solid (4.8 nm inner radius, 5.2 nm shell)
```

**What to observe**:
- How do peak positions shift?
- How do peak intensities change?
- What happens to the separation between peaks?

**Physical model**: 
- **Plasmon hybridization**: Inner and outer surface plasmons couple
- Thin shell: Strong coupling → large splitting between bonding/antibonding
- Thick shell: Weak coupling → modes converge toward solid sphere limit

---

### B. Effect of Metal Type

Available metals with their plasma frequencies (approximate):

```matlab
metal = 'Ag';   % Silver - best plasmonic material (low loss, visible range)
metal = 'Au';   % Gold - stable, good in red/NIR, more lossy than Ag
metal = 'Al';   % Aluminum - UV plasmons, higher loss
metal = 'Cu';   % Copper - similar to gold, less stable
metal = 'K';    % Potassium - alkali metal, lower plasma frequency
```

**What to observe**:
- Peak positions (energy/wavelength)
- Peak widths (narrower = lower damping/loss)
- Relative absorption vs scattering strengths

**Physical explanation**:
- Resonance frequency depends on plasma frequency: ω_p = √(ne²/ε₀m)
- Different metals have different electron densities and damping rates
- Noble metals (Au, Ag) have filled d-bands → better plasmonic properties

---

### C. Optional: Core Index Effects

```matlab
nc = 1.0;    % Air core
nc = 1.5;    % Silica-like core  
nc = 2.0;    % High-index core
nc = 3.5;    % Semiconductor core
```

**What to observe**:
- Red-shift with increasing core index
- Changes in mode character

**Why**: Higher permittivity inside pushes more field into the core → effective longer wavelength

---

## Modified Code Templates

### Template 1: Shell Thickness Study
```matlab
%% Shell thickness parameter scan
a = 0.010;  % 10 nm outer radius
nc = 1;     % air core
nb = 1;     % air background
metal = 'Au';

% Try different shell thicknesses
f_values = [0.1, 0.3, 0.5, 0.7, 0.9];

for i = 1:length(f_values)
    f = f_values(i);
    % ... rest of code runs ...
    % Save or note the peak positions
end
```

### Template 2: Metal Comparison
```matlab
%% Metal comparison
a = 0.010;
f = 0.5;    % Fixed geometry
nc = 1;
nb = 1;

metals = {'Ag', 'Au', 'Al', 'Cu'};

for i = 1:length(metals)
    metal = metals{i};
    % ... rest of code runs ...
end
```

---

## Important Code Sections to Modify

**Lines 16-20**: Main parameters
```matlab
a=0.010;     % Change outer radius (in microns)
f=0.5;       % Change volume fraction (0 to 1)
nc=1;        % Change core index
nb=1;        % Change background index  
metal='Au';  % Change metal: 'Au', 'Ag', 'Al', 'Cu', 'K'
```

**Lines 26-37**: Output options
```matlab
do_fields=1;     % Set to 0 to skip field plots (faster)
slicef='e';      % Field to plot: 'e' for |E|, 'h' for |H|
```

---

## Tips for Success

1. **Run baseline first**: Understand what "normal" looks like
2. **Change one parameter at a time**: Easier to see effects
3. **Document everything**: 
   - Screenshot figures
   - Note peak positions (eV or nm)
   - Note peak heights
   - Observe field patterns
4. **Look for trends**: How parameters systematically affect results
5. **Think physically**: Connect observations to charge distribution models

---

## Physical Models to Reference

### Plasmon Hybridization Model:
- Shell has inner and outer surfaces
- Each has a plasmon mode (like atomic orbitals)
- They hybridize into bonding (symmetric) and antibonding (antisymmetric) modes
- Energy splitting depends on coupling strength (shell thickness)

### Quasi-static Approximation (small particles):
For sphere, resonance when: **ε(ω) = -2ε_background**

For shell, modified by geometry factor depending on aspect ratio.

---

## Common Issues

1. **No peaks visible**: Try different energy range or different metal
2. **Code runs slow**: Set `do_fields=0` for faster calculation
3. **Peaks outside range**: Adjust eV limits or try different metal/size

---

## Expected Observations Summary

| Parameter ↑ | Peak Position | Peak Strength | Peak Width |
|-------------|---------------|---------------|------------|
| Shell thickness | Red-shift initially | Changes | Changes |
| Particle size | Red-shift | Increases | Increases |
| Core index | Red-shift | Changes | ~ |
| Background index | Red-shift | ~ | ~ |

Metal dependence: Ag (visible, narrow) → Au (red-NIR) → Al (UV) → K (IR)

---

## Ready to Start!

1. First run the code as-is to see the baseline
2. Take screenshots and notes
3. Systematically vary parameters
4. Compare and explain your observations using physical models

Good luck with your lab!
