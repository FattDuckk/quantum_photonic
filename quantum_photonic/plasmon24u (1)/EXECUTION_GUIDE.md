# Step-by-Step Execution Guide

## 🎯 Workflow Overview

We'll take this slowly and systematically:

1. **Run MATLAB simulations** → Get numerical data
2. **Analyze in Python notebook** → Visualize and interpret
3. **Compare with theory** → Using your tutorial formulas
4. **Write final report** → Structured format

---

## Step 1: Run MATLAB Simulations

### Option A: Run the comprehensive script (Recommended)

Open MATLAB and execute:

```matlab
cd '/home/fatduck/git/quantum_photonic/quantum_photonic/plasmon24u (1)'
run_simulations_export
```

This will:
- Run all three experiments (basic, thickness study, metal comparison)
- Save results to `simulation_results.json`
- Print summary tables to the console
- Generate figures

**Time: ~5-10 minutes**

### Option B: Run each example individually

If you want to see each simulation separately:

```matlab
% 1. Basic example
example_mie_spectrum_3d_shell_eV

% 2. Thickness study  
example_shell_thickness_study

% 3. Metal comparison
example_metal_comparison
```

Then manually record the data in the notebook.

---

## Step 2: Open and Run the Analysis Notebook

1. Open `plasmon_analysis.ipynb` in VS Code or Jupyter

2. If you ran `run_simulations_export.m`, add this cell at the beginning:

```python
import json

# Load MATLAB results
with open('simulation_results.json', 'r') as f:
    matlab_results = json.load(f)

print("MATLAB data loaded successfully!")
print("Available experiments:", list(matlab_results.keys()))
```

3. Run through each cell systematically:
   - **Part 1**: Basic example analysis
   - **Part 2**: Thickness study with theory comparison  
   - **Part 3**: Metal comparison
   - **Part 4**: Cross-section analysis
   - **Part 5**: Physical interpretation

4. Fill in observations and answer questions in markdown cells

---

## Step 3: Compare with Theory

The notebook already includes:

- **Hybridization model predictions** (S and A parameters from your tutorial)
- **Resonance condition formulas** (from Fröhlich condition)
- **Metal property database** (ωₚ and γ values)
- **Cross-section scaling** (absorption vs scattering)

Compare MATLAB (exact Mie theory) vs Tutorial (approximate quasistatic model):
- Do trends match?
- Where do they agree/disagree?
- Why? (hint: quasistatic assumes particle << wavelength)

---

## Step 4: Write Final Report

Once analysis is complete, I'll help you structure the final report:

### Report Structure:

1. **Aim**
   - What you're investigating
   - Why it matters

2. **Background**  
   - LSP physics from tutorials
   - Plasmon hybridization model
   - Quasistatic theory equations

3. **Methods**
   - Mie theory simulation (exact solution)
   - Parameters used (size, materials, geometry)
   - What was varied (f, metal type)

4. **Results**
   - Spectral data (peak positions, widths)
   - Field distributions
   - Trends with parameters
   - Figures and tables

5. **Result Quality**
   - Convergence (N orders sufficient?)
   - Comparison with theory
   - Systematic trends vs noise

6. **Analysis**
   - Physical interpretation using charge model
   - Connection to hybridization theory
   - Mode identification (dipole/quadrupole)
   - Absorption vs scattering character

7. **Conclusion**
   - Key findings summarized
   - Physical insights
   - Answers to lab questions

8. **References**
   - Your tutorial materials
   - Any additional sources

---

## What to Observe and Record

### From Basic Example:
- [ ] Two distinct peaks? (Yes - bonding & antibonding)
- [ ] Peak 1 energy: _____ eV (lower energy, bonding)
- [ ] Peak 2 energy: _____ eV (higher energy, antibonding)
- [ ] Which peak has more absorption? _____
- [ ] Which peak has more scattering? _____
- [ ] Mode type from angular scattering? (Dipole expected)

### From Thickness Study:
- [ ] As f increases, peak splitting: [ ] Increases [ ] Decreases
- [ ] Thin shell (f=0.1) splitting: _____ eV
- [ ] Thick shell (f=0.9) splitting: _____ eV
- [ ] Does this match theory prediction? [ ] Yes [ ] No
- [ ] Trend: ΔS decreases with f → ΔE decreases

### From Metal Study:
- [ ] Highest energy peaks: _____ (Al expected)
- [ ] Lowest energy peaks: _____ (K expected)
- [ ] Sharpest peaks: _____ (Ag expected - lowest γ)
- [ ] Broadest peaks: _____ (Al expected - highest γ)
- [ ] Do trends match ωₚ and γ values? [ ] Yes [ ] No

---

## Key Physics Questions to Answer

### Q1: What causes the two peaks?
**Answer using:**
- Plasmon hybridization model
- Charge oscillation picture (bonding vs antibonding)
- Mathematical: Two modes with different S values

### Q2: Why does splitting decrease with shell thickness?
**Answer using:**
- Coupling strength ∝ exp(-d/R)
- Thicker shell → surfaces further apart → weaker coupling
- Mathematical: ΔS decreases as f increases

### Q3: Why do different metals have different peak positions?
**Answer using:**
- Plasma frequency ωₚ sets energy scale
- Resonance condition: ε/ε_b ≈ -(1/S - 1)
- Drude model: ε ~ 1 - ωₚ²/ω²
- Higher ωₚ → higher resonance frequency

### Q4: Why Ag has sharpest peaks?
**Answer using:**
- Peak width related to damping γ
- Q-factor ≈ ω/γ
- Ag has lowest γ → highest Q → sharpest peaks

### Q5: Absorption vs scattering - which mode and why?
**Answer using:**
- Bonding mode: symmetric, non-radiative → more absorption
- Antibonding mode: antisymmetric, radiative → more scattering
- Size effects: C_abs ∝ V, C_sca ∝ V²

---

## Timeline

### Session 1 (Now): Run MATLAB
- [ ] Execute `run_simulations_export.m`
- [ ] Verify output files created
- [ ] Save figures
- **Time: 15 minutes**

### Session 2: Analyze in Notebook  
- [ ] Open `plasmon_analysis.ipynb`
- [ ] Load data and run cells
- [ ] Generate comparison plots
- [ ] Fill in observations
- **Time: 1-2 hours**

### Session 3: Write Report
- [ ] Structure sections
- [ ] Write background using tutorial theory
- [ ] Present results with figures
- [ ] Analyze and interpret
- [ ] Conclude
- **Time: 2-3 hours**

---

## Files You'll Have

After completion:

```
📁 Your directory/
├── 📊 MATLAB Results
│   ├── simulation_results.json (numerical data)
│   └── *.fig files (MATLAB figures)
│
├── 📓 Notebook Analysis  
│   ├── plasmon_analysis.ipynb (analysis & plots)
│   └── Figures generated in notebook
│
├── 📝 Documentation
│   ├── LAB_WORKFLOW.md (instructions)
│   ├── PHYSICS_REFERENCE.md (theory)
│   └── RESULTS_TEMPLATE.md (data recording)
│
└── 📄 Final Report
    └── lab_report.pdf (structured document)
```

---

## Ready to Start?

### ✅ Checklist:
- [ ] MATLAB installed and working
- [ ] Python environment with numpy, matplotlib, pandas
- [ ] All .m files in directory
- [ ] Notebook file created
- [ ] Tutorial formulas understood

### 🚀 Next Action:

**Start with:**
```matlab
cd '/home/fatduck/git/quantum_photonic/quantum_photonic/plasmon24u (1)'
run_simulations_export
```

Then let me know when it finishes, and we'll move to notebook analysis!

---

## Need Help?

If you encounter:
- **MATLAB errors**: Check that all functions exist (fmiecoeff, epsAu, etc.)
- **Slow execution**: Reduce number of wavelength points or N orders
- **Missing peaks**: Check energy range and peak finding threshold
- **Python import errors**: Install missing packages

I'm here to help at each step! 🎓
