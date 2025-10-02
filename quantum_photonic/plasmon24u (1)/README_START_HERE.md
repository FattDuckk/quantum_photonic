# Plasmonics Lab - Complete Guide Summary

## 📚 What I've Created For You

I've set up a complete framework to help you with your Plasmonics Computer Lab assignment. Here's what's available:

### 1. **LAB_WORKFLOW.md** ⭐ START HERE
- Complete step-by-step instructions
- What to run and when
- What each output means
- Questions to answer for each part

### 2. **PHYSICS_REFERENCE.md** 📖 THEORY GUIDE
- Quick physics background
- Plasmon hybridization explained
- Charge oscillation pictures
- Drude model and metal properties
- Key equations and concepts
- Common pitfalls to avoid

### 3. **RESULTS_TEMPLATE.md** 📝 RECORD YOUR DATA
- Tables pre-formatted for your results
- Spaces to record observations
- Questions with checkboxes
- Organized by experiment type
- Ready for conversion to lab report

### 4. **Example Scripts Ready to Run:**
- `example_mie_spectrum_3d_shell_eV.m` - Basic example (already existed)
- `example_shell_thickness_study.m` - Systematic shell thickness study (you modified)
- `example_metal_comparison.m` - Compare different metals (exists)

---

## 🎯 Quick Start - What To Do Now

### Step 1: Open MATLAB
Make sure you're in the correct directory:
```matlab
cd '/home/fatduck/git/quantum_photonic/quantum_photonic/plasmon24u (1)'
```

### Step 2: Run Basic Example First
```matlab
example_mie_spectrum_3d_shell_eV
```

**What you'll see:**
- Figure 1: Spectrum with 2 peaks (absorption & scattering)
- Figure 2+: Field distributions at each peak
- Angular scattering patterns

**Take time to understand these outputs** - everything builds on this!

### Step 3: Systematic Studies

**Shell thickness study:**
```matlab
example_shell_thickness_study
```
This will show you how the two peaks shift and split as shell thickness changes.

**Metal comparison:**
```matlab
example_metal_comparison
```
This will show you how different metals (Au, Ag, Al, Cu) behave differently.

### Step 4: Record Results
Use `RESULTS_TEMPLATE.md` to systematically record:
- Peak positions (energy & wavelength)
- Peak widths
- Absorption vs scattering ratios
- Your observations and interpretations

### Step 5: Write Report
Use your recorded results and the physics reference to write up your findings.

---

## 🔑 Key Concepts You Need To Understand

### The Two Peaks Come From Plasmon Hybridization

Think of it like molecular orbitals:
- **Inner surface plasmon** (cavity mode) + **Outer surface plasmon** (sphere mode)
- These couple together to form:
  - **Bonding mode** (lower energy, symmetric) → Peak 1
  - **Antibonding mode** (higher energy, antisymmetric) → Peak 2

### Shell Thickness Controls Coupling

- **Thin shell** → Strong coupling → Large splitting between peaks
- **Thick shell** → Weak coupling → Small splitting, converges to single peak
- This follows: ΔE ∝ exp(-d/R)

### Metal Type Sets Energy Scale

- Different metals have different **plasma frequencies (ωₚ)** and **damping (γ)**
- This shifts peak positions and changes peak widths
- Ag: Low loss, sharp peaks (best for sensing)
- Au: Moderate loss, visible range (good for bio applications)
- Al: High loss, UV range (cheaper but broader peaks)

---

## 📊 What The Lab Asks You To Do

### Task 1: Run and Understand ✅
- [x] Run the example file ✅ **Just run it in MATLAB**
- [x] Understand the figures ✅ **Use LAB_WORKFLOW.md guide**

### Task 2: Interpret ✅
- [ ] What's the difference between two peaks? 
  - **Answer using:** Field plots, absorption vs scattering, charge oscillation picture
- [ ] Which modes are they? 
  - **Answer using:** Angular scattering, field symmetry (dipole or quadrupole)

### Task 3: Explore Parameters ✅
- [ ] Shell thickness effect 
  - **Run:** `example_shell_thickness_study.m`
  - **Observe:** How peaks shift, how splitting changes
- [ ] Metal choice (at least 2) 
  - **Run:** `example_metal_comparison.m`
  - **Observe:** Peak positions, widths, strengths
- [ ] (Optional) Core index effect

### Task 4: Relate to Physics ✅
- [ ] Plasmon hybridization model
- [ ] Charge distribution picture
- [ ] Connection to Drude model

---

## 🔬 Understanding The Code

### What The Code Actually Does:

**You are NOT creating theory** - The code:
1. Takes experimental ε(λ) data for metals (from `epsAu.m`, etc.)
2. Solves Maxwell's equations exactly using Mie theory (`fmiecoeff.m`)
3. Calculates cross-sections from Mie coefficients (`crosssection.m`)
4. Calculates electromagnetic fields (`fmiefields.m`)
5. Plots results for you to interpret

**Your job:** Run it systematically and explain the physics!

### Key Functions:

```
fmiecoeff() ← Solves boundary conditions, returns expansion coefficients
    ↓
crosssection() ← Calculates C_abs, C_sca from coefficients
    ↓
fmiefields() ← Calculates E and H fields in space
    ↓
Your interpretation! ← Physics understanding
```

---

## 💡 Tips for Success

### When Analyzing Spectra:

1. **Count peaks** - You should see 2 for a shell
2. **Identify which is absorption vs scattering** - Check curves
3. **Note peak positions** - Record in eV and nm
4. **Compare peak widths** - Narrower = better Q-factor
5. **Look at trends** - How do peaks move with parameters?

### When Analyzing Fields:

1. **Where are hot spots?** - Inner surface? Outer? Both?
2. **What's the symmetry?** - Symmetric (bonding) or antisymmetric (antibonding)?
3. **Where does power flow?** - Red streamlines show energy
4. **Field enhancement** - How much stronger than incident field?

### When Writing Report:

1. **Don't just describe** - Explain WHY
2. **Connect to theory** - Reference plasmon hybridization, Drude model
3. **Use correct terminology** - Bonding/antibonding, dipole/quadrupole
4. **Show data systematically** - Use tables and clear figures
5. **Be quantitative** - Give actual numbers with units

---

## 🚨 Common Mistakes To Avoid

❌ "The code calculates plasmons"
✅ The code solves Maxwell equations. Plasmons emerge as resonances.

❌ "Higher energy = red-shift"
✅ Higher energy = blue-shift (shorter wavelength, higher frequency)

❌ "Thin shell = high f"
✅ Thin shell = LOW f (f is volume fraction of metal)

❌ "The two peaks are just noise"
✅ The two peaks are fundamental - they're bonding and antibonding modes!

❌ "Metal doesn't matter much"
✅ Metal choice dramatically affects peak position, width, and strength!

---

## 📋 Checklist For Lab Completion

### Experiments:
- [ ] Ran basic example and understood output
- [ ] Identified the two peaks and their character
- [ ] Ran shell thickness study (multiple f values)
- [ ] Ran metal comparison (at least Au and Ag)
- [ ] Saved all relevant figures

### Analysis:
- [ ] Recorded peak positions for all cases
- [ ] Identified mode types (dipole/quadrupole)
- [ ] Explained difference between bonding/antibonding
- [ ] Analyzed field distributions
- [ ] Understood absorption vs scattering behavior

### Physics:
- [ ] Connected to plasmon hybridization model
- [ ] Explained trends using charge oscillation picture
- [ ] Related metal differences to Drude parameters
- [ ] Understood quasi-static approximation

### Report:
- [ ] Clear introduction with physics background
- [ ] Methods section (what you ran)
- [ ] Results with figures and tables
- [ ] Discussion with physical interpretation
- [ ] Conclusions summarizing main findings

---

## 🎓 Learning Objectives

By completing this lab, you should understand:

1. **What Mie theory is** and what it calculates
2. **How nanoshells support multiple plasmon modes** (hybridization)
3. **Physical picture of bonding vs antibonding modes** (charge oscillations)
4. **How geometry affects plasmon resonances** (shell thickness)
5. **How material properties affect plasmons** (metal type)
6. **Difference between absorption and scattering** (size effects)
7. **Where electromagnetic fields are enhanced** (hot spots)

---

## 📚 Document Map

```
START HERE → LAB_WORKFLOW.md (Step-by-step instructions)
                    ↓
            Run experiments in MATLAB
                    ↓
            Record data in RESULTS_TEMPLATE.md
                    ↓
    Reference PHYSICS_REFERENCE.md for interpretation
                    ↓
            Write your lab report!
```

---

## 🆘 If You Get Stuck

### Problem: Don't understand the output
**Solution:** Read LAB_WORKFLOW.md Section 1.2 - explains what each figure shows

### Problem: Don't know what the peaks mean
**Solution:** Read PHYSICS_REFERENCE.md - explains bonding vs antibonding modes

### Problem: Don't know what to conclude
**Solution:** Use RESULTS_TEMPLATE.md questions to guide your thinking

### Problem: Code errors in MATLAB
**Solution:** 
- Check you're in the right directory
- Verify all .m files are present
- Make sure variable names match

### Problem: Physics doesn't make sense
**Solution:** 
- Start with the charge oscillation picture (PHYSICS_REFERENCE.md)
- Think about boundary conditions at interfaces
- Compare to simple sphere first (f→1 limit)

---

## 🎯 Expected Timeline

- **30 min:** Read guides, understand what you need to do
- **1 hour:** Run simulations, generate all figures
- **1 hour:** Analyze results, fill in template
- **1-2 hours:** Write up report with interpretation

**Total: 3.5-4.5 hours** for a thorough job

---

## ✨ Final Encouragement

This lab is about **understanding** not just running code!

The simulations are done FOR you (Mie theory is exact).
Your job is to **interpret the physics** and explain what's happening.

Use the guides I've created - they walk you through everything step by step.

**You've got this!** 🚀

---

## Quick Reference: File Purposes

| File | Purpose | When To Use |
|------|---------|-------------|
| `LAB_WORKFLOW.md` | Step-by-step instructions | Start here, follow along |
| `PHYSICS_REFERENCE.md` | Theory and concepts | When interpreting results |
| `RESULTS_TEMPLATE.md` | Data recording | While running experiments |
| `example_mie_spectrum_3d_shell_eV.m` | Basic demo | First run, understand basics |
| `example_shell_thickness_study.m` | Systematic thickness study | For Task 3 (geometry) |
| `example_metal_comparison.m` | Systematic metal study | For Task 3 (metal choice) |

---

**Now go run those simulations and discover some plasmon physics! 🔬✨**
