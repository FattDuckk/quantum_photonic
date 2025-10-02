# Initial Observations from MATLAB Results

## ✅ Simulations Completed Successfully!

All three experiments ran and data has been exported to `simulation_results.json`

---

## 🔍 Quick Preview of Results

### Experiment 1: Basic Example (K shell)
**Parameters:** K metal, a=20nm, f=0.3, nc=1.5, nb=1.5

**Results:**
- **Peak 1 (bonding):** 0.717 eV (1730 nm) - **Near-IR range**
- **Peak 2 (antibonding):** 3.543 eV (350 nm) - **UV range**
- **Energy splitting:** 2.826 eV - **VERY LARGE!**

**Initial interpretation:**
- K has low plasma frequency (~3.72 eV from tutorial)
- Bonding mode deep in near-infrared
- Large splitting indicates strong coupling (moderate shell thickness)
- This matches tutorial prediction: K → IR/NIR peaks

---

### Experiment 2: Shell Thickness Study (Au shell)
**Fixed:** Au metal, a=10nm, nc=1, nb=1

| f | Shell thickness (nm) | E₁ (eV) | E₂ (eV) | ΔE (eV) | Observation |
|---|---------------------|---------|---------|---------|-------------|
| 0.1 | 0.35 | 1.194 | 1.194 | 0.000 | **Same peak?** |
| 0.3 | 1.12 | 1.861 | 1.861 | 0.000 | **Same peak?** |
| 0.5 | 2.06 | 2.164 | 2.153 | 0.011 | Small splitting |
| 0.7 | 3.31 | 2.321 | 2.299 | 0.022 | Splitting grows |
| 0.9 | 5.36 | 2.416 | 5.578 | 3.162 | **HUGE splitting!** |

**🤔 Interesting observation:**
- This trend is **OPPOSITE** to what hybridization theory predicts!
- Theory says: thin shell → large splitting
- Results show: thick shell → large splitting

**Possible explanations:**
1. **Peak finding issue** - Maybe not finding true peaks for thin shells
2. **Different mode character** - Thin shells might have overlapping modes
3. **Retardation effects** - Particle size not truly quasistatic
4. **Need to look at full spectrum** - Absorption vs scattering peaks separately

**📌 This is interesting for discussion!** Theory vs reality discrepancy.

---

### Experiment 3: Metal Comparison (f=0.5, a=10nm)

| Metal | E₁ (eV) | λ₁ (nm) | E₂ (eV) | λ₂ (nm) | ΔE (eV) | Range |
|-------|---------|---------|---------|---------|---------|-------|
| **Au** | 2.164 | 573 | 2.153 | 576 | 0.011 | Yellow-Green |
| **Ag** | 2.762 | 449 | 2.762 | 449 | 0.000 | Blue |
| **Al** | 5.421 | 229 | 5.451 | 228 | 0.030 | **UV** |
| **Cu** | 2.102 | 590 | 2.072 | 598 | 0.030 | Orange-Red |

**✅ This matches tutorial predictions!**

**Plasma frequency ranking (from tutorial):**
- Al (ω_p ≈ 15 eV) → **Highest energy peaks** ✅ Confirmed!
- Au, Ag, Cu (ω_p ≈ 9 eV) → **Similar energies** ✅ Confirmed!

**Specific observations:**
1. **Al:** 5.4 eV = 229 nm (UV range) - Correct! High ω_p
2. **Ag:** 2.76 eV = 449 nm (blue) - Higher than Au (expected for noble metals)
3. **Au:** 2.16 eV = 573 nm (yellow-green) - Classic gold plasmon range
4. **Cu:** 2.10 eV = 590 nm (orange) - Similar to Au but slightly red

**Peak splitting:**
- Very small for Au, Ag, Cu (~0-0.03 eV)
- This suggests at f=0.5, modes are nearly degenerate
- Or absorption and scattering peaks are very close

---

## 📊 What to Do Next

### 1. Load data in notebook ✅
```python
import json
with open('simulation_results.json', 'r') as f:
    data = json.load(f)
```

### 2. Generate detailed plots
- Spectrum plots for each case
- Thickness study trends
- Metal comparison overlays
- Compare with theoretical predictions

### 3. Investigate anomalies
- **Why does thickness trend reverse?**
  - Look at absorption vs scattering separately
  - Check if we're finding the right peaks
  - Maybe need to look at full spectrum visually

### 4. Compare with theory
- Calculate theoretical S values (already in notebook)
- Predict resonance energies
- Discuss where theory works and where it breaks down

### 5. Physical interpretation
- Use charge oscillation picture
- Explain metal trends (ω_p and γ)
- Discuss quasistatic approximation validity

---

## 🔬 Key Questions to Answer in Analysis

### Q1: Basic Example (K shell)
- **Why such large splitting (2.8 eV)?** 
  - K has lower ω_p → lower energy scale
  - Moderate thickness → significant coupling
  - Both factors contribute

### Q2: Thickness Study Anomaly
- **Why does splitting increase with f instead of decrease?**
  - Need to investigate this carefully
  - Could be fundamental physics or peak-finding artifact
  - **This is important for the report discussion!**

### Q3: Metal Comparison Success
- **Why does this match theory well?**
  - Direct relationship to ω_p
  - Quasistatic prediction: ω_res ∝ ω_p/√(1+2ε_b)
  - All metals at same geometry → fair comparison

---

## 📝 Notes for Report

### Strengths:
✅ Metal comparison trends match theory perfectly  
✅ K shell in IR/NIR as expected (low ω_p)  
✅ Al in UV as expected (high ω_p)  
✅ Au/Ag/Cu cluster together (similar ω_p)

### Anomalies to Discuss:
⚠️ Thickness study shows opposite trend from simple hybridization model  
⚠️ Very small splitting for f=0.1-0.5 in Au shell  
⚠️ Need to examine this more carefully

### Analysis Approach:
1. Present results honestly (including surprises)
2. Compare with theoretical predictions
3. Discuss why theory might differ from Mie simulation
4. Consider limitations of quasistatic approximation
5. Investigate whether we're finding correct peaks

---

## 🎯 Next Action

**Open `plasmon_analysis.ipynb` and run the cells!**

The notebook will:
1. Load the JSON data
2. Create detailed visualizations
3. Calculate theoretical predictions
4. Compare theory vs simulation
5. Help you interpret the physics

**The anomaly in thickness study is actually INTERESTING for your report** - it shows the difference between simple theory and exact solution!

Let me know when you're ready to dive into the notebook analysis! 📊
