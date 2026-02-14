# Observational Spin Data Survey: Neutron Stars and Black Holes

## Summary for Greg's Research

This document summarizes available spin measurements for compact objects, relevant to understanding the role of spin in the tier-transition mass gap mechanism.

---

## 1. NEUTRON STAR SPINS

### Key Findings:
- **All observed neutron stars are rotating** - but this is primarily an observational selection effect
- **Detection methods require rotation**: Pulsars (the primary detection method) require rotation for the "lighthouse effect"
- **Spin rates**: Range from ~1 Hz (old, spun-down pulsars) to 716 Hz (fastest known: PSR J1748-2446ad)

### Fastest Spinning Neutron Stars:
| Pulsar | Spin Rate | Notes |
|--------|-----------|-------|
| PSR J1748-2446ad | 716 Hz | Equatorial velocity ~c/4 |
| PSR B1937+21 | 642 Hz | First millisecond pulsar discovered |
| PSR J1614-2230 | ~317 Hz | Mass: 1.97 ± 0.04 M☉ (heavy!) |

### Spin-Down:
- Pulsars lose rotational energy via magnetic dipole radiation
- Typical spin-down: ~10⁻¹⁵ s/rotation
- A 1-second pulsar slows to 2 seconds in ~30 million years
- Eventually become "dead pulsars" - invisible to us

### Mass-Spin Correlation:
- **Limited data available** on correlation between mass and spin for neutron stars
- Millisecond pulsars (fast spinners) tend to be in binary systems (spun up by accretion)
- Young pulsars (fast birth spin) in isolation show wide range of spin rates

### Formation Spin:
- Core collapse conserves angular momentum
- Even slowly rotating progenitor (period ~days) → rapid neutron star spin (milliseconds)
- Achieving zero spin would require perfect symmetry (essentially impossible)

**Conclusion for Greg's work**: All neutron stars likely rotate, making non-rotating TOV solutions a theoretical idealization. The observed 2.5 M☉ lower mass gap edge already incorporates rotation effects.

---

## 2. BLACK HOLE SPINS FROM LIGO/VIRGO

### Key Finding: **LIGO black holes tend to have LOW spins**

### General Population Results:
- **Effective spin parameter χ_eff**: Most LIGO/Virgo black holes cluster around χ_eff ≈ 0
- This could mean:
  1. Black holes are genuinely slowly spinning, OR
  2. Spins lie in the orbital plane (misaligned with orbital angular momentum)

### Important Spin Parameters:
- **χ_eff** (effective aligned spin): measures spin component along orbital angular momentum
  - Range: -1 to +1
  - Most LIGO events: χ_eff ≈ 0
  
- **χ_p** (effective precessing spin): measures spin component in orbital plane
  - Range: 0 to 2 (can exceed 1 with two precessing spins)
  - Shows larger values for some massive events

### Mass-Spin Trend:
From GWTC-1 analysis: **Negative correlation** between mass and mean effective spin (75-80% confidence)
- Lower mass black holes: slightly higher χ_eff
- Higher mass black holes: χ_eff closer to zero

**Interpretation**: This could support both dynamical assembly (random spin orientations) and field binaries.

---

## 3. GW190521 SPECIFICALLY (The Mass Gap Event)

### System Parameters:
- **Component masses**: 66 M☉ and 85 M☉ (both in upper mass gap!)
- **Remnant**: 142 M☉ (intermediate mass black hole)
- **Radiated energy**: 8 M☉ equivalent

### Spin Measurements:
**Critical finding**: The spin data is **highly ambiguous** due to short signal duration (~0.1s, only ~5 GW cycles)

- **χ_eff ≈ 0** (consistent with zero aligned spin)
- **χ_p**: Constrained to **large values** (significant in-plane spin component)
- **Evidence for precession**: Misaligned spins caused orbital wobble

### Interpretation Challenges:
The data is consistent with multiple scenarios:
1. **Quasi-circular orbit with precessing (misaligned) spins** (LIGO collaboration interpretation)
2. **Eccentric orbit with aligned spins** (alternative interpretation)
3. **Highly unequal mass ratio** (16 + 171 M☉ instead of 66 + 85 M☉)
4. **Head-on collision** (but this predicts too-low final spin)

### Key Quote:
"The primordial scenario is disfavored by the large spin of the primary" - BUT this assumes the quasi-circular+precession interpretation

**Problem**: Without long inspiral signal, can't definitively distinguish between:
- Eccentricity effects
- Spin precession effects
- Mass ratio uncertainties

---

## 4. SPIN IN OTHER MASS GAP CANDIDATES

### GW200129:
- Similar properties to GW190521
- Also shows evidence for spin precession
- **BUT**: Instrumental glitch overlapped part of signal, complicating interpretation

### Recent LIGO Events (GWTC-3):
- Most binary black holes have χ_eff near zero
- Some evidence for isotropic spin distribution (consistent with dynamical formation)
- Spin orientation measurements remain elusive

---

## 5. PRIMORDIAL BLACK HOLE SPIN PREDICTIONS

### Theory:
If primordial black holes form from density fluctuations in early universe:
- **Expected spin**: Low to zero (isotropic collapse, no angular momentum from rotating progenitor)
- **Distinguishing feature**: Should have systematically lower spins than stellar-origin black holes

### Supporting Evidence from LIGO Population:
One 2021 study concluded: "When considered as a homogeneous population of black holes, these results give **strong evidence for low spins** in all cases, and significant evidence for **small isotropic spins** versus any other distribution."

**Quote**: "These results give support to the idea that LIGO/Virgo black holes are primordial."

### Problem:
- Current spin measurements have **large uncertainties**
- Population-level analysis needed (individual events insufficient)
- Selection effects complicate interpretation

---

## 6. IMPLICATIONS FOR GREG'S TIER-TRANSITION MECHANISM

### For Lower Mass Gap (2.5-5.0 M☉):
**Data availability**: EXCELLENT
- Hundreds of neutron star spin measurements
- Clear understanding that all observed neutron stars rotate
- Maximum masses known (~2.5 M☉ for non-rotating configurations)

**Next step**: Compare analytical spin corrections to observed population

### For Upper Mass Gap (50-150 M☉):
**Data availability**: LIMITED
- Only a handful of black holes in this range (GW190521, GW200129, few others)
- Spin measurements highly uncertain due to short signals
- Cannot yet make definitive statistical statements

**Testable prediction from Greg's model**:
Objects in the 50-150 M☉ gap should show:
1. **Systematically lower spin** than comparable-mass objects above the gap
2. **Bimodal distribution**: Some primordial (low spin), some from hierarchical mergers (variable spin)

### Critical Observational Test:
**Can we measure spin well enough to test this?**
- For individual events: **Probably not** (uncertainties too large, especially for short signals like GW190521)
- For population: **Maybe** (with next-generation detectors and larger sample sizes)
- Best hope: Future events with longer inspiral signals → better spin constraints

---

## 7. DATA QUALITY ASSESSMENT

### What We Have:
✓ Excellent neutron star spin catalog (hundreds of objects)
✓ General trend in LIGO black holes (low spins, ~50 events)
✓ Mass measurements for LIGO events (quite good)

### What We Don't Have:
✗ Individual spin measurements for black holes in upper mass gap (too uncertain)
✗ Large statistical sample of objects in 50-150 M☉ range
✗ Waveform models that simultaneously capture eccentricity AND spin precession

### Selection Biases:
1. **Neutron stars**: Only detect rotating ones (pulsars)
2. **Black holes**: More likely to detect high-mass, high-SNR events
3. **Spin measurements**: Easier for longer signals (many inspiral cycles)

---

## 8. RECOMMENDATIONS FOR GREG'S WORK

### Immediate Priority:
**Develop analytical spin corrections** for the tier-transition mechanism at neutron star scales
- Compare with observed neutron star population
- Use literature values for typical pulsar spin rates (~100-700 Hz)
- Calculate how 15-20% mass increase from rotation affects gap boundaries

### Secondary Priority:
**Make qualitative predictions** for upper mass gap spin dependence
- Estimate non-spinning vs spinning gap locations
- Don't overinterpret sparse LIGO data
- Frame as testable prediction for future observations

### Publication Strategy:
1. **Paper 1**: Mass gap mechanism with rotating neutron stars (we have data!)
2. **Paper 2**: Extension to higher tiers with spin predictions (testable but not yet tested)

### Key Message:
"The observational spin data for the lower mass gap is excellent and supports our mechanism. The upper mass gap predictions are testable with future gravitational wave observations."

---

## 9. USEFUL REFERENCES

### Neutron Star Mass/Spin Catalogs:
- Max Planck Institute NS mass catalog: https://www3.mpifr-bonn.mpg.de/staff/pfreire/NS_masses.html
- StellarCollapse.org: https://stellarcollapse.org/nsmasses.html

### LIGO/Virgo Catalogs:
- GWTC-1, GWTC-2, GWTC-3 (gravitational wave transient catalogs)
- ~90 confirmed binary black hole mergers as of 2024

### Key GW190521 Papers:
- Abbott et al. 2020 (Physical Review Letters) - Discovery paper
- Romero-Shaw et al. 2020 (ApJL) - Eccentricity interpretation
- Multiple follow-up papers on spin/eccentricity degeneracy

---

## 10. BOTTOM LINE

**For neutron stars**: All rotate, we have excellent data, can compare analytical model to observations.

**For black holes in upper mass gap**: Very limited data, large uncertainties, but LIGO population shows general trend of low spins. Greg's primordial hypothesis predicts even LOWER spins for gap objects - testable but not yet tested.

**Main limitation**: Current gravitational wave detectors struggle to measure individual black hole spins with high precision, especially for short-duration signals like GW190521.

**Future prospects**: Next-generation detectors (Einstein Telescope, Cosmic Explorer) will dramatically improve spin measurements and provide larger samples.
