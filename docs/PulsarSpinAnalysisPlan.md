# Pulsar Spin-Mass Correlation Analysis - Development Plan

## Overview

This document outlines the development plan for analyzing the correlation between pulsar masses and spin rates to empirically determine the maximum non-rotating neutron star mass (M_max).

## Scientific Objective

**Key Observation:** The four most massive pulsars are also the four fastest-spinning pulsars.

**Hypothesis:** When massive pulsars slow down due to spin-down, they lose centrifugal support and collapse into black holes. This explains why we don't observe slowly-spinning massive pulsars.

**Goal:** Use Hartle-Thorne equations to calculate the equivalent non-rotating mass M_0 from observed (M_obs, P_obs), then extract empirical M_max from the correlation between mass and minimum spin rate P_min(M).

## Data Sources

### Johnston 2023 MeerKAT Census
**Location:** `data/raw/Johnston2023/tables_files/Census_Table1.csv`

**Provides:**
- Pulsar names (PSRJ)
- Spin periods (P) in seconds
- Spin-down energy (logEdot)
- Large sample size (~1000 pulsars)

**Limitation:** No mass measurements (requires binary systems)

### ATNF Pulsar Catalog
**Source:** https://www.atnf.csiro.au/research/pulsar/psrcat/

**Query Parameters:**
- `PSRJ` - Pulsar name
- `P0` - Spin period (seconds)
- `PMASS` - **Pulsar mass (M_sun)** ← Critical field!
- `PMASS_ERR` - Mass measurement error
- `Binary` - Binary companion type (for context)

**Filter:** `PMASS > 0` (only pulsars with mass measurements)

**Expected:** ~50-100 pulsars with mass measurements (subset with sufficient post-Keplerian parameters)

### Known Massive Pulsars (Validation Set)
- J0740+6620: M = 2.08 ± 0.07 M_sun, P = 2.89 ms
- J0348+0432: M = 2.01 ± 0.04 M_sun, P = 39.3 ms
- J1614-2230: M = 1.97 ± 0.04 M_sun, P = 3.15 ms
- J0437-4715: M = 1.44 ± 0.07 M_sun, P = 5.75 ms

## Physics Background

### Hartle-Thorne Formalism

For slowly rotating stars (valid for all observed pulsars):

**Mass correction to first order in Ω:**
```
M(Ω) = M_0 + ΔM_rot(Ω)

where:
ΔM_rot ∝ Ω² × I(M_0)
Ω = 2π/P (angular velocity)
I = moment of inertia
```

**Key insight:** Rotation provides centrifugal support, allowing stars to exceed the non-rotating M_max.

### Spin-Down and Collapse Mechanism

**Observation:** Pulsars spin down via magnetic dipole radiation
- τ_spindown ~ 10^6 - 10^7 years
- P increases over time
- Ω decreases → ΔM_rot decreases

**Critical point:** When M(Ω) drops below M_max for a massive pulsar:
- Star becomes unstable
- Collapses to black hole
- This is why we don't see slowly-spinning massive pulsars!

### Empirical M_max Extraction

**Method:**
1. Bin pulsars by observed mass M_obs
2. For each bin, find P_min(M) - the slowest observed spin
3. Calculate M_0(M_obs, P_min) using Hartle-Thorne
4. Plot M_0 vs M_obs
5. M_0 should converge to M_max for the most massive pulsars

**Expected result:** M_max ~ 2.0-2.3 M_sun (from TOV theory)

**Significance:** If M_max ~ 2.0-2.3 M_sun but no black holes exist below 5 M_sun → **mass amplification required** → supports tier transition mechanism

## Architecture Design

### Code Structure
```
src/bhem/
├── __init__.py              # Package exports
├── eos.py                   # Equation of state (existing)
├── tov_solver.py            # Non-rotating stars (existing)
├── tov_solver_scipy.py      # Alternative solver (existing)
├── tov_tabulated.py         # Tabulated EOS support (existing)
├── hartle_thorne.py         # NEW: Rotating star solver
└── pulsar_data.py           # NEW: Data management

src/
└── analyze_spin_mass_correlation.py  # NEW: Main analysis script

data/
├── raw/
│   ├── Johnston2023/        # Existing: spin period data
│   └── ATNF/                # To be added: mass measurements
└── processed/
    └── pulsar_catalog.csv   # Combined dataset
```

### EOS Compatibility

**Design principle:** Hartle-Thorne solver uses **identical EOS interface** as TOV solver:

```python
def eos(P: float) -> Tuple[float, float]:
    """
    Parameters:
        P: Pressure (dyne/cm^2)
    Returns:
        (rho, epsilon): density and energy density (g/cm^3)
    """
```

**Benefit:** Drop-in compatibility with all existing EOS implementations:
- Analytical polytropes (validated, stiff)
- Tabulated EOS (LS220, APR, SLy4)
- Hybrid neutron/quark models
- Future extensions

## Development Phases

### Phase 1: Hartle-Thorne Solver Implementation

**Goal:** Create rotating star solver with validation

#### Step 1.1: Class Structure and Non-Rotating Baseline
**File:** `src/bhem/hartle_thorne.py`

**Components:**
```python
class HartleThorneSOlver:
    def __init__(self, eos: Callable):
        """Initialize with EOS (same interface as TOVSolver)"""
        
    def solve(self, P_central: float, Omega: float = 0):
        """
        Solve for rotating star configuration
        When Omega=0, must match TOVSolver exactly
        """
```

**Validation:**
- Compare HT(Ω=0) vs TOV for multiple EOS
- Should be **numerically identical** (within integration tolerance)
- Test with: validated neutron, stiff, quark EOS
- Check: M_grav, M_rest, R, pressure profiles

**Success criteria:** 
- Relative difference < 10^-6 for all quantities
- Same profiles r(P), m(r), etc.

#### Step 1.2: First-Order Rotation Corrections
**Implement:**
```python
    def solve_rotating(self, P_central: float, Omega: float):
        """
        Add rotation perturbations:
        - Mass correction ΔM(Ω)
        - Radius correction ΔR(Ω)
        - Moment of inertia I(Ω)
        - Frame dragging ω(r)
        """
```

**Physics to implement:**
- Metric perturbations h_0, h_2, m_0, m_2, k_2, v_2
- Coupled differential equations (6 ODEs + 2 from TOV)
- Boundary conditions at center and surface

**Validation:**
- Test against published M(Ω) curves from literature
- Verify Kepler limit: Ω_max where dM/dΩ → 0
- Check moment of inertia: compare to empirical I-M relations
- Test limiting cases:
  - Low Ω: should approach TOV + quadratic correction
  - High Ω: should show flattening near Kepler limit

**Success criteria:**
- Reproduces published results for canonical 1.4 M_sun NS
- Kepler frequency within 5% of literature values
- Moment of inertia consistent with empirical relations

#### Step 1.3: Inverse Calculation (Critical Function!)
**Implement:**
```python
    def compute_non_rotating_mass(self, M_obs: float, Omega: float, 
                                   tolerance: float = 1e-6):
        """
        Given observed M and Ω, find M_0 such that:
        M(M_0, Ω) = M_obs
        
        Uses root finding on:
        f(M_0) = solve_rotating(M_0, Ω).mass - M_obs = 0
        """
```

**Algorithm:**
- Bracket search to find P_central(M_0)
- For each M_0, solve rotating configuration
- Use bisection/Brent to find M_0 where M(Ω) = M_obs

**Validation - Round-Trip Test:**
```python
# Start with known M_0
M_0_input = 2.0  # M_sun
Omega = 1000     # rad/s

# Forward: compute rotating mass
result = ht_solver.solve_rotating(P_c_for_M0, Omega)
M_obs = result.mass_solar

# Inverse: recover M_0
M_0_recovered = ht_solver.compute_non_rotating_mass(M_obs, Omega)

# Check
assert abs(M_0_recovered - M_0_input) < 1e-4
```

**Success criteria:**
- Round-trip error < 0.01% for range of M_0 and Ω
- Robust convergence (< 20 iterations)
- Handles edge cases (very slow, very fast rotation)

### Phase 2: Data Management

**Goal:** Load, validate, and cross-match pulsar data

#### Step 2.1: Data Structures
**File:** `src/bhem/pulsar_data.py`

**Classes:**
```python
@dataclass
class Pulsar:
    """Single pulsar observation"""
    name: str                    # PSRJ designation
    period: float                # Spin period (seconds)
    period_err: Optional[float]  # Period uncertainty
    mass: Optional[float]        # Pulsar mass (M_sun)
    mass_err: Optional[float]   # Mass uncertainty
    binary_type: Optional[str]   # Companion type
    source: str                  # Data source (Johnston/ATNF/manual)
    
    @property
    def omega(self) -> float:
        """Angular velocity (rad/s)"""
        return 2 * np.pi / self.period
    
    @property
    def has_mass(self) -> bool:
        """Whether mass measurement exists"""
        return self.mass is not None
```

**Validation:**
- Test with manually created instances
- Verify property calculations
- Check error handling for missing data

#### Step 2.2: Johnston 2023 Loader
**Implement:**
```python
class PulsarCatalog:
    def __init__(self):
        self.pulsars: List[Pulsar] = []
    
    def load_johnston2023(self, data_path: Path):
        """
        Load Johnston 2023 Census_Table1.csv
        Provides: PSRJ, P (period)
        Missing: mass measurements
        """
        df = pd.read_csv(data_path / "Census_Table1.csv")
        for _, row in df.iterrows():
            self.pulsars.append(Pulsar(
                name=row['PSRJ'],
                period=row['P'],
                mass=None,  # Not in Johnston 2023
                source='Johnston2023'
            ))
```

**Validation:**
- Load test data
- Check: number of pulsars matches file
- Verify: period values match known pulsars
- Test: proper handling of missing/NaN values

#### Step 2.3: ATNF Catalog Loader
**Implement:**
```python
    def load_atnf(self, data_path: Path):
        """
        Load ATNF catalog
        Provides: PSRJ, P0, PMASS, PMASS_ERR
        Only includes pulsars with mass measurements
        """
```

**Cross-matching logic:**
```python
    def merge_catalogs(self):
        """
        Cross-match Johnston and ATNF by PSRJ name
        Updates Johnston entries with ATNF masses where available
        """
```

**Validation:**
- Compare periods: Johnston vs ATNF (should agree)
- Check known pulsars: J0740+6620, J0348+0432, etc.
- Verify: mass uncertainties are reasonable
- Test: handling of name variations (J2000 format)

#### Step 2.4: Data Quality Filtering
**Implement:**
```python
    def get_pulsars_with_mass(self, 
                              min_mass_snr: float = 5.0,
                              max_period_err: float = 1e-6):
        """
        Filter for high-quality measurements:
        - Mass SNR > threshold (PMASS/PMASS_ERR)
        - Period uncertainty acceptable
        """
```

**Validation:**
- Check filtered sample size
- Verify: only high-quality data retained
- Test: boundary cases (exactly at threshold)

### Phase 3: Analysis Pipeline

**Goal:** Extract empirical M_max from spin-mass correlation

#### Step 3.1: Basic Statistical Analysis
**File:** `src/analyze_spin_mass_correlation.py`

**Implement:**
```python
def analyze_spin_distribution(catalog: PulsarCatalog):
    """
    Basic statistics:
    - Mass distribution (histogram)
    - Period distribution  
    - M vs P scatter plot
    - Identify outliers
    """
```

**Validation:**
- Reproduce known correlation: massive pulsars spin faster
- Check: J0740+6620 is fastest among massive
- Verify: distribution shapes make physical sense

#### Step 3.2: P_min(M) Determination
**Implement:**
```python
def compute_minimum_periods(catalog: PulsarCatalog, 
                           mass_bins: np.ndarray):
    """
    For each mass bin:
    - Find all pulsars in range [M - δM, M + δM]
    - Determine P_min - slowest observed period
    - Calculate statistical uncertainty
    
    Returns: (mass_centers, P_min_values, uncertainties)
    """
```

**Statistical considerations:**
- Bin width: balance resolution vs sample size
- Edge effects: handle bins with few pulsars
- Uncertainty: bootstrap or analytical propagation

**Validation:**
- Test with synthetic data (known P_min(M))
- Check: P_min decreases with increasing M (expected trend)
- Verify: uncertainties scale with sqrt(N_bin)

#### Step 3.3: M_0 Calculation for Each Pulsar
**Implement:**
```python
def calculate_nonrotating_masses(catalog: PulsarCatalog,
                                 ht_solver: HartleThorneSOlver):
    """
    For each pulsar with mass measurement:
    - M_obs, P_obs → Ω_obs
    - Use HT solver: M_0 = f^(-1)(M_obs, Ω_obs)
    - Propagate uncertainties
    
    Returns: DataFrame with M_obs, P, M_0, errors
    """
```

**Uncertainty propagation:**
```python
δM_0 = sqrt[(∂M_0/∂M_obs)² δM_obs² + (∂M_0/∂P)² δP²]
```

**Validation:**
- Check: M_0 < M_obs for all pulsars (rotation adds mass)
- Verify: difference ΔM = M_obs - M_0 increases with faster spin
- Test: error propagation gives reasonable uncertainties

#### Step 3.4: Empirical M_max Extraction
**Implement:**
```python
def extract_empirical_mmax(results: pd.DataFrame):
    """
    Determine M_max from M_0 values:
    
    Method 1: Maximum M_0 among massive pulsars
    Method 2: Asymptotic limit of M_0(M_obs) curve
    Method 3: Weighted average of top N pulsars
    
    Returns: M_max, uncertainty, methodology
    """
```

**Analysis:**
- Plot M_0 vs M_obs - should show convergence
- Statistical tests for plateau/asymptote
- Compare methods: do they agree?

**Validation:**
- Check: M_max ~ 2.0-2.3 M_sun (theory prediction)
- Verify: consistent across different EOS
- Test: sensitivity to outliers

#### Step 3.5: Publication Plots
**Implement:**
```python
def generate_publication_figures():
    """
    Figure 1: M vs P scatter plot with P_min(M) curve
    Figure 2: M_0 vs M_obs showing empirical M_max
    Figure 3: Multi-EOS comparison
    Figure 4: Uncertainty analysis
    """
```

**Style:**
- Publication quality (600 DPI)
- Clear labels, legends
- Error bars visible
- Consistent color scheme

## Validation Strategy

### Level 1: Unit Tests (Each Component)
- Test individual functions in isolation
- Check boundary conditions
- Verify error handling

### Level 2: Integration Tests (Phase Transitions)
- Phase 1 → 2: Does HT solver work with PulsarCatalog?
- Phase 2 → 3: Can analysis pipeline use both?

### Level 3: Physics Validation (Cross-Checks)
- HT(Ω=0) = TOV (exact agreement)
- Round-trip M_0 → M(Ω) → M_0 (< 0.01% error)
- Literature comparison (published M-R curves)
- Known pulsars (J0740+6620 values)

### Level 4: Statistical Validation
- Bootstrap analysis (resample pulsars)
- Cross-validation (leave-one-out)
- Systematic error estimation

### Level 5: Multi-EOS Validation
- Run full pipeline with different EOS
- Check: M_max varies reasonably
- Verify: all give mass amplification conclusion

## Expected Outcomes

### Quantitative Results
1. **Empirical M_max:** 2.0-2.3 M_sun (from pulsar observations)
2. **Spin-mass correlation:** P_min(M) ~ M^(-α), α > 0
3. **Mass amplification:** Gap between M_max and lowest BH (5 M_sun)

### Qualitative Insights
1. **Direct observation:** Massive pulsars collapse when slowing down
2. **Model-independent:** No assumptions about tier mechanism
3. **Statistical significance:** Can be quantified rigorously

### Publication Material
1. **Main result:** Empirical M_max contradicts BH mass gap
2. **Supporting evidence:** Spin-down → collapse mechanism
3. **Figures:** Clear visualization of correlation

## Timeline Estimate

**Phase 1:** Hartle-Thorne Solver
- Step 1.1: 2-3 days (structure + validation)
- Step 1.2: 3-4 days (rotation physics)
- Step 1.3: 2-3 days (inverse calculation)
- **Total:** ~1 week

**Phase 2:** Data Management  
- Step 2.1: 1 day (structures)
- Step 2.2: 1 day (Johnston loader)
- Step 2.3: 1-2 days (ATNF + cross-match)
- Step 2.4: 1 day (filtering)
- **Total:** ~1 week

**Phase 3:** Analysis
- Steps 3.1-3.2: 2 days (statistics)
- Step 3.3: 2-3 days (M_0 calculations)
- Step 3.4: 2-3 days (M_max extraction)
- Step 3.5: 2 days (plots)
- **Total:** ~1.5 weeks

**Grand Total:** ~3.5 weeks for complete implementation + validation

*Note: Timeline assumes incremental development with validation at each step*

## Risk Mitigation

### Technical Risks
**Risk:** HT solver convergence issues
**Mitigation:** Robust root-finding, adaptive step size

**Risk:** Data quality (sparse mass measurements)
**Mitigation:** Multiple validation datasets, quality filters

**Risk:** EOS-dependent results
**Mitigation:** Multi-EOS analysis, sensitivity studies

### Scientific Risks
**Risk:** No clear P_min(M) correlation
**Mitigation:** This would be scientifically interesting! (Why not?)

**Risk:** M_max >> 2.5 M_sun (no gap)
**Mitigation:** Update theory, examine assumptions

**Risk:** Insufficient statistics
**Mitigation:** Bootstrap analysis, conservative error bars

## Success Metrics

### Must Have (Minimum Viable Product)
- ✅ HT solver matches TOV when Ω=0
- ✅ Can calculate M_0 from (M_obs, P_obs)
- ✅ Load pulsar data from ATNF
- ✅ Extract empirical M_max

### Should Have (Full Implementation)
- ✅ Multi-EOS support
- ✅ Uncertainty propagation
- ✅ Publication-quality plots
- ✅ Statistical validation

### Nice to Have (Future Extensions)
- ⭐ Web interface for data exploration
- ⭐ Real-time ATNF catalog updates
- ⭐ Machine learning P_min(M) fitting
- ⭐ 3D visualizations

## References

### Hartle-Thorne Formalism
- Hartle, J. B., & Thorne, K. S. (1968). ApJ, 153, 807
- Hartle, J. B. (1967). ApJ, 150, 1005

### Pulsar Data Sources
- Johnston, S., et al. (2023). MeerKAT Census (arXiv:2305.xxxxx)
- Manchester, R. N., et al. (2005). ATNF Pulsar Catalogue

### Massive Pulsars
- Cromartie, H. T., et al. (2020). J0740+6620 mass (Nature Astronomy)
- Antoniadis, J., et al. (2013). J0348+0432 mass (Science)
- Demorest, P. B., et al. (2010). J1614-2230 mass (Nature)

### Neutron Star Physics
- Lattimer, J. M., & Prakash, M. (2016). Physics Reports
- Özel, F., & Freire, P. (2016). ARA&A

---

**Document Status:** Initial version, to be updated as implementation proceeds
**Last Updated:** 2026-01-18
**Author:** Greg Franklin (with Claude)
