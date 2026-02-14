# CLAUDE.md

## Project Overview

This project develops a paper arguing that the "pair-instability mass gap"
in black hole masses is not a real gap, and that stellar-mass black holes
are quark stars — compact objects with real internal structure, not
singularities. The paper brings together multiple independent lines of
evidence to build a unified case.

## Core Hypothesis

- Black holes are **quark stars**: compact objects where matter has collapsed
  past the neutron degeneracy phase into quark-degenerate matter.
- They have real internal structure, an equation of state, and finite radius
  inside the event horizon.
- **M_g != M_b**: gravitational mass differs from baryonic mass because
  pressure contributes to gravity (standard GR, embedded in TOV equations).
  Spin provides centrifugal support that reduces internal pressure, thereby
  reducing M_g. No claims are made about inertial mass.
- The "pair-instability mass gap" (50-150 M_sun) does not exist. The
  drop-off at ~50 M_sun reflects the **formation boundary** (maximum mass
  of a newly formed BH), not a gap in the population.
- All BHs above ~50 M_sun are products of prior mergers, growing through
  hierarchical generations.

## Key Terminology

- **M_g**: gravitational mass (what LIGO measures, what determines orbits)
- **M_b**: baryonic mass (conserved quantity, total matter content)
- **M_i**: inertial mass (resistance to acceleration)
- **chi**: dimensionless spin parameter, chi = cJ/(GM^2), ranges 0-1
- **chi_eff**: effective spin, mass-weighted projection onto orbital axis
  (what LIGO extracts from waveforms)
- **chi_f**: remnant spin after merger (computed from NR fitting formulas,
  not directly measured)
- **q**: mass ratio m2/m1 (m1 >= m2 by convention), ranges 0-1
- **eta**: symmetric mass ratio, eta = m1*m2/(m1+m2)^2
- **alpha**: [currently unused; reserved for future work]

## Lines of Evidence

1. **TOV equation analysis**: Neutron star → quark star transition via
   Tolman-Oppenheimer-Volkoff equations
2. **Pulsar population analysis**: Most massive neutron stars have high
   spin rates, consistent with spin delaying collapse
3. **Red giant problem**: [connection to be developed]
4. **BH merger spin analysis (GWTC catalog)**:
   - Gaussian remnant spin distribution at chi_f ~ 0.69 (structural limit)
   - Spin efficiency: lopsided mergers deposit more spin per unit m2
   - Excess deficit proportional to spin (spin reduces M_g)
   - Discrete iso-m2 tracks in mass deficit (generational mergers)
   - 50 M_sun formation boundary
   - No mass-to-energy conversion; deficit is spin-induced M_g reduction

## Physical Units

- **Mass**: M_sun (1.989 x 10^30 kg)
- **Distance**: kpc where relevant
- **G = 4.302 x 10^-6 kpc^3 M_sun^-1 Myr^-2** (galactic scale)
- **G = 6.674 x 10^-11 m^3 kg^-1 s^-2** (SI)

## Important Physics

- **Structural spin limit**: When spin exceeds a critical rate, the quark
  star expands (centrifugal force), moment of inertia increases, spin slows.
  Self-regulating feedback loop — a thermostat, not a breaking point.
- **No mass-to-energy conversion**: GW energy comes from rotational kinetic
  energy. The mass deficit is entirely spin-induced M_g reduction. Baryonic
  mass is conserved in mergers.
- **Parsimony**: One GW radiation mechanism throughout inspiral-merger-ringdown
  (quadrupole radiation from physical mass distribution). No switch to
  "spacetime geometry radiating itself."
- **Not opposed to GR**: Opposes singularities, not general relativity.
  Uses GR everywhere observationally validated.
- **Circularity problem**: All LIGO parameters come from Kerr GR template
  fitting. Using GR-derived parameters to test non-GR hypotheses is circular.

## Commands

```bash
# Activate virtual environment (create first if needed)
python -m venv venv
source venv/bin/activate  # or: venv/Scripts/activate on Windows

# Run scripts
PYTHONPATH=src python scripts/<script_name>.py
# On Windows/Cygwin:
PYTHONPATH=src venv/Scripts/python.exe scripts/<script_name>.py

# Build LaTeX paper
cd latex-paper && make

# Dependencies
pip install numpy scipy matplotlib pandas
```

## Data Sources

- `data/gwtc_catalog.csv` — GWTC gravitational wave transient catalog
- Additional data TBD (pulsar catalogs, TOV inputs)

## Project Structure

```
docs/           # Analysis records in markdown
latex-paper/    # LaTeX source and Makefile
src/            # Reusable analysis code
scripts/        # One-off analysis scripts
results/        # Generated figures and tables
data/           # Downloaded observational data
references/     # Downloaded reference papers
```
