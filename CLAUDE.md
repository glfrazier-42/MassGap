# CLAUDE.md

## Project Overview

This project develops a paper arguing that the lower mass gap
(~2.5-5 M_sun) between neutron stars and black holes is the natural
signature of baryonic matter undergoing a phase transition from neutron
degeneracy to quark degeneracy. Black holes are quark stars — compact
objects with real internal structure, not singularities. The paper
brings together multiple independent lines of evidence to build a
unified case.

## Core Hypothesis

- Black holes are **quark stars**: compact objects where matter has collapsed
  past the neutron degeneracy phase into quark-degenerate matter.
- They have real internal structure, an equation of state, and finite radius
  inside the event horizon.
- **M_g != M_b**: gravitational mass differs from baryonic mass because
  pressure contributes to gravity (standard GR, embedded in TOV equations).
  Spin provides centrifugal support that reduces internal pressure, thereby
  reducing M_g. No claims are made about inertial mass.
- The **lower mass gap** (~2.5-5 M_sun) is the phase transition
  signature: matter collapses past neutron degeneracy into quark
  degeneracy, producing a jump in compactness and gravitational mass.

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

## Workflow

- All work is done on the **dev** branch. Periodically merge dev into
  master and push master to GitHub.
- **Never commit other people's work** into the repository. This includes
  third-party papers, datasets with restrictive licenses, etc. The
  `references/` folder will eventually need to be cleared of any such
  files.
- **Full reproducibility** is the goal. A third party should be able to
  clone the repository, download data from public sources, re-run the
  analysis scripts, regenerate all tables and figures, and produce the
  PDF. To that end:
  - `latex-paper/figures/` and `latex-paper/tables/` are in `.gitignore`;
    all tables and figures are built from raw data by scripts.
  - No quantitative results are hard-coded in `mass_gap.tex`.
  - The README (or scripts where feasible) will document how to manually
    download and install data so the paper can be built. Some data sources
    have anti-robot protections requiring manual download.

## Working Directory and Paths

The project root is always the current working directory. Use relative
paths everywhere — in all file tool calls (Read, Write, Edit, Glob).
Never construct or use absolute paths.

## Shell Access

Claude runs without access to a Cygwin shell (no mintty). The Bash tool
is unreliable in this environment. Do not attempt to run commands via the
Bash tool. Instead:
- Use Read, Write, Edit, Glob for all file operations.
- When a command must be run (tests, build, scripts), output it as a
  instruction for the user to execute in their terminal.
- For directory creation, write a file into the target directory; the
  Write tool creates parent directories automatically.

## Commands

- Do not execute 'cd <dir>' before invoking python. Claude will always be
  launched in the MassGap directory.
- Python should always be invoked "PYTHONPATH=src venv/Scripts/python ..."

```bash
# Activate virtual environment (create first if needed)
python -m venv venv
source venv/Scripts/activate  # or: venv/Scripts/activate on Windows

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
