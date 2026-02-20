# data/raw/

This directory holds raw data files that must be downloaded manually before
building the paper. Its contents are excluded from the repository by
`.gitignore`.

## ATNF Pulsar Catalogue

**File:** `psrcat.db`
**Source:** https://www.atnf.csiro.au/research/pulsar/psrcat/
**Instructions:** Download the catalogue database file and place it here.

Used by:
- `scripts/extract_atnf_periods.py` → `data/atnf_pulsar_periods.csv`
- `scripts/extract_atnf_name_mapping.py` → `data/atnf_name_mapping.csv`

## Neutron Star EOS Tables (stellarcollapse.org)

**License:** CC BY-NC-SA — do not redistribute.
**Source:** https://stellarcollapse.org/equationofstate.html

Download the following six HDF5 files and place them in this directory:

| File | EOS | Reference |
|------|-----|-----------|
| `LS180_234r_136t_50y_analmu_20091212_SVNr26.h5` | Lattimer-Swesty K=180 MeV | Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991) |
| `LS220_234r_136t_50y_analmu_20091212_SVNr26.h5` | Lattimer-Swesty K=220 MeV | Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991) |
| `LS375_234r_136t_50y_analmu_20091212_SVNr26.h5` | Lattimer-Swesty K=375 MeV | Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991) |
| `Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5` | SFHo | Steiner et al., ApJ 774, 17 (2013) |
| `Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5` | DD2 | Hempel & Schaffner-Bielich, Nucl. Phys. A 837, 210 (2010) |
| `HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5` | HShen | Shen et al., ApJS 197, 20 (2011) |

After downloading, run:

```bash
PYTHONPATH=src venv/Scripts/python.exe scripts/extract_all_eos.py
```

This extracts cold beta-equilibrium slices and writes text tables to
`data/eos/`, which is what the TOV solver reads.
