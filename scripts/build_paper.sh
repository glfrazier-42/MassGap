#!/usr/bin/env bash
# build_paper.sh -- Full build from a clean checkout.
#
# MANUAL PREREQUISITES (anti-robot protections prevent automation):
#
#   1. ATNF Pulsar Catalogue database
#      Download psrcat.db from:
#        https://www.atnf.csiro.au/research/pulsar/psrcat/
#      Place at: data/raw/psrcat.db
#
#   2. EOS text tables
#      Download HDF5 tables from stellarcollapse.org, extract cold
#      beta-equilibrium slices with scripts/extract_all_eos.py, and
#      place the resulting .txt files in data/eos/.
#      Required files:
#        data/eos/LS180_cold_betaeq.txt
#        data/eos/LS220_cold_betaeq.txt
#        data/eos/LS375_cold_betaeq.txt
#        data/eos/SFHo_cold_betaeq.txt
#        data/eos/DD2_cold_betaeq.txt
#        data/eos/HShen_cold_betaeq.txt
#
#   3. GWTC catalog
#      Download from the Gravitational Wave Open Science Center:
#        https://gwosc.org
#      Place at: data/gwtc_catalog.csv
#
# data/stellar_collapse_ns_masses.csv is committed to the repository
# and does not need to be downloaded.
#
# USAGE:
#   bash scripts/build_paper.sh
#
# Requires Python venv to be set up first:
#   python -m venv venv
#   venv/Scripts/python.exe -m pip install -r requirements.txt

set -e

DIR="$(dirname "$0")"
cd "${DIR}/.."

PYTHON=venv/Scripts/python.exe

# ---------------------------------------------------------------------------
# Step 1: Extract ATNF pulsar data from psrcat.db
# ---------------------------------------------------------------------------
echo "=== Step 1: Extracting ATNF pulsar data ==="
PYTHONPATH=src $PYTHON scripts/extract_atnf_periods.py
PYTHONPATH=src $PYTHON scripts/extract_atnf_name_mapping.py

# ---------------------------------------------------------------------------
# Step 2: Merge pulsar mass and period data
# ---------------------------------------------------------------------------
echo "=== Step 2: Merging pulsar mass and period data ==="
PYTHONPATH=src $PYTHON scripts/merge_mass_period_final.py

# ---------------------------------------------------------------------------
# Step 3: Generate figures and LaTeX macros
# ---------------------------------------------------------------------------
echo "=== Step 3: Generating figures and macros ==="
PYTHONPATH=src $PYTHON scripts/plot_pulsar_masses.py --no-show
PYTHONPATH=src $PYTHON scripts/plot_tov_figures.py --no-show
PYTHONPATH=src $PYTHON scripts/plot_pulsar_population.py --no-show
PYTHONPATH=src $PYTHON scripts/plot_bh_merger_spin.py --no-show
PYTHONPATH=src $PYTHON scripts/investigate_eta_confound.py

# ---------------------------------------------------------------------------
# Step 4: Build the PDF and copy to project root
# ---------------------------------------------------------------------------
echo "=== Step 4: Building PDF ==="
cd latex-paper
make && cp mass_gap.pdf ..
