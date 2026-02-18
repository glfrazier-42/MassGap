#!/usr/bin/env bash

DIR=`dirname $0`
cd ${DIR}/..

PYTHONPATH=src python scripts/plot_pulsar_masses.py --no-show
PYTHONPATH=src python scripts/plot_tov_figures.py --now-show
PYTHONPATH=src python scripts/plot_pulsar_population.py --no-show
PYTHONPATH=src python scripts/plot_bh_merger_spin.py --no-show

cd latex-paper
# force the rebuild
touch mass_gap.tex
make
