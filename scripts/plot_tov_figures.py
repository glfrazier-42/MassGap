"""
TOV Section Figures for the Mass Gap Paper.

Generates:
  Figure 1 (tov_mass_radius.png)  -- M_g vs R for 6 tabulated NS EOS + QS
  Figure 2 (tov_mg_vs_mb.png)     -- M_g vs M_b for the same set

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_tov_figures.py
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_tov_figures.py --no-show
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from bhem import NeutronStarEOS, QuarkStarEOS
from bhem.tov_solver_scipy import TOVSolverScipy

# ---------- configuration ----------

NS_EOS_LIST = {
    'LS180':  {'file': 'LS180_cold_betaeq.txt',  'color': '#1f77b4'},
    'LS220':  {'file': 'LS220_cold_betaeq.txt',  'color': '#2ca02c'},
    'LS375':  {'file': 'LS375_cold_betaeq.txt',  'color': '#ff7f0e'},
    'SFHo':   {'file': 'SFHo_cold_betaeq.txt',   'color': '#d62728'},
    'DD2':    {'file': 'DD2_cold_betaeq.txt',     'color': '#9467bd'},
    'HShen':  {'file': 'HShen_cold_betaeq.txt',   'color': '#8c564b'},
}

QS_B = [60.0, 100.0, 160.0]  # MeV/fm^3
QS_COLOR = ['#000000', '#555555', '#999999']  # black, dark gray, gray

M_SUN = 1.989e33  # g

# pressure scan parameters
NS_P_MIN  = 1e34
NS_P_MAX  = 5e35
NS_NPTS   = 100

QS_P_MIN  = 1e33
QS_P_MAX  = 1e37
QS_NPTS   = 200

# ---------- helpers ----------

def solve_sequence(eos_func, p_min, p_max, n_pts):
    """Scan central pressures and return arrays of (M_g, R, M_b)."""
    solver = TOVSolverScipy(eos_func)
    pressures = np.logspace(np.log10(p_min), np.log10(p_max), n_pts)

    mg_list, r_list, mb_list = [], [], []
    for p_c in pressures:
        res = solver.solve(p_c)
        if res is None:
            continue
        mg_list.append(res.mass_solar)
        r_list.append(res.radius_km)
        mb_list.append(res.rest_mass_solar)

    mg = np.array(mg_list)
    r  = np.array(r_list)
    mb = np.array(mb_list)

    # truncate at maximum gravitational mass (stability boundary)
    if len(mg) > 0:
        i_max = np.argmax(mg)
        mg = mg[:i_max + 1]
        r  = r[:i_max + 1]
        mb = mb[:i_max + 1]

    return mg, r, mb


def run_all(base_dir):
    """Compute TOV sequences for every EOS. Returns dict of results."""
    results = {}

    # neutron star EOS (tabulated)
    for name, info in NS_EOS_LIST.items():
        eos_path = base_dir / 'data' / 'eos' / info['file']
        if not eos_path.exists():
            print("  SKIP %s -- file not found" % name)
            continue
        print("  %s ..." % name)
        ns = NeutronStarEOS(backend='tabulated', table_path=str(eos_path))
        mg, r, mb = solve_sequence(ns.eos_function(), NS_P_MIN, NS_P_MAX,
                                   NS_NPTS)
        results[name] = {'mg': mg, 'r': r, 'mb': mb,
                         'color': info['color'], 'type': 'ns'}

    # quark star (MIT bag model)
    for qs_b, qs_color in zip(QS_B, QS_COLOR):
        print("  QS (B=%g) ..." % qs_b)
        qs = QuarkStarEOS(backend='analytical', B=qs_b)
        mg, r, mb = solve_sequence(qs.eos_function(), QS_P_MIN, QS_P_MAX,
                                   QS_NPTS)
        results['QS B=%g' % qs_b] = {'mg': mg, 'r': r, 'mb': mb,
                                     'color': qs_color, 'type': 'qs'}

    return results


# ---------- plotting ----------

def plot_mass_radius(results, outpath):
    """Figure 1: gravitational mass vs radius."""
    fig, ax = plt.subplots(figsize=(4.0, 3.0))

    for name, d in results.items():
        if len(d['mg']) == 0:
            continue
        ls = '-' if d['type'] == 'ns' else '--'
        lw = 1.4 if d['type'] == 'ns' else 2.0
        ax.plot(d['r'], d['mg'], ls, color=d['color'], lw=lw, label=name)

        # mark M_max
        i_max = np.argmax(d['mg'])
        ax.plot(d['r'][i_max], d['mg'][i_max], 'o', color=d['color'],
                ms=5, zorder=5)

    ax.set_xlabel(r'Radius (km)', fontsize=9)
    ax.set_ylabel(r'Gravitational mass $M_\text{g}\;(\text{M}_\odot)$', fontsize=9)
    ax.legend(fontsize=8, ncol=1)
    ax.set_xlim(6, 22)
    ax.set_ylim(0, 3.5)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=300)
    print("  Saved %s" % outpath)
    return fig


def plot_mg_vs_mb(results, outpath):
    """Figure 2: gravitational mass vs baryonic mass."""
    fig, ax = plt.subplots(figsize=(3.0, 3.0))

    # reference line M_g = M_b
    mb_ref = np.linspace(0, 4, 200)
    ax.plot(mb_ref, mb_ref, 'k:', lw=0.8, label=r'$M_g = M_b$')

    for name, d in results.items():
        if len(d['mg']) == 0:
            continue
        ls = '-' if d['type'] == 'ns' else '--'
        lw = 1.4 if d['type'] == 'ns' else 2.0
        ax.plot(d['mb'], d['mg'], ls, color=d['color'], lw=lw, label=name)

        # mark M_max
        i_max = np.argmax(d['mg'])
        ax.plot(d['mb'][i_max], d['mg'][i_max], 'o', color=d['color'],
                ms=5, zorder=5)

    ax.set_xlabel(r'Baryonic mass $M_\text{b}\;(\text{M}_\odot)$', fontsize=9)
    ax.set_ylabel(r'Gravitational mass $M_\text{g}\;(\text{M}_\odot)$', fontsize=9)
    # ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(0, 3.5)
    ax.set_ylim(0, 3.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(str(outpath), dpi=300)
    print("  Saved %s" % outpath)
    return fig


# ---------- main ----------

def main(args):
    base_dir = Path(__file__).resolve().parent.parent
    fig_dir  = base_dir / 'latex-paper' / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("Computing TOV sequences ...")
    results = run_all(base_dir)

    # summary
    print("\nSummary (stable branch up to M_max):")
    print("  %-12s  %8s  %8s  %8s" % ("EOS", "M_max", "R(M_max)", "M_b(M_max)"))
    for name, d in results.items():
        if len(d['mg']) == 0:
            continue
        i = np.argmax(d['mg'])
        print("  %-12s  %8.3f  %8.2f  %8.3f"
              % (name, d['mg'][i], d['r'][i], d['mb'][i]))

    print("\nGenerating figures ...")
    fig1 = plot_mass_radius(results, fig_dir / 'tov_mass_radius.png')
    fig2 = plot_mg_vs_mb(results, fig_dir / 'tov_mg_vs_mb.png')

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate TOV section figures for the mass gap paper')
    parser.add_argument('--no-show', action='store_true',
                        help='Skip interactive figure display')
    args = parser.parse_args()
    main(args)
