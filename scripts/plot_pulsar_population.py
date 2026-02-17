"""
Pulsar Population Survival Figures for the Mass Gap Paper.

Generates:
  Figure 1 (pulsar_population_logbins.png):
    Log-decade density bar chart with non-decreasing reference line.
    Shows pulsars vanishing at long periods.

  Figure 2 (pulsar_population_linear.png):
    Linear-bin observed counts vs steady-state dipole spin-down model.
    Growing gap between prediction and reality at long periods.

  LaTeX macros (population_macros.tex):
    Key statistics for use in the paper.

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_pulsar_population.py
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_pulsar_population.py --no-show
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats


# ============================================================
# Data loading
# ============================================================

def load_normal_pulsars(base_dir):
    """Load ATNF catalog, return periods in ms for normal pulsars (P >= 100 ms)."""
    path = base_dir / 'data' / 'atnf_pulsar_periods.csv'
    atnf = pd.read_csv(path)
    periods_ms = atnf['period'].values * 1000.0  # seconds -> ms
    total = len(periods_ms)
    normal = periods_ms[periods_ms >= 100.0]
    return normal, total, len(normal)


# ============================================================
# Figure 1: Log-decade density
# ============================================================

LOG_BIN_EDGES = [100, 300, 1000, 3000, 10000]  # ms


def compute_log_bins(periods_ms):
    """Bin normal pulsars into log-spaced bins, compute densities and deficits."""
    edges = LOG_BIN_EDGES
    n_bins = len(edges) - 1

    counts = []
    for i in range(n_bins):
        c = int(np.sum((periods_ms >= edges[i]) & (periods_ms < edges[i + 1])))
        counts.append(c)

    log_widths = [np.log10(edges[i + 1] / edges[i]) for i in range(n_bins)]
    densities = [counts[i] / log_widths[i] for i in range(n_bins)]

    # Non-decreasing reference: each value = max of all previous densities
    reference = []
    running_max = 0.0
    for d in densities:
        running_max = max(running_max, d)
        reference.append(running_max)

    # Poisson deficit test (adjacent bins)
    deficits = []
    sigmas = []
    for i in range(n_bins - 1):
        expected_min = counts[i] * (log_widths[i + 1] / log_widths[i])
        if counts[i + 1] < expected_min:
            deficit = expected_min - counts[i + 1]
            p_val = stats.poisson.cdf(counts[i + 1], expected_min)
            if 0 < p_val < 0.5:
                sigma = -stats.norm.ppf(p_val)
            elif p_val == 0:
                sigma = (expected_min - counts[i + 1]) / np.sqrt(expected_min)
            else:
                sigma = 0.0
            deficits.append(deficit)
            sigmas.append(sigma)
        else:
            deficits.append(0.0)
            sigmas.append(0.0)

    return {
        'edges': edges,
        'counts': counts,
        'log_widths': log_widths,
        'densities': densities,
        'reference': reference,
        'deficits': deficits,
        'sigmas': sigmas,
    }


def plot_log_bins(ax, info):
    """Plot log-decade density bar chart."""
    edges = info['edges']
    densities = info['densities']
    reference = info['reference']
    n_bins = len(densities)

    # Bar positions: use log-scale centers
    centers = [np.sqrt(edges[i] * edges[i + 1]) for i in range(n_bins)]
    widths = [edges[i + 1] - edges[i] for i in range(n_bins)]

    # Color bars: blue if at/above peak so far, red-orange if below
    colors = []
    for i in range(n_bins):
        if i == 0:
            colors.append('#4878CF')
        elif densities[i] < reference[i] * 0.99:
            colors.append('#C44E52')  # deficit
        else:
            colors.append('#4878CF')

    ax.bar(centers, densities, width=[w * 0.7 for w in widths],
           color=colors, edgecolor='k', linewidth=0.6, zorder=3)

    ax.set_xscale('log')
    ax.set_xlabel('Spin period (ms)')
    ax.set_ylabel('Pulsars per log-decade')
    ax.grid(True, alpha=0.2, which='both')

    # Annotate deficit bars with sigma
    for i in range(n_bins):
        if colors[i] == '#C44E52' and i > 0:
            sig = info['sigmas'][i - 1]
            if sig > 0:
                ax.text(centers[i], densities[i] + 30,
                        r'%.1f$\sigma$ deficit' % sig,
                        ha='center', va='bottom', fontsize=8, color='#C44E52')


# ============================================================
# Figure 2: Linear-bin observed vs model
# ============================================================

# LINEAR_BIN_WIDTH = 500.0   # ms
LINEAR_BIN_WIDTH = 100.0   # ms
LINEAR_BIN_LOW = 100.0     # ms
# LINEAR_BIN_HIGH = 7600.0   # ms  (15 bins)
LINEAR_BIN_HIGH = 1600.0   # ms  (15 bins)
NORM_BIN_IDX = 1           # normalise model to this bin (0-indexed; bin 1 = 200-300 ms)


def compute_linear_bins(periods_ms):
    """Bin normal pulsars into fixed-width linear bins."""
    edges = np.arange(LINEAR_BIN_LOW, LINEAR_BIN_HIGH + 1, LINEAR_BIN_WIDTH)
    n_bins = len(edges) - 1
    centers = (edges[:-1] + edges[1:]) / 2.0

    counts = np.array([
        int(np.sum((periods_ms >= edges[i]) & (periods_ms < edges[i + 1])))
        for i in range(n_bins)
    ])

    # Model: expected count ∝ P_c (dwell time ∝ P under dipole spin-down)
    model_raw = centers.copy()

    # Normalise model to match observed at one anchor bin
    scale = counts[NORM_BIN_IDX] / model_raw[NORM_BIN_IDX]
    model = model_raw * scale

    return {
        'edges': edges,
        'centers': centers,
        'counts': counts,
        'model': model,
    }


def plot_linear_bins(ax, info):
    """Plot observed counts vs dipole spin-down model."""
    centers = info['centers']
    counts = info['counts']
    model = info['model']
    width = LINEAR_BIN_WIDTH

    ax.bar(centers, counts, width=width * 0.85,
           color='#4878CF', edgecolor='k', linewidth=0.5,
           label='Observed', zorder=3, alpha=0.85)

    ax.step(np.append(info['edges'][:-1], info['edges'][-1]),
            np.append(model, model[-1]),
            where='post', color='#C44E52', lw=2.0,
            label=r'Steady-state model ($\propto P$)', zorder=4)

    ax.set_yscale('log')
    ax.set_xlabel('Spin period (ms)')
    ax.set_ylabel('Pulsars per bin')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.2, which='major')
    ax.grid(True, alpha=0.1, which='minor')
    ax.set_xlim(0, LINEAR_BIN_HIGH + 100)
    ax.set_ylim(bottom=1)


# ============================================================
# LaTeX macros
# ============================================================

def write_macros(total_all, total_normal, log_info, linear_info, outpath):
    """Write population analysis macros to LaTeX file."""
    macros = {}
    macros['popTotal'] = str(total_all)
    macros['popNormal'] = str(total_normal)

    edges = log_info['edges']
    bin_keys = ['A', 'B', 'C', 'D']
    for i, key in enumerate(bin_keys):
        macros['popBinRange%s' % key] = (
            '%d--%d\\,ms' % (edges[i], edges[i + 1])
        )
        macros['popCount%s' % key] = str(log_info['counts'][i])
        macros['popDensity%s' % key] = '%.0f' % log_info['densities'][i]

    # Deficit stats (between adjacent bins)
    trans_keys = ['AB', 'BC', 'CD']
    for j, key in enumerate(trans_keys):
        macros['popDeficit%s' % key] = '%.0f' % log_info['deficits'][j]
        macros['popSigma%s' % key] = '%.1f' % log_info['sigmas'][j]

    # Linear-bin stats
    macros['popLinearBinWidth'] = '%.0f' % LINEAR_BIN_WIDTH
    macros['popLinearNBins'] = str(len(linear_info['counts']))

    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, 'w') as f:
        f.write('%% Auto-generated by plot_pulsar_population.py -- do not edit\n')
        for name, val in macros.items():
            f.write('\\newcommand{\\%s}{%s}\n' % (name, val))

    print("Wrote %d macros to %s" % (len(macros), outpath))


# ============================================================
# Main
# ============================================================

def main(args):
    base_dir = Path(__file__).resolve().parent.parent
    fig_dir = base_dir / 'latex-paper' / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir = base_dir / 'latex-paper' / 'tables'
    table_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    periods_ms, total_all, total_normal = load_normal_pulsars(base_dir)
    print("Total pulsars in catalog: %d" % total_all)
    print("Normal pulsars (P >= 100 ms): %d" % total_normal)

    # Compute bins
    log_info = compute_log_bins(periods_ms)
    linear_info = compute_linear_bins(periods_ms)

    # Print summary
    print("\nLog-decade bins:")
    for i in range(len(log_info['counts'])):
        print("  %5d-%5d ms: count=%4d, density=%.0f pulsars/decade"
              % (log_info['edges'][i], log_info['edges'][i + 1],
                 log_info['counts'][i], log_info['densities'][i]))

    print("\nLinear bins (%.0f ms wide):" % LINEAR_BIN_WIDTH)
    for i in range(len(linear_info['counts'])):
        print("  %5.0f-%5.0f ms: observed=%4d, model=%.0f"
              % (linear_info['edges'][i], linear_info['edges'][i + 1],
                 linear_info['counts'][i], linear_info['model'][i]))

    # Figure 1: Log-decade density
    fig1, ax1 = plt.subplots(figsize=(7, 5))
    plot_log_bins(ax1, log_info)
    fig1.tight_layout()
    out1 = fig_dir / 'pulsar_population_logbins.png'
    fig1.savefig(str(out1), dpi=300)
    print("\nSaved %s" % out1)

    # Figure 2: Linear observed vs model
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    plot_linear_bins(ax2, linear_info)
    fig2.tight_layout()
    out2 = fig_dir / 'pulsar_population_linear.png'
    fig2.savefig(str(out2), dpi=300)
    print("Saved %s" % out2)

    # LaTeX macros
    write_macros(total_all, total_normal, log_info, linear_info,
                 table_dir / 'population_macros.tex')

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pulsar population survival figures for the mass gap paper')
    parser.add_argument('--no-show', action='store_true',
                        help='Skip interactive figure display')
    args = parser.parse_args()
    main(args)
