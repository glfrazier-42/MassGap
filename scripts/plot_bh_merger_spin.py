"""
BH Merger Spin Analysis Figures for the Mass Gap Paper.

Generates four figures for comparison (user picks which go in the paper):

  Figure 1a (bh_merger_excess_vs_spin.png):
    Excess fractional deficit vs |chi_eff|, with QS model curves.

  Figure 1b (bh_merger_excess_vs_chif.png):
    Excess fractional deficit vs computed chi_f, with alpha*chi_f^2 fit.
    This is the physically cleaner x-axis (remnant spin, not progenitor spin).

  Figure 2a (bh_merger_deficit_per_m2.png):
    Deficit per unit m2 vs q, with GR and QS prediction curves.

  Figure 2b (bh_merger_excess_per_m2.png):
    Excess deficit per unit m2 vs q, with QS prediction curve.
    Spearman here isolates the spin-efficiency effect from geometric factors.

  LaTeX macros (merger_macros.tex):
    Key statistics for use in the paper.

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_bh_merger_spin.py
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_bh_merger_spin.py --no-show
"""

import argparse
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr


# ============================================================
# Data loading
# ============================================================

def safe_float(val):
    """Convert to float, return NaN for empty/invalid."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def load_gwtc(filepath):
    """Load GWTC catalog, return list of event dicts."""
    events = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            events.append(row)
    return events


def filter_bbh(events):
    """Filter for confident BBH events."""
    bbh = []
    for e in events:
        m1 = safe_float(e.get('mass_1_source', ''))
        m2 = safe_float(e.get('mass_2_source', ''))
        chi = safe_float(e.get('chi_eff', ''))
        p_astro = safe_float(e.get('p_astro', ''))
        m_final = safe_float(e.get('final_mass_source', ''))
        m_total = safe_float(e.get('total_mass_source', ''))

        if np.isnan(m1) or np.isnan(m2) or m1 < 3 or m2 < 3:
            continue
        if np.isnan(p_astro) or p_astro < 0.5:
            continue

        name = e.get('commonName', e.get('id', ''))
        bbh.append({
            'name': name,
            'm1': m1, 'm2': m2,
            'm_total': m_total if not np.isnan(m_total) else m1 + m2,
            'm_final': m_final,
            'chi_eff': chi,
        })

    return bbh


# ============================================================
# Physics computations
# ============================================================

def gr_nonspinning_frac(eta):
    """GR non-spinning fractional mass radiated (Buonanno, Kidder & Lehner 2008)."""
    return 0.0572 * eta + 0.498 * eta**2


def L_orb(eta):
    """Orbital angular momentum contribution to remnant spin (NR fit)."""
    return 2 * np.sqrt(3) * eta - 3.5171 * eta**2 + 2.5763 * eta**3


def compute_chi_f(q, eta, chi_eff):
    """Compute remnant spin from NR fitting formula."""
    l = L_orb(eta)
    S_over_M2 = chi_eff * (1 + q**2) / (1 + q)**2
    chi_f = l + S_over_M2
    return np.clip(chi_f, 0.0, 1.0)


def compute_event_quantities(bbh):
    """Compute derived quantities for BBH events with valid final mass + chi_eff."""
    results = []
    for e in bbh:
        if np.isnan(e['m_final']) or e['m_final'] <= 0:
            continue
        if np.isnan(e['chi_eff']):
            continue

        q = e['m2'] / e['m1']
        eta = e['m1'] * e['m2'] / e['m_total']**2
        frac_obs = (e['m_total'] - e['m_final']) / e['m_total']
        frac_gr = gr_nonspinning_frac(eta)
        excess = frac_obs - frac_gr
        deficit_abs = e['m_total'] - e['m_final']
        deficit_per_m2 = deficit_abs / e['m2']
        gr_deficit_per_m2 = (frac_gr * e['m_total']) / e['m2']
        chi_f = compute_chi_f(q, eta, e['chi_eff'])

        results.append({
            'name': e['name'],
            'm1': e['m1'], 'm2': e['m2'], 'm_total': e['m_total'],
            'm_final': e['m_final'],
            'q': q, 'eta': eta,
            'chi_eff': e['chi_eff'],
            'chi_f': chi_f,
            'frac_obs': frac_obs,
            'frac_gr': frac_gr,
            'excess': excess,
            'deficit_abs': deficit_abs,
            'deficit_per_m2': deficit_per_m2,
            'gr_deficit_per_m2': gr_deficit_per_m2,
            'excess_per_m2': deficit_per_m2 - gr_deficit_per_m2,
        })

    return results


def fit_qs_alpha(evts):
    """Fit the one-parameter QS model: excess = alpha * chi_f^2.

    Least squares through origin: alpha = sum(y*x) / sum(x^2)
    where y = excess, x = chi_f^2.
    """
    chi_f2 = np.array([e['chi_f']**2 for e in evts])
    excess = np.array([e['excess'] for e in evts])
    alpha = np.dot(excess, chi_f2) / np.dot(chi_f2, chi_f2)
    return alpha


# ============================================================
# Formatting helpers
# ============================================================

def fmt_p(p):
    r"""Format p-value as LaTeX scientific notation, e.g. 1.5 \times 10^{-11}."""
    if p < 0.001:
        coeff, exponent = ('%.1e' % p).split('e')
        return r'%s \times 10^{%d}' % (coeff, int(exponent))
    return '%.4f' % p


def spearman_annotation(rho, p):
    """Return formatted annotation string for Spearman correlation."""
    return (r'Spearman $\rho = %+.3f$' '\n' r'$p = %s$') % (rho, fmt_p(p))


def m1_color(m):
    """Color a point by primary mass bin."""
    if m >= 50:
        return '#d62728'    # red
    elif m >= 35:
        return '#ff7f0e'    # orange
    elif m >= 20:
        return '#1f77b4'    # blue
    return '#2ca02c'        # green


def add_m1_legend(ax):
    """Add mass-bin color legend to axis."""
    for label, color in [('$m_1 < 20$', '#2ca02c'),
                          ('$20 \\leq m_1 < 35$', '#1f77b4'),
                          ('$35 \\leq m_1 < 50$', '#ff7f0e'),
                          ('$m_1 \\geq 50$', '#d62728')]:
        ax.scatter([], [], c=color, s=30, edgecolors='k', linewidths=0.3,
                   label=label)


Q_BIN_EDGES = [0.15, 0.35, 0.50, 0.60, 0.70, 0.80, 0.95]


def binned_medians(x, y, bin_edges, min_count=3):
    """Compute binned medians. Returns (centers, medians) arrays."""
    centers, medians = [], []
    for i in range(len(bin_edges) - 1):
        mask = (x >= bin_edges[i]) & (x < bin_edges[i + 1])
        if np.sum(mask) >= min_count:
            centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            medians.append(np.median(y[mask]))
    return np.array(centers), np.array(medians)


# ============================================================
# Figure 1a: Excess vs |chi_eff| with QS model curves
# ============================================================

def plot_fig1a(ax, evts, alpha_qs):
    chi_abs = np.array([abs(e['chi_eff']) for e in evts])
    excess = np.array([e['excess'] for e in evts])
    m1 = np.array([e['m1'] for e in evts])
    colors = [m1_color(m) for m in m1]

    ax.scatter(chi_abs, excess, c=colors, s=22, alpha=0.7,
               edgecolors='k', linewidths=0.3, zorder=5)
    ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=2)

    # Binned medians
    bc, bm = binned_medians(chi_abs, excess,
                            np.linspace(0, chi_abs.max() * 1.01, 9))
    ax.plot(bc, bm, 's-', color='#d62728', linewidth=2, markersize=7,
            markeredgecolor='k', markeredgewidth=0.5, zorder=10,
            label='Binned median')

    # QS model prediction curves for different q
    chi_eff_fine = np.linspace(0, 0.65, 100)
    for q_val, ls, lbl in [(1.0, '-', '$q = 1.0$'),
                            (0.7, '--', '$q = 0.7$'),
                            (0.4, ':', '$q = 0.4$')]:
        eta_val = q_val / (1 + q_val)**2
        l_orb_val = L_orb(eta_val)
        s_factor = (1 + q_val**2) / (1 + q_val)**2
        chi_f_curve = np.clip(l_orb_val + chi_eff_fine * s_factor, 0, 1)
        excess_curve = alpha_qs * chi_f_curve**2
        ax.plot(chi_eff_fine, excess_curve, ls, color='#9467bd',
                linewidth=2, label='QS model, %s' % lbl, zorder=7)

    rho, p = spearmanr(chi_abs, excess)
    ax.text(0.97, 0.97, spearman_annotation(rho, p),
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

    add_m1_legend(ax)
    ax.legend(fontsize=7, loc='upper left')
    ax.set_xlabel(r'$|\chi_{\rm eff}|$', fontsize=11)
    ax.set_ylabel('Excess fractional deficit over non-spinning GR', fontsize=10)
    ax.grid(True, alpha=0.2)
    return rho, p


# ============================================================
# Figure 1b: Excess vs chi_f with QS model parabola
# ============================================================

def plot_fig1b(ax, evts, alpha_qs):
    chi_f = np.array([e['chi_f'] for e in evts])
    excess = np.array([e['excess'] for e in evts])
    m1 = np.array([e['m1'] for e in evts])
    colors = [m1_color(m) for m in m1]

    ax.scatter(chi_f, excess, c=colors, s=22, alpha=0.7,
               edgecolors='k', linewidths=0.3, zorder=5)
    ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=2)

    # Binned medians
    bc, bm = binned_medians(chi_f, excess,
                            np.linspace(chi_f.min() * 0.95, chi_f.max() * 1.02, 9))
    ax.plot(bc, bm, 's-', color='#d62728', linewidth=2, markersize=7,
            markeredgecolor='k', markeredgewidth=0.5, zorder=10,
            label='Binned median')

    rho, p = spearmanr(chi_f, excess)
    ax.text(0.97, 0.03, spearman_annotation(rho, p),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

    add_m1_legend(ax)
    ax.legend(fontsize=7, loc='upper left')
    ax.set_xlabel(r'Remnant spin $\chi_f$ (computed from NR fit)', fontsize=11)
    ax.set_ylabel('Excess fractional deficit over non-spinning GR', fontsize=10)
    ax.grid(True, alpha=0.2)
    return rho, p


# ============================================================
# Figure 2a: Deficit/m2 vs q with GR + QS curves
# ============================================================

def plot_fig2a(ax, evts, alpha_qs):
    q_arr = np.array([e['q'] for e in evts])
    deficit_pm2 = np.array([e['deficit_per_m2'] for e in evts])

    ax.scatter(q_arr, deficit_pm2, c='#1f77b4', s=22, alpha=0.6,
               edgecolors='k', linewidths=0.3, zorder=5, label='Observed')

    # Theory curves
    q_fine = np.linspace(0.10, 1.0, 200)
    eta_fine = q_fine / (1 + q_fine)**2
    frac_gr = gr_nonspinning_frac(eta_fine)
    l_orb_fine = L_orb(eta_fine)

    # GR non-spinning
    gr_dm2 = frac_gr * (1 + q_fine) / q_fine
    ax.plot(q_fine, gr_dm2, 'r--', linewidth=2.5,
            label='Singularity model (non-spinning GR)', zorder=8)

    # QS model: deficit includes spin-induced M_g reduction
    frac_qs = frac_gr + alpha_qs * l_orb_fine**2
    qs_dm2 = frac_qs * (1 + q_fine) / q_fine
    ax.plot(q_fine, qs_dm2, '-', color='#9467bd', linewidth=2.5,
            label='Quark star model', zorder=8)

    # Binned medians
    bc, bm = binned_medians(q_arr, deficit_pm2, Q_BIN_EDGES)
    ax.plot(bc, bm, 's-', color='#2ca02c', linewidth=2, markersize=7,
            markeredgecolor='k', markeredgewidth=0.5, zorder=10,
            label='Binned median')

    rho, p = spearmanr(q_arr, deficit_pm2)
    ax.text(0.97, 0.97, spearman_annotation(rho, p),
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

    ax.legend(fontsize=7.5, loc='upper left')
    ax.set_xlabel(r'Mass ratio $q = m_2/m_1$', fontsize=11)
    ax.set_ylabel(r'Mass deficit per unit $m_2$ ($M_\odot / M_\odot$)',
                  fontsize=10)
    ax.grid(True, alpha=0.2)
    return rho, p


# ============================================================
# Figure 2b: Excess deficit/m2 vs q with QS curve
# ============================================================

def plot_fig2b(ax, evts, alpha_qs):
    q_arr = np.array([e['q'] for e in evts])
    excess_pm2 = np.array([e['excess_per_m2'] for e in evts])

    ax.scatter(q_arr, excess_pm2, c='#1f77b4', s=22, alpha=0.6,
               edgecolors='k', linewidths=0.3, zorder=5, label='Observed')
    ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=2)

    # QS prediction: excess_per_m2 = alpha * L_orb(eta)^2 * (1+q)/q
    q_fine = np.linspace(0.10, 1.0, 200)
    eta_fine = q_fine / (1 + q_fine)**2
    l_orb_fine = L_orb(eta_fine)
    qs_excess_dm2 = alpha_qs * l_orb_fine**2 * (1 + q_fine) / q_fine
    ax.plot(q_fine, qs_excess_dm2, '-', color='#9467bd', linewidth=2.5,
            label='Quark star model', zorder=8)

    # Binned medians
    bc, bm = binned_medians(q_arr, excess_pm2, Q_BIN_EDGES)
    ax.plot(bc, bm, 's-', color='#2ca02c', linewidth=2, markersize=7,
            markeredgecolor='k', markeredgewidth=0.5, zorder=10,
            label='Binned median')

    rho, p = spearmanr(q_arr, excess_pm2)
    ax.text(0.97, 0.97, spearman_annotation(rho, p),
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlabel(r'Mass ratio $q = m_2/m_1$', fontsize=11)
    ax.set_ylabel(r'Excess deficit per unit $m_2$ (over non-spinning GR)',
                  fontsize=10)
    ax.grid(True, alpha=0.2)
    return rho, p


# ============================================================
# LaTeX macros
# ============================================================

def write_macros(bbh, evts, stats, alpha_qs, outpath):
    """Write merger analysis macros to LaTeX file."""
    macros = {}

    macros['gwNbbh'] = str(len(bbh))
    macros['gwNdeficit'] = str(len(evts))

    excess_arr = np.array([e['excess'] for e in evts])
    macros['gwMedianExcess'] = '%.4f' % np.median(excess_arr)
    frac_pos = np.sum(excess_arr > 0) / len(excess_arr)
    macros['gwFracPositive'] = '%.0f\\%%' % (frac_pos * 100)

    frac_deficit = np.array([e['frac_obs'] for e in evts])
    macros['gwMedianFracDeficit'] = '%.3f' % np.median(frac_deficit)

    macros['gwAlphaQS'] = '%.4f' % alpha_qs

    # All Spearman results
    for key, (rho, p) in stats.items():
        macros['gwRho%s' % key] = '%+.3f' % rho
        macros['gwP%s' % key] = fmt_p(p)

    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, 'w') as f:
        f.write('%% Auto-generated by plot_bh_merger_spin.py -- do not edit\n')
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

    # Load and filter data
    data_path = base_dir / 'data' / 'gwtc_catalog.csv'
    events = load_gwtc(str(data_path))
    bbh = filter_bbh(events)
    print("Confident BBH events (m1>3, m2>3, p_astro>0.5): %d" % len(bbh))

    evts = compute_event_quantities(bbh)
    print("Events with final mass and chi_eff: %d" % len(evts))

    if len(evts) < 5:
        print("ERROR: Too few events to analyse. Check data file.")
        return

    # Fit QS model
    alpha_qs = fit_qs_alpha(evts)
    print("\nQS model fit: excess = %.4f * chi_f^2" % alpha_qs)

    excess_arr = np.array([e['excess'] for e in evts])
    frac_pos = np.sum(excess_arr > 0) / len(excess_arr)
    print("Median excess: %+.4f" % np.median(excess_arr))
    print("Fraction positive: %.0f%%" % (frac_pos * 100))

    # Generate all four figures
    stats = {}

    fig1a, ax1a = plt.subplots(figsize=(7, 5))
    rho, p = plot_fig1a(ax1a, evts, alpha_qs)
    stats['ExcessChi'] = (rho, p)
    fig1a.tight_layout()
    fig1a.savefig(str(fig_dir / 'bh_merger_excess_vs_spin.png'), dpi=300)
    print("\nFig 1a saved (excess vs |chi_eff|)")

    fig1b, ax1b = plt.subplots(figsize=(7, 5))
    rho, p = plot_fig1b(ax1b, evts, alpha_qs)
    stats['ExcessChif'] = (rho, p)
    fig1b.tight_layout()
    fig1b.savefig(str(fig_dir / 'bh_merger_excess_vs_chif.png'), dpi=300)
    print("Fig 1b saved (excess vs chi_f)")

    fig2a, ax2a = plt.subplots(figsize=(7, 5))
    rho, p = plot_fig2a(ax2a, evts, alpha_qs)
    stats['DeficitQ'] = (rho, p)
    fig2a.tight_layout()
    fig2a.savefig(str(fig_dir / 'bh_merger_deficit_per_m2.png'), dpi=300)
    print("Fig 2a saved (deficit/m2 vs q)")

    fig2b, ax2b = plt.subplots(figsize=(7, 5))
    rho, p = plot_fig2b(ax2b, evts, alpha_qs)
    stats['ExcessPerMtwoQ'] = (rho, p)
    fig2b.tight_layout()
    fig2b.savefig(str(fig_dir / 'bh_merger_excess_per_m2.png'), dpi=300)
    print("Fig 2b saved (excess deficit/m2 vs q)")

    # Also compute Spearman for excess_frac vs q (cleanest geometric test)
    q_arr = np.array([e['q'] for e in evts])
    rho_eq, p_eq = spearmanr(q_arr, excess_arr)
    stats['ExcessQ'] = (rho_eq, p_eq)

    # Macros
    write_macros(bbh, evts, stats, alpha_qs,
                 table_dir / 'merger_macros.tex')

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("  BBH events:              %d" % len(bbh))
    print("  With deficit + chi_eff:  %d" % len(evts))
    print("  QS model alpha:          %.4f" % alpha_qs)
    print()
    for key, (rho, p) in stats.items():
        print("  %-20s  rho=%+.3f  p=%s" % (key, rho, fmt_p(p)))

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='BH merger spin analysis figures for the mass gap paper')
    parser.add_argument('--no-show', action='store_true',
                        help='Skip interactive figure display')
    args = parser.parse_args()
    main(args)
