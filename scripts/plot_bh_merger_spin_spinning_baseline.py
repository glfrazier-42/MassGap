"""BH Merger Spin Analysis: Spinning GR Baseline Comparison.

Tests whether the excess mass deficit observed in GWTC binary black hole
mergers survives after accounting for the spin-dependent GR prediction.

The original analysis (plot_bh_merger_spin.py) compares observed deficit
against the non-spinning GR formula (Buonanno, Kidder & Lehner 2008):

    f_GR_ns = 0.0572 * eta + 0.498 * eta^2

A referee objection: aligned-spin progenitors should radiate MORE energy
even in standard GR (spin-orbit coupling shrinks the effective ISCO),
so the "excess" over the non-spinning baseline might simply be the
known spinning-GR effect rather than a quark-star signature.

This script computes both baselines:
    f_GR_ns  -- non-spinning (existing formula)
    f_GR_sp  -- spinning extension via Kerr ISCO binding energy (new)

and repeats the two paper tests against both, allowing direct comparison.

Spinning baseline (BKL Kerr-ISCO extension):
    f_GR_sp(eta, chi_eff) = e_isco_kerr(chi_eff) * eta + 0.498 * eta^2

where e_isco_kerr(a) = 1 - E_ISCO(a)/mu is the analytic Kerr ISCO
binding energy (Bardeen, Press & Teukolsky 1972).  At chi_eff = 0 this
reduces identically to f_GR_ns.

Output figures (latex-paper/figures/):
    bh_merger_spinning_test1.png  -- excess vs chi_f, both baselines
    bh_merger_spinning_test2.png  -- deficit/m2 vs q, three model curves

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_bh_merger_spin_spinning_baseline.py
    PYTHONPATH=src venv/Scripts/python.exe scripts/plot_bh_merger_spin_spinning_baseline.py --no-show
"""

import argparse
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import spearmanr

from bhem.gr_baselines import f_gr_nonspin, f_gr_spin


# ---------------------------------------------------------------------------
# Data loading  (identical filter to plot_bh_merger_spin.py)
# ---------------------------------------------------------------------------

def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def load_gwtc(filepath):
    events = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            events.append(row)
    return events


def filter_bbh(events):
    bbh = []
    for e in events:
        m1      = safe_float(e.get('mass_1_source', ''))
        m2      = safe_float(e.get('mass_2_source', ''))
        chi     = safe_float(e.get('chi_eff', ''))
        p_astro = safe_float(e.get('p_astro', ''))
        m_final = safe_float(e.get('final_mass_source', ''))
        m_total = safe_float(e.get('total_mass_source', ''))

        if np.isnan(m1) or np.isnan(m2) or m1 < 3 or m2 < 3:
            continue
        if np.isnan(p_astro) or p_astro < 0.5:
            continue

        name = e.get('commonName', e.get('id', ''))
        bbh.append({
            'name':    name,
            'm1':      m1,
            'm2':      m2,
            'm_total': m_total if not np.isnan(m_total) else m1 + m2,
            'm_final': m_final,
            'chi_eff': chi,
        })
    return bbh


# ---------------------------------------------------------------------------
# Physics  (NR fitting formula for chi_f matches plot_bh_merger_spin.py)
# ---------------------------------------------------------------------------

def L_orb(eta):
    """Orbital angular momentum contribution to remnant spin."""
    return 2 * np.sqrt(3) * eta - 3.5171 * eta ** 2 + 2.5763 * eta ** 3


def compute_chi_f(q, eta, chi_eff):
    l = L_orb(eta)
    S = chi_eff * (1 + q ** 2) / (1 + q) ** 2
    return np.clip(l + S, 0.0, 1.0)


def fit_alpha(excess_arr, chi_f_arr):
    """Least-squares fit: excess = alpha * chi_f^2 (through origin)."""
    x = chi_f_arr ** 2
    return np.dot(excess_arr, x) / np.dot(x, x)


# ---------------------------------------------------------------------------
# Event computation
# ---------------------------------------------------------------------------

def compute_events(bbh):
    """Return per-event dict with both baselines for events that have
    final_mass and chi_eff."""
    results = []
    for e in bbh:
        if np.isnan(e['m_final']) or e['m_final'] <= 0:
            continue
        if np.isnan(e['chi_eff']):
            continue

        q        = e['m2'] / e['m1']
        eta      = e['m1'] * e['m2'] / e['m_total'] ** 2
        frac_obs = (e['m_total'] - e['m_final']) / e['m_total']
        chi_f    = compute_chi_f(q, eta, e['chi_eff'])

        f_ns = f_gr_nonspin(eta)
        f_sp = f_gr_spin(eta, e['chi_eff'])

        results.append({
            'name':         e['name'],
            'm1':           e['m1'],
            'm2':           e['m2'],
            'm_total':      e['m_total'],
            'm_final':      e['m_final'],
            'q':            q,
            'eta':          eta,
            'chi_eff':      e['chi_eff'],
            'chi_f':        chi_f,
            'frac_obs':     frac_obs,
            'f_gr_ns':      f_ns,
            'f_gr_sp':      f_sp,
            'excess_ns':    frac_obs - f_ns,   # vs non-spinning GR
            'excess_sp':    frac_obs - f_sp,   # vs spinning GR
            'deficit_pm2':  (e['m_total'] - e['m_final']) / e['m2'],
            'gr_ns_pm2':    f_ns * e['m_total'] / e['m2'],
            'gr_sp_pm2':    f_sp * e['m_total'] / e['m2'],
            'spin_corr':    f_sp - f_ns,       # magnitude of spin correction
        })
    return results


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def fmt_p(p):
    """Format p-value for use inside a matplotlib mathtext $...$ block."""
    if p < 0.001:
        coeff, exp = ('%.1e' % p).split('e')
        return r'%s \times 10^{%d}' % (coeff, int(exp))
    return '%.4f' % p


def spearman_box(ax, rho, p, loc='upper right'):
    corners = {
        'upper right': (0.97, 0.97, 'right', 'top'),
        'upper left':  (0.03, 0.97, 'left',  'top'),
        'lower right': (0.97, 0.03, 'right', 'bottom'),
    }
    x, y, ha, va = corners[loc]
    txt = r'Spearman $\rho = %+.3f$' '\n' r'$p = %s$' % (rho, fmt_p(p))
    ax.text(x, y, txt, transform=ax.transAxes, ha=ha, va=va, fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))


def binned_medians(x, y, bin_edges, min_count=3):
    centers, medians = [], []
    for i in range(len(bin_edges) - 1):
        mask = (x >= bin_edges[i]) & (x < bin_edges[i + 1])
        if np.sum(mask) >= min_count:
            centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            medians.append(np.median(y[mask]))
    return np.array(centers), np.array(medians)


# ---------------------------------------------------------------------------
# Figure 1: Test 1 side-by-side — excess vs chi_f for both baselines
# ---------------------------------------------------------------------------

def plot_test1(evts, alpha_ns, alpha_sp, fig_dir):
    chi_f      = np.array([e['chi_f']    for e in evts])
    excess_ns  = np.array([e['excess_ns'] for e in evts])
    excess_sp  = np.array([e['excess_sp'] for e in evts])
    m1         = np.array([e['m1']        for e in evts])

    rho_ns, p_ns = spearmanr(chi_f, excess_ns)
    rho_sp, p_sp = spearmanr(chi_f, excess_sp)

    chif_fine = np.linspace(0.55, 0.95, 200)
    qs_ns_curve = alpha_ns * chif_fine ** 2
    qs_sp_curve = alpha_sp * chif_fine ** 2

    fig, (ax_ns, ax_sp) = plt.subplots(1, 2, figsize=(13, 5),
                                        sharey=False)
    fig.suptitle('Test 1: Excess Deficit vs Remnant Spin\n'
                 '(left: vs non-spinning GR baseline; '
                 'right: vs spinning GR baseline)',
                 fontsize=11, fontweight='bold')

    for ax, excess, rho, p, alpha, label, color in [
        (ax_ns, excess_ns, rho_ns, p_ns, alpha_ns,
         'Non-spinning GR baseline', '#d62728'),
        (ax_sp, excess_sp, rho_sp, p_sp, alpha_sp,
         'Spinning GR baseline (Kerr ISCO)', '#e377c2'),
    ]:
        sc = ax.scatter(chi_f, excess, c=m1, cmap='viridis',
                        s=22, alpha=0.7, edgecolors='k', linewidths=0.3,
                        zorder=5, vmin=5, vmax=80)
        ax.axhline(0, color='gray', linewidth=0.8, linestyle='--', zorder=2)

        # Binned medians
        bc, bm = binned_medians(chi_f, excess,
                                np.linspace(chi_f.min() * 0.95,
                                            chi_f.max() * 1.02, 9))
        ax.plot(bc, bm, 's-', color='#d62728', linewidth=2, markersize=7,
                markeredgecolor='k', markeredgewidth=0.5, zorder=10,
                label='Binned median')

        # QS model parabola fitted on this baseline
        qs_curve = alpha * chif_fine ** 2
        ax.plot(chif_fine, qs_curve, '-', color='#9467bd', linewidth=2,
                label=r'QS model ($\alpha \chi_f^2$, $\alpha=%.4f$)' % alpha,
                zorder=7)

        spearman_box(ax, rho, p, loc='upper left')
        ax.legend(fontsize=7.5, loc='lower right')
        ax.set_xlabel(r'Remnant spin $\chi_f$', fontsize=11)
        ax.set_ylabel('Excess fractional deficit', fontsize=10)
        ax.set_title(label, fontsize=10)
        ax.grid(True, alpha=0.2)

        plt.colorbar(sc, ax=ax, label='Primary mass $m_1$ ($M_\odot$)',
                     pad=0.02).ax.tick_params(labelsize=7)

    fig.tight_layout()
    outpath = fig_dir / 'bh_merger_spinning_test1.png'
    fig.savefig(str(outpath), dpi=300)
    print("Saved %s" % outpath)
    return rho_ns, p_ns, rho_sp, p_sp


# ---------------------------------------------------------------------------
# Figure 2: Test 2 — deficit/m2 vs q with three model curves
# ---------------------------------------------------------------------------

def plot_test2(evts, alpha_ns, alpha_sp, chi_eff_med, fig_dir):
    q_arr      = np.array([e['q']          for e in evts])
    def_pm2    = np.array([e['deficit_pm2'] for e in evts])

    rho_ns, p_ns = spearmanr(q_arr, def_pm2)   # same as original Test 2

    q_fine   = np.linspace(0.10, 1.0, 300)
    eta_fine = q_fine / (1 + q_fine) ** 2

    # Non-spinning GR curve
    gr_ns_dm2 = f_gr_nonspin(eta_fine) * (1 + q_fine) / q_fine

    # Spinning GR curve (evaluated at sample median chi_eff)
    gr_sp_dm2 = f_gr_spin(eta_fine, chi_eff_med) * (1 + q_fine) / q_fine

    # QS model curves (fitted on respective baselines)
    l_orb_fine = L_orb(eta_fine)
    qs_ns_dm2 = (f_gr_nonspin(eta_fine) + alpha_ns * l_orb_fine ** 2) \
                * (1 + q_fine) / q_fine
    qs_sp_dm2 = (f_gr_spin(eta_fine, chi_eff_med) + alpha_sp * l_orb_fine ** 2) \
                * (1 + q_fine) / q_fine

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle('Test 2: Mass Deficit per Unit $m_2$ vs Mass Ratio\n'
                 '(three model curves: non-spinning GR, spinning GR, '
                 'quark star)',
                 fontsize=10, fontweight='bold')

    ax.scatter(q_arr, def_pm2, c='#1f77b4', s=22, alpha=0.6,
               edgecolors='k', linewidths=0.3, zorder=5, label='Observed')

    ax.plot(q_fine, gr_ns_dm2, 'r--', linewidth=2,
            label='Singularity model, non-spinning GR', zorder=7)
    ax.plot(q_fine, gr_sp_dm2, 'b-.', linewidth=2,
            label=r'Singularity model, spinning GR ($\bar\chi_{\rm eff}=%.2f$)'
                  % chi_eff_med,
            zorder=7)
    ax.plot(q_fine, qs_ns_dm2, '-', color='#9467bd', linewidth=2,
            label='Quark star model (fit vs non-spinning GR)', zorder=8)
    ax.plot(q_fine, qs_sp_dm2, color='#17becf', linewidth=1.5,
            linestyle=(0, (3, 1)),
            label='Quark star model (fit vs spinning GR)', zorder=8)

    # Binned medians
    Q_BINS = [0.15, 0.35, 0.50, 0.60, 0.70, 0.80, 0.95]
    bc, bm = binned_medians(q_arr, def_pm2, Q_BINS)
    ax.plot(bc, bm, 's-', color='#2ca02c', linewidth=2, markersize=7,
            markeredgecolor='k', markeredgewidth=0.5, zorder=10,
            label='Binned median')

    spearman_box(ax, rho_ns, p_ns, loc='upper right')
    ax.legend(fontsize=7.5, loc='upper left')
    ax.set_xlabel(r'Mass ratio $q = m_2/m_1$', fontsize=11)
    ax.set_ylabel(r'Mass deficit per unit $m_2$ ($M_\odot/M_\odot$)',
                  fontsize=10)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()

    outpath = fig_dir / 'bh_merger_spinning_test2.png'
    fig.savefig(str(outpath), dpi=300)
    print("Saved %s" % outpath)
    return rho_ns, p_ns


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args):
    matplotlib.use('Agg' if args.no_show else 'TkAgg')

    base_dir = Path(__file__).resolve().parent.parent
    fig_dir  = base_dir / 'latex-paper' / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)

    data_path = base_dir / 'data' / 'gwtc_catalog.csv'
    events = load_gwtc(str(data_path))
    bbh    = filter_bbh(events)
    print("Confident BBH events (m1>3, m2>3, p_astro>0.5): %d" % len(bbh))

    evts = compute_events(bbh)
    print("Events with final_mass and chi_eff: %d" % len(evts))

    if len(evts) < 5:
        print("ERROR: Too few events. Check data file.")
        return

    excess_ns  = np.array([e['excess_ns']  for e in evts])
    excess_sp  = np.array([e['excess_sp']  for e in evts])
    chi_f      = np.array([e['chi_f']      for e in evts])
    spin_corr  = np.array([e['spin_corr']  for e in evts])
    chi_eff    = np.array([e['chi_eff']    for e in evts])
    chi_eff_med = float(np.median(chi_eff))

    # Fit QS model on both baselines
    alpha_ns = fit_alpha(excess_ns, chi_f)
    alpha_sp = fit_alpha(excess_sp, chi_f)

    # ── Summary statistics ──────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("SPINNING GR BASELINE: COMPARISON SUMMARY")
    print("=" * 70)
    print()
    print("Sample: %d events  |  median chi_eff = %.3f" % (len(evts), chi_eff_med))
    print()
    print("  %-30s  %8s  %8s" % ("Statistic", "Non-spin", "Spinning"))
    print("  " + "-" * 50)

    def row(label, ns_val, sp_val, fmt="%.4f"):
        print("  %-30s  %8s  %8s" % (label, fmt % ns_val, fmt % sp_val))

    row("Median excess deficit",
        np.median(excess_ns), np.median(excess_sp))
    row("Mean excess deficit",
        np.mean(excess_ns),   np.mean(excess_sp))
    row("Fraction positive (%)",
        100 * np.mean(excess_ns > 0), 100 * np.mean(excess_sp > 0), fmt="%.1f%%")
    row("QS model alpha",
        alpha_ns, alpha_sp)

    rho_ns, p_ns = spearmanr(chi_f, excess_ns)
    rho_sp, p_sp = spearmanr(chi_f, excess_sp)
    print()
    print("  Test 1 — Spearman(excess vs chi_f):")
    print("    Non-spinning baseline:  rho=%+.3f  p=%.2e" % (rho_ns, p_ns))
    print("    Spinning baseline:      rho=%+.3f  p=%.2e" % (rho_sp, p_sp))

    print()
    print("  Spin correction magnitude (f_GR_spin - f_GR_nonspin):")
    print("    Median: %.4f  (%.1f%% of non-spinning baseline)"
          % (np.median(spin_corr),
             100 * np.median(spin_corr) / np.median([e['f_gr_ns'] for e in evts])))
    print("    Max:    %.4f" % np.max(spin_corr))
    print("    Min:    %.4f" % np.min(spin_corr))
    print()
    print("  Interpretation:")
    if np.median(excess_sp) > 0 and p_sp < 0.05:
        print("    Excess deficit SURVIVES the spinning GR correction.")
        print("    The quark star model explains residual structure that")
        print("    spinning-progenitor GR does not.")
    elif np.median(excess_sp) <= 0:
        print("    WARNING: Spinning GR baseline eliminates the excess.")
        print("    The quark star model claim requires reassessment.")
    else:
        print("    Excess persists but correlation weakens; review figures.")
    print("=" * 70)

    # ── Figures ─────────────────────────────────────────────────────────────
    plot_test1(evts, alpha_ns, alpha_sp, fig_dir)
    plot_test2(evts, alpha_ns, alpha_sp, chi_eff_med, fig_dir)

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Spinning GR baseline comparison for BH merger analysis')
    parser.add_argument('--no-show', action='store_true',
                        help='Skip interactive figure display (save only)')
    main(parser.parse_args())
