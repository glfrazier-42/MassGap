"""Investigate the eta-confound in Test 1 (excess deficit vs chi_f).

The concern: chi_f is dominated by L_orb(eta), and the excess also has
eta-dependence through f_GR(eta).  A spurious correlation between excess
and chi_f could arise simply because both track eta, independently of any
real spin-excess relationship.

Strategy
--------
Decompose chi_f into its two components:

    chi_f_orb  = L_orb(eta)                       -- pure eta, no spin
    chi_f_spin = chi_eff * (1+q^2) / (1+q)^2      -- progenitor spin, no eta

Then run the following Spearman tests on (excess, X):

    1. X = chi_f          -- full remnant spin (paper's Test 1)
    2. X = chi_f_orb      -- orbital part only  (pure eta proxy)
    3. X = chi_f_spin     -- spin part only      (pure progenitor spin proxy)
    4. X = eta            -- symmetric mass ratio directly
    5. X = chi_eff        -- raw progenitor spin

If the confound is driving the result:
    rho(excess, chi_f_orb) ≈ rho(excess, chi_f)    (eta dominates)
    rho(excess, chi_f_spin) is weak or absent

If the result is genuine spin physics:
    rho(excess, chi_f_spin) is significant
    the partial correlation controlling for eta remains significant

Partial Spearman correlation (excess vs chi_f | eta):
    Rank all three variables.
    Regress rank(excess) on rank(eta) -> residuals_e
    Regress rank(chi_f)  on rank(eta) -> residuals_c
    Pearson r of residuals_e and residuals_c, with its own p-value.

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/investigate_eta_confound.py
"""

import csv
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr, pearsonr

from bhem.gr_baselines import f_gr_spin


def fmt_rho(r):
    return '%+.3f' % r


def fmt_p(p):
    if p < 0.001:
        coeff, exp = ('%.1e' % p).split('e')
        return r'%s \times 10^{%d}' % (coeff, int(exp))
    return '%.4f' % p


# ---------------------------------------------------------------------------
# Data loading (identical to plot_bh_merger_spin.py)
# ---------------------------------------------------------------------------

def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def load_events(filepath):
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
            'name': name,
            'm1': m1, 'm2': m2,
            'm_total': m_total if not np.isnan(m_total) else m1 + m2,
            'm_final': m_final,
            'chi_eff': chi,
        })
    return bbh


# ---------------------------------------------------------------------------
# Physics (identical to plot_bh_merger_spin.py)
# ---------------------------------------------------------------------------

def L_orb(eta):
    return 2 * np.sqrt(3) * eta - 3.5171 * eta**2 + 2.5763 * eta**3


def compute_quantities(bbh):
    rows = []
    for e in bbh:
        if np.isnan(e['m_final']) or e['m_final'] <= 0:
            continue
        if np.isnan(e['chi_eff']):
            continue

        q        = e['m2'] / e['m1']
        eta      = e['m1'] * e['m2'] / e['m_total']**2
        frac_obs = (e['m_total'] - e['m_final']) / e['m_total']
        frac_gr  = f_gr_spin(eta, e['chi_eff'])
        excess   = frac_obs - frac_gr

        l_orb    = L_orb(eta)
        s_spin   = e['chi_eff'] * (1 + q**2) / (1 + q)**2
        chi_f    = np.clip(l_orb + s_spin, 0.0, 1.0)

        rows.append({
            'name':       e['name'],
            'eta':        eta,
            'chi_eff':    e['chi_eff'],
            'chi_f':      chi_f,
            'chi_f_orb':  l_orb,          # L_orb(eta) only
            'chi_f_spin': s_spin,          # progenitor spin only
            'excess':     excess,
        })
    return rows


# ---------------------------------------------------------------------------
# Partial Spearman correlation: corr(X, Y | Z)
# ---------------------------------------------------------------------------

def partial_spearman(x, y, z):
    """Partial rank correlation of x and y controlling for z.

    Method: rank all three, regress rank(x) and rank(y) separately on
    rank(z), then Pearson-correlate the residuals.  Returns (r, p).
    """
    n = len(x)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rz = np.argsort(np.argsort(z)).astype(float)

    # OLS residuals of rx on rz
    b_xz = np.dot(rz - rz.mean(), rx - rx.mean()) / np.dot(rz - rz.mean(), rz - rz.mean())
    res_x = rx - (rz * b_xz + (rx.mean() - b_xz * rz.mean()))

    # OLS residuals of ry on rz
    b_yz = np.dot(rz - rz.mean(), ry - ry.mean()) / np.dot(rz - rz.mean(), rz - rz.mean())
    res_y = ry - (rz * b_yz + (ry.mean() - b_yz * rz.mean()))

    r, p = pearsonr(res_x, res_y)
    return r, p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    base_dir  = Path(__file__).resolve().parent.parent
    data_path = base_dir / 'data' / 'gwtc_catalog.csv'

    events = load_events(str(data_path))
    bbh    = filter_bbh(events)
    rows   = compute_quantities(bbh)
    print("Events used: %d" % len(rows))

    excess      = np.array([r['excess']     for r in rows])
    eta         = np.array([r['eta']        for r in rows])
    chi_eff     = np.array([r['chi_eff']    for r in rows])
    chi_f       = np.array([r['chi_f']      for r in rows])
    chi_f_orb   = np.array([r['chi_f_orb']  for r in rows])
    chi_f_spin  = np.array([r['chi_f_spin'] for r in rows])

    print()
    print("=" * 65)
    print("DECOMPOSITION: chi_f variance explained by each component")
    print("=" * 65)
    rho_orb_f,  _ = spearmanr(chi_f_orb,  chi_f)
    rho_spin_f, _ = spearmanr(chi_f_spin, chi_f)
    print("  Spearman(chi_f_orb,  chi_f) = %+.3f  (how much chi_f tracks eta)"
          % rho_orb_f)
    print("  Spearman(chi_f_spin, chi_f) = %+.3f  (how much chi_f tracks chi_eff)"
          % rho_spin_f)

    print()
    print("=" * 65)
    print("BIVARIATE SPEARMAN:  excess vs various predictors")
    print("=" * 65)
    tests = [
        ("chi_f       (paper Test 1)",  chi_f),
        ("chi_f_orb   (L_orb only)",    chi_f_orb),
        ("chi_f_spin  (progenitor spin only)", chi_f_spin),
        ("eta         (mass ratio direct)",    eta),
        ("chi_eff     (raw progenitor spin)",  chi_eff),
    ]
    for label, x in tests:
        rho, p = spearmanr(x, excess)
        stars = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
        print("  %-40s  rho=%+.3f  p=%.2e  %s" % (label, rho, p, stars))

    print()
    print("=" * 65)
    print("PARTIAL SPEARMAN:  excess vs chi_f | eta")
    print("(controls for shared eta-dependence)")
    print("=" * 65)
    r_partial, p_partial = partial_spearman(excess, chi_f, eta)
    print("  partial rho(excess, chi_f | eta) = %+.3f  p=%.2e" % (r_partial, p_partial))

    print()
    print("=" * 65)
    print("WITHIN-ETA-BIN SPEARMAN")
    print("(does the excess-chi_f correlation persist at fixed eta?)")
    print("=" * 65)
    # Split into terciles of eta
    terciles = np.percentile(eta, [33.3, 66.7])
    bins = [
        ("low eta    (q < %.2f)" % terciles[0],
         eta < terciles[0]),
        ("mid eta    (%.2f - %.2f)" % (terciles[0], terciles[1]),
         (eta >= terciles[0]) & (eta < terciles[1])),
        ("high eta   (q > %.2f)" % terciles[1],
         eta >= terciles[1]),
    ]
    for label, mask in bins:
        n = np.sum(mask)
        if n < 5:
            print("  %-35s  n=%d  (too few)" % (label, n))
            continue
        rho, p = spearmanr(chi_f[mask], excess[mask])
        # also spin-only component
        rho_sp, p_sp = spearmanr(chi_f_spin[mask], excess[mask])
        print("  %-35s  n=%2d  "
              "rho(chi_f,excess)=%+.3f p=%.2e  |  "
              "rho(chi_f_spin,excess)=%+.3f p=%.2e"
              % (label, n, rho, p, rho_sp, p_sp))

    # -----------------------------------------------------------------------
    # Write LaTeX macros
    # -----------------------------------------------------------------------
    rho_chif_orb_chif,  _ = spearmanr(chi_f_orb,  chi_f)
    rho_chif_spin_chif, _ = spearmanr(chi_f_spin, chi_f)
    rho_excess_chif_spin, p_excess_chif_spin = spearmanr(chi_f_spin, excess)
    rho_excess_eta, p_excess_eta = spearmanr(eta, excess)
    r_part, p_part = partial_spearman(excess, chi_f, eta)

    macros = {
        'ecRhoChifOrbChif':   fmt_rho(rho_chif_orb_chif),
        'ecRhoChifSpinChif':  fmt_rho(rho_chif_spin_chif),
        'ecRhoExcessChifSpin': fmt_rho(rho_excess_chif_spin),
        'ecPExcessChifSpin':  fmt_p(p_excess_chif_spin),
        'ecRhoExcessEta':     fmt_rho(rho_excess_eta),
        'ecPExcessEta':       fmt_p(p_excess_eta),
        'ecRhoPartialChif':   fmt_rho(r_part),
        'ecPPartialChif':     fmt_p(p_part),
    }
    table_dir = base_dir / 'latex-paper' / 'tables'
    write_macros(macros, table_dir / 'eta_confound_macros.tex')

    print()
    print("=" * 65)
    print("INTERPRETATION GUIDE")
    print("=" * 65)
    print("""
  Confound diagnosis:
    If rho(excess, chi_f_orb) ≈ rho(excess, chi_f):
        -> correlation is driven by eta, not spin  [confound present]
    If partial rho(excess, chi_f | eta) << rho(excess, chi_f):
        -> eta is the true driver                  [confound present]
    If within-eta-bin correlations are weak/insignificant:
        -> same conclusion                         [confound present]

  Clean signal diagnosis:
    If rho(excess, chi_f_spin) is significant:
        -> progenitor spin (not eta) is contributing
    If partial rho remains significant:
        -> spin effect survives eta control        [result is genuine]
    If within-bin correlations remain positive:
        -> same conclusion                         [result is genuine]
""")


def write_macros(stats, outpath):
    """Write LaTeX macros for the eta-confound robustness paragraph."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, 'w') as f:
        f.write('%% Auto-generated by investigate_eta_confound.py -- do not edit\n')
        for name, val in stats.items():
            f.write('\\newcommand{\\%s}{%s}\n' % (name, val))
    print("Wrote %d macros to %s" % (len(stats), outpath))


if __name__ == '__main__':
    main()
