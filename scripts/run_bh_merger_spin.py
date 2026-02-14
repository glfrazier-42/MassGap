"""
BH Merger Spin Analysis: Testing the Spin-Collapse Mass Gap Hypothesis.

Hypothesis: The pair-instability mass gap (50-150 Msun) is not a formation
gap but a gap in M_g caused by the quark -> preon collapse transition.
Slow-spinning progenitors collapse earlier and land in the gap; fast-spinning
progenitors reach higher mass before collapsing above the gap.

Predictions:
  1. BHs in the mass gap should have low chi_eff (from slow-spinning progenitors)
  2. Mass deficit (m1+m2-m_final) should correlate with remnant spin
  3. The fractional energy radiated should vary systematically with mass

Data: GWTC catalog from LIGO/Virgo/KAGRA

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/run_bh_merger_spin.py
"""

import os
import csv
import numpy as np


def load_gwtc(filepath="data/gwtc_catalog.csv"):
    """Load GWTC catalog, return list of event dicts."""
    events = []
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            events.append(row)
    return events


def safe_float(val):
    """Convert to float, return NaN for empty/invalid."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def main():
    events = load_gwtc()

    # ── Filter for confident BBH events with chi_eff ──
    bbh = []
    for e in events:
        m1 = safe_float(e.get('mass_1_source', ''))
        m2 = safe_float(e.get('mass_2_source', ''))
        chi = safe_float(e.get('chi_eff', ''))
        p_astro = safe_float(e.get('p_astro', ''))
        m_final = safe_float(e.get('final_mass_source', ''))
        m_total = safe_float(e.get('total_mass_source', ''))
        chi_lo = safe_float(e.get('chi_eff_lower', ''))
        chi_hi = safe_float(e.get('chi_eff_upper', ''))
        m_final_lo = safe_float(e.get('final_mass_source_lower', ''))
        m_final_hi = safe_float(e.get('final_mass_source_upper', ''))

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
            'chi_eff_lower': chi_lo,
            'chi_eff_upper': chi_hi,
            'm_final_lower': m_final_lo,
            'm_final_upper': m_final_hi,
        })

    print("=" * 78)
    print("BH MERGER SPIN ANALYSIS: Mass Gap Hypothesis")
    print("=" * 78)
    print(f"\nConfident BBH events: {len(bbh)}")
    with_chi = [e for e in bbh if not np.isnan(e['chi_eff'])]
    with_final = [e for e in bbh if not np.isnan(e['m_final'])]
    print(f"  With chi_eff: {len(with_chi)}")
    print(f"  With final_mass: {len(with_final)}")

    # ==================================================================
    # ANALYSIS 1: chi_eff vs primary mass
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 1: chi_eff vs Primary Mass (m1)")
    print(f"{'=' * 78}")
    print()
    print("Prediction: BHs in the mass gap (m1 > 50 Msun) should have")
    print("lower chi_eff than sub-gap BHs (m1 < 50 Msun), because gap")
    print("objects originated from slow-spinning progenitors.")
    print()

    # Split into mass bins
    bins = [
        ("m1 < 20", lambda e: e['m1'] < 20),
        ("20 <= m1 < 35", lambda e: 20 <= e['m1'] < 35),
        ("35 <= m1 < 50", lambda e: 35 <= e['m1'] < 50),
        ("50 <= m1 < 65", lambda e: 50 <= e['m1'] < 65),
        ("m1 >= 65 (deep gap)", lambda e: e['m1'] >= 65),
    ]

    print(f"  {'Bin':>22s}  {'N':>4s}  {'mean':>7s}  {'median':>7s}"
          f"  {'std':>6s}  {'|chi|<0.1':>9s}")
    print(f"  {'-' * 65}")

    bin_data = {}
    for label, filt in bins:
        subset = [e['chi_eff'] for e in with_chi if filt(e)]
        if len(subset) == 0:
            print(f"  {label:>22s}  {0:4d}  {'--':>7s}  {'--':>7s}"
                  f"  {'--':>6s}  {'--':>9s}")
            continue
        arr = np.array(subset)
        frac_low = np.sum(np.abs(arr) < 0.1) / len(arr)
        print(f"  {label:>22s}  {len(arr):4d}  {arr.mean():+7.3f}"
              f"  {np.median(arr):+7.3f}  {arr.std():6.3f}"
              f"  {frac_low:9.1%}")
        bin_data[label] = arr

    # Sub-gap vs gap comparison
    sub_gap = np.array([e['chi_eff'] for e in with_chi if e['m1'] < 50])
    in_gap = np.array([e['chi_eff'] for e in with_chi if e['m1'] >= 50])

    print()
    print(f"  Sub-gap (m1 < 50):  N={len(sub_gap)}, "
          f"mean={sub_gap.mean():+.3f}, median={np.median(sub_gap):+.3f}")
    print(f"  In-gap  (m1 >= 50): N={len(in_gap)}, "
          f"mean={in_gap.mean():+.3f}, median={np.median(in_gap):+.3f}")

    # Statistical test
    from scipy.stats import mannwhitneyu, spearmanr
    if len(in_gap) >= 3:
        U, p_mw = mannwhitneyu(sub_gap, in_gap, alternative='two-sided')
        print(f"\n  Mann-Whitney U test (sub-gap vs gap): U={U:.0f}, p={p_mw:.4f}")
        if p_mw < 0.05:
            print("  --> SIGNIFICANT difference in chi_eff distributions.")
        else:
            print("  --> No significant difference.")

    # Spearman: chi_eff vs m1 across all events
    m1_arr = np.array([e['m1'] for e in with_chi])
    chi_arr = np.array([e['chi_eff'] for e in with_chi])
    rho, p_sp = spearmanr(m1_arr, chi_arr)
    print(f"\n  Spearman (chi_eff vs m1, all events): rho={rho:+.3f}, p={p_sp:.4f}")

    # Also test |chi_eff| vs m1 (magnitude, ignoring sign)
    rho_abs, p_abs = spearmanr(m1_arr, np.abs(chi_arr))
    print(f"  Spearman (|chi_eff| vs m1): rho={rho_abs:+.3f}, p={p_abs:.4f}")

    # Individual events in/above the gap
    print(f"\n  Events with m1 >= 50 (individual chi_eff):")
    print(f"  {'Name':>25s}  {'m1':>6s}  {'m2':>6s}  {'chi_eff':>8s}"
          f"  {'chi_lo':>7s}  {'chi_hi':>7s}")
    print(f"  {'-' * 65}")
    gap_events = sorted([e for e in with_chi if e['m1'] >= 50],
                        key=lambda e: e['m1'], reverse=True)
    for e in gap_events:
        lo = f"{e['chi_eff_lower']:+.2f}" if not np.isnan(e['chi_eff_lower']) else "  --"
        hi = f"{e['chi_eff_upper']:+.2f}" if not np.isnan(e['chi_eff_upper']) else "  --"
        print(f"  {e['name']:>25s}  {e['m1']:6.1f}  {e['m2']:6.1f}"
              f"  {e['chi_eff']:+8.2f}  {lo:>7s}  {hi:>7s}")

    # ==================================================================
    # ANALYSIS 1b: Angular Momentum Correction
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 1b: Spin-Up Correction (Angular Momentum)")
    print(f"{'=' * 78}")
    print()
    print("The naive prediction assumed chi_progenitor ~ chi_BH. But collapse")
    print("conserves angular momentum J while changing M_g:")
    print()
    print("  chi_BH = c * J / (G * M_g^2)")
    print()
    print("J is set by the progenitor: J = I_qs * omega_qs")
    print("After collapse, M_g changes but J does not.")
    print()
    print("For two BHs from slow vs fast progenitors:")
    print("  chi_slow / chi_fast = (J_slow / J_fast) * (M_g,fast / M_g,slow)^2")
    print()
    print("The first factor is < 1 (slow spinner has less J).")
    print("The second factor is > 1 and SQUARED (above-gap mass >> gap mass).")
    print("The prediction FLIPS when M_g,fast/M_g,slow > sqrt(J_fast/J_slow).")
    print()

    # Parametric exploration
    print("  Parametric: when does the prediction flip?")
    print(f"  {'J_fast/J_slow':>14s}  {'sqrt':>6s}  {'M_g ratio needed':>18s}"
          f"  {'e.g. 80 -> ?':>14s}")
    print(f"  {'-' * 58}")
    for j_ratio in [2, 3, 5, 10, 20, 50]:
        threshold = np.sqrt(j_ratio)
        example = threshold * 80
        flip = "FLIPS" if example < 300 else "unlikely"
        print(f"  {j_ratio:14d}  {threshold:6.2f}  {threshold:18.2f}"
              f"  {example:10.0f} Msun  {flip}")

    print()
    print("  For J_fast/J_slow ~ 5-10 (reasonable), the prediction flips")
    print("  if M_g above the gap is ~2-3x the gap mass. This is plausible")
    print("  (gap ~ 80 Msun, above-gap > 150 Msun -> ratio ~ 2).")
    print()

    # Compute J_eff = |chi_eff| * m1^2 as proxy for progenitor J
    print("  ANGULAR MOMENTUM PROXY: J_eff = |chi_eff| * m1^2")
    print("  (proportional to primary BH angular momentum)")
    print()

    j_sub = np.array([abs(e['chi_eff']) * e['m1']**2
                       for e in with_chi if e['m1'] < 50])
    j_gap = np.array([abs(e['chi_eff']) * e['m1']**2
                       for e in with_chi if e['m1'] >= 50])

    print(f"  Sub-gap (m1 < 50):  N={len(j_sub)}, "
          f"median J_eff = {np.median(j_sub):.0f}")
    print(f"  In-gap  (m1 >= 50): N={len(j_gap)}, "
          f"median J_eff = {np.median(j_gap):.0f}")

    if len(j_gap) >= 3:
        U_j, p_j = mannwhitneyu(j_sub, j_gap, alternative='two-sided')
        print(f"\n  Mann-Whitney U (J_eff sub-gap vs gap): U={U_j:.0f}, p={p_j:.4f}")

    # J_eff by mass bin
    print(f"\n  {'Bin':>22s}  {'N':>4s}  {'median J_eff':>12s}"
          f"  {'mean J_eff':>12s}")
    print(f"  {'-' * 55}")
    for label, filt in bins:
        subset = [abs(e['chi_eff']) * e['m1']**2
                  for e in with_chi if filt(e)]
        if len(subset) == 0:
            continue
        arr = np.array(subset)
        print(f"  {label:>22s}  {len(arr):4d}  {np.median(arr):12.0f}"
              f"  {arr.mean():12.0f}")

    # Spearman on J_eff vs m1
    j_all = np.array([abs(e['chi_eff']) * e['m1']**2 for e in with_chi])
    rho_j, p_j_sp = spearmanr(m1_arr, j_all)
    print(f"\n  Spearman (J_eff vs m1): rho={rho_j:+.3f}, p={p_j_sp:.4f}")

    print()
    print("  INTERPRETATION:")
    print("  If gap BHs have higher chi but LOWER J than above-gap BHs,")
    print("  that is consistent with the spin-up model: low progenitor J")
    print("  is amplified by the 1/M_g^2 factor in the chi formula.")
    print("  If gap BHs have higher J, that contradicts slow-spinning")
    print("  progenitors regardless of the spin-up correction.")

    # ==================================================================
    # ANALYSIS 2: Mass deficit vs total mass
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 2: Fractional Mass Deficit")
    print(f"{'=' * 78}")
    print()
    print("If spin reduces M_g, the apparent 'mass radiated as GWs' includes")
    print("a spin-induced M_g reduction. The fractional deficit should correlate")
    print("with mass (more massive mergers -> more angular momentum -> more spin).")
    print()

    deficit_events = [e for e in with_final
                      if not np.isnan(e['m_final']) and e['m_final'] > 0]

    if len(deficit_events) > 0:
        m_tot_arr = np.array([e['m_total'] for e in deficit_events])
        m_fin_arr = np.array([e['m_final'] for e in deficit_events])
        deficit_arr = m_tot_arr - m_fin_arr
        frac_deficit = deficit_arr / m_tot_arr

        print(f"  Events with final mass measured: {len(deficit_events)}")
        print(f"  Mass deficit (m1+m2 - m_final):")
        print(f"    Range: {deficit_arr.min():.1f} - {deficit_arr.max():.1f} Msun")
        print(f"    Median: {np.median(deficit_arr):.1f} Msun")
        print(f"  Fractional deficit (deficit / m_total):")
        print(f"    Range: {frac_deficit.min():.3f} - {frac_deficit.max():.3f}")
        print(f"    Median: {np.median(frac_deficit):.3f}")
        print(f"    Mean: {np.mean(frac_deficit):.3f}")

        # Correlation: fractional deficit vs total mass
        rho_d, p_d = spearmanr(m_tot_arr, frac_deficit)
        print(f"\n  Spearman (frac_deficit vs m_total): rho={rho_d:+.3f}, p={p_d:.4f}")

        # Correlation: fractional deficit vs chi_eff
        chi_for_deficit = np.array([e['chi_eff'] for e in deficit_events
                                     if not np.isnan(e['chi_eff'])])
        frac_for_chi = np.array([
            (e['m_total'] - e['m_final']) / e['m_total']
            for e in deficit_events if not np.isnan(e['chi_eff'])])
        if len(chi_for_deficit) > 3:
            rho_dc, p_dc = spearmanr(chi_for_deficit, frac_for_chi)
            print(f"  Spearman (frac_deficit vs chi_eff): rho={rho_dc:+.3f}, p={p_dc:.4f}")

        # Binned by total mass
        mass_bins = [
            ("M_total < 30", lambda e: e['m_total'] < 30),
            ("30 <= M < 50", lambda e: 30 <= e['m_total'] < 50),
            ("50 <= M < 80", lambda e: 50 <= e['m_total'] < 80),
            ("80 <= M < 120", lambda e: 80 <= e['m_total'] < 120),
            ("M_total >= 120", lambda e: e['m_total'] >= 120),
        ]

        print(f"\n  {'Bin':>18s}  {'N':>4s}  {'med deficit':>12s}"
              f"  {'med frac':>9s}  {'med chi_eff':>11s}")
        print(f"  {'-' * 60}")

        for label, filt in mass_bins:
            subset = [e for e in deficit_events if filt(e)]
            if len(subset) == 0:
                continue
            deficits = np.array([e['m_total'] - e['m_final'] for e in subset])
            fracs = deficits / np.array([e['m_total'] for e in subset])
            chis = np.array([e['chi_eff'] for e in subset
                             if not np.isnan(e['chi_eff'])])
            chi_str = f"{np.median(chis):+.3f}" if len(chis) > 0 else "--"
            print(f"  {label:>18s}  {len(subset):4d}  {np.median(deficits):12.1f}"
                  f"  {np.median(fracs):9.3f}  {chi_str:>11s}")

        # Show individual high-mass events
        print(f"\n  High-mass events (M_total >= 100):")
        print(f"  {'Name':>25s}  {'m1':>6s}  {'m2':>6s}  {'M_tot':>6s}"
              f"  {'M_fin':>6s}  {'deficit':>7s}  {'frac':>6s}  {'chi':>6s}")
        print(f"  {'-' * 78}")
        high_mass = sorted([e for e in deficit_events if e['m_total'] >= 100],
                           key=lambda e: e['m_total'], reverse=True)
        for e in high_mass:
            deficit = e['m_total'] - e['m_final']
            frac = deficit / e['m_total']
            chi_str = f"{e['chi_eff']:+.2f}" if not np.isnan(e['chi_eff']) else "--"
            print(f"  {e['name']:>25s}  {e['m1']:6.1f}  {e['m2']:6.1f}"
                  f"  {e['m_total']:6.1f}  {e['m_final']:6.1f}"
                  f"  {deficit:7.1f}  {frac:6.3f}  {chi_str:>6s}")

    # ==================================================================
    # ANALYSIS 3: GR prediction for mass deficit
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 3: Observed vs GR-Predicted Mass Deficit")
    print(f"{'=' * 78}")
    print()
    print("In standard GR, the fractional mass radiated depends primarily on")
    print("the symmetric mass ratio eta = m1*m2/(m1+m2)^2 and the spins.")
    print("For non-spinning, equal-mass mergers: ~5% radiated.")
    print("For spinning mergers, this can reach ~10%.")
    print()
    print("Under the spin-collapse hypothesis, the apparent deficit includes")
    print("BOTH GW radiation AND spin-induced M_g reduction. So observed")
    print("deficits should systematically EXCEED the GR prediction for")
    print("non-spinning mergers.")
    print()

    if len(deficit_events) > 0:
        # GR prediction for non-spinning mergers (Buonanno, Kidder & Lehner 2008)
        # E_rad/M = (1 - sqrt(8/9)) * eta + 0.498 * eta^2
        # = 0.0572 * eta + 0.498 * eta^2
        # Gives ~4.5% at eta=0.25 (equal mass), ~3% at eta=0.16 (q=0.5)

        print(f"  {'Name':>25s}  {'q':>5s}  {'eta':>5s}  {'frac_obs':>9s}"
              f"  {'frac_GR':>8s}  {'excess':>7s}  {'chi':>6s}")
        print(f"  {'-' * 72}")

        excess_list = []
        chi_for_excess = []
        m1_for_excess = []

        for e in sorted(deficit_events, key=lambda x: x['m_total'], reverse=True):
            q = e['m2'] / e['m1']
            eta = e['m1'] * e['m2'] / e['m_total']**2
            frac_obs = (e['m_total'] - e['m_final']) / e['m_total']

            # Non-spinning GR prediction (Buonanno, Kidder & Lehner 2008)
            frac_gr = 0.0572 * eta + 0.498 * eta**2

            excess = frac_obs - frac_gr
            chi_str = f"{e['chi_eff']:+.2f}" if not np.isnan(e['chi_eff']) else "--"

            if e['m_total'] >= 80 or excess > 0.03:
                print(f"  {e['name']:>25s}  {q:5.2f}  {eta:5.3f}"
                      f"  {frac_obs:9.3f}  {frac_gr:8.3f}"
                      f"  {excess:+7.3f}  {chi_str:>6s}")

            excess_list.append(excess)
            if not np.isnan(e['chi_eff']):
                chi_for_excess.append(e['chi_eff'])
                m1_for_excess.append(e['m1'])

        excess_arr = np.array(excess_list)
        print(f"\n  All events:")
        print(f"    Median excess over non-spinning GR: {np.median(excess_arr):+.3f}")
        print(f"    Mean excess: {np.mean(excess_arr):+.3f}")
        print(f"    Fraction with positive excess: "
              f"{np.sum(excess_arr > 0)}/{len(excess_arr)}"
              f" ({np.sum(excess_arr > 0)/len(excess_arr):.1%})")

        if len(chi_for_excess) > 3:
            rho_ex, p_ex = spearmanr(np.abs(chi_for_excess), excess_arr[:len(chi_for_excess)])
            print(f"    Spearman (excess vs |chi_eff|): rho={rho_ex:+.3f}, p={p_ex:.4f}")

            rho_em, p_em = spearmanr(m1_for_excess, excess_arr[:len(m1_for_excess)])
            print(f"    Spearman (excess vs m1): rho={rho_em:+.3f}, p={p_em:.4f}")

    # ==================================================================
    # ANALYSIS 4: Remnant Spin (chi_f) Distribution
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 4: Computed Remnant Spin (chi_f)")
    print(f"{'=' * 78}")
    print()
    print("chi_f is not directly measured by LIGO. It is computed from m1, m2,")
    print("and chi_eff using NR fitting formulas. We use:")
    print()
    print("  chi_f = L_orb(eta) + S/M^2")
    print("  L_orb(eta) = 2*sqrt(3)*eta - 3.5171*eta^2 + 2.5763*eta^3")
    print("  S/M^2 = chi_eff * (1 + q^2) / (1 + q)^2")
    print("  (assuming a1 = a2 = chi_eff)")
    print()
    print("Under the preon star hypothesis, chi_f should be capped at ~0.7")
    print("by the structural spin limit (centrifugal expansion feedback).")
    print("Under standard GR, chi_f can reach ~0.95 for aligned high-spin")
    print("progenitors.")
    print()

    chi_f_list = []
    chi_f_events = []
    for e in bbh:
        if np.isnan(e['chi_eff']):
            continue
        q = e['m2'] / e['m1']
        eta = e['m1'] * e['m2'] / e['m_total']**2

        # Orbital angular momentum contribution (non-spinning)
        L_orb = 2 * np.sqrt(3) * eta - 3.5171 * eta**2 + 2.5763 * eta**3

        # Spin angular momentum contribution
        # Assuming a1 = a2 = chi_eff (equal aligned spins)
        S_over_M2 = e['chi_eff'] * (1 + q**2) / (1 + q)**2

        chi_f = L_orb + S_over_M2
        chi_f = min(max(chi_f, 0.0), 1.0)  # physical bounds

        chi_f_list.append(chi_f)
        chi_f_events.append({**e, 'chi_f': chi_f, 'q': q, 'eta': eta})

    chi_f_arr = np.array(chi_f_list)

    print(f"  Events with computed chi_f: {len(chi_f_arr)}")
    print(f"  chi_f distribution:")
    print(f"    Min:    {chi_f_arr.min():.3f}")
    print(f"    Max:    {chi_f_arr.max():.3f}")
    print(f"    Mean:   {chi_f_arr.mean():.3f}")
    print(f"    Median: {np.median(chi_f_arr):.3f}")
    print(f"    Std:    {chi_f_arr.std():.3f}")

    # Fraction above various thresholds
    print()
    for threshold in [0.65, 0.70, 0.75, 0.80, 0.85, 0.90]:
        n_above = np.sum(chi_f_arr > threshold)
        print(f"    chi_f > {threshold:.2f}: {n_above:3d} / {len(chi_f_arr)}"
              f" ({n_above/len(chi_f_arr):.1%})")

    # Symmetry test: compare upper and lower tails around the median
    med_f = np.median(chi_f_arr)
    print(f"\n  SYMMETRY TEST (is the upper tail just measurement error?)")
    print(f"  If the true chi_f is capped at ~{med_f:.2f} and the upper tail")
    print(f"  is from error, the lower tail should mirror it.")
    print()
    print(f"  {'Distance from median':>22s}  {'N above':>8s}  {'N below':>8s}"
          f"  {'Ratio':>6s}")
    print(f"  {'-' * 50}")
    for delta in [0.05, 0.10, 0.15, 0.20, 0.25]:
        n_above = np.sum(chi_f_arr > med_f + delta)
        n_below = np.sum(chi_f_arr < med_f - delta)
        ratio = n_above / n_below if n_below > 0 else float('inf')
        print(f"  {delta:22.2f}  {n_above:8d}  {n_below:8d}  {ratio:6.2f}")

    # Skewness
    from scipy.stats import skew, skewtest
    sk = skew(chi_f_arr)
    try:
        sk_stat, sk_p = skewtest(chi_f_arr)
        print(f"\n  Skewness: {sk:+.3f}")
        print(f"  Skewness test: z={sk_stat:.2f}, p={sk_p:.4f}")
        if sk_p < 0.05:
            if sk > 0:
                print("  --> Significantly right-skewed (upper tail is heavier)")
            else:
                print("  --> Significantly left-skewed (lower tail is heavier)")
        else:
            print("  --> Not significantly skewed (symmetric)")
    except Exception:
        print(f"\n  Skewness: {sk:+.3f}")

    # Events with highest chi_f
    top_events = sorted(chi_f_events, key=lambda e: e['chi_f'], reverse=True)
    print(f"\n  Top 15 events by chi_f:")
    print(f"  {'Name':>25s}  {'m1':>6s}  {'m2':>6s}  {'q':>5s}"
          f"  {'chi_eff':>8s}  {'chi_f':>6s}")
    print(f"  {'-' * 62}")
    for e in top_events[:15]:
        print(f"  {e['name']:>25s}  {e['m1']:6.1f}  {e['m2']:6.1f}"
              f"  {e['q']:5.2f}  {e['chi_eff']:+8.3f}  {e['chi_f']:6.3f}")

    # Events with lowest chi_f
    print(f"\n  Bottom 10 events by chi_f:")
    print(f"  {'Name':>25s}  {'m1':>6s}  {'m2':>6s}  {'q':>5s}"
          f"  {'chi_eff':>8s}  {'chi_f':>6s}")
    print(f"  {'-' * 62}")
    for e in top_events[-10:]:
        print(f"  {e['name']:>25s}  {e['m1']:6.1f}  {e['m2']:6.1f}"
              f"  {e['q']:5.2f}  {e['chi_eff']:+8.3f}  {e['chi_f']:6.3f}")

    # chi_f vs chi_eff correlation
    chi_eff_for_f = np.array([e['chi_eff'] for e in chi_f_events])
    rho_ff, p_ff = spearmanr(chi_eff_for_f, chi_f_arr)
    print(f"\n  Spearman (chi_f vs chi_eff): rho={rho_ff:+.3f}, p={p_ff:.4f}")

    # ==================================================================
    # ANALYSIS 5: Mass Ratio Distribution (Generational Structure?)
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 5: Mass Ratio Distribution")
    print(f"{'=' * 78}")
    print()
    print("If BHs grow through generational mergers (Gen1 ~5-10 Msun from NS")
    print("collapse, Gen2 ~10-20 from Gen1 mergers, etc.), mass ratios should")
    print("cluster at discrete values: q~1 (same-generation mergers) and")
    print("q~0.5 (adjacent-generation mergers).")
    print()

    q_all = np.array([e['m2'] / e['m1'] for e in bbh])
    print(f"  Mass ratio (q = m2/m1) distribution:")
    print(f"    N = {len(q_all)}")
    print(f"    Min: {q_all.min():.3f}")
    print(f"    Max: {q_all.max():.3f}")
    print(f"    Mean: {q_all.mean():.3f}")
    print(f"    Median: {np.median(q_all):.3f}")
    print(f"    Std: {q_all.std():.3f}")

    # Binned distribution
    q_edges = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
    print(f"\n  {'q range':>12s}  {'N':>4s}  {'fraction':>8s}  {'bar':s}")
    print(f"  {'-' * 45}")
    for i in range(len(q_edges) - 1):
        n = np.sum((q_all >= q_edges[i]) & (q_all < q_edges[i+1]))
        frac = n / len(q_all)
        bar = '#' * int(frac * 100)
        print(f"  {q_edges[i]:.1f} - {q_edges[i+1]:.2f}  {n:4d}  {frac:8.1%}  {bar}")

    # m2 distribution: does the secondary mass cluster at preferred values?
    m2_all = np.array([e['m2'] for e in bbh])
    print(f"\n  Secondary mass (m2) distribution:")
    print(f"    Min: {m2_all.min():.1f}")
    print(f"    Max: {m2_all.max():.1f}")
    print(f"    Mean: {m2_all.mean():.1f}")
    print(f"    Median: {np.median(m2_all):.1f}")

    m2_edges = [0, 8, 12, 16, 20, 25, 30, 40, 50, 110]
    print(f"\n  {'m2 range':>14s}  {'N':>4s}  {'fraction':>8s}  {'bar':s}")
    print(f"  {'-' * 50}")
    for i in range(len(m2_edges) - 1):
        n = np.sum((m2_all >= m2_edges[i]) & (m2_all < m2_edges[i+1]))
        frac = n / len(m2_all)
        bar = '#' * int(frac * 100)
        print(f"  {m2_edges[i]:3.0f} - {m2_edges[i+1]:3.0f} Msun"
              f"  {n:4d}  {frac:8.1%}  {bar}")

    # For deficit events: which mass ratios produce the observed deficit bands?
    if len(deficit_events) > 0:
        print(f"\n  Fractional deficit by mass ratio band:")
        print(f"  {'q range':>14s}  {'N':>4s}  {'GR pred':>8s}"
              f"  {'med obs':>8s}  {'med excess':>11s}")
        print(f"  {'-' * 52}")
        q_bands = [(0.7, 1.01, "0.70-1.00"),
                   (0.5, 0.7, "0.50-0.70"),
                   (0.3, 0.5, "0.30-0.50"),
                   (0.1, 0.3, "0.10-0.30")]
        for q_lo, q_hi, label in q_bands:
            subset = [e for e in deficit_events
                      if q_lo <= e['m2']/e['m1'] < q_hi]
            if len(subset) == 0:
                continue
            fracs = np.array([(e['m_total'] - e['m_final'])/e['m_total']
                              for e in subset])
            q_mid = (q_lo + q_hi) / 2
            eta_mid = q_mid / (1 + q_mid)**2
            gr_pred = 0.0572 * eta_mid + 0.498 * eta_mid**2
            excess = np.median(fracs) - gr_pred
            print(f"  {label:>14s}  {len(subset):4d}  {gr_pred:8.4f}"
                  f"  {np.median(fracs):8.4f}  {excess:+11.4f}")

    # ==================================================================
    # ANALYSIS 6: Spin Efficiency — Deficit per Unit m2 vs q
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("ANALYSIS 6: Spin Efficiency — Does Low q Deposit More Spin per m2?")
    print(f"{'=' * 78}")
    print()
    print("Physical prediction: in a lopsided merger (low q), the small m2")
    print("orbits around and plunges into the large m1, depositing angular")
    print("momentum efficiently. In an equal-mass merger (q~1), the symmetric")
    print("collision wastes energy on compression. So the EXCESS mass deficit")
    print("per unit m2 should increase as q decreases.")
    print()

    if len(deficit_events) > 0:
        # Compute per-event quantities
        spin_eff_events = []
        for e in deficit_events:
            q = e['m2'] / e['m1']
            eta = e['m1'] * e['m2'] / e['m_total']**2
            frac_obs = (e['m_total'] - e['m_final']) / e['m_total']
            frac_gr = 0.0572 * eta + 0.498 * eta**2
            excess = frac_obs - frac_gr
            deficit_abs = e['m_total'] - e['m_final']
            # Deficit per unit m2 (Msun of deficit per Msun of m2)
            deficit_per_m2 = deficit_abs / e['m2']
            # GR-predicted deficit per m2
            gr_deficit_per_m2 = (frac_gr * e['m_total']) / e['m2']
            # Excess deficit per m2
            excess_per_m2 = deficit_per_m2 - gr_deficit_per_m2

            spin_eff_events.append({
                'name': e['name'], 'q': q, 'eta': eta,
                'm1': e['m1'], 'm2': e['m2'], 'm_total': e['m_total'],
                'frac_obs': frac_obs, 'frac_gr': frac_gr,
                'excess': excess,
                'deficit_per_m2': deficit_per_m2,
                'gr_deficit_per_m2': gr_deficit_per_m2,
                'excess_per_m2': excess_per_m2,
                'chi_eff': e['chi_eff'],
            })

        q_arr = np.array([e['q'] for e in spin_eff_events])
        excess_arr = np.array([e['excess'] for e in spin_eff_events])
        deficit_per_m2_arr = np.array([e['deficit_per_m2']
                                       for e in spin_eff_events])
        excess_per_m2_arr = np.array([e['excess_per_m2']
                                      for e in spin_eff_events])

        # Correlations
        rho_eq, p_eq = spearmanr(q_arr, excess_arr)
        rho_dm, p_dm = spearmanr(q_arr, deficit_per_m2_arr)
        rho_em, p_em = spearmanr(q_arr, excess_per_m2_arr)

        print(f"  Spearman correlations (all {len(spin_eff_events)} events):")
        print(f"    Excess frac deficit vs q:  rho={rho_eq:+.3f}, p={p_eq:.4f}")
        print(f"    Deficit/m2 vs q:           rho={rho_dm:+.3f}, p={p_dm:.4f}")
        print(f"    Excess deficit/m2 vs q:    rho={rho_em:+.3f}, p={p_em:.4f}")
        print()

        if rho_eq < 0 and p_eq < 0.05:
            print("  --> CONFIRMED: Lower q produces larger excess deficit.")
        elif rho_eq < 0:
            print("  --> Trend in expected direction but not significant.")
        else:
            print("  --> No trend detected.")

        # Check confound: does chi_eff correlate with q?
        chi_q = [(e['q'], e['chi_eff']) for e in spin_eff_events
                 if not np.isnan(e['chi_eff'])]
        if len(chi_q) > 3:
            q_c = np.array([x[0] for x in chi_q])
            chi_c = np.array([x[1] for x in chi_q])
            rho_cq, p_cq = spearmanr(q_c, chi_c)
            print(f"\n  Confound check: chi_eff vs q:")
            print(f"    rho={rho_cq:+.3f}, p={p_cq:.4f}")
            if p_cq < 0.05:
                print("    WARNING: chi_eff correlates with q.")
                print("    The excess could be from progenitor spin, not geometry.")
            else:
                print("    No correlation — the excess is NOT from progenitor spin.")

        # Binned analysis with finer bins
        q_bin_edges = [0.15, 0.35, 0.50, 0.60, 0.70, 0.80, 0.95]
        print(f"\n  {'q bin':>12s}  {'N':>3s}  {'med q':>6s}  {'med frac':>9s}"
              f"  {'GR pred':>8s}  {'excess':>7s}  {'def/m2':>7s}"
              f"  {'GR/m2':>6s}  {'exc/m2':>7s}")
        print(f"  {'-' * 78}")

        for i in range(len(q_bin_edges) - 1):
            q_lo, q_hi = q_bin_edges[i], q_bin_edges[i + 1]
            subset = [e for e in spin_eff_events if q_lo <= e['q'] < q_hi]
            if len(subset) < 2:
                continue
            mq = np.median([e['q'] for e in subset])
            mf = np.median([e['frac_obs'] for e in subset])
            mg = np.median([e['frac_gr'] for e in subset])
            me = np.median([e['excess'] for e in subset])
            md = np.median([e['deficit_per_m2'] for e in subset])
            mgm = np.median([e['gr_deficit_per_m2'] for e in subset])
            mem = np.median([e['excess_per_m2'] for e in subset])
            label = f"{q_lo:.2f}-{q_hi:.2f}"
            print(f"  {label:>12s}  {len(subset):3d}  {mq:6.3f}  {mf:9.4f}"
                  f"  {mg:8.4f}  {me:+7.4f}  {md:7.3f}"
                  f"  {mgm:6.3f}  {mem:+7.4f}")

        # Theoretical prediction: chi_f per unit (m2/M_total)
        print(f"\n  THEORETICAL: chi_f per unit (m2/M_total) from NR formula:")
        print(f"  {'q':>6s}  {'eta':>6s}  {'chi_f':>6s}  {'m2/M':>6s}"
              f"  {'chi_f / (m2/M)':>14s}")
        print(f"  {'-' * 44}")
        for q_val in [0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]:
            eta_val = q_val / (1 + q_val)**2
            L_orb = (2 * np.sqrt(3) * eta_val - 3.5171 * eta_val**2
                     + 2.5763 * eta_val**3)
            m2_frac = q_val / (1 + q_val)
            eff = L_orb / m2_frac
            print(f"  {q_val:6.2f}  {eta_val:6.4f}  {L_orb:6.3f}  {m2_frac:6.3f}"
                  f"  {eff:14.3f}")

    # ==================================================================
    # PLOT
    # ==================================================================
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n  matplotlib not available -- skipping plot.")
        return

    fig, axes = plt.subplots(2, 3, figsize=(19, 11))

    # Panel 1: chi_eff vs m1
    ax = axes[0, 0]
    m1_plot = np.array([e['m1'] for e in with_chi])
    chi_plot = np.array([e['chi_eff'] for e in with_chi])
    chi_lo_plot = np.array([abs(safe_float(e.get('chi_eff_lower', np.nan)))
                            for e in with_chi
                            if not np.isnan(e['chi_eff'])])
    chi_hi_plot = np.array([abs(safe_float(e.get('chi_eff_upper', np.nan)))
                            for e in with_chi
                            if not np.isnan(e['chi_eff'])])

    # Color by mass bin
    colors = np.where(m1_plot >= 65, 'red',
             np.where(m1_plot >= 50, 'orange',
             np.where(m1_plot >= 35, 'steelblue', 'navy')))

    ax.scatter(m1_plot, chi_plot, c=colors, s=20, alpha=0.7, edgecolors='k',
               linewidths=0.3, zorder=5)
    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.axvline(50, color='red', linewidth=1, linestyle=':', alpha=0.5,
               label='Gap boundary (50 Msun)')
    ax.set_xlabel("Primary mass m1 (Msun)", fontsize=11)
    ax.set_ylabel("chi_eff", fontsize=11)
    ax.set_title(f"chi_eff vs m1 (rho={rho:+.2f}, p={p_sp:.3f})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Running median
    sort_idx = np.argsort(m1_plot)
    m1_sorted = m1_plot[sort_idx]
    chi_sorted = chi_plot[sort_idx]
    window = max(len(m1_sorted) // 8, 5)
    running_med_m1 = []
    running_med_chi = []
    for i in range(window, len(m1_sorted) - window):
        running_med_m1.append(np.median(m1_sorted[i-window:i+window]))
        running_med_chi.append(np.median(chi_sorted[i-window:i+window]))
    ax.plot(running_med_m1, running_med_chi, 'r-', linewidth=2, alpha=0.7,
            label='Running median')
    ax.legend(fontsize=8)

    # Panel 2: chi_f histogram
    ax = axes[0, 1]
    ax.hist(chi_f_arr, bins=30, range=(0.3, 1.0), color='steelblue',
            edgecolor='k', linewidth=0.5, alpha=0.8)
    ax.axvline(0.69, color='green', linewidth=2, linestyle='--',
               label='Equal-mass non-spinning (0.69)')
    ax.axvline(np.median(chi_f_arr), color='red', linewidth=2,
               linestyle='-', label=f'Median ({np.median(chi_f_arr):.3f})')
    ax.axvline(0.7, color='orange', linewidth=1.5, linestyle=':',
               label='Structural spin limit (~0.7)?')
    ax.set_xlabel(r"Computed $\chi_f$ (remnant spin)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(f"Remnant Spin Distribution (N={len(chi_f_arr)})")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    # Panel 3: Fractional deficit vs total mass, colored by m2
    ax = axes[1, 0]
    if len(deficit_events) > 0:
        m_tot_plot = np.array([e['m_total'] for e in deficit_events])
        frac_plot = np.array([(e['m_total'] - e['m_final']) / e['m_total']
                              for e in deficit_events])
        m2_plot = np.array([e['m2'] for e in deficit_events])

        # Color by m2
        sc = ax.scatter(m_tot_plot, frac_plot, c=m2_plot, cmap='viridis',
                        s=25, alpha=0.8, edgecolors='k', linewidths=0.3,
                        zorder=5, vmin=5, vmax=60)
        cbar = plt.colorbar(sc, ax=ax, label='Secondary mass m2 (Msun)',
                             pad=0.02)
        cbar.ax.tick_params(labelsize=7)

        # Draw fixed-m2 curves (GR non-spinning prediction)
        m2_curves = [7, 10, 15, 20, 25, 30, 40]
        curve_colors = ['#440154', '#443983', '#31688e', '#21918c',
                        '#35b779', '#90d743', '#fde725']
        for m2_val, cc in zip(m2_curves, curve_colors):
            m_total_range = np.linspace(2 * m2_val, 250, 200)
            m1_range = m_total_range - m2_val
            eta_range = m1_range * m2_val / m_total_range**2
            frac_range = 0.0572 * eta_range + 0.498 * eta_range**2
            ax.plot(m_total_range, frac_range, color=cc, linewidth=1.2,
                    alpha=0.7, label=f'm2={m2_val}')

        ax.set_xlabel("Total mass (Msun)", fontsize=11)
        ax.set_ylabel("Fractional mass deficit", fontsize=11)
        ax.set_title("Mass Deficit vs Total Mass (colored by m2)")
        ax.legend(fontsize=5.5, ncol=2, loc='upper right')
        ax.set_xlim(10, 250)
        ax.set_ylim(0.01, 0.065)
        ax.grid(True, alpha=0.3)

    # Panel 4: Excess deficit over GR prediction vs chi_eff
    ax = axes[1, 1]
    if len(deficit_events) > 0:
        excess_plot = []
        chi_exc_plot = []
        m1_exc_plot = []
        for e in deficit_events:
            if np.isnan(e['chi_eff']):
                continue
            q = e['m2'] / e['m1']
            eta = e['m1'] * e['m2'] / e['m_total']**2
            frac_obs = (e['m_total'] - e['m_final']) / e['m_total']
            frac_gr = 0.0572 * eta + 0.498 * eta**2
            excess_plot.append(frac_obs - frac_gr)
            chi_exc_plot.append(e['chi_eff'])
            m1_exc_plot.append(e['m1'])

        excess_plot = np.array(excess_plot)
        chi_exc_plot = np.array(chi_exc_plot)
        m1_exc_plot = np.array(m1_exc_plot)

        colors_exc = np.where(m1_exc_plot >= 65, 'red',
                     np.where(m1_exc_plot >= 50, 'orange',
                     np.where(m1_exc_plot >= 35, 'steelblue', 'navy')))
        ax.scatter(chi_exc_plot, excess_plot, c=colors_exc, s=20, alpha=0.7,
                   edgecolors='k', linewidths=0.3, zorder=5)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.set_xlabel("chi_eff", fontsize=11)
        ax.set_ylabel("Excess deficit over non-spinning GR", fontsize=11)
        ax.set_title("Excess Mass Deficit vs Spin")
        ax.grid(True, alpha=0.3)

        # Legend for colors
        for label, color in [('m1 < 35', 'navy'), ('35-50', 'steelblue'),
                              ('50-65 (gap)', 'orange'), ('65+ (deep gap)', 'red')]:
            ax.scatter([], [], c=color, s=30, edgecolors='k', linewidths=0.3,
                       label=label)
        ax.legend(fontsize=7, loc='upper left')

    # Panel 5: Excess deficit vs q (spin efficiency)
    ax = axes[0, 2]
    if len(deficit_events) > 0 and len(spin_eff_events) > 0:
        q_se = np.array([e['q'] for e in spin_eff_events])
        excess_se = np.array([e['excess'] for e in spin_eff_events])

        ax.scatter(q_se, excess_se, c='steelblue', s=20, alpha=0.6,
                   edgecolors='k', linewidths=0.3, zorder=5)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')

        # Binned medians
        q_centers = []
        excess_medians = []
        for i in range(len(q_bin_edges) - 1):
            q_lo, q_hi = q_bin_edges[i], q_bin_edges[i + 1]
            sub = [e['excess'] for e in spin_eff_events
                   if q_lo <= e['q'] < q_hi]
            if len(sub) >= 2:
                q_centers.append((q_lo + q_hi) / 2)
                excess_medians.append(np.median(sub))
        ax.plot(q_centers, excess_medians, 'ro-', linewidth=2, markersize=8,
                markeredgecolor='k', markeredgewidth=0.5, zorder=10,
                label='Binned median')

        ax.set_xlabel("Mass ratio q = m2/m1", fontsize=11)
        ax.set_ylabel("Excess fractional deficit over GR", fontsize=11)
        ax.set_title(f"Spin Efficiency: Excess vs q (rho={rho_eq:+.2f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Panel 6: Deficit per unit m2 vs q
    ax = axes[1, 2]
    if len(deficit_events) > 0 and len(spin_eff_events) > 0:
        deficit_pm2 = np.array([e['deficit_per_m2'] for e in spin_eff_events])
        gr_pm2 = np.array([e['gr_deficit_per_m2'] for e in spin_eff_events])

        ax.scatter(q_se, deficit_pm2, c='steelblue', s=20, alpha=0.6,
                   edgecolors='k', linewidths=0.3, zorder=5,
                   label='Observed deficit/m2')

        # GR prediction curve
        q_fine = np.linspace(0.15, 1.0, 200)
        eta_fine = q_fine / (1 + q_fine)**2
        frac_fine = 0.0572 * eta_fine + 0.498 * eta_fine**2
        # deficit/m2 = frac * M_total / m2 = frac * (1+q)/q
        gr_dm2_fine = frac_fine * (1 + q_fine) / q_fine
        ax.plot(q_fine, gr_dm2_fine, 'r--', linewidth=2,
                label='GR non-spinning')

        # Binned medians
        q_centers2 = []
        dm2_medians = []
        for i in range(len(q_bin_edges) - 1):
            q_lo, q_hi = q_bin_edges[i], q_bin_edges[i + 1]
            sub = [e['deficit_per_m2'] for e in spin_eff_events
                   if q_lo <= e['q'] < q_hi]
            if len(sub) >= 2:
                q_centers2.append((q_lo + q_hi) / 2)
                dm2_medians.append(np.median(sub))
        ax.plot(q_centers2, dm2_medians, 'bo-', linewidth=2, markersize=8,
                markeredgecolor='k', markeredgewidth=0.5, zorder=10,
                label='Binned median (observed)')

        ax.set_xlabel("Mass ratio q = m2/m1", fontsize=11)
        ax.set_ylabel("Mass deficit per unit m2 (Msun/Msun)", fontsize=11)
        ax.set_title("Deficit per m2 vs Mass Ratio")
        ax.legend(fontsize=7.5)
        ax.grid(True, alpha=0.3)

    plt.suptitle("BH Merger Spin Analysis: Quark Star Hypothesis",
                 fontsize=13, fontweight='bold')
    plt.tight_layout()

    os.makedirs("results", exist_ok=True)
    outpath = "results/bh_merger_spin_analysis.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {outpath}")
    plt.close()

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print(f"\n{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}")
    print()
    print(f"  1. chi_eff vs m1:")
    print(f"     Spearman rho = {rho:+.3f}, p = {p_sp:.4f}")
    print(f"     |chi_eff| vs m1: rho = {rho_abs:+.3f}, p = {p_abs:.4f}")
    if len(in_gap) >= 3:
        print(f"     Sub-gap mean chi_eff: {sub_gap.mean():+.3f}"
              f"  vs  Gap mean: {in_gap.mean():+.3f}")
        print(f"     Mann-Whitney p = {p_mw:.4f}")
    print()
    if len(deficit_events) > 0:
        print(f"  2. Fractional mass deficit:")
        print(f"     Median: {np.median(frac_deficit):.3f}")
        print(f"     Correlation with m_total: rho={rho_d:+.3f}, p={p_d:.4f}")
    print()
    print(f"  3. Excess over non-spinning GR:")
    if len(excess_list) > 0:
        print(f"     Median excess: {np.median(excess_arr):+.3f}")
        print(f"     {np.sum(excess_arr > 0)}/{len(excess_arr)}"
              f" ({np.sum(excess_arr > 0)/len(excess_arr):.0%}) events"
              f" exceed non-spinning GR prediction")
    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
