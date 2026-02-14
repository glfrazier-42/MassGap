"""
Population Survival Analysis v2 (no-matplotlib version)

Filters out recycled millisecond pulsars (P < 100 ms). Focuses on normal
pulsars where spindown predicts accumulation at longer periods.

Usage:
  python population_survival_v2.py            # hardcoded counts
  python population_survival_v2.py --from-csv  # regenerate from CSV

Requires: numpy, scipy
(Add --from-csv to also require pandas)
"""

import numpy as np
from scipy import stats
import argparse
from pathlib import Path


def load_from_csv(data_dir):
    """Load ATNF catalog and bin normal pulsars."""
    import pandas as pd
    atnf = pd.read_csv(data_dir / "atnf_pulsar_periods.csv")
    periods_ms = atnf['period'].values * 1000
    total = len(periods_ms)

    # Filter to normal pulsars
    normal = periods_ms[periods_ms >= 100.0]

    bin_edges = [100, 300, 1000, 3000, 10000]
    counts = []
    for i in range(len(bin_edges) - 1):
        c = int(np.sum((normal >= bin_edges[i]) & (normal < bin_edges[i+1])))
        counts.append(c)
    return bin_edges, counts, total, len(normal)


def get_hardcoded():
    """ATNF catalog counts for normal pulsars (P >= 100 ms)."""
    # From previous analysis of 4268 total pulsars
    bin_edges = [100, 300, 1000, 3000, 10000]
    counts = [545, 1660, 929, 153]
    return bin_edges, counts, 4268, sum(counts)


def run_analysis(bin_edges, counts, total_all, total_normal):
    bin_names = [f"{bin_edges[i]}-{bin_edges[i+1]} ms"
                 for i in range(len(counts))]
    log_widths = [np.log10(bin_edges[i+1] / bin_edges[i])
                  for i in range(len(counts))]
    densities = [counts[i] / log_widths[i] for i in range(len(counts))]

    print(f"Total pulsars in catalog: {total_all}")
    print(f"Normal pulsars (P >= 100 ms): {total_normal}")
    print()

    # --------------------------------------------------------
    # Density per decade
    # --------------------------------------------------------
    print("=" * 65)
    print("DENSITY PER DECADE (normal pulsars, P >= 100 ms)")
    print("=" * 65)
    print()
    print("Spindown predicts this should INCREASE with period.")
    print("At minimum it should stay FLAT. A decrease = pulsars vanishing.")
    print()
    print(f"{'Bin':>15s}  {'Count':>6s}  {'Pulsars/decade':>14s}  {'Trend':>12s}")
    print("-" * 55)

    for i in range(len(counts)):
        if i == 0:
            trend = ""
        else:
            change = densities[i] / densities[i-1]
            if change > 1.05:
                trend = f"UP x{change:.1f}"
            elif change < 0.95:
                trend = f"DOWN x{change:.2f}"
            else:
                trend = "~flat"
        print(f"{bin_names[i]:>15s}  {counts[i]:6d}  {densities[i]:14.0f}"
              f"  {trend:>12s}")

    # --------------------------------------------------------
    # Adjacent-bin deficit test
    # --------------------------------------------------------
    print()
    print("=" * 65)
    print("ADJACENT-BIN DEFICIT TEST")
    print("=" * 65)
    print()
    print("Null: density per decade non-decreasing.")
    print("Min expected = prev_count * (next_log_width / prev_log_width)")
    print()

    print(f"{'Transition':>35s}  {'Obs':>5s}  {'MinExp':>6s}  "
          f"{'Deficit':>7s}  {'Signif':>10s}")
    print("-" * 70)

    for i in range(len(counts) - 1):
        expected_min = counts[i] * (log_widths[i+1] / log_widths[i])

        if counts[i+1] < expected_min:
            deficit = expected_min - counts[i+1]
            p_val = stats.poisson.cdf(counts[i+1], expected_min)
            if p_val > 0 and p_val < 0.5:
                sigma = -stats.norm.ppf(p_val)
                sig_str = f"{sigma:.1f}sigma"
            elif p_val == 0:
                z = (expected_min - counts[i+1]) / np.sqrt(expected_min)
                sig_str = f">{z:.0f}sigma"
            else:
                sig_str = "n/s"
            print(f"{bin_names[i]:>15s} -> {bin_names[i+1]:>15s}  "
                  f"{counts[i+1]:5d}  {expected_min:6.0f}  {deficit:+7.0f}  "
                  f"{sig_str:>10s}")
        else:
            excess = counts[i+1] - expected_min
            print(f"{bin_names[i]:>15s} -> {bin_names[i+1]:>15s}  "
                  f"{counts[i+1]:5d}  {expected_min:6.0f}  +{int(excess):6d}  "
                  f"{'OK':>10s}")

    # --------------------------------------------------------
    # Sensitivity to detection bias
    # --------------------------------------------------------
    print()
    print("=" * 65)
    print("SENSITIVITY TO DETECTION BIAS")
    print("=" * 65)
    print()

    for i in range(len(counts) - 1):
        expected_min = counts[i] * (log_widths[i+1] / log_widths[i])
        if counts[i+1] < expected_min:
            n_obs = counts[i+1]
            print(f"Transition: {bin_names[i]} -> {bin_names[i+1]}")
            print(f"Observed: {n_obs}, Minimum expected: {expected_min:.0f}")
            print()
            for bias_pct in [0, 10, 20, 30, 40, 50, 60]:
                adj = expected_min * (1 - bias_pct / 100)
                if n_obs < adj:
                    p = stats.poisson.cdf(n_obs, adj)
                    if p > 0 and p < 0.5:
                        s = -stats.norm.ppf(p)
                        print(f"  {bias_pct:2d}% bias: exp={adj:.0f}, "
                              f"deficit={adj - n_obs:.0f}, {s:.1f}sigma")
                    elif p == 0:
                        z = (adj - n_obs) / np.sqrt(adj)
                        print(f"  {bias_pct:2d}% bias: exp={adj:.0f}, "
                              f"deficit={adj - n_obs:.0f}, >{z:.0f}sigma")
                else:
                    print(f"  {bias_pct:2d}% bias: exp={adj:.0f}, "
                          f"NO deficit")
            print()

    # --------------------------------------------------------
    # Birth period note
    # --------------------------------------------------------
    print("=" * 65)
    print("NOTE: Birth period distribution makes deficit MORE significant")
    print("=" * 65)
    print()
    print("Many normal pulsars are born with P0 ~ 10-100 ms and enter")
    print("the slow bins directly. This adds to their expected population,")
    print("making the flat-density null hypothesis CONSERVATIVE.")
    print("The true deficit is larger than shown above.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Population survival analysis (normal pulsars only)')
    parser.add_argument('--from-csv', action='store_true',
                        help='Regenerate from atnf_pulsar_periods.csv')
    args = parser.parse_args()

    if args.from_csv:
        data_dir = Path("data")
        edges, cts, tot, norm = load_from_csv(data_dir)
    else:
        edges, cts, tot, norm = get_hardcoded()

    run_analysis(edges, cts, tot, norm)
