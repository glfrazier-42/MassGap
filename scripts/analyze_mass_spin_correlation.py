#!/usr/bin/env python3
"""
Analyze mass-period correlation to test for "forbidden zone"

This tests the hypothesis that massive pulsars can only exist at fast spins,
creating an asymmetric distribution in M-P space.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def quadrant_analysis(masses, periods):
    """
    Analyze distribution in quadrants of M-P space
    
    Returns statistics for each quadrant
    """
    mass_median = np.median(masses)
    period_median = np.median(periods)
    
    # Define quadrants
    # Q1: Low mass, fast spin (low P)
    # Q2: High mass, fast spin (low P)  
    # Q3: High mass, slow spin (high P) <- FORBIDDEN ZONE
    # Q4: Low mass, slow spin (high P)
    
    q1 = (masses < mass_median) & (periods < period_median)
    q2 = (masses >= mass_median) & (periods < period_median)
    q3 = (masses >= mass_median) & (periods >= period_median)  # Forbidden?
    q4 = (masses < mass_median) & (periods >= period_median)
    
    return {
        'mass_median': mass_median,
        'period_median': period_median,
        'Q1_low_mass_fast': np.sum(q1),
        'Q2_high_mass_fast': np.sum(q2),
        'Q3_high_mass_slow': np.sum(q3),  # Should be depleted!
        'Q4_low_mass_slow': np.sum(q4),
    }

def main():
    # Load data
    data_csv = Path("data/processed/pulsar_mass_period_combined.csv")
    df = pd.read_csv(data_csv)
    
    masses = df['mass'].values
    periods = df['period'].values * 1000  # Convert to ms
    freqs = df['freq_hz'].values
    
    print("=" * 70)
    print("MASS-SPIN CORRELATION ANALYSIS")
    print("=" * 70)
    print(f"Dataset: {len(df)} pulsars")
    print()
    
    # Quadrant analysis
    quad_stats = quadrant_analysis(masses, periods)
    
    print("QUADRANT ANALYSIS (using medians as boundaries):")
    print(f"  Mass median:   {quad_stats['mass_median']:.3f} M☉")
    print(f"  Period median: {quad_stats['period_median']:.2f} ms")
    print()
    print("Distribution:")
    print(f"  Q1 (Low M,  Fast spin): {quad_stats['Q1_low_mass_fast']:2d} pulsars")
    print(f"  Q2 (High M, Fast spin): {quad_stats['Q2_high_mass_fast']:2d} pulsars")
    print(f"  Q3 (High M, Slow spin): {quad_stats['Q3_high_mass_slow']:2d} pulsars  <- FORBIDDEN ZONE?")
    print(f"  Q4 (Low M,  Slow spin): {quad_stats['Q4_low_mass_slow']:2d} pulsars")
    print()
    
    # Expected for independent distributions
    expected_per_quadrant = len(df) / 4
    print(f"Expected for independent normal distributions: {expected_per_quadrant:.1f} per quadrant")
    print()
    
    # Chi-square test for independence
    observed = np.array([
        quad_stats['Q1_low_mass_fast'],
        quad_stats['Q2_high_mass_fast'],
        quad_stats['Q3_high_mass_slow'],
        quad_stats['Q4_low_mass_slow']
    ])
    expected = np.full(4, expected_per_quadrant)
    chi_square = np.sum((observed - expected)**2 / expected)
    
    print(f"Chi-square statistic: {chi_square:.2f}")
    print(f"(Higher values indicate stronger deviation from independence)")
    print()
    
    # Check the specific prediction: Q3 should be depleted
    q3_ratio = quad_stats['Q3_high_mass_slow'] / expected_per_quadrant
    print("FORBIDDEN ZONE TEST:")
    print(f"  Q3 observed: {quad_stats['Q3_high_mass_slow']}")
    print(f"  Q3 expected: {expected_per_quadrant:.1f}")
    print(f"  Ratio (obs/exp): {q3_ratio:.2f}")
    if q3_ratio < 0.5:
        print("  ✓ STRONG DEPLETION - Forbidden zone confirmed!")
    elif q3_ratio < 0.75:
        print("  ~ Moderate depletion - Some evidence for forbidden zone")
    else:
        print("  ✗ No significant depletion")
    print()
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Mass vs Period
    ax1 = axes[0]
    scatter = ax1.scatter(periods, masses, c=freqs, s=100, alpha=0.6, 
                         cmap='viridis', edgecolors='black', linewidth=0.5)
    
    # Add median lines
    ax1.axhline(quad_stats['mass_median'], color='red', linestyle='--', 
                alpha=0.5, label=f"Median M = {quad_stats['mass_median']:.2f} M☉")
    ax1.axvline(quad_stats['period_median'], color='red', linestyle='--', 
                alpha=0.5, label=f"Median P = {quad_stats['period_median']:.1f} ms")
    
    # Highlight forbidden zone
    ax1.axhspan(quad_stats['mass_median'], masses.max()+0.1, 
                xmin=0.5, xmax=1.0, alpha=0.1, color='red', 
                label='Predicted Forbidden Zone')
    
    # Annotate key pulsars
    key_pulsars = ['J0740+6620', 'J0348+0432', 'J1614-2230', 'J0437-4715']
    for kp in key_pulsars:
        match = df[df['PSRJ'] == kp]
        if len(match) > 0:
            row = match.iloc[0]
            ax1.annotate(kp, xy=(row['period']*1000, row['mass']), 
                        xytext=(10, 10), textcoords='offset points',
                        fontsize=8, alpha=0.7,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))
    
    ax1.set_xlabel('Period (ms)', fontsize=12)
    ax1.set_ylabel('Mass (M☉)', fontsize=12)
    ax1.set_title('Mass vs Period - Testing for Forbidden Zone', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    cbar1 = plt.colorbar(scatter, ax=ax1)
    cbar1.set_label('Spin Frequency (Hz)', fontsize=10)
    
    # Plot 2: Log scale to see structure better
    ax2 = axes[1]
    scatter2 = ax2.scatter(periods, masses, c=freqs, s=100, alpha=0.6,
                          cmap='viridis', edgecolors='black', linewidth=0.5)
    
    ax2.set_xscale('log')
    ax2.axhline(quad_stats['mass_median'], color='red', linestyle='--', alpha=0.5)
    ax2.axvline(quad_stats['period_median'], color='red', linestyle='--', alpha=0.5)
    
    ax2.set_xlabel('Period (ms) [log scale]', fontsize=12)
    ax2.set_ylabel('Mass (M☉)', fontsize=12)
    ax2.set_title('Mass vs Period (Log Scale)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')
    
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label('Spin Frequency (Hz)', fontsize=10)
    
    plt.tight_layout()
    
    # Save figure
    output_fig = Path("C:/cygwin64/home/glfra/BlackHoleExplosionMechanism/figures/mass_period_forbidden_zone.png")
    output_fig.parent.mkdir(exist_ok=True, parents=True)
    plt.savefig(output_fig, dpi=300, bbox_inches='tight')
    print(f"✓ Saved figure to: {output_fig}")
    print()
    
    # Additional analysis: Slowest pulsar at each mass range
    print("SLOWEST PULSAR AT EACH MASS RANGE:")
    mass_bins = [(1.0, 1.5), (1.5, 2.0), (2.0, 2.5), (2.5, 3.0)]
    for m_low, m_high in mass_bins:
        in_bin = df[(df['mass'] >= m_low) & (df['mass'] < m_high)]
        if len(in_bin) > 0:
            slowest = in_bin.loc[in_bin['period'].idxmax()]
            fastest = in_bin.loc[in_bin['period'].idxmin()]
            print(f"  {m_low:.1f}-{m_high:.1f} M☉: Slowest = {slowest['PSRJ']}, "
                  f"P = {slowest['period']*1000:.1f} ms, M = {slowest['mass']:.3f} M☉"
                  f"  Fastest = {fastest['PSRJ']}, "
                  f"P = {fastest['period']*1000:.1f} ms, M = {fastest['mass']:.3f} M☉")
        else:
            print(f"  {m_low:.1f}-{m_high:.1f} M☉: No pulsars")
    print()
    
    print("=" * 70)
    print("INTERPRETATION:")
    print("=" * 70)
    if q3_ratio < 0.5:
        print("✓ Strong evidence for forbidden zone!")
        print("  Massive pulsars are found ONLY at fast spins.")
        print("  This supports the collapse mechanism:")
        print("  → Massive pulsars slow down via magnetic dipole radiation")
        print("  → When they slow below critical Ω, they collapse to BH")
        print("  → We don't observe slow, massive pulsars")
    else:
        print("✗ Forbidden zone not clearly evident")
        print("  Distribution appears more symmetric")
        print("  May need larger sample or different boundaries")
    
    plt.show()

if __name__ == "__main__":
    main()
