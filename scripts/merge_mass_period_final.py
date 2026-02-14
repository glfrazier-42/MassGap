#!/usr/bin/env python3
"""
Merge stellar collapse masses with ATNF periods

Uses multiple matching strategies:
1. Direct J name match
2. B name to J name mapping
3. Manual additions for known aliases
"""

import pandas as pd
import numpy as np
from pathlib import Path

def normalize_name(name):
    """Normalize pulsar names for matching"""
    name = str(name).strip().upper()
    # Remove spaces
    name = name.replace(' ', '')
    return name

def main():
    print("=" * 70)
    print("MERGING STELLAR COLLAPSE MASSES WITH ATNF PERIODS")
    print("=" * 70)
    print()
    
    # Load data
    masses_csv = Path("C:/cygwin64/home/glfra/BlackHoleExplosionMechanism/data/processed/stellar_collapse_ns_masses.csv")
    periods_csv = Path("C:/cygwin64/home/glfra/BlackHoleExplosionMechanism/data/processed/atnf_pulsar_periods.csv")
    b_to_j_csv = Path("C:/cygwin64/home/glfra/BlackHoleExplosionMechanism/data/processed/atnf_name_mapping.csv")
    
    masses_df = pd.read_csv(masses_csv)
    periods_df = pd.read_csv(periods_csv)
    b_to_j_df = pd.read_csv(b_to_j_csv)
    
    # Remove companions
    masses_df = masses_df[~masses_df['name'].str.contains('companion', case=False, na=False)]
    
    print(f"Stellar collapse masses: {len(masses_df)}")
    print(f"ATNF periods:            {len(periods_df)}")
    print(f"B→J mappings:            {len(b_to_j_df)}")
    print()
    
    # Create B→J lookup
    b_to_j = dict(zip(b_to_j_df['B_name'], b_to_j_df['J_name']))
    
    # Strategy 1: Direct J name match
    print("Strategy 1: Direct J name matching...")
    masses_df['name_norm'] = masses_df['name'].apply(normalize_name)
    periods_df['PSRJ_norm'] = periods_df['PSRJ'].apply(normalize_name)
    
    direct_match = masses_df.merge(
        periods_df,
        left_on='name_norm',
        right_on='PSRJ_norm',
        how='inner'
    )
    print(f"  ✓ Matched {len(direct_match)} pulsars by J name")
    
    # Strategy 2: B name to J name match
    print("Strategy 2: B name → J name matching...")
    
    # Add J name column to masses for B names
    masses_df['J_name_from_B'] = masses_df['name'].map(b_to_j)
    masses_with_j = masses_df[masses_df['J_name_from_B'].notna()].copy()
    masses_with_j['J_name_norm'] = masses_with_j['J_name_from_B'].apply(normalize_name)
    
    b_match = masses_with_j.merge(
        periods_df,
        left_on='J_name_norm',
        right_on='PSRJ_norm',
        how='inner'
    )
    print(f"  ✓ Matched {len(b_match)} pulsars by B→J mapping")
    
    # Combine matches (avoid duplicates)
    all_matches = pd.concat([direct_match, b_match], ignore_index=True)
    all_matches = all_matches.drop_duplicates(subset=['name'])
    
    print()
    print(f"TOTAL MATCHED: {len(all_matches)} pulsars")
    print()
    
    if len(all_matches) == 0:
        print("ERROR: No matches found!")
        return
    
    # Calculate derived quantities
    all_matches['omega'] = 2 * np.pi / all_matches['period']
    all_matches['freq_hz'] = 1.0 / all_matches['period']
    
    # Select output columns
    output_df = all_matches[[
        'name', 'PSRJ', 'mass', 'error_plus', 'error_minus', 'error_symmetric',
        'period', 'omega', 'freq_hz', 'category', 'reference'
    ]].copy()
    
    # Sort by mass
    output_df = output_df.sort_values('mass', ascending=False)
    
    # Display results
    print("=" * 70)
    print("MATCHED PULSARS:")
    print("=" * 70)
    print()
    print(f"{'Original Name':<20s} {'PSRJ':<15s} {'Mass(M☉)':<12s} {'Period(ms)':<12s} {'Category':<12s}")
    print("-" * 70)
    
    for _, row in output_df.iterrows():
        if pd.notna(row['error_symmetric']):
            err_str = f"±{row['error_symmetric']:.3f}"
        else:
            err_str = f"+{row['error_plus']:.3f}/-{row['error_minus']:.3f}"
        mass_str = f"{row['mass']:.3f}{err_str}"
        period_ms = row['period'] * 1000
        print(f"{row['name']:<20s} {row['PSRJ']:<15s} {mass_str:<12s} {period_ms:<12.2f} {row['category']:<12s}")
    
    print()
    
    # Check key pulsars
    known_massive = ['J0740+6620', 'J0348+0432', 'J1614-2230', 'J0437-4715']
    print("Key massive pulsars:")
    for km in known_massive:
        match = output_df[output_df['PSRJ'].str.contains(km, case=False, na=False)]
        if len(match) > 0:
            row = match.iloc[0]
            print(f"  ✓ {row['PSRJ']:<15s} M = {row['mass']:.3f} M☉, P = {row['period']*1000:.2f} ms")
        else:
            print(f"  ✗ {km:<15s} NOT FOUND")
    print()
    
    # Statistics
    masses = output_df['mass'].values
    periods = output_df['period'].values
    
    print("=" * 70)
    print("STATISTICS:")
    print("=" * 70)
    print(f"Total matched: {len(output_df)}")
    print(f"Mass range:   {masses.min():.3f} - {masses.max():.3f} M☉")
    print(f"Period range: {periods.min()*1000:.2f} - {periods.max()*1000:.2f} ms")
    print(f"M > 2.0 M☉:   {sum(masses > 2.0)}")
    print(f"P < 10 ms:    {sum(periods < 0.01)}")
    print()
    
    # By category
    print("By category:")
    for cat in output_df['category'].unique():
        count = sum(output_df['category'] == cat)
        print(f"  {cat:<20s}: {count:3d} pulsars")
    print()
    
    # Save
    output_csv = Path("C:/cygwin64/home/glfra/BlackHoleExplosionMechanism/data/processed/pulsar_mass_period_combined.csv")
    output_df.to_csv(output_csv, index=False)
    print(f"✓ Saved to: {output_csv}")
    print()
    
    # Show what we couldn't match
    print("=" * 70)
    print("UNMATCHED PULSARS (from stellar_collapse):")
    print("=" * 70)
    matched_names = set(output_df['name'])
    unmatched = masses_df[~masses_df['name'].isin(matched_names)]
    print(f"Could not match {len(unmatched)} pulsars:")
    for _, row in unmatched.head(20).iterrows():
        print(f"  {row['name']:<30s} M = {row['mass']:.3f} M☉  [{row['category']}]")
    if len(unmatched) > 20:
        print(f"  ... and {len(unmatched) - 20} more")
    print()
    print("(These likely have X-ray or common names, not J/B names)")

if __name__ == "__main__":
    main()
