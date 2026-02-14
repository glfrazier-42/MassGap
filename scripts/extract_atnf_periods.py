"""
Extract pulsar periods from ATNF psrcat.db

This extracts PSRJ, PSRB, P0 (or F0 to calculate period) from the
full ATNF catalog to get periods for all pulsars.
"""

import pandas as pd
from pathlib import Path

def parse_atnf_periods(psrcat_db_path):
    """
    Parse psrcat.db to extract pulsar names and periods
    
    Returns DataFrame with: PSRJ, PSRB, period
    """
    
    pulsars = []
    current_entry = {}
    
    with open(psrcat_db_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            # Check for pulsar separator
            if line.startswith('@'):
                # Save current entry if it has a name and period
                if 'PSRJ' in current_entry and 'period' in current_entry:
                    pulsars.append(current_entry.copy())
                # Reset for next pulsar
                current_entry = {}
                continue
            
            # Parse parameter lines
            parts = line.split()
            if len(parts) >= 2:
                param_name = parts[0]
                param_value = parts[1]
                
                if param_name == 'PSRJ':
                    current_entry['PSRJ'] = param_value
                elif param_name == 'PSRB':
                    current_entry['PSRB'] = param_value
                elif param_name == 'P0':
                    # Period in seconds
                    try:
                        current_entry['period'] = float(param_value)
                    except ValueError:
                        pass
                elif param_name == 'F0' and 'period' not in current_entry:
                    # Frequency in Hz, convert to period
                    try:
                        freq = float(param_value)
                        if freq > 0:
                            current_entry['period'] = 1.0 / freq
                    except ValueError:
                        pass
    
    # Convert to DataFrame
    df = pd.DataFrame(pulsars)
    
    # Fill missing PSRB with empty string
    if 'PSRB' not in df.columns:
        df['PSRB'] = ''
    else:
        df['PSRB'] = df['PSRB'].fillna('')
    
    return df

def main():
    print("=" * 70)
    print("EXTRACTING PERIODS FROM ATNF CATALOG")
    print("=" * 70)
    print()
    
    # Path to ATNF catalog
    psrcat_db = Path("data/raw/psrcat.db")
    
    if not psrcat_db.exists():
        print(f"ERROR: {psrcat_db} not found")
        return
    
    print("Parsing psrcat.db for periods...")
    periods_df = parse_atnf_periods(psrcat_db)
    
    print(f"OK Found {len(periods_df)} pulsars with periods")
    print()
    
    # Show some examples
    print("First 15 pulsars:")
    print(f"{'PSRJ':<20s} {'PSRB':<15s} {'Period(s)':<15s} {'Period(ms)':<15s}")
    print("-" * 70)
    for _, row in periods_df.head(15).iterrows():
        period_s = row['period']
        period_ms = period_s * 1000
        psrb = row.get('PSRB', '')
        print(f"{row['PSRJ']:<20s} {psrb:<15s} {period_s:<15.6f} {period_ms:<15.2f}")
    print()
    
    # Check for our key pulsars
    known_pulsars = ['J0437-4715', 'J1614-2230', 'J0348+0432', 'J0740+6620']
    print("Key massive pulsars:")
    for kp in known_pulsars:
        match = periods_df[periods_df['PSRJ'] == kp]
        if len(match) > 0:
            period_ms = match.iloc[0]['period'] * 1000
            print(f"  OK {kp:<15s} P = {period_ms:8.3f} ms")
        else:
            print(f"  MISS {kp:<15s} NOT FOUND")
    print()
    
    # Statistics
    periods = periods_df['period'].values
    print(f"Period range: {periods.min()*1000:.3f} - {periods.max()*1000:.1f} ms")
    print(f"Millisecond pulsars (P < 10 ms): {sum(periods < 0.01)}")
    print()
    
    # Save to CSV
    output_csv = Path("data/atnf_pulsar_periods.csv")
    periods_df.to_csv(output_csv, index=False)
    
    print(f"OK Saved to: {output_csv}")
    print()
    print("Next: Merge with stellar_collapse masses using both PSRJ and PSRB")

if __name__ == "__main__":
    main()
