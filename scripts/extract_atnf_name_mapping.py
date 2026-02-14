"""
Extract B name to J name mapping from ATNF psrcat.db

This creates a lookup table to convert old-style pulsar names
(B names, X-ray names, etc.) to modern J names for matching.
"""

import pandas as pd
from pathlib import Path

def parse_psrcat_name_mapping(psrcat_db_path):
    """
    Parse psrcat.db to extract PSRB -> PSRJ mappings
    
    Returns dict: {B_name: J_name}
    """
    
    name_mapping = {}
    current_psrj = None
    current_psrb = None
    
    with open(psrcat_db_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            # Check for pulsar separator
            if line.startswith('@'):
                # Save the mapping if we have both names
                if current_psrb and current_psrj:
                    name_mapping[current_psrb] = current_psrj
                # Reset for next pulsar
                current_psrj = None
                current_psrb = None
                continue
            
            # Parse parameter lines
            parts = line.split()
            if len(parts) >= 2:
                param_name = parts[0]
                param_value = parts[1]
                
                if param_name == 'PSRJ':
                    current_psrj = param_value
                elif param_name == 'PSRB':
                    current_psrb = param_value
    
    return name_mapping

def main():
    print("=" * 70)
    print("EXTRACTING B->J NAME MAPPING FROM ATNF CATALOG")
    print("=" * 70)
    print()
    
    # Path to ATNF catalog
    psrcat_db = Path("data/raw/psrcat.db")
    
    if not psrcat_db.exists():
        print(f"ERROR: {psrcat_db} not found")
        return
    
    print("Parsing psrcat.db...")
    name_mapping = parse_psrcat_name_mapping(psrcat_db)
    
    print(f"OK Found {len(name_mapping)} B->J name mappings")
    print()
    
    # Show some examples
    print("Example mappings:")
    for i, (b_name, j_name) in enumerate(list(name_mapping.items())[:15], 1):
        print(f"  {i:2}. {b_name:15s} -> {j_name}")
    print()
    
    # Check for our known pulsars with B names
    known_b_names = ['B1913+16', 'B1534+12', 'B1957+20', 'B2127+11C']
    print("Known pulsars with B names:")
    for b_name in known_b_names:
        if b_name in name_mapping:
            print(f"  OK {b_name:15s} -> {name_mapping[b_name]}")
        else:
            # Try without suffix
            base_name = b_name.rstrip('ABC')
            if base_name in name_mapping:
                print(f"  ~ {b_name:15s} -> {name_mapping[base_name]} (base name)")
            else:
                print(f"  MISS {b_name:15s} NOT FOUND")
    print()
    
    # Save to CSV
    output_csv = Path("data/atnf_name_mapping.csv")
    
    df = pd.DataFrame({
        'B_name': list(name_mapping.keys()),
        'J_name': list(name_mapping.values())
    })
    df.to_csv(output_csv, index=False)
    
    print(f"OK Saved mapping to: {output_csv}")
    print()
    print("Next: Use this mapping to add J names to stellar_collapse data")

if __name__ == "__main__":
    main()
