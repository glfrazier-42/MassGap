"""
Extract cold beta-equilibrium slices from multiple stellarcollapse.org EOS tables

This script processes all downloaded EOS tables and extracts the T=0,
beta-equilibrium slice needed for TOV integration.

Usage:
    python extract_all_eos.py
"""

import os
import h5py
import numpy as np
from pathlib import Path
from extract_stellarcollapse_eos import extract_cold_eos

# EOS tables to process
EOS_TABLES = {
    'LS180': {
        'file': 'LS180_234r_136t_50y_analmu_20091212_SVNr26.h5',
        'description': 'Lattimer-Swesty, K=180 MeV (SOFT)',
        'reference': 'Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991)'
    },
    'LS220': {
        'file': 'LS220_234r_136t_50y_analmu_20091212_SVNr26.h5',
        'description': 'Lattimer-Swesty, K=220 MeV (MEDIUM)',
        'reference': 'Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991)'
    },
    'LS375': {
        'file': 'LS375_234r_136t_50y_analmu_20091212_SVNr26.h5',
        'description': 'Lattimer-Swesty, K=375 MeV (VERY STIFF)',
        'reference': 'Lattimer & Swesty, Nucl. Phys. A 535, 331 (1991)'
    },
    'SFHo': {
        'file': 'Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5',
        'description': 'Steiner-Fischer-Hempel-Oertel (STIFF)',
        'reference': 'Steiner et al., ApJ 774, 17 (2013)'
    },
    'DD2': {
        'file': 'Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5',
        'description': 'Density-Dependent RMF (STIFF)',
        'reference': 'Hempel & Schaffner-Bielich, Nucl. Phys. A 837, 210 (2010)'
    },
    'HShen': {
        'file': 'HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5',
        'description': 'H. Shen et al. (2011)',
        'reference': 'Shen et al., ApJS 197, 20 (2011)'
    }
}

def extract_cold_slice(h5_file, output_file):
    """
    Extract T=0, beta-equilibrium slice from stellarcollapse.org HDF5 table.
    
    Parameters:
        h5_file: Path to input HDF5 file
        output_file: Path to output text file
    
    Returns:
        Number of points extracted
    """
    print(f"  Reading: {h5_file.name}")
    
    with h5py.File(h5_file, 'r') as f:
        # Load density grid
        logrho = np.array(f['logrho'])  # log10(density in g/cm³)
        rho = 10**logrho  # Convert to linear
        
        # Load temperature grid  
        logtemp = np.array(f['logtemp'])  # log10(temperature in MeV)
        
        # Load electron fraction grid
        ye = np.array(f['ye'])
        
        # Find indices for T≈0 and beta equilibrium
        # T=0 is typically the first temperature point
        itemp = 0  
        
        # Load thermodynamic quantities (3D arrays: [rho, temp, ye])
        logpress = np.array(f['logpress'])  # log10(pressure in dyne/cm²)
        energy = np.array(f['energy'])  # Energy per baryon (MeV)
        
        # For beta equilibrium, we want Ye that minimizes energy at each density
        # (at T=0)
        print(f"  Finding beta equilibrium...")
        
        # Arrays to store result
        rho_eq = []
        press_eq = []
        energy_eq = []
        ye_eq = []
        
        for irho in range(len(logrho)):
            # Energy at this density for all Ye (at T=0)
            energy_slice = energy[irho, itemp, :]
            
            # Find Ye that minimizes energy (beta equilibrium)
            iye_min = np.argmin(energy_slice)
            
            rho_eq.append(rho[irho])
            press_eq.append(10**logpress[irho, itemp, iye_min])
            energy_eq.append(energy_slice[iye_min])
            ye_eq.append(ye[iye_min])
        
        # Convert to arrays
        rho_eq = np.array(rho_eq)
        press_eq = np.array(press_eq)
        energy_eq = np.array(energy_eq)
        ye_eq = np.array(ye_eq)
        
        # Convert energy per baryon (MeV) to energy density (erg/cm³)
        # ε = ρ * (1 + E/c²) where E is energy per baryon
        # E is in MeV, need to convert to erg
        MeV_to_erg = 1.60218e-6  # MeV to erg
        c_light = 2.99792458e10  # cm/s
        m_baryon = 1.66054e-24  # g (average nucleon mass)
        
        # energy_eq is in MeV per baryon
        # Convert to dimensionless: E/(m_b c²)
        # m_b c² ≈ 931.5 MeV
        m_b_c2_MeV = 931.5
        epsilon_eq = rho_eq * c_light**2 * (1.0 + energy_eq / m_b_c2_MeV)
        
        # Write output file
        print(f"  Writing: {output_file.name}")
        header = (
            f"Cold (T=0) beta-equilibrium slice\n"
            f"Extracted from: {h5_file.name}\n"
            f"Columns: density(g/cm³), pressure(dyne/cm²), "
            f"energy_density(erg/cm³), Ye\n"
            f"Number of points: {len(rho_eq)}\n"
        )
        
        data = np.column_stack([rho_eq, press_eq, epsilon_eq, ye_eq])
        np.savetxt(output_file, data, 
                   header=header,
                   fmt='%.10e',
                   comments='# ')
        
        return len(rho_eq)


def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    eos_dir = base_dir / 'data' / 'eos_tables'
    output_dir = eos_dir
    
    print("="*70)
    print("EXTRACTING COLD BETA-EQUILIBRIUM SLICES")
    print("="*70)
    print()
    
    # Process each EOS
    results = {}
    
    for name, info in EOS_TABLES.items():
        print(f"{name}: {info['description']}")
        print(f"  Reference: {info['reference']}")
        
        # Check if input file exists
        input_file = eos_dir / info['file']
        if not input_file.exists():
            print(f"  WARN Skipping - file not found: {info['file']}")
            print()
            continue
        
        # Output file name
        output_file = output_dir / f"{name}_cold_betaeq.txt"
        
        # Check if already extracted
        if output_file.exists():
            print(f"  OK Already extracted: {output_file.name}")
            # Count lines to get number of points
            with open(output_file) as f:
                npoints = sum(1 for line in f if not line.startswith('#'))
            results[name] = {
                'status': 'exists',
                'points': npoints,
                'file': output_file
            }
            print()
            continue
        
        # Extract
        try:
            # npoints = extract_cold_slice(input_file, output_file)
            rho_clean, press_clean, energy_clean = extract_cold_eos(input_file, output_file)
            npoints = len(rho_clean)
            print(f"  OK Extracted {npoints} points")
            results[name] = {
                'status': 'extracted',
                'points': npoints,
                'file': output_file
            }
        except Exception as e:
            print(f"  FAIL Error: {e}")
            results[name] = {
                'status': 'error',
                'error': str(e)
            }
        
        print()
    
    # Summary
    print("="*70)
    print("EXTRACTION SUMMARY")
    print("="*70)
    print()
    
    for name, info in EOS_TABLES.items():
        if name in results:
            status = results[name]
            if status['status'] in ['extracted', 'exists']:
                print(f"OK {name}: {status['points']} points -> {status['file'].name}")
            else:
                print(f"FAIL {name}: {status.get('error', 'Unknown error')}")
        else:
            print(f"WARN {name}: Not processed")
    
    print()
    print("Next step:")
    print("  python tests/test_multi_eos_robustness.py")
    print()


if __name__ == '__main__':
    main()
