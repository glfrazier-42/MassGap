"""
Extract cold, beta-equilibrium EOS table from stellarcollapse.org HDF5 files.

This script:
1. Reads the HDF5 file structure
2. Extracts the T=0, beta-equilibrium slice
3. Converts to simple 3-column format: rho, pressure, energy_density
4. Outputs in CGS units suitable for TOV solver
"""

import h5py
import numpy as np
import sys
import os

def examine_hdf5_structure(filename):
    """Print the structure of the HDF5 file to understand its contents."""
    print(f"\n{'='*60}")
    print(f"Examining HDF5 file: {filename}")
    print(f"{'='*60}\n")
    
    with h5py.File(filename, 'r') as f:
        print("Top-level datasets and groups:")
        print("-" * 40)
        for key in f.keys():
            item = f[key]
            if isinstance(item, h5py.Dataset):
                print(f"  Dataset: {key}")
                print(f"    Shape: {item.shape}")
                print(f"    Dtype: {item.dtype}")
                if item.shape == ():  # scalar
                    print(f"    Value: {item[()]}")
                elif len(item.shape) == 1 and item.shape[0] < 10:
                    print(f"    Values: {item[:]}")
                print()
            elif isinstance(item, h5py.Group):
                print(f"  Group: {key}")
                print()

def extract_cold_eos(filename, output_file="eos_table.dat"):
    """
    Extract cold (T=0), beta-equilibrium EOS from stellarcollapse HDF5 table.
    
    Parameters:
    -----------
    filename : str
        Path to the HDF5 file
    output_file : str
        Output filename for the extracted table
    """
    
    print(f"\nExtracting EOS from {filename}...")
    
    with h5py.File(filename, 'r') as f:
        # Read the grid dimensions
        # These tables typically have:
        # - logrho: log10(density in g/cm³)
        # - logtemp: log10(temperature in MeV)
        # - ye: electron fraction
        
        logrho = f['logrho'][:]  # log10(rho)
        logtemp = f['logtemp'][:]  # log10(T)
        ye = f['ye'][:]  # electron fraction
        
        print(f"\nGrid dimensions:")
        print(f"  Density points: {len(logrho)} (log10(rho) from {logrho[0]:.2f} to {logrho[-1]:.2f})")
        print(f"  Temperature points: {len(logtemp)} (log10(T/MeV) from {logtemp[0]:.2f} to {logtemp[-1]:.2f})")
        print(f"  Ye points: {len(ye)} (from {ye[0]:.3f} to {ye[-1]:.3f})")
        
        # Read thermodynamic quantities
        # These are typically 3D arrays: [nrho, ntemp, nye]
        # We need to find the right keys - let's check what's available
        
        print(f"\nAvailable thermodynamic quantities:")
        thermo_keys = [k for k in f.keys() if k not in ['logrho', 'logtemp', 'ye']]
        for key in thermo_keys[:10]:  # Show first 10
            print(f"  {key}: shape {f[key].shape}")
        
        # Typically: 'logpress', 'logenergy', or 'press', 'energy'
        # Let's try to find them
        
        if 'logpress' in f:
            logpress = f['logpress'][:]  # log10(pressure in dyne/cm²)
        elif 'press' in f:
            press_raw = f['press'][:]
            logpress = np.log10(press_raw)
        else:
            raise KeyError("Could not find pressure data (logpress or press)")
            
        if 'logenergy' in f:
            logenergy = f['logenergy'][:]  # log10(energy density in erg/cm³)
        elif 'energy' in f:
            energy_raw = f['energy'][:]
            logenergy = np.log10(energy_raw)
        else:
            raise KeyError("Could not find energy density data (logenergy or energy)")
        
        # For cold beta-equilibrium:
        # - Take coldest temperature (index 0, usually T ~ 0.01 MeV)
        # - Find Ye corresponding to beta-equilibrium (need to check if provided)
        
        itemp = 0  # Coldest temperature
        
        print(f"\nExtracting data at T = {10**logtemp[itemp]:.4f} MeV")
        print(f"Array shapes: logpress={logpress.shape}, logenergy={logenergy.shape}")
        
        # Figure out the actual array structure
        # Could be [nrho, ntemp, nye] or [ntemp, nrho, nye] or other
        rho = 10**logrho  # Convert to linear g/cm³
        
        # Try to determine dimension order by checking shapes
        nrho, ntemp, nye = len(logrho), len(logtemp), len(ye)
        
        if logpress.shape == (nrho, ntemp, nye):
            # [nrho, ntemp, nye]
            press_2d = 10**logpress[:, itemp, :]  # [nrho, nye]
            energy_2d = 10**logenergy[:, itemp, :]  # [nrho, nye]
            dim_order = "[nrho, ntemp, nye]"
        elif logpress.shape == (ntemp, nrho, nye):
            # [ntemp, nrho, nye]
            press_2d = 10**logpress[itemp, :, :]  # [nrho, nye]
            energy_2d = 10**logenergy[itemp, :, :]  # [nrho, nye]
            dim_order = "[ntemp, nrho, nye]"
        elif logpress.shape == (nye, ntemp, nrho):
            # [nye, ntemp, nrho]
            press_2d = 10**logpress[:, itemp, :].T  # [nrho, nye]
            energy_2d = 10**logenergy[:, itemp, :].T  # [nrho, nye]
            dim_order = "[nye, ntemp, nrho]"
        else:
            raise ValueError(f"Unexpected array shape: {logpress.shape} vs expected ({nrho}, {ntemp}, {nye})")
        
        print(f"Detected dimension order: {dim_order}")
        print(f"2D slice shape: {press_2d.shape}")
        
        # Find Ye closest to 0.1 (typical for cold beta-equilibrium neutron star matter)
        target_ye = 0.1
        iye = np.argmin(np.abs(ye - target_ye))
        actual_ye = ye[iye]
        
        print(f"Using Ye = {actual_ye:.4f} (closest to {target_ye})")
        
        # Extract 1D arrays at this Ye
        press_1d = press_2d[:, iye]
        energy_1d = energy_2d[:, iye]
        
        print(f"1D arrays shape: rho={rho.shape}, press={press_1d.shape}, energy={energy_1d.shape}")
        
        # Filter out any bad values (NaN, inf, or negative)
        valid = np.isfinite(rho) & np.isfinite(press_1d) & np.isfinite(energy_1d)
        valid &= (press_1d > 0) & (energy_1d > 0)
        
        rho_clean = rho[valid]
        press_clean = press_1d[valid]
        energy_clean = energy_1d[valid]
        
        print(f"\nExtracted {len(rho_clean)} valid data points")
        print(f"  Density range: {rho_clean[0]:.3e} to {rho_clean[-1]:.3e} g/cm^3")
        print(f"  Pressure range: {press_clean[0]:.3e} to {press_clean[-1]:.3e} dyne/cm^2")
        print(f"  Energy density range: {energy_clean[0]:.3e} to {energy_clean[-1]:.3e} erg/cm^3")
        
        # Save to file
        header = f"""# EOS Table extracted from {os.path.basename(filename)}
# Source: stellarcollapse.org
# Temperature: T = {10**logtemp[itemp]:.4f} MeV (cold)
# Electron fraction: Ye = {actual_ye:.4f}
# Units: CGS (g/cm³, dyne/cm², erg/cm³)
# Columns: rho (g/cm³), pressure (dyne/cm²), energy_density (erg/cm³)
"""
        
        with open(output_file, 'w') as out:
            out.write(header)
            for r, p, e in zip(rho_clean, press_clean, energy_clean):
                out.write(f"{r:.6e}  {p:.6e}  {e:.6e}\n")
        
        print(f"\nOK Saved to: {output_file}")
        
        return rho_clean, press_clean, energy_clean


def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_stellarcollapse_eos.py <hdf5_file> [output_file]")
        print("\nExample:")
        print("  python extract_stellarcollapse_eos.py LS220_234r_136t_50y_analmu_20091212_SVNr26.h5")
        sys.exit(1)
    
    hdf5_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "eos_table.dat"
    
    if not os.path.exists(hdf5_file):
        print(f"Error: File not found: {hdf5_file}")
        sys.exit(1)
    
    # First examine the file structure
    examine_hdf5_structure(hdf5_file)
    
    # Then extract the EOS
    try:
        extract_cold_eos(hdf5_file, output_file)
    except Exception as e:
        print(f"\nERROR Error extracting EOS: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
