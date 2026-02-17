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

    Beta-equilibrium is found by minimizing ``logenergy`` (specific
    internal energy) over the Ye grid at each density, at the lowest
    temperature slice.  The stellarcollapse.org ``logenergy`` key stores
    ``log10(eps_specific + energy_shift)`` in erg/g; the constant shift
    does not affect which Ye is the minimum.  We subtract the shift to
    recover clean eps_specific.

    Output columns (CGS):
        rho          rest-mass density         g/cm^3
        pressure     total pressure            dyne/cm^2
        eps_specific specific internal energy   erg/g

    The TOV consumer (NeutronStarEOS.eos_function) reconstructs total
    energy density as  rho_total = rho * (1 + eps_specific / c^2).

    Parameters
    ----------
    filename : str or Path
        Path to the HDF5 file.
    output_file : str or Path
        Output filename for the extracted table.
    """

    print(f"\nExtracting EOS from {filename}...")

    with h5py.File(str(filename), 'r') as f:
        # ---- grid axes ----
        logrho  = f['logrho'][:]    # log10(rho / [g/cm^3])
        logtemp = f['logtemp'][:]   # log10(T   / [MeV])
        ye      = f['ye'][:]        # electron fraction

        nrho, ntemp, nye = len(logrho), len(logtemp), len(ye)

        print(f"\nGrid dimensions:")
        print(f"  Density points: {nrho} (log10 rho from {logrho[0]:.2f} to {logrho[-1]:.2f})")
        print(f"  Temperature points: {ntemp} (log10 T/MeV from {logtemp[0]:.2f} to {logtemp[-1]:.2f})")
        print(f"  Ye points: {nye} (from {ye[0]:.3f} to {ye[-1]:.3f})")

        # ---- thermodynamic arrays ----
        logpress_3d  = f['logpress'][:]       # log10(P / [dyne/cm^2])
        logenergy_3d = f['logenergy'][:]      # log10((eps_specific + shift) / [erg/g])
        energy_shift = float(f['energy_shift'][()])  # erg/g

        print(f"  energy_shift = {energy_shift:.6e} erg/g")

        # ---- detect dimension ordering ----
        # All stellarcollapse.org files use [nye, ntemp, nrho]
        itemp = 0  # coldest temperature slice

        if logpress_3d.shape == (nrho, ntemp, nye):
            dim_order = "[nrho, ntemp, nye]"
            logpress_2d  = logpress_3d[:, itemp, :]    # [nrho, nye]
            logenergy_2d = logenergy_3d[:, itemp, :]
        elif logpress_3d.shape == (ntemp, nrho, nye):
            dim_order = "[ntemp, nrho, nye]"
            logpress_2d  = logpress_3d[itemp, :, :]
            logenergy_2d = logenergy_3d[itemp, :, :]
        elif logpress_3d.shape == (nye, ntemp, nrho):
            dim_order = "[nye, ntemp, nrho]"
            logpress_2d  = logpress_3d[:, itemp, :].T  # [nrho, nye]
            logenergy_2d = logenergy_3d[:, itemp, :].T
        else:
            raise ValueError(
                f"Unexpected array shape: {logpress_3d.shape} "
                f"vs expected ({nrho}, {ntemp}, {nye})")

        print(f"\nExtracting at T = {10**logtemp[itemp]:.4f} MeV")
        print(f"Detected dimension order: {dim_order}")

        # ---- find beta-equilibrium Ye at each density ----
        # Minimise logenergy over Ye.  Since log10 is monotonic and
        # energy_shift is a constant, argmin is identical to minimising
        # the unshifted specific energy.
        rho = 10**logrho                      # g/cm^3

        rho_out  = np.empty(nrho)
        P_out    = np.empty(nrho)
        eps_out  = np.empty(nrho)
        ye_out   = np.empty(nrho)

        for irho in range(nrho):
            e_slice = logenergy_2d[irho, :]    # log10(eps+shift)(Ye)
            iye_min = np.argmin(e_slice)        # beta-equilibrium Ye

            rho_out[irho] = rho[irho]
            P_out[irho]   = 10**logpress_2d[irho, iye_min]   # dyne/cm^2
            ye_out[irho]  = ye[iye_min]

            # Recover specific internal energy (erg/g), removing the shift
            eps_out[irho] = 10**e_slice[iye_min] - energy_shift

        # ---- filter bad values ----
        # eps_specific can be negative (bound matter), so only require
        # finite values and positive pressure.
        valid = np.isfinite(P_out) & np.isfinite(eps_out) & (P_out > 0)
        rho_out  = rho_out[valid]
        P_out    = P_out[valid]
        eps_out  = eps_out[valid]
        ye_out   = ye_out[valid]

        # Enforce strict pressure monotonicity (required by the TOV
        # interpolator).  Some EOS tables (e.g. LS180) have wiggles
        # near the nuclear saturation transition; keep only the
        # running-maximum envelope.
        mono = np.ones(len(P_out), dtype=bool)
        P_max_so_far = P_out[0]
        for i in range(1, len(P_out)):
            if P_out[i] <= P_max_so_far:
                mono[i] = False
            else:
                P_max_so_far = P_out[i]
        if not mono.all():
            n_drop = (~mono).sum()
            print(f"  Dropped {n_drop} non-monotonic pressure points")
            rho_out = rho_out[mono]
            P_out   = P_out[mono]
            eps_out = eps_out[mono]
            ye_out  = ye_out[mono]

        print(f"\nExtracted {len(rho_out)} valid data points")
        print(f"  Density range: {rho_out[0]:.3e} to {rho_out[-1]:.3e} g/cm^3")
        print(f"  Pressure range: {P_out[0]:.3e} to {P_out[-1]:.3e} dyne/cm^2")
        print(f"  Ye range (beta-eq): {ye_out.min():.4f} to {ye_out.max():.4f}")

        # ---- write output ----
        header = (
            f"# EOS Table extracted from {os.path.basename(str(filename))}\n"
            f"# Source: stellarcollapse.org\n"
            f"# Temperature: T = {10**logtemp[itemp]:.4f} MeV (cold)\n"
            f"# Electron fraction: Ye = beta-equilibrium (minimised energy per baryon)\n"
            f"# Units: CGS\n"
            f"# Columns: rho (g/cm^3), pressure (dyne/cm^2), eps_specific (erg/g)\n"
        )

        with open(str(output_file), 'w') as out:
            out.write(header)
            for r, p, e in zip(rho_out, P_out, eps_out):
                out.write(f"{r:.6e}  {p:.6e}  {e:.6e}\n")

        print(f"\nOK Saved to: {output_file}")

        return rho_out, P_out, eps_out


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
