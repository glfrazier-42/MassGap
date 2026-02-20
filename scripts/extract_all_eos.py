"""
Extract cold beta-equilibrium slices from multiple stellarcollapse.org EOS tables.

Reads HDF5 files from data/raw/ and writes 3-column text tables
(rho, P, eps_specific) to data/eos/.

See data/raw/README.md for the list of HDF5 files that must be downloaded
from stellarcollapse.org before running this script.

Usage:
    PYTHONPATH=src venv/Scripts/python.exe scripts/extract_all_eos.py
"""

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


def main():
    base_dir   = Path(__file__).resolve().parent.parent
    h5_dir     = base_dir / 'data' / 'raw'
    output_dir = base_dir / 'data' / 'eos'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("EXTRACTING COLD BETA-EQUILIBRIUM SLICES")
    print("=" * 70)
    print(f"\nHDF5 source dir: {h5_dir}")
    print(f"Output dir:      {output_dir}\n")

    results = {}

    for name, info in EOS_TABLES.items():
        print(f"{name}: {info['description']}")
        print(f"  Reference: {info['reference']}")

        input_file  = h5_dir / info['file']
        output_file = output_dir / f"{name}_cold_betaeq.txt"

        if not input_file.exists():
            print(f"  WARN Skipping - file not found: {info['file']}")
            print()
            continue

        try:
            rho, P, eps = extract_cold_eos(str(input_file), str(output_file))
            npoints = len(rho)
            print(f"  OK Extracted {npoints} points")
            results[name] = {'status': 'extracted', 'points': npoints,
                             'file': output_file}
        except Exception as e:
            print(f"  FAIL Error: {e}")
            results[name] = {'status': 'error', 'error': str(e)}

        print()

    # Summary
    print("=" * 70)
    print("EXTRACTION SUMMARY")
    print("=" * 70)
    print()

    for name in EOS_TABLES:
        if name in results:
            s = results[name]
            if s['status'] == 'extracted':
                print(f"  OK   {name}: {s['points']} points -> {s['file'].name}")
            else:
                print(f"  FAIL {name}: {s.get('error', 'Unknown error')}")
        else:
            print(f"  WARN {name}: Not processed (HDF5 not found)")

    print()


if __name__ == '__main__':
    main()
