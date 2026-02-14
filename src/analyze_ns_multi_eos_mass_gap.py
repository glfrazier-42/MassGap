"""
Multi-EOS Mass Gap Robustness Testing - REFACTORED WITH MODE SUPPORT

Analyzes maximum masses across multiple equations of state with three modes:

Modes:
------
  neutron: Pure neutron star analysis with linear extrapolation (original)
  hybrid:  Hybrid star analysis showing direct M_max from NS+quark core
  both:    Run both modes and create side-by-side comparison

Usage:
------
  python analyze_ns_multi_eos_mass_gap.py --mode neutron    # Original behavior
  python analyze_ns_multi_eos_mass_gap.py --mode hybrid     # Hybrid stars
  python analyze_ns_multi_eos_mass_gap.py --mode both       # Comparison
  python analyze_ns_multi_eos_mass_gap.py --mode hybrid --P-transition 5e34
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from bhem import NeutronStarEOS, QuarkStarEOS
from bhem.eos import HybridEOS
from bhem.tov_solver_scipy import TOVSolverScipy

# EOS to test
NS_EOS_LIST = {
    'LS180': {
        'file': 'LS180_cold_betaeq.txt',
        'label': 'LS180 (K=180, SOFT)',
        'color': 'blue'
    },
    'LS220': {
        'file': 'LS220_cold_betaeq.txt', 
        'label': 'LS220 (K=220, MEDIUM)',
        'color': 'green'
    },
    'LS375': {
        'file': 'LS375_cold_betaeq.txt',
        'label': 'LS375 (K=375, STIFF)',
        'color': 'orange'
    },
    'SFHo': {
        'file': 'SFHo_cold_betaeq.txt',
        'label': 'SFHo (STIFF)',
        'color': 'red'
    },
    'DD2': {
        'file': 'DD2_cold_betaeq.txt',
        'label': 'DD2 (STIFF)',
        'color': 'purple'
    },
    'HShen': {
        'file': 'HShen_cold_betaeq.txt',
        'label': 'HShen',
        'color': 'brown'
    }
}

# B values to test (MeV/fm³) - only for neutron mode
B_VALUES = [50, 60, 70]

# Constants
M_SUN = 1.989e33  # g
C_LIGHT = 2.99792458e10  # cm/s
G_CONST = 6.67430e-8  # cm³/g/s²


def compute_rest_mass(result):
    """Compute baryonic (rest) mass from TOV result"""
    if hasattr(result, 'rest_mass'):
        return result.rest_mass
    
    r = result.r_array
    m = result.m_array
    rho = result.rho_array
    
    if len(r) == 0:
        return 0.0
    
    # M_rest = ∫ ρ √(g_rr) 4πr² dr
    metric_factor = 1.0 / np.sqrt(1.0 - 2*G_CONST*m / (r*C_LIGHT**2))
    integrand = rho * metric_factor * 4 * np.pi * r**2
    M_rest = np.trapz(integrand, r)
    
    return M_rest


def analyze_star(name, eos_file, base_dir, mode='neutron', P_transition=1.0e35):
    """
    Analyze compact star using EOS table.
    
    Parameters:
    -----------
    name : str
        EOS name
    eos_file : str
        Path to EOS table file
    base_dir : Path
        Project base directory
    mode : str
        'neutron' for pure NS, 'hybrid' for NS with quark core
    P_transition : float
        Phase transition pressure for hybrid mode (dyne/cm²)
    
    Returns dict with:
        - M_max: maximum gravitational mass (Msun)
        - M_rest_max: rest mass at maximum (Msun)
        - masses, radii, rest_masses: full solution arrays
        - mode: which mode was used
    """
    mode_label = 'HYBRID STAR' if mode == 'hybrid' else 'NEUTRON STAR'
    print(f"\n{'='*70}")
    print(f"ANALYZING {mode_label}: {name}")
    print(f"{'='*70}")
    
    eos_path = base_dir / 'data' / 'eos_tables' / eos_file
    if not eos_path.exists():
        print(f"[X] File not found: {eos_file}")
        return None
    
    # Load EOS based on mode
    print(f"Loading EOS: {eos_file} (mode: {mode})")
    try:
        neutron_eos_obj = NeutronStarEOS(backend='tabulated', table_path=str(eos_path))
        
        if mode == 'hybrid':
            # Create hybrid EOS with quark core
            quark_eos_obj = QuarkStarEOS(backend='analytical', B=60.0)
            eos_obj = HybridEOS(neutron_eos_obj, quark_eos_obj, P_transition)
            print(f"  Using hybrid EOS with P_transition = {P_transition:.2e} dyne/cm^2")
        else:
            eos_obj = neutron_eos_obj
            
    except Exception as e:
        print(f"x Error loading EOS: {e}")
        return None
    
    eos_func = eos_obj.eos_function()
    solver = TOVSolverScipy(eos_func)
    
    # Scan central pressures
    config_label = 'hybrid star' if mode == 'hybrid' else 'neutron star'
    print(f"Scanning {config_label} configurations...")
    p_min = 1e34  # dyne/cm²
    p_max = 5e35
    n_points = 100
    
    pressures = np.logspace(np.log10(p_min), np.log10(p_max), n_points)
    
    masses = []
    radii = []
    rest_masses = []
    
    for p_c in pressures:
        result = solver.solve(p_c)
        if result is not None:
            masses.append(result.mass_solar)
            radii.append(result.radius_km)
            M_rest = compute_rest_mass(result)
            rest_masses.append(M_rest / M_SUN)
    
    masses = np.array(masses)
    radii = np.array(radii)
    rest_masses = np.array(rest_masses)
    
    if len(masses) == 0:
        print("x No solutions found!")
        return None
    
    # Find maximum REST MASS (stability criterion: dM_rest/dρ_c > 0)
    i_max = np.argmax(rest_masses)
    M_max = masses[i_max]
    R_max = radii[i_max]
    M_rest_max = rest_masses[i_max]
    
    max_label = 'Hybrid Star' if mode == 'hybrid' else 'Neutron Star'
    print(f"\n{max_label} Maximum:")
    print(f"  M_max = {M_max:.4f} Msun")
    print(f"  M_rest_max = {M_rest_max:.4f} Msun")
    print(f"  R(M_max) = {R_max:.2f} km")
    
    return {
        'name': name,
        'mode': mode,
        'M_max': M_max,
        'R_max': R_max,
        'M_rest_max': M_rest_max,
        'masses': masses[:i_max+1],
        'radii': radii[:i_max+1],
        'rest_masses': rest_masses[:i_max+1]
    }


def analyze_quark_star(B, verbose=True):
    """Analyze quark star with given bag constant B (MeV/fm³)"""
    if verbose:
        print(f"\n{'-'*70}")
        print(f"ANALYZING QUARK STAR: B = {B} MeV/fm³")
        print(f"{'-'*70}")
    
    qs_eos = QuarkStarEOS(backend='analytical', B=B)
    qs_func = qs_eos.eos_function()
    solver = TOVSolverScipy(qs_func)
    
    if verbose:
        print("Scanning quark star configurations...")
    
    p_min = 1e33
    p_max = 1e37
    n_points = 200
    
    pressures = np.logspace(np.log10(p_min), np.log10(p_max), n_points)
    
    masses = []
    rest_masses = []
    
    for p_c in pressures:
        result = solver.solve(p_c)
        if result is not None:
            masses.append(result.mass_solar)
            M_rest = compute_rest_mass(result)
            rest_masses.append(M_rest / M_SUN)
    
    masses = np.array(masses)
    rest_masses = np.array(rest_masses)
    
    if len(masses) == 0:
        if verbose:
            print("x No quark star solutions found!")
        return None
    
    # Find maximum REST MASS
    i_max = np.argmax(rest_masses)
    M_max = masses[i_max]
    M_rest_max = rest_masses[i_max]
    
    # Get slope at maximum for extrapolation
    i_fit_start = max(0, i_max - 5)
    slope, intercept = np.polyfit(
        rest_masses[i_fit_start:i_max+1],
        masses[i_fit_start:i_max+1],
        1
    )
    
    if verbose:
        print(f"\nQuark Star Maximum:")
        print(f"  M_max = {M_max:.4f} Msun")
        print(f"  M_rest_max = {M_rest_max:.4f} Msun")
        print(f"  dM_grav/dM_rest at max = {slope:.3f}")
    
    return {
        'B': B,
        'M_max': M_max,
        'M_rest_max': M_rest_max,
        'masses': masses[:i_max+1],
        'rest_masses': rest_masses[:i_max+1],
        'slope_at_max': slope,
        'intercept': intercept
    }


def predict_mass_gap(ns_result, qs_result, verbose=True):
    """Predict mass gap using NS-to-QS transition with linear extrapolation"""
    M_rest_target = ns_result['M_rest_max']
    M_max_QS = qs_result['M_max']
    M_rest_max_QS = qs_result['M_rest_max']
    slope = qs_result['slope_at_max']
    
    delta_M_rest = M_rest_target - M_rest_max_QS
    
    if delta_M_rest < 0 and verbose:
        print(f"! Warning: NS M_rest ({M_rest_target:.3f}) < QS M_rest_max ({M_rest_max_QS:.3f})")
    
    M_grav_predicted = M_max_QS + slope * delta_M_rest
    M_max_NS = ns_result['M_max']
    gap_width = M_grav_predicted - M_max_NS
    
    if verbose:
        print(f"\nMass Gap Prediction:")
        print(f"  NS M_rest_max = {M_rest_target:.4f} Msun")
        print(f"  Predicted M_grav = {M_grav_predicted:.3f} Msun")
        print(f"  NS M_max = {M_max_NS:.3f} Msun")
        print(f"  >>> GAP WIDTH = {gap_width:.3f} Msun <<<")
    
    return {
        'M_grav_predicted': M_grav_predicted,
        'gap_width': gap_width,
        'delta_M_rest': delta_M_rest,
        'extrapolation_valid': delta_M_rest > 0
    }


def run_full_analysis(base_dir, mode='neutron', P_transition=1.0e35):
    """Run complete multi-EOS analysis"""
    
    mode_title = {
        'neutron': 'PURE NEUTRON STAR',
        'hybrid': 'HYBRID STAR (NS + QUARK CORE)'
    }[mode]
    
    print("="*70)
    print(f"MULTI-EOS ANALYSIS: {mode_title}")
    print("="*70)
    print(f"\nTesting across:")
    print(f"  - {len(NS_EOS_LIST)} neutron star equations of state")
    if mode == 'hybrid':
        print(f"  - P_transition = {P_transition:.2e} dyne/cm^2")
    if mode == 'neutron':
        print(f"  - {len(B_VALUES)} quark matter bag constants")
    print()
    
    # Step 1: Analyze all stars
    step1_label = 'HYBRID STAR' if mode == 'hybrid' else 'NEUTRON STAR'
    print("\n" + "="*70)
    print(f"STEP 1: ANALYZING {step1_label} EOS")
    print("="*70)
    
    ns_results = {}
    for name, info in NS_EOS_LIST.items():
        result = analyze_star(name, info['file'], base_dir, mode=mode, P_transition=P_transition)
        if result:
            ns_results[name] = result
    
    if len(ns_results) == 0:
        print(f"\nx No {step1_label.lower()} solutions found!")
        return None
    
    # For neutron mode only: analyze quark stars and predict gaps
    qs_results = {}
    gap_predictions = {}
    
    if mode == 'neutron':
        print("\n" + "="*70)
        print("STEP 2: ANALYZING QUARK STAR EOS")
        print("="*70)
        
        for B in B_VALUES:
            result = analyze_quark_star(B, verbose=True)
            if result:
                qs_results[B] = result
        
        if len(qs_results) == 0:
            print("\nx No quark star solutions found!")
            return None
        
        print("\n" + "="*70)
        print("STEP 3: PREDICTING MASS GAPS")
        print("="*70)
        
        for ns_name in ns_results:
            gap_predictions[ns_name] = {}
            print(f"\n{'='*70}")
            print(f"NS EOS: {ns_name}")
            print(f"{'='*70}")
            
            for B in qs_results:
                print(f"\n  Testing B = {B} MeV/fm³:")
                print(f"  {'-'*66}")
                gap_result = predict_mass_gap(ns_results[ns_name], qs_results[B], verbose=True)
                gap_predictions[ns_name][B] = gap_result
    
    # Summary table
    print("\n" + "="*70)
    if mode == 'neutron':
        print("MASS GAP PREDICTIONS - SUMMARY TABLE")
    else:
        print("MAXIMUM MASSES - SUMMARY TABLE")
    print("="*70)
    print()
    
    if mode == 'neutron':
        # Original table with gap predictions
        header = f"{'EOS':<10} {'NS M_max':<10}"
        for B in B_VALUES:
            header += f" | {'B='+str(B):<10}"
        print(header)
        print("=" * 70)
        
        all_gaps = {B: [] for B in B_VALUES}
        
        for ns_name in sorted(ns_results.keys()):
            M_max_NS = ns_results[ns_name]['M_max']
            row = f"{ns_name:<10} {M_max_NS:<10.3f}"
            
            for B in B_VALUES:
                gap = gap_predictions[ns_name][B]['gap_width']
                valid = gap_predictions[ns_name][B]['extrapolation_valid']
                
                if valid:
                    row += f" | {gap:<10.3f}"
                    all_gaps[B].append(gap)
                else:
                    row += f" | {'! ' + f'{gap:.2f}':<10}"
            
            print(row)
        
        print("-" * 70)
        row = f"{'RANGE':<10} {'':<10}"
        for B in B_VALUES:
            if len(all_gaps[B]) > 0:
                gap_min = min(all_gaps[B])
                gap_max = max(all_gaps[B])
                row += f" | {gap_min:.2f}-{gap_max:.2f}   "
            else:
                row += f" | {'N/A':<10}"
        print(row)
        
    else:
        # Hybrid mode: show maximum masses and gap position
        header = f"{'EOS':<15} {'M_max':<12} {'R_max (km)':<12} {'vs Observed Gap':<20}"
        print(header)
        print("=" * 70)
        
        for ns_name in sorted(ns_results.keys()):
            M_max = ns_results[ns_name]['M_max']
            R_max = ns_results[ns_name]['R_max']
            
            # Check position relative to observed gap (2.5-5.0 Msun)
            if M_max < 2.5:
                gap_pos = "Below gap"
            elif M_max > 5.0:
                gap_pos = "Above gap"
            else:
                gap_pos = "[OK] IN GAP (2.5-5)"
            
            print(f"{ns_name:<15} {M_max:<12.3f} {R_max:<12.2f} {gap_pos:<20}")
    
    print("="*70)
    
    # Physical interpretation
    print("\n" + "="*70)
    print("PHYSICAL INTERPRETATION")
    print("="*70)
    
    if mode == 'neutron':
        for B in B_VALUES:
            if len(all_gaps[B]) > 0:
                gap_min = min(all_gaps[B])
                gap_max = max(all_gaps[B])
                spread = gap_max - gap_min
                
                print(f"\nB = {B} MeV/fm³:")
                print(f"  Gap range: {gap_min:.2f} - {gap_max:.2f} Msun")
                print(f"  Spread: {spread:.2f} Msun")
                
                if gap_min >= 2.5 and gap_max <= 5.0:
                    print(f"  Match: [OK][OK] All within observed gap (2.5-5 Msun)")
                elif gap_min < 5.0 and gap_max > 2.5:
                    print(f"  Match: [OK] Overlaps observed gap")
                else:
                    print(f"  Match: x Outside observed range")
    else:
        print("\nHybrid star maximum masses:")
        M_values = [ns_results[name]['M_max'] for name in ns_results]
        print(f"  Range: {min(M_values):.3f} - {max(M_values):.3f} Msun")
        print(f"  Mean: {np.mean(M_values):.3f} Msun")
        
        in_gap = sum(1 for M in M_values if 2.5 <= M <= 5.0)
        print(f"  Count in gap (2.5-5.0): {in_gap}/{len(M_values)}")
        
        if in_gap == 0:
            print("  Result: [OK] No hybrid stars in observed gap")
        else:
            print("  Result: [X] Some hybrid stars fall in gap")
    
    return {
        'ns_results': ns_results,
        'qs_results': qs_results,
        'gap_predictions': gap_predictions if mode == 'neutron' else {}
    }


def main():
    parser = argparse.ArgumentParser(
        description='Analyze maximum masses across multiple EOS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  neutron : Pure NS with extrapolation (original)
  hybrid  : Hybrid NS+quark core (direct M_max)
  both    : Run both modes for comparison

Example:
  python analyze_ns_multi_eos_mass_gap.py --mode hybrid --P-transition 1e35
        """
    )
    parser.add_argument('--mode', type=str, default='neutron',
                       choices=['neutron', 'hybrid', 'both'],
                       help='Analysis mode')
    parser.add_argument('--P-transition', type=float, default=1.0e35,
                       help='Phase transition pressure (dyne/cm²)')
    
    args = parser.parse_args()
    base_dir = Path(__file__).parent.parent
    
    if args.mode == 'both':
        print("\n" + "="*70)
        print("RUNNING BOTH MODES FOR COMPARISON")
        print("="*70)
        
        print("\n" + "#"*70)
        print("# MODE 1: PURE NEUTRON STARS")
        print("#"*70)
        results_neutron = run_full_analysis(base_dir, mode='neutron')
        
        print("\n\n" + "#"*70)
        print("# MODE 2: HYBRID STARS (NEUTRON + QUARK CORE)")
        print("#"*70)
        results_hybrid = run_full_analysis(base_dir, mode='hybrid', P_transition=args.P_transition)
        
        # Comparison
        print("\n" + "="*70)
        print("COMPARISON: PURE NS vs HYBRID")
        print("="*70)
        print(f"\n{'EOS':<15} {'NS M_max':<12} {'Hybrid M_max':<12} {'Difference':<12}")
        print("-"*70)
        
        if results_neutron and results_hybrid:
            for name in sorted(results_neutron['ns_results'].keys()):
                if name in results_hybrid['ns_results']:
                    M_ns = results_neutron['ns_results'][name]['M_max']
                    M_hybrid = results_hybrid['ns_results'][name]['M_max']
                    diff = M_hybrid - M_ns
                    print(f"{name:<15} {M_ns:<12.3f} {M_hybrid:<12.3f} {diff:+12.3f}")
        
        results = {'neutron': results_neutron, 'hybrid': results_hybrid}
    else:
        results = run_full_analysis(base_dir, mode=args.mode, P_transition=args.P_transition)
    
    if results:
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE [OK]")
        print("="*70)
    else:
        print("\nx Analysis failed")


if __name__ == '__main__':
    main()
