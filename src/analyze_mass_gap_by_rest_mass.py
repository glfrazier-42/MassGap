"""
Improved Mass Gap Analysis - Matching REST Masses

KEY INSIGHT: The gap appears when we look at configurations with the SAME rest mass
but different tiers. The same rest mass produces different gravitational masses
because pressure contributions differ between tiers.

This script:
1. Finds neutron star configurations across a range of masses
2. For each neutron star, finds a quark star with the SAME rest mass
3. Compares their gravitational masses
4. Maps out where the gap appears
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from bhem.tov_solver import TOVSolver, Constants
from bhem.eos import STANDARD_EOS


def find_quark_star_at_rest_mass(target_rest_mass, quark_solutions, tolerance=0.05):
    """
    Find quark star configuration with rest mass closest to target
    
    Parameters:
    -----------
    target_rest_mass : float
        Target rest mass in M_sun
    quark_solutions : list
        List of TOVResult objects for quark stars
    tolerance : float
        Maximum fractional difference in rest mass
        
    Returns:
    --------
    TOVResult or None
    """
    if not quark_solutions:
        return None
    
    # Find closest match
    best_match = min(quark_solutions, 
                     key=lambda s: abs(s.rest_mass_solar - target_rest_mass))
    
    # Check if within tolerance
    frac_diff = abs(best_match.rest_mass_solar - target_rest_mass) / target_rest_mass
    
    if frac_diff > tolerance:
        return None
    
    return best_match


def analyze_gap_by_rest_mass(neutron_solutions, quark_solutions, verbose=True):
    """
    Analyze mass gap by matching rest masses
    
    For each neutron star configuration, find the quark star with matching
    rest mass and calculate the gravitational mass difference.
    """
    
    print("=" * 70)
    print("MASS GAP ANALYSIS - MATCHING REST MASSES")
    print("=" * 70)
    
    results = []
    
    print(f"\n{'M_rest':<12} {'NS: M_grav':<12} {'QS: M_grav':<12} {'Gap (M_grav)':<15} {'Gap Type':<10}")
    print("-" * 70)
    
    for ns in neutron_solutions:
        # Find matching quark star
        qs = find_quark_star_at_rest_mass(ns.rest_mass_solar, quark_solutions, 
                                           tolerance=0.05)
        
        if qs is not None:
            gap = qs.mass_solar - ns.mass_solar
            gap_type = "QS > NS" if gap > 0 else "NS > QS"
            
            results.append({
                'M_rest': ns.rest_mass_solar,
                'M_grav_NS': ns.mass_solar,
                'M_grav_QS': qs.mass_solar,
                'gap': gap,
                'ns_config': ns,
                'qs_config': qs
            })
            
            print(f"{ns.rest_mass_solar:<12.3f} {ns.mass_solar:<12.3f} "
                  f"{qs.mass_solar:<12.3f} {gap:<15.3f} {gap_type:<10}")
    
    if not results:
        print("\nNo matching rest masses found!")
        return None
    
    # Find maximum gap
    max_gap_result = max(results, key=lambda r: abs(r['gap']))
    
    print("\n" + "=" * 70)
    print("MAXIMUM GAP CONFIGURATION")
    print("=" * 70)
    
    print(f"\nRest Mass (same for both): {max_gap_result['M_rest']:.3f} M_sun")
    print(f"\nNeutron Star:")
    print(f"  M_grav = {max_gap_result['M_grav_NS']:.3f} M_sun")
    print(f"  M_pressure = {max_gap_result['ns_config'].pressure_mass_solar:.4f} M_sun")
    print(f"  Fraction = {max_gap_result['ns_config'].pressure_mass_solar/max_gap_result['M_grav_NS']:.2%}")
    
    print(f"\nQuark Star:")
    print(f"  M_grav = {max_gap_result['M_grav_QS']:.3f} M_sun")
    print(f"  M_pressure = {max_gap_result['qs_config'].pressure_mass_solar:.4f} M_sun")
    print(f"  Fraction = {max_gap_result['qs_config'].pressure_mass_solar/max_gap_result['M_grav_QS']:.2%}")
    
    print(f"\n{'GRAVITATIONAL MASS GAP:':-^70}")
    print(f"  Delta M_grav = {max_gap_result['gap']:.3f} M_sun")
    print(f"  Range: {max_gap_result['M_grav_NS']:.2f} to {max_gap_result['M_grav_QS']:.2f} M_sun")
    
    pressure_diff = (max_gap_result['qs_config'].pressure_mass_solar - 
                     max_gap_result['ns_config'].pressure_mass_solar)
    print(f"\n{'PRESSURE CONTRIBUTION DIFFERENCE:':-^70}")
    print(f"  Delta M_pressure = {pressure_diff:.4f} M_sun")
    print(f"  This explains {abs(pressure_diff/max_gap_result['gap'])*100:.1f}% of the gap")
    
    return results


def plot_gap_analysis(results, neutron_solutions, quark_solutions):
    """
    Visualize the mass gap as a function of rest mass
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Extract data
    M_rest = [r['M_rest'] for r in results]
    M_grav_NS = [r['M_grav_NS'] for r in results]
    M_grav_QS = [r['M_grav_QS'] for r in results]
    gaps = [r['gap'] for r in results]
    
    # 1. Gravitational mass vs Rest mass
    ax1 = axes[0, 0]
    ax1.plot(M_rest, M_grav_NS, 'bo-', label='Neutron Star', linewidth=2, markersize=6)
    ax1.plot(M_rest, M_grav_QS, 'ro-', label='Quark Star', linewidth=2, markersize=6)
    ax1.fill_between(M_rest, M_grav_NS, M_grav_QS, alpha=0.3, color='yellow', 
                     label='Mass Gap')
    ax1.set_xlabel('Rest Mass (M_sun)', fontsize=12)
    ax1.set_ylabel('Gravitational Mass (M_sun)', fontsize=12)
    ax1.set_title('Same Rest Mass -> Different Gravitational Masses', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 2. Gap size vs Rest mass
    ax2 = axes[0, 1]
    ax2.plot(M_rest, gaps, 'g^-', linewidth=2, markersize=8)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Rest Mass (M_sun)', fontsize=12)
    ax2.set_ylabel('Gap = M_grav(QS) - M_grav(NS) (M_sun)', fontsize=12)
    ax2.set_title('Gravitational Mass Gap vs Rest Mass', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Highlight maximum gap
    max_idx = np.argmax(np.abs(gaps))
    ax2.plot(M_rest[max_idx], gaps[max_idx], 'r*', markersize=20, 
             label=f'Max gap: {gaps[max_idx]:.3f} M_sun')
    ax2.legend(fontsize=10)
    
    # 3. Pressure contributions
    ax3 = axes[1, 0]
    M_pressure_NS = [r['ns_config'].pressure_mass_solar for r in results]
    M_pressure_QS = [r['qs_config'].pressure_mass_solar for r in results]
    
    ax3.plot(M_rest, M_pressure_NS, 'bo-', label='Neutron Star', linewidth=2)
    ax3.plot(M_rest, M_pressure_QS, 'ro-', label='Quark Star', linewidth=2)
    ax3.set_xlabel('Rest Mass (M_sun)', fontsize=12)
    ax3.set_ylabel('Pressure Contribution to Mass (M_sun)', fontsize=12)
    ax3.set_title('Pressure Contributions at Same Rest Mass', 
                  fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # 4. Pressure fraction
    ax4 = axes[1, 1]
    frac_NS = [r['ns_config'].pressure_mass_solar/r['M_grav_NS'] for r in results]
    frac_QS = [r['qs_config'].pressure_mass_solar/r['M_grav_QS'] for r in results]
    
    ax4.plot(M_rest, [f*100 for f in frac_NS], 'bo-', label='Neutron Star', linewidth=2)
    ax4.plot(M_rest, [f*100 for f in frac_QS], 'ro-', label='Quark Star', linewidth=2)
    ax4.set_xlabel('Rest Mass (M_sun)', fontsize=12)
    ax4.set_ylabel('Pressure Contribution (%)', fontsize=12)
    ax4.set_title('Pressure as Percentage of Gravitational Mass', 
                  fontsize=13, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    output_path = Path(__file__).parent.parent / 'docs' / 'images' / 'mass_gap_by_rest_mass.png'
    output_path.parent.mkdir(exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")
    
    return fig


def main():
    """
    Main analysis comparing configurations at matched rest masses
    """
    print("MASS GAP ANALYSIS - REST MASS MATCHING")
    print("Key Insight: Same rest mass -> different gravitational mass -> gap!")
    print()
    
    # Calculate neutron star solutions
    print("Calculating neutron star solutions...")
    eos_ns = STANDARD_EOS['neutron_validated']
    solver_ns = TOVSolver(eos_ns)
    
    P_range_ns = np.logspace(33, 36, 100)
    neutron_solutions = []
    for P_c in P_range_ns:
        result = solver_ns.solve(P_c, dr=1e3)
        if result is not None:
            neutron_solutions.append(result)
    
    print(f"  Found {len(neutron_solutions)} valid neutron star configurations")
    
    # Calculate quark star solutions
    print("Calculating quark star solutions...")
    eos_qs = STANDARD_EOS['quark_standard']
    solver_qs = TOVSolver(eos_qs)
    
    P_range_qs = np.logspace(33, 37, 100)
    quark_solutions = []
    for P_c in P_range_qs:
        result = solver_qs.solve(P_c, dr=1e3)
        if result is not None:
            quark_solutions.append(result)
    
    print(f"  Found {len(quark_solutions)} valid quark star configurations")
    
    # Analyze gap by matching rest masses
    results = analyze_gap_by_rest_mass(neutron_solutions, quark_solutions)
    
    if results:
        print("\nGenerating plots...")
        plot_gap_analysis(results, neutron_solutions, quark_solutions)
        
        print("\n" + "=" * 70)
        print("ANALYSIS COMPLETE")
        print("=" * 70)
        print("\nKey Finding:")
        print("  At the same rest mass, different tiers produce different gravitational")
        print("  masses due to different pressure contributions. This creates an")
        print("  observable gap in the gravitational mass spectrum.")


if __name__ == '__main__':
    main()
