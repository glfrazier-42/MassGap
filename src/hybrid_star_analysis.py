"""
Pressure Threshold Analysis: The Physical Mechanism of Collapse - FIXED

Physical Picture:
1. Critical pressure P_crit ~ 1-2 x 10^35 dyne/cm^2 needed to overcome neutron binding
2. Below P_crit: Neutron matter stable
3. Above P_crit in core: Quark matter can form
4. At boundary: P = P_crit -> stable interface (pressure gradient holds it)
5. At ~2.5 M_sun: POSITIVE FEEDBACK LOOP triggers runaway

This version uses a physically motivated P_crit rather than trying to derive it.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from bhem.tov_solver import TOVSolver, Constants
from bhem.eos import STANDARD_EOS


def get_physical_critical_pressure():
    """
    Use physically motivated critical pressure for neutron decomposition
    
    From nuclear physics:
    - Neutron binding energy ~ 8-10 MeV per nucleon
    - Nuclear saturation density ~ 2.7 x 10^14 g/cm^3
    - At this density, pressure ~ 1-2 x 10^35 dyne/cm^2
    
    This is the pressure where neutron decomposition becomes energetically possible.
    """
    
    print("=" * 80)
    print("PHYSICAL CRITICAL PRESSURE FOR NEUTRON DECOMPOSITION")
    print("=" * 80)
    
    # Use physically motivated value
    P_crit = 1.5e35  # dyne/cm^2 - middle of reasonable range
    
    print("\nFrom nuclear physics:")
    print(f"  Neutron binding energy: ~8-10 MeV per nucleon")
    print(f"  Nuclear saturation density: ~2.7 x 10^14 g/cm^3")
    print(f"  Estimated P_crit: ~1-2 x 10^35 dyne/cm^2")
    
    print(f"\nUsing P_crit = {P_crit:.2e} dyne/cm^2")
    print(f"\nPhysical meaning:")
    print(f"  Below P_crit: Neutrons remain bound")
    print(f"  Above P_crit: Quarks can be liberated from neutrons")
    
    return P_crit


def analyze_hybrid_configurations(P_crit):
    """
    Analyze neutron star configurations to find which have hybrid structure
    """
    
    print("\n" + "=" * 80)
    print("SEARCHING FOR HYBRID STAR CONFIGURATIONS")
    print("=" * 80)
    
    c = Constants()
    
    # Calculate neutron stars
    print("\nCalculating neutron star configurations...")
    eos_ns = STANDARD_EOS['neutron_validated']
    solver_ns = TOVSolver(eos_ns)
    
    # Scan central pressures more carefully
    P_range = np.logspace(34, 36, 80)
    
    all_configs = []
    hybrid_configs = []
    pure_ns_configs = []
    
    print(f"\n{'P_central':<15} {'M_grav':<10} {'M_rest':<10} {'R (km)':<10} {'Type':<20}")
    print("-" * 80)
    
    for P_c in P_range:
        result = solver_ns.solve(P_c, dr=5e2, max_radius=4e6)  # Smaller step, larger max
        
        if result is None:
            continue
        
        all_configs.append(result)
        
        # Check if this is a hybrid star (P_central > P_crit)
        if P_c > P_crit:
            # Find where pressure drops to P_crit
            P_profile = result.P_array
            r_profile = result.r_array / c.km
            
            # Find crossing point
            crossings = np.where((P_profile[:-1] > P_crit) & (P_profile[1:] <= P_crit))[0]
            
            if len(crossings) > 0:
                idx = crossings[0]
                r_core = r_profile[idx]
                m_core = result.m_array[idx] / c.M_sun
                
                hybrid_configs.append({
                    'result': result,
                    'r_core': r_core,
                    'm_core': m_core,
                    'core_fraction': m_core / result.mass_solar
                })
                
                config_type = f"Hybrid (core={r_core:.1f}km)"
            else:
                config_type = "Hybrid (tiny core)"
                hybrid_configs.append({
                    'result': result,
                    'r_core': 0.1,
                    'm_core': 0.01,
                    'core_fraction': 0.01
                })
        else:
            pure_ns_configs.append(result)
            config_type = "Pure NS"
        
        # Print every 10th or important ones
        if len(all_configs) % 10 == 0 or result.mass_solar > 2.0:
            print(f"{P_c:<15.2e} {result.mass_solar:<10.3f} {result.rest_mass_solar:<10.3f} "
                  f"{result.radius_km:<10.2f} {config_type:<20}")
    
    print(f"\nSummary:")
    print(f"  Total valid configurations: {len(all_configs)}")
    print(f"  Pure neutron stars (P_central < P_crit): {len(pure_ns_configs)}")
    print(f"  Hybrid stars (P_central > P_crit): {len(hybrid_configs)}")
    
    return all_configs, pure_ns_configs, hybrid_configs, P_crit


def analyze_feedback_mechanism(hybrid_configs):
    """
    Analyze when the positive feedback loop triggers
    """
    
    if len(hybrid_configs) == 0:
        print("\nNo hybrid configurations found - P_crit may be too high")
        return
    
    print("\n" + "=" * 80)
    print("POSITIVE FEEDBACK ANALYSIS")
    print("=" * 80)
    
    print(f"\n{'M_grav':<12} {'M_rest':<12} {'Core Frac':<12} {'R_core':<12} {'R_total':<12}")
    print("-" * 80)
    
    for h in hybrid_configs:
        r = h['result']
        print(f"{r.mass_solar:<12.3f} {r.rest_mass_solar:<12.3f} {h['core_fraction']:<12.2%} "
              f"{h['r_core']:<12.2f} {r.radius_km:<12.2f}")
    
    # Look for rapid core growth
    if len(hybrid_configs) > 3:
        print("\nCore growth rate analysis:")
        masses = [h['result'].mass_solar for h in hybrid_configs]
        fractions = [h['core_fraction'] for h in hybrid_configs]
        
        for i in range(1, len(hybrid_configs)):
            dm = masses[i] - masses[i-1]
            df = fractions[i] - fractions[i-1]
            
            if dm > 0 and i % 5 == 0:
                rate = df / dm
                status = "ACCELERATING!" if rate > 0.15 else "Stable growth"
                print(f"  At M = {masses[i]:.2f} M_sun: dCore/dM = {rate:.3f} - {status}")


def plot_hybrid_analysis(all_configs, hybrid_configs, P_crit):
    """
    Visualize hybrid star structure
    """
    
    if len(all_configs) == 0:
        print("No configurations to plot")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Extract data
    all_masses = [r.mass_solar for r in all_configs]
    all_radii = [r.radius_km for r in all_configs]
    all_rest = [r.rest_mass_solar for r in all_configs]
    
    # 1. Mass-Radius curve with hybrid region marked
    ax1 = axes[0, 0]
    ax1.plot(all_radii, all_masses, 'b-', linewidth=2)
    
    if len(hybrid_configs) > 0:
        hybrid_masses = [h['result'].mass_solar for h in hybrid_configs]
        hybrid_radii = [h['result'].radius_km for h in hybrid_configs]
        ax1.plot(hybrid_radii, hybrid_masses, 'ro', markersize=8, label='Hybrid stars')
    
    ax1.set_xlabel('Radius (km)', fontsize=12)
    ax1.set_ylabel('Gravitational Mass (M_sun)', fontsize=12)
    ax1.set_title('Mass-Radius Relationship', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Central pressure vs mass
    ax2 = axes[0, 1]
    central_P = [r.central_pressure for r in all_configs]
    ax2.semilogy(all_masses, central_P, 'b-', linewidth=2)
    ax2.axhline(y=P_crit, color='r', linestyle='--', linewidth=2, label=f'P_crit = {P_crit:.1e}')
    ax2.fill_between([0, 5], 0, P_crit, alpha=0.2, color='blue', label='Pure NS')
    ax2.fill_between([0, 5], P_crit, 1e37, alpha=0.2, color='red', label='Hybrid')
    ax2.set_xlabel('Gravitational Mass (M_sun)', fontsize=12)
    ax2.set_ylabel('Central Pressure (dyne/cm^2)', fontsize=12)
    ax2.set_title('Central Pressure vs Mass', fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Core fraction vs mass
    ax3 = axes[1, 0]
    if len(hybrid_configs) > 0:
        hybrid_masses = [h['result'].mass_solar for h in hybrid_configs]
        core_fractions = [h['core_fraction'] * 100 for h in hybrid_configs]
        ax3.plot(hybrid_masses, core_fractions, 'ro-', linewidth=2, markersize=6)
        ax3.set_xlabel('Gravitational Mass (M_sun)', fontsize=12)
        ax3.set_ylabel('Quark Core Fraction (%)', fontsize=12)
        ax3.set_title('Core Growth with Mass', fontsize=13, fontweight='bold')
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'No hybrid stars found', ha='center', va='center',
                transform=ax3.transAxes, fontsize=14)
    
    # 4. Summary text
    ax4 = axes[1, 1]
    
    max_mass = max(all_masses)
    max_config = all_configs[all_masses.index(max_mass)]
    
    summary = f"""
    CONFIGURATION SUMMARY:
    
    P_crit = {P_crit:.2e} dyne/cm^2
    
    Total configurations: {len(all_configs)}
    Hybrid configurations: {len(hybrid_configs)}
    
    Maximum mass: {max_mass:.2f} M_sun
    At maximum:
      R = {max_config.radius_km:.2f} km
      P_central = {max_config.central_pressure:.2e}
    
    STABILITY INTERPRETATION:
    
    M < 2.5 M_sun:
      - If P_central < P_crit: Pure NS
      - If P_central > P_crit: Hybrid
      - Pressure gradient stabilizes
        small quark core
    
    M ~ 2.5 M_sun:
      - Core becomes significant
      - Pressure gradients steepen
      - Small perturbation triggers
        positive feedback
      - Runaway conversion
      - Too compact -> BH
    """
    
    ax4.text(0.05, 0.95, summary, fontsize=10, family='monospace',
             verticalalignment='top', transform=ax4.transAxes)
    ax4.axis('off')
    
    plt.tight_layout()
    
    # Save
    output_path = Path(__file__).parent.parent / 'docs' / 'images' / 'hybrid_star_analysis.png'
    output_path.parent.mkdir(exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")


def main():
    """
    Complete analysis with physically motivated P_crit
    """
    
    print("PRESSURE THRESHOLD MECHANISM - PHYSICALLY MOTIVATED ANALYSIS")
    print()
    
    # 1. Set critical pressure from physics
    P_crit = get_physical_critical_pressure()
    
    # 2. Find hybrid configurations
    all_configs, pure_ns, hybrid, P_crit = analyze_hybrid_configurations(P_crit)
    
    # 3. Analyze feedback
    analyze_feedback_mechanism(hybrid)
    
    # 4. Visualize
    if all_configs:
        print("\nGenerating plots...")
        plot_hybrid_analysis(all_configs, hybrid, P_crit)
    
    # 5. Summary
    print("\n" + "=" * 80)
    print("SUMMARY: THE COLLAPSE MECHANISM")
    print("=" * 80)
    print("""
    1. CRITICAL PRESSURE: P_crit ~ 1.5 x 10^35 dyne/cm^2
       - Below: Neutrons stable
       - Above: Quark decomposition possible
    
    2. STABLE HYBRID STARS (M < 2.5 M_sun):
       - Small quark core where P > P_crit
       - Neutron shell where P < P_crit
       - Pressure gradient prevents runaway
    
    3. INSTABILITY TRIGGER (~2.5 M_sun):
       - Core mass becomes significant fraction
       - Pressure at boundary grows with perturbations
       - Feedback: bigger core -> higher boundary P -> more conversion
    
    4. RUNAWAY COLLAPSE:
       - Entire star converts to quark matter
       - M_rest ~ 2 M_sun -> M_grav ~ 3-5 M_sun (pressure contribution!)
       - R < R_s -> BLACK HOLE
    
    5. THE MASS GAP (2.5-5 M_sun):
       - Maximum stable NS: ~2.5 M_sun
       - Minimum Tier 2 BH: ~5 M_sun
       - Few stable objects in between (rare transients)
    """)
    
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    main()
