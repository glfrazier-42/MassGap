"""
Equations of State (EOS) for Compact Stars - WITH TABULATED BACKEND SUPPORT

Implements various EOS models for:
- Neutron stars (Tier 1)
- Quark stars (Tier 2)
- Transition regions

Each EOS can use either:
- Analytical backend: Direct calculation from formulas
- Tabulated backend: Interpolation from data tables

Backend is chosen at initialization time.

VALIDATED EOS DATA SOURCES:
============================

For proper validation, tabulated EOS should use REAL data from validated sources:

1. CompOSE Database (Recommended):
   URL: https://compose.obspm.fr/
   
   Recommended tables:
   - APR (Akmal-Pandharipande-Ravenhall): https://compose.obspm.fr/eos/68
     * M_max = 2.19 M☉, R_1.4 = 11.37 km
     * Reference: Akmal et al., Phys. Rev. C 58, 1804 (1998)
   
   - SLy4: https://compose.obspm.fr/eos/67
     * M_max = 2.05 M☉, R_1.4 = 11.67 km  
     * Reference: Douchin & Haensel, A&A 380, 151 (2001)
   
   Citation for CompOSE:
   Typel, S., Oertel, M., Klähn, T., Phys. Part. Nucl. 46, 633 (2015)

2. Installation:
   - Download eos.zip from CompOSE
   - Extract to data/eos_tables/<EOS_NAME>/
   - Key file: eos.thermo (contains P, ε vs density)

3. Using 'generate_from_analytical':
   This mode generates tables from analytical formulas for TESTING ONLY.
   It defeats the purpose of validation (comparing something to itself).
   Use real tables for proper validation!

See EOS_TABLE_INSTRUCTIONS.md for detailed download instructions.
"""

import numpy as np
from typing import Tuple, Optional, Union
from pathlib import Path
from scipy.interpolate import CubicSpline, interp1d


# Physical constants
c = 2.99792458e10  # cm/s
c2 = c**2
MeV_to_erg = 1.60218e-6  # MeV to erg conversion
fm_to_cm = 1e-13  # fm to cm
MeV_fm3_to_dyne_cm2 = MeV_to_erg / fm_to_cm**3  # Pressure unit conversion


class NeutronStarEOS:
    """
    Neutron star equations of state with analytical or tabulated backend
    
    Parameters
    ----------
    backend : str
        'analytical' or 'tabulated'
    eos_type : str
        For analytical: 'validated' or 'stiff'
        For tabulated: name of table file or 'generate_from_analytical'
    table_path : str or Path, optional
        Path to directory containing EOS tables
    **kwargs : 
        Additional parameters (e.g., K, Gamma for analytical)
    """
    
    def __init__(self, backend='analytical', eos_type='validated', 
                 table_path=None, **kwargs):
        self.backend = backend
        self.eos_type = eos_type
        
        if backend == 'analytical':
            self._init_analytical(eos_type, **kwargs)
        elif backend == 'tabulated':
            self._init_tabulated(eos_type, table_path, **kwargs)
        else:
            raise ValueError(f"Unknown backend: {backend}. Use 'analytical' or 'tabulated'")
    
    def _init_analytical(self, eos_type, **kwargs):
        """Initialize analytical backend"""
        if eos_type == 'validated':
            # Validated polytrope: gives M_max ~ 2.5 M_sun
            self.Gamma = kwargs.get('Gamma', 2.0)
            self.K = kwargs.get('K', 2.34e5)  # CGS units
        elif eos_type == 'stiff':
            # Stiffer EOS for testing
            self.Gamma = kwargs.get('Gamma', 2.5)
            self.K = kwargs.get('K', 3.0e4)
        else:
            raise ValueError(f"Unknown analytical EOS type: {eos_type}")
    
    def _init_tabulated(self, eos_type, table_path, **kwargs):
        """Initialize tabulated backend"""
        if eos_type == 'generate_from_analytical':
            # Generate table from current analytical EOS
            print("WARNING: Generating table from analytical EOS for TESTING ONLY")
            print("This is NOT proper validation - use real tables from CompOSE!")
            self._generate_table_from_analytical(**kwargs)
        else:
            # Load from file
            self._load_table(table_path)
    
    def _generate_table_from_analytical(self, **kwargs):
        """Generate interpolation table from analytical EOS"""
        # Temporarily switch to analytical to generate table
        original_backend = self.backend
        self.backend = 'analytical'
        # Default to 'validated' if not specified
        analytical_type = kwargs.get('analytical_type', 'validated')
        self._init_analytical(analytical_type, **kwargs)
        
        # Generate table over relevant density range
        # Nuclear density: ~2.8e14 g/cm^3
        # Need to extend to very low densities to cover surface region (P -> 0)
        # For Gamma=2, P=K*rho^2, so rho_min ~ 1e7 gives P ~ 2e19 (near zero)
        # This ensures the table covers the full range from center to surface
        rho_min = 1e7   # g/cm^3 - low enough to cover stellar surface
        rho_max = 1e16  # g/cm^3 - high enough for maximum central density
        n_points = 2000  # More points for better accuracy over wider range
        
        self.table_rho = np.logspace(np.log10(rho_min), np.log10(rho_max), n_points)
        self.table_P = np.zeros_like(self.table_rho)
        self.table_eps = np.zeros_like(self.table_rho)
        
        for i, rho in enumerate(self.table_rho):
            P = self._pressure_analytical(rho)
            eps = self._energy_density_analytical(rho, P)
            self.table_P[i] = P
            self.table_eps[i] = eps
        
        # Create interpolators
        self._setup_interpolators()
        
        # Switch back to tabulated backend
        self.backend = original_backend
        print(f"Generated table with {n_points} points")
        print(f"Density range: {rho_min:.2e} to {rho_max:.2e} g/cm^3")
        print(f"Pressure range: {self.table_P[0]:.2e} to {self.table_P[-1]:.2e} dyne/cm^2")
    
    def _load_table(self, table_path):
        """Load EOS table from file"""
        if table_path is None:
            raise ValueError("table_path must be provided when loading from file")
        
        table_file = Path(table_path)
        
        if not table_file.exists():
            raise FileNotFoundError(f"Table file not found: {table_file}")
        
        # Load table
        print(f"Loading EOS table from {table_file}")
        data = np.loadtxt(table_file)
        
        if data.shape[1] < 2:
            raise ValueError("Table must have at least 2 columns: rho, P")
        
        self.table_rho = data[:, 0]
        self.table_P = data[:, 1]
        
        if data.shape[1] >= 3:
            self.table_eps = data[:, 2]
        else:
            # Calculate energy density assuming relativistic correction
            self.table_eps = self.table_rho * (1.0 + self.table_P / (self.table_rho * c2))
        
        # Create interpolators
        self._setup_interpolators()
        
        print(f"Loaded table with {len(self.table_rho)} points")
        print(f"Density range: {self.table_rho[0]:.2e} to {self.table_rho[-1]:.2e} g/cm^3")
        print(f"Pressure range: {self.table_P[0]:.2e} to {self.table_P[-1]:.2e} dyne/cm^2")
    
    def _setup_interpolators(self):
        """Setup interpolation functions"""
        # Use cubic spline for smooth interpolation
        # P(rho), eps(rho)
        self.interp_P = CubicSpline(self.table_rho, self.table_P, extrapolate=False)
        self.interp_eps = CubicSpline(self.table_rho, self.table_eps, extrapolate=False)
        
        # Inverse: rho(P) - use linear interpolation as P is monotonic
        self.interp_rho = interp1d(self.table_P, self.table_rho, 
                                    kind='cubic', 
                                    bounds_error=False,
                                    fill_value=(self.table_rho[0], self.table_rho[-1]))
    
    def _pressure_analytical(self, rho):
        """Calculate pressure analytically (polytropic)"""
        if rho <= 0:
            return 0.0
        return self.K * rho**self.Gamma
    
    def _energy_density_analytical(self, rho, P=None):
        """Calculate energy density analytically"""
        if rho <= 0:
            return 0.0
        if P is None:
            P = self._pressure_analytical(rho)
        # Include relativistic correction
        return rho * (1.0 + P / (rho * c2))
    
    def pressure(self, rho):
        """
        Get pressure for given density
        
        Parameters
        ----------
        rho : float
            Mass density in g/cm^3
            
        Returns
        -------
        P : float
            Pressure in dyne/cm^2
        """
        if self.backend == 'analytical':
            return self._pressure_analytical(rho)
        elif self.backend == 'tabulated':
            if rho < self.table_rho[0] or rho > self.table_rho[-1]:
                return 0.0  # Out of table range
            return float(self.interp_P(rho))
    
    def energy_density(self, rho):
        """
        Get energy density for given density
        
        Parameters
        ----------
        rho : float
            Mass density in g/cm^3
            
        Returns
        -------
        eps : float
            Energy density in g/cm^3
        """
        if self.backend == 'analytical':
            return self._energy_density_analytical(rho)
        elif self.backend == 'tabulated':
            if rho < self.table_rho[0] or rho > self.table_rho[-1]:
                return 0.0
            return float(self.interp_eps(rho))
    
    def density_from_pressure(self, P):
        """
        Get density for given pressure (inverse function)
        
        Parameters
        ----------
        P : float
            Pressure in dyne/cm^2
            
        Returns
        -------
        rho : float
            Mass density in g/cm^3
        """
        if P <= 0:
            return 0.0
        
        if self.backend == 'analytical':
            # For polytrope: rho = (P/K)^(1/Gamma)
            return (P / self.K)**(1.0 / self.Gamma)
        elif self.backend == 'tabulated':
            if P < self.table_P[0] or P > self.table_P[-1]:
                return 0.0
            return float(self.interp_rho(P))
    
    def eos_function(self):
        """
        Return EOS as callable for the TOV solver.

        Returns
        -------
        eos : callable
            Function  P (dyne/cm^2) -> (rho_total, rho_rest)
              rho_total : float  — total energy density / c^2  (g/cm^3)
              rho_rest  : float  — baryon rest-mass density    (g/cm^3)
        """
        def eos(P):
            if P <= 0:
                return 0.0, 0.0
            rho_rest = self.density_from_pressure(P)
            if self.backend == 'tabulated':
                # table column 3 is specific internal energy (erg/g),
                # not volumetric energy density — see extraction script.
                eps_specific = self.energy_density(rho_rest)   # erg/g
                rho_total = rho_rest * (1.0 + eps_specific / c2)
            else:
                # analytical _energy_density already returns epsilon/c^2
                rho_total = self.energy_density(rho_rest)      # g/cm^3
            return rho_total, rho_rest
        return eos


class QuarkStarEOS:
    """
    Quark star equations of state with analytical or tabulated backend
    
    Uses MIT bag model for analytical backend
    """
    
    def __init__(self, backend='analytical', B=60.0, table_path=None, **kwargs):
        self.backend = backend
        self.B = B  # Bag constant in MeV/fm^3
        
        if backend == 'analytical':
            self._init_analytical(B)
        elif backend == 'tabulated':
            if kwargs.get('generate_from_analytical', False):
                self._generate_table_from_analytical(B)
            else:
                self._load_table(table_path)
        else:
            raise ValueError(f"Unknown backend: {backend}")
    
    def _init_analytical(self, B):
        """Initialize analytical MIT bag model"""
        # Convert B to CGS units
        self.B_cgs = B * MeV_fm3_to_dyne_cm2
    
    def _generate_table_from_analytical(self, B):
        """Generate table from analytical MIT bag model"""
        print(f"Generating quark star table from MIT bag model (B={B} MeV/fm^3)...")
        
        # Switch to analytical temporarily
        original_backend = self.backend
        self.backend = 'analytical'
        self._init_analytical(B)
        
        # Generate table
        P_min = 1e34  # dyne/cm^2
        P_max = 1e37
        n_points = 1000
        
        self.table_P = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
        self.table_rho = np.zeros_like(self.table_P)
        self.table_eps = np.zeros_like(self.table_P)
        
        for i, P in enumerate(self.table_P):
            eps = 3.0 * P + 4.0 * self.B_cgs
            rho = eps / c2
            self.table_rho[i] = rho
            self.table_eps[i] = rho  # eps/c^2 in mass density units
        
        # Setup interpolators
        self._setup_interpolators()
        
        self.backend = original_backend
        print(f"Generated quark star table with {n_points} points")
    
    def _load_table(self, table_path):
        """Load table from file"""
        # Similar to NeutronStarEOS._load_table
        # Implementation would go here
        raise NotImplementedError("Loading quark star tables from file not yet implemented")
    
    def _setup_interpolators(self):
        """Setup interpolation functions"""
        # rho(P), eps(P)
        self.interp_rho = CubicSpline(self.table_P, self.table_rho, extrapolate=False)
        self.interp_eps = CubicSpline(self.table_P, self.table_eps, extrapolate=False)
        
        # Inverse: P(rho)
        self.interp_P = interp1d(self.table_rho, self.table_P,
                                 kind='cubic',
                                 bounds_error=False,
                                 fill_value=(self.table_P[0], self.table_P[-1]))
    
    def pressure(self, rho):
        """Get pressure for given density"""
        if self.backend == 'analytical':
            # For MIT bag: P = (1/3)(eps - 4B)
            # eps = rho * c^2, so P = (1/3)(rho*c^2 - 4B)
            if rho <= 0:
                return 0.0
            return (rho * c2 - 4.0 * self.B_cgs) / 3.0
        elif self.backend == 'tabulated':
            if rho < self.table_rho[0] or rho > self.table_rho[-1]:
                return 0.0
            return float(self.interp_P(rho))
    
    def energy_density(self, rho):
        """Get energy density"""
        if self.backend == 'analytical':
            if rho <= 0:
                return 0.0
            # For ultra-relativistic quarks
            return rho
        elif self.backend == 'tabulated':
            if rho < self.table_rho[0] or rho > self.table_rho[-1]:
                return 0.0
            return float(self.interp_eps(rho))
    
    def density_from_pressure(self, P):
        """Get density from pressure (inverse)"""
        if P <= 0:
            return 0.0
        
        if self.backend == 'analytical':
            # From P = (1/3)(rho*c^2 - 4B): rho = (3P + 4B)/c^2
            eps = 3.0 * P + 4.0 * self.B_cgs
            return eps / c2
        elif self.backend == 'tabulated':
            if P < self.table_P[0] or P > self.table_P[-1]:
                return 0.0
            return float(self.interp_rho(P))
    
    def eos_function(self):
        """
        Return EOS as callable for the TOV solver.

        Returns
        -------
        eos : callable
            Function  P (dyne/cm^2) -> (rho_total, rho_rest)
              rho_total : float  — epsilon/c^2, total energy density  (g/cm^3)
              rho_rest  : float  — 4B/c^2, bag vacuum energy density  (g/cm^3)

        For the MIT bag model, epsilon = 3P + 4B:
          - 4B   = QCD vacuum (bag) energy — structural ground-state mass
          - 3P   = ultrarelativistic quark kinetic energy — generates pressure
          rho_total = (3P + 4B) / c^2
          rho_rest  = 4B / c^2  (bag energy only, no kinetic contribution)

        This parallels the NS convention where rho_rest = baryon rest-mass
        density (no kinetic/thermal energy).  Both strip the internal energy
        that generates pressure and keep only the ground-state structural
        energy.  The difference M_g - M_b then measures pressure's
        contribution to gravitating mass for both matter types.
        """
        rho_rest_const = 4.0 * self.B_cgs / c2   # bag energy only

        def eos(P):
            if P <= 0:
                return 0.0, 0.0
            rho_total = self.density_from_pressure(P)   # (3P + 4B)/c^2
            return rho_total, rho_rest_const
        return eos


class HybridEOS:
    """
    Hybrid neutron/quark star EOS with phase transition
    
    This EOS automatically switches between neutron matter and quark matter
    phases based on pressure. It's designed for modeling hybrid stars with
    quark cores surrounded by neutron envelopes.
    
    The phase transition occurs at P_transition:
    - P > P_transition: quark matter phase (core)
    - P < P_transition: neutron matter phase (envelope)
    
    Parameters
    ----------
    neutron_eos : NeutronStarEOS
        EOS for the neutron matter phase
    quark_eos : QuarkStarEOS
        EOS for the quark matter phase  
    P_transition : float
        Transition pressure in dyne/cm^2
        
    Notes
    -----
    The TOV solver integrates from center outward using P -> (rho, eps).
    This class provides a clean interface where the phase is automatically
    determined by the local pressure, with no modifications needed to the
    TOV solver.
    
    For TOV integration:
    1. Start with high central pressure (quark phase)
    2. Integrate outward, pressure decreases
    3. When P drops below P_transition, automatically switch to neutron phase
    4. Continue to surface
    
    This produces hybrid stars naturally: quark core + neutron envelope.
    """
    
    def __init__(self, neutron_eos, quark_eos, P_transition):
        self.neutron_eos = neutron_eos
        self.quark_eos = quark_eos
        self.P_transition = P_transition
        
        # Cache transition point properties for diagnostics
        self.rho_transition_neutron = neutron_eos.density_from_pressure(P_transition)
        self.rho_transition_quark = quark_eos.density_from_pressure(P_transition)
        
        print(f"\nHybrid EOS initialized:")
        print(f"  P_transition = {P_transition:.3e} dyne/cm^2")
        print(f"  rho (neutron side) = {self.rho_transition_neutron:.3e} g/cm^3")
        print(f"  rho (quark side) = {self.rho_transition_quark:.3e} g/cm^3")
        
        # Check for density discontinuity
        if abs(self.rho_transition_neutron - self.rho_transition_quark) > 1e10:
            delta_rho = abs(self.rho_transition_neutron - self.rho_transition_quark)
            print(f"  WARNING: Density discontinuity of {delta_rho:.3e} g/cm^3 at transition")
    
    def density_from_pressure(self, P):
        """
        Get density for given pressure
        
        This is the critical method for TOV integration.
        
        Parameters
        ----------
        P : float
            Pressure in dyne/cm^2
            
        Returns
        -------
        rho : float
            Mass density in g/cm^3
        """
        if P <= 0:
            return 0.0
            
        if P > self.P_transition:
            # Quark phase (high pressure, core)
            return self.quark_eos.density_from_pressure(P)
        else:
            # Neutron phase (lower pressure, envelope)
            return self.neutron_eos.density_from_pressure(P)
    
    def energy_density(self, rho):
        """
        Get energy density for given density
        
        Note: At the phase transition, there may be a density discontinuity.
        This method determines phase based on which side of the transition
        the density corresponds to.
        
        Parameters
        ----------
        rho : float
            Mass density in g/cm^3
            
        Returns  
        -------
        eps : float
            Energy density in g/cm^3
        """
        if rho <= 0:
            return 0.0
        
        # Determine which phase we're in by checking pressure
        P_neutron = self.neutron_eos.pressure(rho)
        P_quark = self.quark_eos.pressure(rho)
        
        # Use the phase that's actually active at this density
        # (the one that gives pressure consistent with being on its side of transition)
        if P_quark > self.P_transition and rho > self.rho_transition_neutron:
            # High density, high pressure -> quark phase
            return self.quark_eos.energy_density(rho)
        else:
            # Lower density -> neutron phase
            return self.neutron_eos.energy_density(rho)
    
    def pressure(self, rho):
        """
        Get pressure for given density
        
        Note: This is the inverse of density_from_pressure, which may not be
        uniquely defined at phase transitions. This implementation chooses the
        phase based on density thresholds.
        
        Parameters
        ----------
        rho : float
            Mass density in g/cm^3
            
        Returns
        -------
        P : float
            Pressure in dyne/cm^2
        """
        if rho <= 0:
            return 0.0
        
        # Use same logic as energy_density for consistency
        P_neutron = self.neutron_eos.pressure(rho)
        P_quark = self.quark_eos.pressure(rho)
        
        if P_quark > self.P_transition and rho > self.rho_transition_neutron:
            return P_quark
        else:
            return P_neutron
    
    def eos_function(self):
        """
        Return EOS as callable for the TOV solver.

        Delegates to whichever phase is active at the given pressure.

        Returns
        -------
        eos : callable
            Function  P (dyne/cm^2) -> (rho_total, rho_rest)
        """
        ns_eos = self.neutron_eos.eos_function()
        qs_eos = self.quark_eos.eos_function()

        def eos(P):
            if P <= 0:
                return 0.0, 0.0
            if P > self.P_transition:
                return qs_eos(P)
            else:
                return ns_eos(P)
        return eos

    def get_phase_at_pressure(self, P):
        """
        Diagnostic: determine which phase is active at given pressure
        
        Parameters
        ----------
        P : float
            Pressure in dyne/cm^2
            
        Returns
        -------
        phase : str
            'quark' or 'neutron'
        """
        return 'quark' if P > self.P_transition else 'neutron'


# Backward compatibility: factory functions
def get_neutron_eos(backend='analytical', eos_type='validated', **kwargs):
    """Get neutron star EOS"""
    return NeutronStarEOS(backend=backend, eos_type=eos_type, **kwargs)


def get_quark_eos(backend='analytical', B=60.0, **kwargs):
    """Get quark star EOS"""
    return QuarkStarEOS(backend=backend, B=B, **kwargs)


def get_hybrid_eos(neutron_eos=None, quark_eos=None, P_transition=None, **kwargs):
    """
    Get hybrid neutron/quark star EOS
    
    Parameters
    ----------
    neutron_eos : NeutronStarEOS, optional
        If None, creates validated neutron EOS
    quark_eos : QuarkStarEOS, optional  
        If None, creates standard quark EOS (B=60 MeV/fm^3)
    P_transition : float, optional
        Transition pressure in dyne/cm^2
        If None, estimates from literature (~2-5 times nuclear saturation)
    """
    if neutron_eos is None:
        neutron_eos = get_neutron_eos(backend='analytical', eos_type='validated')
    
    if quark_eos is None:
        quark_eos = get_quark_eos(backend='analytical', B=60.0)
    
    if P_transition is None:
        # Estimate transition pressure
        # Nuclear saturation: rho_0 ~ 2.8e14 g/cm^3
        # Typical P at rho_0 for polytrope: P ~ K * rho^Gamma ~ 2.34e5 * (2.8e14)^2 ~ 1.8e34
        # Deconfinement typically occurs at 2-5 times nuclear density
        # So P_transition ~ 1e35 dyne/cm^2 is reasonable
        P_transition = 1.0e35  # dyne/cm^2
        print(f"Using default P_transition = {P_transition:.3e} dyne/cm^2")
    
    return HybridEOS(neutron_eos, quark_eos, P_transition)


def test_backends():
    """Test both analytical and tabulated backends"""
    import matplotlib.pyplot as plt
    
    print("\n=== Testing NeutronStarEOS backends ===\n")
    
    # Create EOS with both backends
    eos_analytical = NeutronStarEOS(backend='analytical', eos_type='validated')
    eos_tabulated = NeutronStarEOS(backend='tabulated', eos_type='generate_from_analytical',
                                   Gamma=2.0, K=2.34e5)
    
    # Test pressure range
    rho_range = np.logspace(14, 15.5, 100)  # g/cm^3
    
    P_analytical = [eos_analytical.pressure(rho) for rho in rho_range]
    P_tabulated = [eos_tabulated.pressure(rho) for rho in rho_range]
    
    # Plot comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.loglog(rho_range, P_analytical, 'b-', label='Analytical', linewidth=2)
    ax1.loglog(rho_range, P_tabulated, 'r--', label='Tabulated', linewidth=2)
    ax1.set_xlabel('Density (g/cm³)')
    ax1.set_ylabel('Pressure (dyne/cm²)')
    ax1.set_title('Pressure vs Density')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Relative difference
    rel_diff = [(p_t - p_a) / p_a * 100 for p_a, p_t in zip(P_analytical, P_tabulated)]
    ax2.semilogx(rho_range, rel_diff, 'g-', linewidth=2)
    ax2.set_xlabel('Density (g/cm³)')
    ax2.set_ylabel('Relative Difference (%)')
    ax2.set_title('Tabulated vs Analytical')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('eos_backend_comparison.png', dpi=150)
    print(f"\nComparison plot saved to eos_backend_comparison.png")
    print(f"Maximum relative difference: {max(abs(d) for d in rel_diff):.6f}%")


# =============================================================================
# BACKWARD COMPATIBILITY: STANDARD_EOS dictionary
# =============================================================================
# For scripts that expect the old STANDARD_EOS dictionary interface

def _quark_standard_eos():
    """Standard quark star EOS with B=60 MeV/fm^3"""
    qs = QuarkStarEOS(backend='analytical', B=60.0)
    return qs.eos_function()

def _neutron_validated_eos():
    """Validated neutron star EOS (polytrope with K=2.34e5, Gamma=2.0)"""
    ns = NeutronStarEOS(backend='analytical', eos_type='validated')
    return ns.eos_function()

def _neutron_stiff_eos():
    """Stiff neutron star EOS (polytrope with K=3.0e4, Gamma=2.5)"""
    ns = NeutronStarEOS(backend='analytical', eos_type='stiff')
    return ns.eos_function()

# Dictionary of standard EOS for backward compatibility
STANDARD_EOS = {
    'quark_standard': _quark_standard_eos(),
    'neutron_validated': _neutron_validated_eos(),
    'neutron_stiff': _neutron_stiff_eos()
}


if __name__ == '__main__':
    test_backends()
