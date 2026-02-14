"""
TOV Solver with Tabulated EOS Support

This module provides a simple interface for using tabulated equations of state
(from stellarcollapse.org or other sources) with the TOV solver.

Usage:
    from tov_tabulated import TabulatedEOS, TOVSolver
    
    # Load EOS from text file
    eos = TabulatedEOS('path/to/eos_file.txt')
    
    # Create solver
    solver = TOVSolver(eos)
    
    # Solve for star at given central pressure
    result = solver.solve(P_central=1e35)
    
    print(f"Mass: {result['M'] / M_SUN:.3f} M_sun")
    print(f"Radius: {result['R'] / 1e5:.2f} km")
"""

import sys
import os

# Add parent directory to path to import from bhem package
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
from scipy.interpolate import CubicSpline, interp1d

# Constants
M_SUN = 1.989e33  # g
C_LIGHT = 2.99792458e10  # cm/s
G_CONST = 6.67430e-8  # cm³/g/s²
KM = 1e5  # cm


class TabulatedEOS:
    """
    Equation of state loaded from tabulated data
    
    Expected file format (text file with 3+ columns):
    # Comments start with #
    # Column 1: density (g/cm³)
    # Column 2: pressure (dyne/cm²)
    # Column 3: energy density (erg/cm³)
    density1  pressure1  energy_density1
    density2  pressure2  energy_density2
    ...
    
    Parameters
    ----------
    filepath : str
        Path to EOS table file
    """
    
    def __init__(self, filepath):
        self.filepath = filepath
        self._load_table(filepath)
        self._setup_interpolators()
    
    def _load_table(self, filepath):
        """Load EOS table from file"""
        print(f"Loading EOS: {filepath}")
        
        # Load data (skip comment lines)
        data = np.loadtxt(filepath)
        
        if data.shape[1] < 2:
            raise ValueError(f"EOS file must have at least 2 columns (rho, P), got {data.shape[1]}")
        
        self.rho_table = data[:, 0]  # g/cm³
        self.P_table = data[:, 1]    # dyne/cm²
        
        # Energy density
        if data.shape[1] >= 3:
            self.eps_table = data[:, 2]  # erg/cm³
        else:
            # Estimate: ε ≈ ρ c² (non-relativistic approximation)
            print("  Warning: No energy density column, using eps = rho*c^2")
            self.eps_table = self.rho_table * C_LIGHT**2
        
        # Additional columns (e.g., Ye) can be stored but aren't used
        
        # Validate table
        if not np.all(np.diff(self.rho_table) > 0):
            raise ValueError("Density must be strictly increasing")
        if not np.all(np.diff(self.P_table) > 0):
            raise ValueError("Pressure must be strictly increasing")
        
        print(f"  Points: {len(self.rho_table)}")
        print(f"  rho range: {self.rho_table[0]:.2e} - {self.rho_table[-1]:.2e} g/cm^3")
        print(f"  P range: {self.P_table[0]:.2e} - {self.P_table[-1]:.2e} dyne/cm^2")
    
    def _setup_interpolators(self):
        """Create interpolation functions"""
        # Forward functions: P(ρ), ε(ρ)
        self.P_of_rho = CubicSpline(self.rho_table, self.P_table, extrapolate=False)
        self.eps_of_rho = CubicSpline(self.rho_table, self.eps_table, extrapolate=False)
        
        # Inverse function: ρ(P)
        # Use linear interpolation for stability (pressure is monotonic)
        self.rho_of_P = interp1d(
            self.P_table, 
            self.rho_table,
            kind='cubic',
            bounds_error=False,
            fill_value=(self.rho_table[0], self.rho_table[-1])
        )
    
    def pressure(self, rho):
        """
        Get pressure at given density
        
        Parameters
        ----------
        rho : float
            Mass density (g/cm³)
        
        Returns
        -------
        P : float
            Pressure (dyne/cm²)
        """
        if rho <= 0:
            return 0.0
        
        # Check bounds
        if rho < self.rho_table[0] or rho > self.rho_table[-1]:
            return 0.0
        
        return float(self.P_of_rho(rho))
    
    def energy_density(self, rho):
        """
        Get energy density at given density
        
        Parameters
        ----------
        rho : float
            Mass density (g/cm³)
        
        Returns
        -------
        eps : float
            Energy density (erg/cm³)
        """
        if rho <= 0:
            return 0.0
        
        # Check bounds
        if rho < self.rho_table[0] or rho > self.rho_table[-1]:
            return 0.0
        
        return float(self.eps_of_rho(rho))
    
    def density_from_pressure(self, P):
        """
        Get density for given pressure (inverse)
        
        Parameters
        ----------
        P : float
            Pressure (dyne/cm²)
        
        Returns
        -------
        rho : float
            Mass density (g/cm³)
        """
        if P <= 0:
            return 0.0
        
        # Check bounds
        if P < self.P_table[0] or P > self.P_table[-1]:
            return 0.0
        
        return float(self.rho_of_P(P))
    
    def __call__(self, P):
        """
        EOS callable: returns (rho, epsilon) for given pressure
        
        This is the interface expected by TOV solvers.
        
        Parameters
        ----------
        P : float
            Pressure (dyne/cm²)
        
        Returns
        -------
        rho : float
            Mass density (g/cm³)
        eps : float
            Energy density (g/cm³) - note: converted from erg/cm³
        """
        if P <= 0:
            return 0.0, 0.0
        
        rho = self.density_from_pressure(P)
        eps_erg = self.energy_density(rho)
        
        # Convert energy density from erg/cm³ to g/cm³ (divide by c²)
        eps = eps_erg / C_LIGHT**2
        
        return rho, eps


class TOVSolver:
    """
    Tolman-Oppenheimer-Volkoff equation solver using scipy
    
    Solves the general relativistic stellar structure equations to find
    mass, radius, and internal structure of compact stars.
    
    Parameters
    ----------
    eos : TabulatedEOS or callable
        Equation of state: eos(P) -> (rho, eps)
    """
    
    def __init__(self, eos):
        self.eos = eos
    
    def tov_derivatives(self, r, y):
        """
        TOV equation derivatives
        
        Parameters
        ----------
        r : float
            Radius (cm)
        y : array
            [m, P] where m is enclosed mass (g), P is pressure (dyne/cm²)
        
        Returns
        -------
        dy/dr : array
            [dm/dr, dP/dr]
        """
        m, P = y
        
        if P <= 0 or r <= 0:
            return [0.0, 0.0]
        
        # Get density and energy density from EOS
        rho, eps = self.eos(P)
        
        if rho <= 0:
            return [0.0, 0.0]
        
        # dm/dr = 4π r² ρ
        dmdr = 4.0 * np.pi * r**2 * rho
        
        # dP/dr = -[G(ρ + P/c²)(m + 4πr³P/c²)] / [r(r - 2Gm/c²)]
        numerator = G_CONST * (rho + P/C_LIGHT**2) * (m + 4.0*np.pi*r**3 * P/C_LIGHT**2)
        denominator = r * (r - 2.0*G_CONST*m/C_LIGHT**2)
        
        if denominator <= 0:
            # Approaching Schwarzschild radius
            return [0.0, 0.0]
        
        dPdr = -numerator / denominator
        
        return [dmdr, dPdr]
    
    def solve(self, P_central, method='DOP853', rtol=1e-8, atol=1e-10):
        """
        Solve TOV equations for given central pressure
        
        Parameters
        ----------
        P_central : float
            Central pressure (dyne/cm²)
        method : str
            Integration method (default: DOP853)
        rtol : float
            Relative tolerance
        atol : float
            Absolute tolerance
        
        Returns
        -------
        result : dict
            Dictionary with keys:
            - 'converged': bool - whether solution converged
            - 'M': float - total mass (g)
            - 'R': float - radius (cm)
            - 'm': array - enclosed mass profile (g)
            - 'r': array - radius grid (cm)
            - 'P': array - pressure profile (dyne/cm²)
            - 'rho': array - density profile (g/cm³)
        """
        from scipy.integrate import solve_ivp
        
        # Get central density
        rho_central, _ = self.eos(P_central)
        
        if rho_central <= 0:
            return {'converged': False}
        
        # Initial conditions at r = 1 cm (numerically stable)
        r_start = 1.0  # cm
        m_start = (4.0/3.0) * np.pi * r_start**3 * rho_central
        
        y0 = [m_start, P_central]
        
        # Integration range
        r_max = 3e6  # 30 km
        
        # Event: detect surface (P drops to nearly zero)
        def surface_event(r, y):
            return y[1] - 1e28  # P < 1e28 dyne/cm² (essentially zero)
        
        surface_event.terminal = True
        surface_event.direction = -1
        
        # Solve
        try:
            sol = solve_ivp(
                fun=self.tov_derivatives,
                t_span=(r_start, r_max),
                y0=y0,
                method=method,
                events=surface_event,
                dense_output=True,
                rtol=rtol,
                atol=atol,
                max_step=1e4
            )
            
            if not sol.success:
                return {'converged': False}
            
            # Get surface radius
            if len(sol.t_events[0]) > 0:
                R = sol.t_events[0][0]
            else:
                # Use last point if event didn't trigger
                R = sol.t[-1]
                if sol.y[1, -1] > 1e30:
                    # Pressure still high - didn't reach surface
                    return {'converged': False}
            
            # Get final mass
            M = sol.sol(R)[0]
            
            # Get profiles on uniform grid
            n_points = 500
            r_grid = np.linspace(0, R, n_points)
            y_grid = sol.sol(r_grid)
            
            m_grid = y_grid[0]
            P_grid = y_grid[1]
            
            # Get density profile
            rho_grid = np.array([self.eos(max(P, 0))[0] for P in P_grid])
            
            return {
                'converged': True,
                'M': M,
                'R': R,
                'r': r_grid,
                'm': m_grid,
                'P': P_grid,
                'rho': rho_grid
            }
            
        except Exception as e:
            print(f"Integration error: {e}")
            return {'converged': False}


def test_tabulated_eos():
    """Test the tabulated EOS and solver"""
    print("\n" + "="*70)
    print("TESTING TABULATED EOS AND TOV SOLVER")
    print("="*70 + "\n")
    
    # Example: load LS220 EOS (if it exists)
    import os
    from pathlib import Path
    
    base_dir = Path(__file__).parent.parent.parent
    eos_file = base_dir / 'data' / 'eos' / 'LS220_cold_betaeq.txt'
    
    if not eos_file.exists():
        print(f"EOS file not found: {eos_file}")
        print("Run extract_stellarcollapse_eos.py first!")
        return
    
    # Load EOS
    print("Loading EOS...")
    eos = TabulatedEOS(str(eos_file))
    print()
    
    # Create solver
    print("Creating TOV solver...")
    solver = TOVSolver(eos)
    print()
    
    # Solve for a star
    print("Solving for neutron star structure...")
    print("Central pressure: 1.0e35 dyne/cm^2")
    print()
    
    result = solver.solve(P_central=1e35)
    
    if result['converged']:
        M_solar = result['M'] / M_SUN
        R_km = result['R'] / KM
        
        print(f"OK Solution converged!")
        print(f"  Mass: {M_solar:.4f} M_sun")
        print(f"  Radius: {R_km:.2f} km")
        print(f"  Compactness: {G_CONST*result['M']/(result['R']*C_LIGHT**2):.3f}")
    else:
        print(f"FAIL Solution failed to converge")
    
    print()
    print("="*70)
    print("TEST COMPLETE")
    print("="*70)


if __name__ == '__main__':
    test_tabulated_eos()
