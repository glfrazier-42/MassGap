"""
TOV Solver using scipy - FIXED VERSION

Corrected numerical issues:
- Better starting radius (1 cm instead of 1e-10 cm)
- Relaxed tolerances for stability
- Improved event detection
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Callable, Tuple, Optional

try:
    from .tov_solver import Constants, TOVResult
except ImportError:
    from tov_solver import Constants, TOVResult


class TOVSolverScipy:
    """scipy-based TOV solver with fixed numerical stability"""
    
    def __init__(self, eos: Callable[[float], Tuple[float, float]], method: str = 'DOP853'):
        self.eos = eos
        self.method = method
        self.c = Constants()
        self.last_sol = None
    
    def tov_derivatives(self, r: float, y: np.ndarray) -> np.ndarray:
        """TOV equations"""
        if r < 0.1:  # Very close to origin
            return np.array([0.0, 0.0, 0.0])
        
        m, P, m_rest = y
        
        if P <= 0:
            return np.array([0.0, 0.0, 0.0])
        
        rho, epsilon = self.eos(P)
        
        # dm/dr
        dmdr = 4.0 * np.pi * r**2 * rho
        
        # dP/dr
        numerator = ((rho + P/self.c.c2) * (m + 4.0*np.pi*r**3 * P/self.c.c2) * self.c.G)
        denominator = r * (r - 2.0*self.c.G*m/self.c.c2)
        
        if denominator <= 0:
            return np.array([0.0, 0.0, 0.0])
        
        dPdr = -numerator / denominator
        
        # dm_rest/dr
        metric_factor = np.sqrt(abs(1.0 - 2.0*self.c.G*m/(r*self.c.c2)))
        dm_rest_dr = 4.0 * np.pi * r**2 * rho * metric_factor
        
        return np.array([dmdr, dPdr, dm_rest_dr])
    
    def surface_event(self, r: float, y: np.ndarray) -> float:
        """Event: detect when P drops below small threshold"""
        return y[1] - 1e28  # Stop when P < 1e28 (essentially zero)
    
    surface_event.terminal = True
    surface_event.direction = -1
    
    def solve(self, P_central: float, r_max: float = 3e6) -> Optional[TOVResult]:
        """Solve TOV equations"""
        
        # Get central density
        rho_central, _ = self.eos(P_central)
        
        # Start at 1 cm (numerically stable, essentially at center)
        r_start = 1.0  # cm
        m_start = (4.0/3.0) * np.pi * r_start**3 * rho_central
        P_start = P_central
        m_rest_start = m_start
        
        y0 = np.array([m_start, P_start, m_rest_start])
        r_span = (r_start, r_max)
        
        try:
            sol = solve_ivp(
                fun=self.tov_derivatives,
                t_span=r_span,
                y0=y0,
                method=self.method,
                events=self.surface_event,
                dense_output=True,
                rtol=1e-8,   # Relaxed from 1e-10
                atol=1e-10,  # Relaxed from 1e-12
                max_step=1e4
            )
            
            self.last_sol = sol
            
            if not sol.success:
                print(f"Integration failed: {sol.message}")
                return None
            
            # Get surface radius
            if len(sol.t_events[0]) > 0:
                R = sol.t_events[0][0]
            else:
                # Use last point if event not triggered
                R = sol.t[-1]
                if sol.y[1, -1] > 1e30:  # Pressure still high
                    print(f"Warning: Surface not found")
                    return None
            
            # Evaluate at surface
            y_surface = sol.sol(R)
            M = y_surface[0]
            M_rest = y_surface[2]
            
            # Get profiles
            n_points = min(1000, len(sol.t) * 5)
            r_grid = np.linspace(0, R, n_points)
            y_grid = sol.sol(r_grid)
            
            r_array = r_grid
            m_array = y_grid[0]
            P_array = y_grid[1]
            m_rest_array = y_grid[2]
            
            rho_array = np.array([self.eos(max(P, 0))[0] for P in P_array])
            
            M_pressure = M - M_rest
            
            return TOVResult(
                radius=R,
                mass=M,
                mass_solar=M / self.c.M_sun,
                rest_mass=M_rest,
                rest_mass_solar=M_rest / self.c.M_sun,
                pressure_mass=M_pressure,
                pressure_mass_solar=M_pressure / self.c.M_sun,
                central_pressure=P_central,
                central_density=rho_central,
                radius_km=R / self.c.km,
                r_array=r_array,
                m_array=m_array,
                P_array=P_array,
                rho_array=rho_array
            )
            
        except Exception as e:
            print(f"Integration error: {e}")
            return None
    
    def find_maximum_mass(self, P_min: float = 1e33, P_max: float = 1e37, n_points: int = 50):
        """Find maximum mass"""
        P_range = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
        solutions = []
        
        print(f"Scanning with {self.method}...")
        
        for i, P_c in enumerate(P_range):
            if (i+1) % 10 == 0:
                print(f"  {i+1}/{n_points}")
            
            result = self.solve(P_c)
            if result is not None:
                solutions.append(result)
        
        if not solutions:
            raise ValueError("No valid solutions")
        
        print(f"Found {len(solutions)} solutions")
        
        masses = [s.mass_solar for s in solutions]
        max_idx = np.argmax(masses)
        
        return solutions[max_idx], solutions


def test_scipy_solver():
    """Test the solver"""
    try:
        from .eos import NeutronStarEOS
    except ImportError:
        from eos import NeutronStarEOS
    
    print("\n=== Testing FIXED scipy TOV Solver ===\n")
    
    eos_obj = NeutronStarEOS(backend='analytical', eos_type='validated')
    eos = eos_obj.eos_function()
    
    P_central = 5e34
    
    print(f"Solving with DOP853...")
    solver = TOVSolverScipy(eos, method='DOP853')
    result = solver.solve(P_central)
    
    if result is not None:
        print("\n" + str(result))
        print(f"\nPressure fraction: {result.pressure_mass_solar/result.mass_solar*100:.1f}%")
    else:
        print("Failed!")


if __name__ == '__main__':
    test_scipy_solver()
