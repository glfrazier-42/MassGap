
"""
TOV (Tolman-Oppenheimer-Volkoff) Equation Solver

Solves the relativistic structure equations for compact stars:
- Neutron stars (Tier 1)
- Quark stars (Tier 2)
- Higher tier configurations

Key focus: Calculating pressure contribution to gravitational mass
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from dataclasses import dataclass
from typing import Callable, Tuple, Optional


# Physical constants (CGS units)
class Constants:
    """Physical constants in CGS units"""
    G = 6.67430e-8          # Gravitational constant (cm^3 g^-1 s^-2)
    c = 2.99792458e10       # Speed of light (cm/s)
    c2 = c**2               # c^2
    M_sun = 1.98847e33      # Solar mass (g)
    km = 1e5                # km to cm conversion


@dataclass
class TOVResult:
    """Results from TOV integration"""
    radius: float           # Stellar radius (cm)
    mass: float            # Gravitational mass (g)
    mass_solar: float      # Gravitational mass (M_sun)
    rest_mass: float       # Baryonic rest mass (g)
    rest_mass_solar: float # Rest mass (M_sun)
    pressure_mass: float   # Pressure contribution to mass (g)
    pressure_mass_solar: float  # Pressure mass (M_sun)
    central_pressure: float     # Central pressure (dyne/cm^2)
    central_density: float      # Central density (g/cm^3)
    radius_km: float           # Radius in km
    r_array: np.ndarray        # Radial coordinates
    m_array: np.ndarray        # Mass profile m(r)
    P_array: np.ndarray        # Pressure profile P(r)
    rho_array: np.ndarray      # Density profile rho(r)
    
    def __str__(self):
        return (f"TOV Solution:\n"
                f"  R = {self.radius_km:.2f} km\n"
                f"  M_grav = {self.mass_solar:.3f} M_sun\n"
                f"  M_rest = {self.rest_mass_solar:.3f} M_sun\n"
                f"  M_pressure = {self.pressure_mass_solar:.3f} M_sun\n"
                f"  P_central = {self.central_pressure:.3e} dyne/cm^2\n"
                f"  rho_central = {self.central_density:.3e} g/cm^3")


class TOVSolver:
    """
    Solves TOV equations for compact star structure
    
    The TOV equations in geometric units with G=c=1:
    dm/dr = 4πr²ρ(r)
    dP/dr = -(ρ + P)(m + 4πr³P) / (r² - 2m/r)
    
    With restoration of G and c:
    dm/dr = 4πr²ρ(r)
    dP/dr = -(ρ + P/c²)(m + 4πr³P/c²) * G / (r(r - 2Gm/c²))
    """
    
    def __init__(self, eos: Callable[[float], Tuple[float, float]]):
        """
        Initialize TOV solver with equation of state
        
        Parameters:
        -----------
        eos : callable
            Function that takes pressure P (dyne/cm^2) and returns:
            (density rho in g/cm^3, energy_density epsilon in g/cm^3)
        """
        self.eos = eos
        self.c = Constants()
        
    def tov_equations(self, r: float, y: np.ndarray) -> np.ndarray:
        """
        TOV differential equations
        
        Parameters:
        -----------
        r : float
            Radial coordinate (cm)
        y : array [m, P, m_rest]
            m: enclosed gravitational mass (g)
            P: pressure (dyne/cm^2)
            m_rest: enclosed rest mass (g)
            
        Returns:
        --------
        dydr : array [dm/dr, dP/dr, dm_rest/dr]
        """
        if r < 1e-10:  # Avoid singularity at origin
            return np.array([0.0, 0.0, 0.0])
        
        m, P, m_rest = y
        
        # Check for surface (P <= 0)
        if P <= 0:
            return np.array([0.0, 0.0, 0.0])
        
        # Get density from EOS
        rho, epsilon = self.eos(P)
        
        # TOV equations with proper units
        # dm/dr = 4πr²ρ
        dmdr = 4.0 * np.pi * r**2 * rho
        
        # dP/dr = -(ρ + P/c²)(m + 4πr³P/c²) * G / (r(r - 2Gm/c²))
        numerator = (rho + P/self.c.c2) * (m + 4.0*np.pi*r**3 * P/self.c.c2) * self.c.G
        denominator = r * (r - 2.0*self.c.G*m/self.c.c2)
        
        # Avoid singularity near Schwarzschild radius
        if denominator <= 0:
            return np.array([0.0, 0.0, 0.0])
        
        dPdr = -numerator / denominator
        
        # Rest mass: dm_rest/dr = 4πr²ρ * sqrt(1 - 2Gm/(rc²))
        # This accounts for gravitational redshift
        metric_factor = np.sqrt(abs(1.0 - 2.0*self.c.G*m/(r*self.c.c2)))
        dm_rest_dr = 4.0 * np.pi * r**2 * rho * metric_factor
        
        return np.array([dmdr, dPdr, dm_rest_dr])
    
    def solve(self, P_central: float, dr: float = 1e3, max_radius: float = 3e6) -> Optional[TOVResult]:
        """
        Solve TOV equations for given central pressure
        
        Parameters:
        -----------
        P_central : float
            Central pressure (dyne/cm^2)
        dr : float
            Integration step size (cm), default 1 km
        max_radius : float
            Maximum integration radius (cm), default 30 km
            
        Returns:
        --------
        TOVResult or None if integration fails
        """
        
        # Initial conditions at center
        m0 = 0.0
        P0 = P_central
        m_rest0 = 0.0
        y0 = np.array([m0, P0, m_rest0])
        
        # Set up radial grid
        r_grid = np.arange(dr, max_radius, dr)
        
        # Arrays to store results
        r_arr = [0.0]
        m_arr = [m0]
        P_arr = [P0]
        m_rest_arr = [m_rest0]
        rho_central, _ = self.eos(P0)
        rho_arr = [rho_central]
        
        # Integrate outward
        y = y0
        for r in r_grid:
            # One step of integration
            dydr = self.tov_equations(r, y)
            y_new = y + dydr * dr
            
            # Check for surface (P <= 0)
            if y_new[1] <= 0:
                # Found surface - do final step to P=0
                if dydr[1] != 0:
                    dr_final = -y[1] / dydr[1]
                    if dr_final > 0 and dr_final < dr:
                        r_final = r + dr_final
                        y_final = y + dydr * dr_final
                    else:
                        r_final = r
                        y_final = y
                        y_final[1] = 0  # Set pressure to exactly zero
                else:
                    r_final = r
                    y_final = y
                    y_final[1] = 0
                
                # Store final point
                r_arr.append(r_final)
                m_arr.append(y_final[0])
                P_arr.append(0.0)
                m_rest_arr.append(y_final[2])
                rho_arr.append(0.0)
                break
            
            # Store current point
            r_arr.append(r)
            m_arr.append(y_new[0])
            P_arr.append(y_new[1])
            m_rest_arr.append(y_new[2])
            rho, _ = self.eos(y_new[1])
            rho_arr.append(rho)
            
            y = y_new
        else:
            # Didn't find surface within max_radius
            print(f"Warning: Surface not found within {max_radius/1e5:.1f} km")
            return None
        
        # Convert to numpy arrays
        r_arr = np.array(r_arr)
        m_arr = np.array(m_arr)
        P_arr = np.array(P_arr)
        m_rest_arr = np.array(m_rest_arr)
        rho_arr = np.array(rho_arr)
        
        # Final values
        R = r_arr[-1]
        M = m_arr[-1]
        M_rest = m_rest_arr[-1]
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
            r_array=r_arr,
            m_array=m_arr,
            P_array=P_arr,
            rho_array=rho_arr
        )
    
    def find_maximum_mass(self, P_min: float = 1e33, P_max: float = 1e37, 
                         n_points: int = 50) -> Tuple[TOVResult, np.ndarray]:
        """
        Find maximum stable mass by scanning central pressures
        
        Returns the configuration at maximum mass and array of all solutions
        """
        P_range = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
        solutions = []
        
        for P_c in P_range:
            result = self.solve(P_c)
            if result is not None:
                solutions.append(result)
        
        if not solutions:
            raise ValueError("No valid solutions found")
        
        # Find maximum mass
        masses = [s.mass_solar for s in solutions]
        max_idx = np.argmax(masses)
        
        return solutions[max_idx], solutions


def calculate_pressure_contribution(result: TOVResult) -> dict:
    """
    Detailed analysis of pressure contribution to gravitational mass
    
    Returns dictionary with breakdown of mass contributions
    """
    c = Constants()
    
    # Integrate pressure contribution: ∫ 4πr²(P/c²) dr with metric factors
    r = result.r_array
    P = result.P_array
    m = result.m_array
    
    # Metric factor: sqrt(1 - 2Gm/(rc²))
    metric = np.sqrt(np.maximum(0, 1.0 - 2.0*c.G*m/(r*c.c2 + 1e-100)))
    
    # Pressure contribution density
    P_contribution = P / c.c2
    
    # Integrate
    dr = np.diff(r)
    r_mid = 0.5 * (r[1:] + r[:-1])
    P_mid = 0.5 * (P_contribution[1:] + P_contribution[:-1])
    m_mid = 0.5 * (m[1:] + m[:-1])
    metric_mid = 0.5 * (metric[1:] + metric[:-1])
    
    dM_pressure = 4.0 * np.pi * r_mid**2 * P_mid * metric_mid * dr
    M_pressure_integrated = np.sum(dM_pressure)
    
    return {
        'M_pressure_total': result.pressure_mass,
        'M_pressure_integrated': M_pressure_integrated,
        'M_pressure_solar': result.pressure_mass_solar,
        'fraction_of_total': result.pressure_mass / result.mass,
        'pressure_profile': P,
        'radius_profile': r
    }
