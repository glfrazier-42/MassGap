"""
BHEM (Black Hole Explosion Mechanism) Package

Numerical tools for analyzing compact star structure and mass gaps
"""

from .tov_solver import TOVSolver, TOVResult, calculate_pressure_contribution
from .eos import NeutronStarEOS, QuarkStarEOS, get_neutron_eos, get_quark_eos

__version__ = '0.1.0'
__all__ = [
    'TOVSolver',
    'TOVResult', 
    'calculate_pressure_contribution',
    'NeutronStarEOS',
    'QuarkStarEOS',
    'get_neutron_eos',
    'get_quark_eos',
]
