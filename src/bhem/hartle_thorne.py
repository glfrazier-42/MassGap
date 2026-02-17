"""
Hartle-Thorne Slow-Rotation Formalism
======================================

Computes the stability boundary M_max(P_spin) for rotating neutron
stars, combining the TOV maximum mass with the quasi-universal
rotational enhancement from Breu & Rezzolla (2016).

Also provides a moment-of-inertia calculation via the Hartle (1967)
frame-dragging ODE for validation.

References
----------
Hartle, J. B. 1967, ApJ, 150, 1005
Breu, C. & Rezzolla, L. 2016, MNRAS, 459, 646
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

# Physical constants (CGS)
G = 6.67430e-8          # cm^3 g^-1 s^-2
c_light = 2.99792458e10 # cm s^-1
c2 = c_light**2
M_sun = 1.98847e33      # g
km = 1e5                # cm


# ---------------------------------------------------------------
#  TOV scan for M_TOV and R_TOV from a tabulated EOS
# ---------------------------------------------------------------

def _find_max_mass(eos_path, P_min=1e34, P_max=5e36, n_scan=80):
    """
    Find the maximum-mass TOV configuration for a tabulated EOS.

    Uses NeutronStarEOS (tabulated backend) + TOVSolverScipy — the same
    validated code path used by plot_tov_figures.py.

    Parameters
    ----------
    eos_path : str or Path
        Three-column table (rho, P, eps_specific) in CGS.
    P_min, P_max : float
        Central-pressure scan range (dyne/cm^2).
    n_scan : int
        Number of trial pressures.

    Returns
    -------
    M_tov : float   Maximum gravitational mass (g).
    R_tov : float   Radius at M_TOV (cm).
    """
    from .eos import NeutronStarEOS
    from .tov_solver_scipy import TOVSolverScipy

    ns = NeutronStarEOS(backend='tabulated', table_path=str(eos_path))
    solver = TOVSolverScipy(ns.eos_function())

    best_M, best_R = 0.0, 0.0
    for Pc in np.logspace(np.log10(P_min), np.log10(P_max), n_scan):
        res = solver.solve(Pc)
        if res is not None and res.mass > best_M:
            best_M = res.mass
            best_R = res.radius

    if best_M == 0:
        raise RuntimeError("No converged TOV solution for %s" % eos_path)

    return best_M, best_R


# ---------------------------------------------------------------
#  Kepler (mass-shedding) period
# ---------------------------------------------------------------

def kepler_period_ms(M_g, R_cm):
    """
    Kepler period in milliseconds, including an empirical GR
    correction (Haensel, Zdunik & Bejger 2009).

        P_K ~ 0.77 ms (M/M_sun)^{-1/2} (R/10 km)^{3/2}
    """
    return 0.77 * (M_g / M_sun)**(-0.5) * (R_cm / (10 * km))**1.5


# ---------------------------------------------------------------
#  Public API: stability boundary
# ---------------------------------------------------------------

def stability_boundary(eos_path, periods_ms):
    """
    Maximum stable NS mass as a function of spin period.

    Uses the quasi-universal relation (Breu & Rezzolla 2016)::

        M_max(f) = M_TOV [1 + a2 (f/f_K)^2 + a4 (f/f_K)^4]

    where f = 1/P is the spin frequency and f_K the Kepler
    (mass-shedding) frequency of the maximum-mass static star.

    Parameters
    ----------
    eos_path : str or Path
        Tabulated EOS file (3-column CGS).
    periods_ms : array-like
        Spin periods in milliseconds.

    Returns
    -------
    M_max : ndarray
        Maximum stable mass (M_sun) at each period.
    info : dict
        M_tov (M_sun), R_tov_km, P_K_ms.
    """
    M_g, R_cm = _find_max_mass(eos_path)
    M_tov = M_g / M_sun
    R_km = R_cm / km
    P_K = kepler_period_ms(M_g, R_cm)

    # Breu & Rezzolla (2016) quasi-universal coefficients
    a2, a4 = 0.132, 0.071

    P = np.asarray(periods_ms, dtype=float)
    chi = P_K / P                    # = f/f_K = Omega/Omega_K
    M_max = M_tov * (1 + a2 * chi**2 + a4 * chi**4)
    M_max[chi > 1] = np.nan         # past mass-shedding

    return M_max, dict(M_tov=M_tov, R_tov_km=R_km, P_K_ms=P_K)


# ---------------------------------------------------------------
#  Moment of inertia via Hartle (1967) frame-dragging ODE
# ---------------------------------------------------------------

def compute_moment_of_inertia(tov_result):
    """
    Moment of inertia from the Hartle (1967) frame-dragging ODE.

    Solves::

        d/dr [r^4 j(r) d omega_tilde/dr]
            + 4 r^3 (dj/dr) omega_tilde = 0

    where j(r) = exp(-Phi) sqrt(1 - 2Gm/rc^2), and extracts I
    from exterior matching at the stellar surface.

    Parameters
    ----------
    tov_result : TOVResult
        Solved TOV configuration (from TOVSolver or TOVSolverScipy).

    Returns
    -------
    I : float or None
        Moment of inertia in g cm^2, or None on failure.
    """
    r = tov_result.r_array
    m = tov_result.m_array
    P = tov_result.P_array
    R = tov_result.radius
    M = tov_result.mass

    # ---- 1. metric potential Phi(r) ----
    n = len(r)
    Phi = np.zeros(n)
    for i in range(1, n):
        ri, mi, Pi = r[i], m[i], P[i]
        if ri < 1.0:
            continue
        denom = ri * (ri * c2 - 2 * G * mi)
        if denom <= 0:
            Phi[i] = Phi[i - 1]
            continue
        dPhi = G * (mi + 4 * np.pi * ri**3 * Pi / c2) / denom
        Phi[i] = Phi[i - 1] + dPhi * (r[i] - r[i - 1])

    # Fix surface boundary: Phi(R) = 1/2 ln(1 - 2GM/Rc^2)
    compactness = 2 * G * M / (R * c2)
    if compactness >= 1:
        return None
    Phi += 0.5 * np.log(1 - compactness) - Phi[-1]

    # ---- 2. j(r) = exp(-Phi) sqrt(1 - 2Gm/rc^2) ----
    mask = r > 100.0
    r_g = r[mask]
    m_g = m[mask]
    Phi_g = Phi[mask]
    f = np.maximum(1 - 2 * G * m_g / (r_g * c2), 1e-30)
    j_arr = np.exp(-Phi_g) * np.sqrt(f)
    j_spl = CubicSpline(r_g, j_arr)

    # ---- 3. solve frame-dragging ODE ----
    #   d omega_tilde / dr = xi
    #   d xi / dr = -(4/r + j'/j) xi - 4 j'/(r j) omega_tilde
    def rhs(rv, y):
        ob, xi = y
        jv = float(j_spl(rv))
        djv = float(j_spl(rv, 1))
        if abs(jv) < 1e-40:
            return [0.0, 0.0]
        return [xi,
                -(4.0 / rv + djv / jv) * xi
                - 4.0 * djv / (rv * jv) * ob]

    r0, r1 = r_g[0], R * 0.999
    sol = solve_ivp(rhs, (r0, r1), [1.0, 0.0],
                    method='DOP853', rtol=1e-10, atol=1e-14,
                    dense_output=True, max_step=1e4)
    if not sol.success:
        return None

    # ---- 4. extract I from exterior matching ----
    #   omega_tilde(R) = Omega - 2GJ/(c^2 R^3)
    #   d omega_tilde/dr|_R = 6GJ/(c^2 R^4)
    #   => I = c^2 R^4 xi_R / (6G (omega_tilde_R + R xi_R / 3))
    om_R = sol.sol(r1)[0]
    xi_R = sol.sol(r1)[1]
    d = om_R + r1 * xi_R / 3.0
    if abs(d) < 1e-30:
        return None

    return c2 * r1**4 * xi_R / (6.0 * G * d)
