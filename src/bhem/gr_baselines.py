"""GR baseline predictions for binary black hole merger mass deficit.

Provides the non-spinning Buonanno, Kidder & Lehner (2008) formula and a
spinning extension based on the Kerr ISCO binding energy.

The non-spinning formula:
    f_GR = 0.0572 * eta + 0.498 * eta^2

is the leading-order (linear in eta) prediction from the ISCO binding
energy of a Schwarzschild black hole, plus a fitted quadratic correction.
The spinning extension replaces the Schwarzschild ISCO binding energy
with the analytic Kerr ISCO binding energy evaluated at the effective
progenitor spin chi_eff.  At chi_eff = 0 the two formulas are identical
to numerical precision.

Reference:
    Buonanno, Kidder & Lehner, Phys. Rev. D 77, 026004 (2008),
    arXiv:0709.3839
    Bardeen, Press & Teukolsky, ApJ 178, 347 (1972) [ISCO formula]
"""

import numpy as np


# ---------------------------------------------------------------------------
# Kerr ISCO geometry
# ---------------------------------------------------------------------------

def r_isco_kerr(a):
    """ISCO radius for a Kerr black hole with dimensionless spin *a*.

    Parameters
    ----------
    a : float or array_like
        Dimensionless spin, clipped to (-0.9999, 0.9999) to avoid the
        coordinate singularity at |a| = 1.  Positive = prograde orbit.

    Returns
    -------
    float or ndarray
        ISCO radius in units of M (gravitational radii, M=1 convention).

    Notes
    -----
    From Bardeen, Press & Teukolsky (1972), valid for prograde (a > 0)
    and retrograde (a < 0) orbits.

    Limiting cases:
        a =  0    ->  r_ISCO = 6 M   (Schwarzschild)
        a -> +1   ->  r_ISCO -> 1 M  (extremal prograde)
        a -> -1   ->  r_ISCO -> 9 M  (extremal retrograde)
    """
    a = np.asarray(a, dtype=float)
    a = np.clip(a, -0.9999, 0.9999)
    a2 = a ** 2
    abs_a = np.abs(a)
    Z1 = (1
          + (1 - a2) ** (1 / 3)
          * ((1 + abs_a) ** (1 / 3) + (1 - abs_a) ** (1 / 3)))
    Z2 = np.sqrt(3 * a2 + Z1 ** 2)
    return 3 + Z2 - np.sign(a) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))


def e_isco_kerr(a):
    """Fractional binding energy at the Kerr ISCO: 1 - E_ISCO/mu.

    Parameters
    ----------
    a : float or array_like
        Dimensionless spin, clipped to (-0.9999, 0.9999).

    Returns
    -------
    float or ndarray
        Binding energy fraction.  Equals 1 - sqrt(8/9) ≈ 0.05719 at a = 0.
        Increases for prograde spin (a > 0), decreases for retrograde (a < 0).
    """
    a = np.asarray(a, dtype=float)
    r = r_isco_kerr(a)
    sqrt_r = np.sqrt(r)
    num = r ** 2 - 2 * r + a * sqrt_r
    denom = r * np.sqrt(np.maximum(r ** 2 - 3 * r + 2 * a * sqrt_r, 0.0))
    e_isco = np.where(denom > 0, num / denom, 1.0)   # E_ISCO / mu
    return 1.0 - e_isco                               # binding energy fraction


# ---------------------------------------------------------------------------
# GR baseline predictions
# ---------------------------------------------------------------------------

def f_gr_nonspin(eta):
    """Non-spinning GR fractional mass deficit (Buonanno, Kidder & Lehner 2008).

    f = (1 - sqrt(8/9)) * eta + 0.498 * eta^2

    The linear coefficient is the exact Schwarzschild ISCO binding energy
    1 - sqrt(8/9) ≈ 0.05719, which BKL round to 0.0572 in their paper.
    Using the exact value keeps f_gr_spin(eta, 0) == f_gr_nonspin(eta)
    to numerical precision.

    Parameters
    ----------
    eta : float or array_like
        Symmetric mass ratio m1*m2/(m1+m2)^2, in (0, 0.25].
    """
    eta = np.asarray(eta, dtype=float)
    return e_isco_kerr(0.0) * eta + 0.498 * eta ** 2


def f_gr_spin(eta, chi_eff):
    """Spinning extension of the BKL (2008) fractional mass deficit.

    Replaces the Schwarzschild ISCO binding energy (0.0572 per unit eta)
    with the Kerr ISCO binding energy evaluated at chi_eff, while keeping
    the quadratic finite-mass-ratio correction unchanged:

        f = e_isco_kerr(chi_eff) * eta + 0.498 * eta^2

    At chi_eff = 0 this reduces exactly to f_gr_nonspin(eta).
    For chi_eff > 0 (prograde / aligned spins) the baseline is higher,
    meaning standard GR already predicts more deficit.
    For chi_eff < 0 (retrograde / anti-aligned) the baseline is lower.

    Parameters
    ----------
    eta : float or array_like
        Symmetric mass ratio.
    chi_eff : float or array_like
        Effective spin parameter (mass-weighted projection onto orbital axis).
    """
    eta = np.asarray(eta, dtype=float)
    chi_eff = np.asarray(chi_eff, dtype=float)
    return e_isco_kerr(chi_eff) * eta + 0.498 * eta ** 2
