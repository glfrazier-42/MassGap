"""Unit tests for src/bhem/gr_baselines.py.

Run with:
    PYTHONPATH=src venv/Scripts/python.exe -m pytest tests/test_gr_baselines.py -v
or standalone (no pytest required):
    PYTHONPATH=src venv/Scripts/python.exe tests/test_gr_baselines.py
"""

import math
import sys
import numpy as np

from bhem.gr_baselines import (
    r_isco_kerr,
    e_isco_kerr,
    f_gr_nonspin,
    f_gr_spin,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def assert_close(actual, expected, rtol=1e-4, label=""):
    diff = abs(actual - expected)
    scale = max(abs(expected), 1e-12)
    if diff / scale > rtol:
        raise AssertionError(
            "%s: expected %.8f, got %.8f (rel err %.2e, tol %.2e)"
            % (label, expected, actual, diff / scale, rtol)
        )


# ---------------------------------------------------------------------------
# r_isco_kerr
# ---------------------------------------------------------------------------

def test_isco_radius_schwarzschild():
    """Schwarzschild ISCO is exactly 6 M."""
    assert_close(r_isco_kerr(0.0), 6.0, rtol=1e-6,
                 label="r_ISCO(a=0)")


def test_isco_radius_extremal_prograde():
    """Extremal prograde Kerr ISCO approaches 1 M."""
    r = r_isco_kerr(0.9999)
    assert 1.0 < r < 1.1, \
        "r_ISCO(a->1) should be just above 1 M, got %.4f" % r


def test_isco_radius_extremal_retrograde():
    """Extremal retrograde Kerr ISCO approaches 9 M."""
    r = r_isco_kerr(-0.9999)
    assert 8.9 < r < 9.1, \
        "r_ISCO(a->-1) should be near 9 M, got %.4f" % r


def test_isco_radius_prograde_decreases():
    """ISCO radius decreases monotonically with prograde spin."""
    spins = np.linspace(0.0, 0.99, 50)
    r_vals = r_isco_kerr(spins)
    assert np.all(np.diff(r_vals) < 0), \
        "r_ISCO should decrease monotonically for a in [0, 1)"


def test_isco_radius_retrograde_increases():
    """ISCO radius increases monotonically with retrograde spin magnitude."""
    spins = np.linspace(0.0, -0.99, 50)
    r_vals = r_isco_kerr(spins)
    assert np.all(np.diff(r_vals) > 0), \
        "r_ISCO should increase monotonically for a in [0, -1)"


# ---------------------------------------------------------------------------
# e_isco_kerr
# ---------------------------------------------------------------------------

def test_isco_energy_schwarzschild():
    """Schwarzschild ISCO binding energy = 1 - sqrt(8/9) ≈ 0.05719."""
    expected = 1.0 - math.sqrt(8.0 / 9.0)
    assert_close(e_isco_kerr(0.0), expected, rtol=1e-4,
                 label="e_ISCO(a=0)")


def test_isco_energy_prograde_increases():
    """Binding energy increases monotonically with prograde spin."""
    spins = np.linspace(0.0, 0.99, 50)
    e_vals = e_isco_kerr(spins)
    assert np.all(np.diff(e_vals) > 0), \
        "Binding energy should increase for prograde spin"


def test_isco_energy_retrograde_decreases():
    """Retrograde spin reduces binding energy below Schwarzschild."""
    assert e_isco_kerr(-0.5) < e_isco_kerr(0.0)
    assert e_isco_kerr(-0.9) < e_isco_kerr(-0.5)


# ---------------------------------------------------------------------------
# f_gr_nonspin
# ---------------------------------------------------------------------------

def test_nonspin_formula_coefficients():
    """f_gr_nonspin matches the exact BKL formula using 1-sqrt(8/9)."""
    import math
    linear = 1.0 - math.sqrt(8.0 / 9.0)   # exact Schwarzschild ISCO binding energy
    for eta in [0.05, 0.10, 0.15, 0.20, 0.25]:
        expected = linear * eta + 0.498 * eta ** 2
        assert_close(f_gr_nonspin(eta), expected, rtol=1e-6,
                     label="f_gr_nonspin(eta=%.2f)" % eta)


def test_nonspin_zero_eta():
    """f_gr_nonspin(0) = 0."""
    assert f_gr_nonspin(0.0) == 0.0


def test_nonspin_equal_mass():
    """For eta=0.25 (equal mass), f_GR ≈ 4.54%."""
    f = f_gr_nonspin(0.25)
    assert 0.040 < f < 0.055, \
        "f_gr_nonspin(eta=0.25) = %.4f, expected ~0.045" % f


# ---------------------------------------------------------------------------
# KEY TEST: spinning baseline reduces to non-spinning at chi_eff = 0
# ---------------------------------------------------------------------------

def test_spin_zero_equals_nonspin():
    """f_gr_spin(eta, chi_eff=0) must equal f_gr_nonspin(eta) for all eta.

    This is the primary consistency check: at zero spin the two baselines
    must be numerically identical.
    """
    etas = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
    for eta in etas:
        ns = f_gr_nonspin(eta)
        sp = f_gr_spin(eta, 0.0)
        assert_close(sp, ns, rtol=1e-4,
                     label="f_gr_spin(eta=%.2f, chi_eff=0)" % eta)


# ---------------------------------------------------------------------------
# Monotonicity and direction of spin correction
# ---------------------------------------------------------------------------

def test_spin_positive_exceeds_nonspin():
    """Prograde spin (chi_eff > 0) predicts more energy radiated."""
    eta = 0.25
    ns = f_gr_nonspin(eta)
    for chi in [0.1, 0.2, 0.3, 0.5, 0.7]:
        sp = f_gr_spin(eta, chi)
        assert sp > ns, (
            "f_gr_spin(eta=0.25, chi_eff=%.1f)=%.6f should exceed "
            "f_gr_nonspin=%.6f" % (chi, sp, ns)
        )


def test_spin_negative_below_nonspin():
    """Retrograde spin (chi_eff < 0) predicts less energy radiated."""
    eta = 0.25
    ns = f_gr_nonspin(eta)
    for chi in [-0.1, -0.2, -0.3, -0.5]:
        sp = f_gr_spin(eta, chi)
        assert sp < ns, (
            "f_gr_spin(eta=0.25, chi_eff=%.1f)=%.6f should be below "
            "f_gr_nonspin=%.6f" % (chi, sp, ns)
        )


def test_spin_monotone_in_chi_eff():
    """f_gr_spin increases monotonically with chi_eff."""
    eta = 0.25
    chis = np.linspace(-0.8, 0.8, 80)
    vals = np.array([f_gr_spin(eta, c) for c in chis])
    assert np.all(np.diff(vals) > 0), \
        "f_gr_spin should increase monotonically with chi_eff"


def test_spin_correction_magnitude():
    """Spin correction is small but non-negligible at typical GWTC spins.

    For chi_eff = 0.2 and eta = 0.25 the correction should be in the
    range 0.001 -- 0.005 (i.e. 0.1 -- 0.5 percentage points).
    """
    delta = f_gr_spin(0.25, 0.2) - f_gr_nonspin(0.25)
    assert 0.001 < delta < 0.010, \
        "Spin correction at chi_eff=0.2 is %.4f, expected 0.001--0.010" % delta


# ---------------------------------------------------------------------------
# Standalone runner (no pytest needed)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        test_isco_radius_schwarzschild,
        test_isco_radius_extremal_prograde,
        test_isco_radius_extremal_retrograde,
        test_isco_radius_prograde_decreases,
        test_isco_radius_retrograde_increases,
        test_isco_energy_schwarzschild,
        test_isco_energy_prograde_increases,
        test_isco_energy_retrograde_decreases,
        test_nonspin_formula_coefficients,
        test_nonspin_zero_eta,
        test_nonspin_equal_mass,
        test_spin_zero_equals_nonspin,
        test_spin_positive_exceeds_nonspin,
        test_spin_negative_below_nonspin,
        test_spin_monotone_in_chi_eff,
        test_spin_correction_magnitude,
    ]

    failed = 0
    for t in tests:
        try:
            t()
            print("  PASS  %s" % t.__name__)
        except AssertionError as e:
            print("  FAIL  %s: %s" % (t.__name__, e))
            failed += 1

    print("\n%d/%d tests passed." % (len(tests) - failed, len(tests)))
    sys.exit(failed)
