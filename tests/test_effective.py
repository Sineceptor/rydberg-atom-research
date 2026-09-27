"""Effective transverse-field Ising couplings of the one-excitation-per-pair sector."""

import numpy as np
import pytest

from rydberg_chain import effective, geometry


@pytest.mark.parametrize("eps,g", [(0.0, 8.0), (0.3, 11.39), (0.6, 11.39), (0.45, 3.0)])
def test_couplings_match_direct_pair_energies(eps, g):
    """Fit U(tau_k, tau_{k+m}) = A + B (tau_k - tau_{k+m}) + J tau_k tau_{k+m} to 1/r^6 energies."""
    mmax = 5
    out = effective.tfim_couplings(eps, g, 2.5, mmax=mmax)
    x = geometry.chain_positions(2 * (mmax + 1), eps)
    atom = {+1: 0, -1: 1}                        # tau = +1: left atom excited
    design = np.array([[1, tk - tl, tk * tl] for tk in (1, -1) for tl in (1, -1)], float)
    for m in range(1, mmax + 1):
        energies = [g / abs(x[2 * m + atom[tl]] - x[atom[tk]]) ** 6 for tk in (1, -1) for tl in (1, -1)]
        (A, B, J), *_ = np.linalg.lstsq(design, energies, rcond=None)
        assert np.allclose(design @ [A, B, J], energies, rtol=1e-12)     # exact decomposition
        assert out["J_m"][m - 1] == pytest.approx(J, rel=1e-10, abs=1e-15)
        assert out["B_m"][m - 1] == pytest.approx(B, rel=1e-10, abs=1e-15)
        assert J <= 0                            # convexity of 1/r^6: ferromagnetic pseudo-spins


def test_nearest_pair_closed_form():
    eps, g = 0.6, 11.39
    s = 1 - eps
    out = effective.tfim_couplings(eps, g, 2.5, mmax=3)
    V_s, V_w = geometry.strong_weak_couplings(eps, g)
    assert out["V_s"] == pytest.approx(V_s) and out["V_w"] == pytest.approx(V_w)
    assert out["u_minus"][0] == pytest.approx(V_w)
    assert abs(out["J1"]) == pytest.approx((V_w + g / (2 + s) ** 6 - 2 * g / 64) / 4)
    assert out["m"] == [1, 2, 3] and len(out["J_m"]) == len(out["B_m"]) == 3


@pytest.mark.parametrize("eps,g,delta", [(0.3, 11.39, 2.5), (0.6, 11.39, 2.5), (0.2, 8.0, 3.0)])
def test_transverse_field_matches_single_pair(eps, g, delta):
    """Half the splitting of the one-excitation doublet of one pair is h + O(Omega^4)."""
    omega = 0.02
    V_s = g / (1 - eps) ** 6
    w = omega / 2
    H = np.array([[0, w, w, 0], [w, -delta, 0, w], [w, 0, -delta, w], [0, w, w, V_s - 2 * delta]])
    E = np.linalg.eigvalsh(H)
    h = effective.transverse_field(eps, g, delta, omega)
    assert h == pytest.approx(effective.tfim_couplings(eps, g, delta, omega)["h"])
    assert h == pytest.approx(omega ** 2 / 4 * (1 / delta + 1 / (V_s - delta)))
    # For 0 < delta < V_s the two lowest levels are the |10>, |01> doublet.
    assert (E[1] - E[0]) / 2 == pytest.approx(h, rel=1e-3)


def test_transverse_field_outside_window_is_nan():
    assert np.isnan(effective.transverse_field(0.3, 11.39, 0.0))
    assert np.isnan(effective.transverse_field(0.3, 11.39, -1.0))
    assert np.isnan(effective.transverse_field(0.0, 2.0, 3.0))       # delta > V_s


def test_pfeuty_order_parameter():
    op = effective.tfim_order_parameter
    assert op(1.0, 0.0) == 1.0
    assert op(-2.0, 1.0) == pytest.approx(0.75 ** 0.125)
    assert op(1.0, 1.0) == 0.0 and op(1.0, 3.0) == 0.0 and op(0.0, 0.1) == 0.0
    assert np.isnan(op(1.0, float("nan")))
    values = op(np.array([1.0, 1.0, 1.0]), np.array([0.2, 0.99, 1.5]))
    assert values.shape == (3,) and values[0] > values[1] > values[2] == 0.0


def test_critical_eps():
    eps_c = [effective.critical_eps(11.39, d) for d in (2.0, 2.5, 3.0)]
    for d, e in zip((2.0, 2.5, 3.0), eps_c):
        J1 = effective.tfim_couplings(e, 11.39, d, mmax=1)["J1"]
        assert abs(J1) == pytest.approx(effective.transverse_field(e, 11.39, d), rel=1e-9)
    assert eps_c[0] < eps_c[1] < eps_c[2]          # h ~ 1/Delta: larger Delta melts later
    assert np.isnan(effective.critical_eps(200.0, 2.5, eps_max=0.3))
