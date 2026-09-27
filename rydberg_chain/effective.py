"""Effective transverse-field Ising model for one excitation per strong pair.

Pair k holds atoms L = 2k and R = 2k + 1 at 2k and 2k + s, with s = 1 - eps.
The pseudo-spin is tau^z = +1 when L is excited and -1 when R is excited.
Two excitations in pairs k and k + m have energy

    U_m = A_m + B_m (tau_k - tau_{k+m}) + J_m tau_k tau_{k+m},
    J_m = (u0 - (u_plus + u_minus)/2)/2,  B_m = (u_plus - u_minus)/4,

with u0 = g/(2m)^6 and u_plus/u_minus = g/(2m + s)^6, g/(2m - s)^6.  Second-order
perturbation theory in Omega through |00> and |11> gives the transverse field

    h = (Omega^2/4) (1/Delta + 1/(V_s - Delta)),  V_s = g/s^6.

See docs/DIMERIZED_THEORY.md for the derivation.
"""

import numpy as np
from scipy.optimize import brentq


def tfim_couplings(eps, g, delta, omega=1.0, mmax=8):
    """Leading-order couplings of the effective pseudo-spin chain.

    Returns a dict with keys
      'm'        pair separations 1..mmax,
      'J_m'      Ising couplings (negative = ferromagnetic between pseudo-spins),
      'B_m'      field-like terms (telescope to boundary fields in the bulk),
      'u0', 'u_plus', 'u_minus'  the underlying pair energies,
      'h'        transverse field (nan unless 0 < delta < V_s, where the
                 perturbative expression is not defined),
      'J1'       nearest-pair coupling J_m[0],
      's', 'V_s', 'V_w'  strong gap and nearest-neighbour energies.
    Lists are Python floats so the dict can be written to JSON directly.
    """
    s = 1.0 - eps
    m = np.arange(1, mmax + 1, dtype=float)
    u0 = g / (2 * m) ** 6
    u_plus = g / (2 * m + s) ** 6
    u_minus = g / (2 * m - s) ** 6
    J = 0.5 * (u0 - 0.5 * (u_plus + u_minus))
    B = 0.25 * (u_plus - u_minus)
    V_s, V_w = g / s ** 6, g / (1.0 + eps) ** 6
    h = transverse_field(eps, g, delta, omega)
    return {"m": [int(v) for v in m], "J_m": J.tolist(), "B_m": B.tolist(),
            "u0": u0.tolist(), "u_plus": u_plus.tolist(), "u_minus": u_minus.tolist(),
            "h": h, "J1": float(J[0]), "s": s, "V_s": V_s, "V_w": V_w}


def transverse_field(eps, g, delta, omega=1.0):
    """h = (omega^2/4)(1/delta + 1/(V_s - delta)); nan outside 0 < delta < V_s."""
    V_s = g / (1.0 - eps) ** 6
    if not 0.0 < delta < V_s:
        return float("nan")
    return float(0.25 * omega ** 2 * (1.0 / delta + 1.0 / (V_s - delta)))


def tfim_order_parameter(J, h):
    """Pfeuty's <tau^z> = (1 - (h/J)^2)^(1/8) of the nearest-neighbour TFIM; 0 for |h| >= |J|."""
    J, h = np.abs(np.asarray(J, dtype=float)), np.abs(np.asarray(h, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        x = 1.0 - (h / J) ** 2
    out = np.where(h < J, np.maximum(x, 0.0) ** 0.125, 0.0)
    out = np.where(np.isnan(h) | np.isnan(J), np.nan, out)
    return float(out) if out.ndim == 0 else out


def critical_eps(g, delta, omega=1.0, eps_min=0.0, eps_max=0.95):
    """Dimerization where |J_1| = h (nearest-neighbour TFIM estimate of Z2 melting).

    Returns nan when |J_1| - h does not change sign on [eps_min, eps_max].
    """
    def f(eps):
        return abs(tfim_couplings(eps, g, delta, omega, mmax=1)["J1"]) - transverse_field(eps, g, delta, omega)

    a, b = f(eps_min), f(eps_max)
    if not (np.isfinite(a) and np.isfinite(b)) or a * b > 0:
        return float("nan")
    return float(brentq(f, eps_min, eps_max, xtol=1e-12))
