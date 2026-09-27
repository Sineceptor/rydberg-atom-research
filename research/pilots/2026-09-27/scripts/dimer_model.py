"""Pilot: dimerized (alternating-spacing) Rydberg chain for ED and DMRG.

H = (Omega/2) sum_i X_i - Delta sum_i n_i + sum_{i<j} C6/|x_i - x_j|^6 n_i n_j
Units: Omega = 1, mean spacing a = 1, g = C6/(Omega a^6) = (R_b/a)^6.
Bond (2k, 2k+1) has length 1 - eps (strong), bond (2k+1, 2k+2) has 1 + eps (weak).
"""
import sys
import numpy as np

sys.path.insert(0, "/home/user/rydberg-atom-research/study")


def positions(number, eps):
    gaps = [(1 - eps) if i % 2 == 0 else (1 + eps) for i in range(number - 1)]
    return np.concatenate([[0.0], np.cumsum(gaps)])


def ed_ground(number, eps, g, delta, k=2):
    from rydberg_model import interactions, hamiltonian
    from scipy.sparse.linalg import eigsh
    x = positions(number, eps)
    coupling = interactions(np.column_stack([x, np.zeros(number)]), g)
    H = hamiltonian(coupling, omega=1.0, delta=delta)
    v0 = np.random.default_rng(1).normal(size=H.shape[0])
    energies, vectors = eigsh(H, k=k, which="SA", tol=1e-12, v0=v0)
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def build_finite(number, eps, g, delta, rmax=None):
    from tenpy.networks.site import BosonSite
    from tenpy.models.lattice import Chain
    from tenpy.models.model import CouplingModel, MPOModel

    site = BosonSite(Nmax=1, conserve=None)
    lat = Chain(number, site, bc="open", bc_MPS="finite")
    model = CouplingModel(lat)
    x = positions(number, eps)
    rmax = number - 1 if rmax is None else rmax
    for i in range(number):
        model.add_onsite_term(-delta, i, "N")
        model.add_onsite_term(0.5, i, "B")
        model.add_onsite_term(0.5, i, "Bd")
    for i in range(number):
        for j in range(i + 1, min(number, i + rmax + 1)):
            model.add_coupling_term(g / abs(x[j] - x[i]) ** 6, i, j, "N", "N")
    return lat, MPOModel(lat, model.calc_H_MPO())


def dmrg_ground(number, eps, g, delta, chi=64, rmax=10, excited=False, seed_state=None):
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
    lat, M = build_finite(number, eps, g, delta, rmax)
    # Start from a Z2 product state; DMRG mixer lets it relax to any state.
    init = seed_state or ["1" if i % 2 == 0 else "0" for i in range(number)]
    psi = MPS.from_product_state(lat.mps_sites(), init, bc="finite")
    params = {"mixer": True, "max_E_err": 1e-11, "max_S_err": 1e-8,
              "trunc_params": {"chi_max": chi, "svd_min": 1e-12},
              "max_sweeps": 60, "min_sweeps": 6}
    info = dmrg.run(psi, M, params)
    return info["E"], psi, M


def x_of(j, eps):
    """Position of site j in the infinite dimerized chain (unit cell length 2)."""
    return 2.0 * (j // 2) + (1 - eps) * (j % 2)


def build_infinite(eps, g, delta, rmax=12, cell=2):
    from tenpy.networks.site import BosonSite
    from tenpy.models.lattice import Chain
    from tenpy.models.model import CouplingModel, MPOModel
    site = BosonSite(Nmax=1, conserve=None)
    lat = Chain(cell, site, bc="periodic", bc_MPS="infinite")
    model = CouplingModel(lat)
    for i in range(cell):
        model.add_onsite_term(-delta, i, "N")
        model.add_onsite_term(0.5, i, "B")
        model.add_onsite_term(0.5, i, "Bd")
        for j in range(i + 1, i + rmax + 1):
            model.add_coupling_term(g / (x_of(j, eps) - x_of(i, eps)) ** 6, i, j, "N", "N")
    return lat, MPOModel(lat, model.calc_H_MPO())


def idmrg_point(eps, g, delta, chi=64, rmax=12, cell=2, init=None):
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
    lat, M = build_infinite(eps, g, delta, rmax, cell)
    init = init or (["1", "0"] * (cell // 2))
    psi = MPS.from_product_state(lat.mps_sites(), init, bc="infinite")
    params = {"mixer": True, "max_E_err": 1e-10, "max_S_err": 1e-7,
              "trunc_params": {"chi_max": chi, "svd_min": 1e-12},
              "max_sweeps": 80, "min_sweeps": 10}
    info = dmrg.run(psi, M, params)
    n = psi.expectation_value("N")
    return {"eps": eps, "g": g, "delta": delta, "chi": chi, "cell": cell,
            "E_per_site": float(info["E"]), "n": [float(v) for v in n],
            "order": float(abs(n[0] - n[1])), "density": float(np.mean(n)),
            "xi": float(psi.correlation_length()), "S_bond": [float(s) for s in psi.entanglement_entropy()]}


def tfim_estimate(eps, g, delta):
    """Leading-order effective Ising couplings in the one-excitation-per-dimer sector."""
    Vs, Vw = g / (1 - eps) ** 6, g / (1 + eps) ** 6
    J = 0.0
    for m in range(1, 8):
        u0, up, um = g / (2 * m) ** 6, g / (2 * m + 1 - eps) ** 6, g / (2 * m - 1 + eps) ** 6
        J += ((up + um) / 2 - u0) / 2
    h = 0.25 / delta + 0.25 / (Vs + Vw - delta)
    return J, h


def dimer_site():
    """One strong pair as a single 4-state site: |nL nR> in order 00, 10, 01, 11."""
    import tenpy.linalg.np_conserved as npc
    from tenpy.networks.site import Site
    nL = np.diag([0., 1., 0., 1.])
    nR = np.diag([0., 0., 1., 1.])
    XL = np.zeros((4, 4)); XR = np.zeros((4, 4))
    for a, b in [(0, 1), (2, 3)]:
        XL[a, b] = XL[b, a] = 1.0     # flips the left atom
    for a, b in [(0, 2), (1, 3)]:
        XR[a, b] = XR[b, a] = 1.0     # flips the right atom
    leg = npc.LegCharge.from_trivial(4)
    return Site(leg, ["00", "10", "01", "11"], nL=nL, nR=nR, XL=XL, XR=XR,
                tau=nL - nR, nLnR=nL @ nR)


def _dimer_terms(model, cells, eps, g, delta, mmax, finite_size=None):
    """Exactly the atom Hamiltonian, regrouped: strong bond inside a site, rest between sites."""
    Vs = g / (1 - eps) ** 6
    for a in range(cells):
        model.add_onsite_term(0.5, a, "XL")
        model.add_onsite_term(0.5, a, "XR")
        model.add_onsite_term(-delta, a, "nL")
        model.add_onsite_term(-delta, a, "nR")
        model.add_onsite_term(Vs, a, "nLnR")
        for m in range(1, mmax + 1):
            b = a + m
            if finite_size is not None and b >= finite_size:
                break
            for op_a, off_a in [("nL", 0.0), ("nR", 1 - eps)]:
                for op_b, off_b in [("nL", 0.0), ("nR", 1 - eps)]:
                    d = 2 * m + off_b - off_a
                    model.add_coupling_term(g / d ** 6, a, b, op_a, op_b)


def build_dimer(eps, g, delta, cells, infinite, mmax=6):
    from tenpy.models.lattice import Chain
    from tenpy.models.model import CouplingModel, MPOModel
    site = dimer_site()
    if infinite:
        lat = Chain(cells, site, bc="periodic", bc_MPS="infinite")
        model = CouplingModel(lat)
        _dimer_terms(model, cells, eps, g, delta, mmax)
    else:
        lat = Chain(cells, site, bc="open", bc_MPS="finite")
        model = CouplingModel(lat)
        _dimer_terms(model, cells, eps, g, delta, mmax, finite_size=cells)
    return lat, MPOModel(lat, model.calc_H_MPO())


def dimer_dmrg(eps, g, delta, cells, infinite, chi=64, mmax=6, init=None, mixer=True):
    from tenpy.networks.mps import MPS
    from tenpy.algorithms import dmrg
    lat, M = build_dimer(eps, g, delta, cells, infinite, mmax)
    init = init or ["10"] * cells
    psi = MPS.from_product_state(lat.mps_sites(), init, bc="infinite" if infinite else "finite")
    params = {"mixer": mixer, "max_E_err": 1e-11, "max_S_err": 1e-8,
              "trunc_params": {"chi_max": chi, "svd_min": 1e-12},
              "max_sweeps": 100, "min_sweeps": 10}
    info = dmrg.run(psi, M, params)
    return info["E"], psi, M
