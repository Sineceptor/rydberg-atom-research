"""Finite DMRG with TeNPy, either one atom per site or one strong pair per site.

The pair representation groups atoms (2k, 2k+1) into a 4-state site with basis
|n_L n_R> = 00, 10, 01, 11.  The strong-bond energy V_s n_L n_R becomes on-site
and only inter-pair couplings enter the MPO.  It is the same Hamiltonian, not an
approximation (unless mmax truncates the inter-pair range).

Both ground-state functions return (energy, psi, obs); entropies in obs are in bits.
"""

import numpy as np
import tenpy.linalg.np_conserved as npc
from tenpy.algorithms import dmrg as tenpy_dmrg
from tenpy.models.lattice import Chain
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPS
from tenpy.networks.site import BosonSite, Site

from .geometry import coupling_matrix

LN2 = np.log(2.0)


def dmrg_options(chi=64, max_sweeps=60, min_sweeps=6, max_E_err=1e-12, max_S_err=1e-8,
                 svd_min=1e-12, mixer=True, mixer_sweeps=10):
    """Options dict for tenpy.algorithms.dmrg.run (two-site finite DMRG).

    Convergence requires |dE| < max_E_err and |dS| < max_S_err between sweeps
    (S = mean bond entropy).  States converge more slowly than energies (error
    ~ sqrt(dE/gap)); for observables at the 1e-7 level in small or gapped chains
    pass max_S_err=1e-10.  In long even chains the soft domain-wall mode makes
    such tight settings run to max_sweeps.  The density-matrix mixer is switched
    off after mixer_sweeps sweeps, or earlier once the criteria are met.
    """
    options = {"mixer": mixer, "max_E_err": max_E_err, "max_S_err": max_S_err,
               "trunc_params": {"chi_max": chi, "svd_min": svd_min},
               "max_sweeps": max_sweeps, "min_sweeps": min_sweeps, "combine": True}
    if mixer:
        options["mixer_params"] = {"amplitude": 1e-5, "decay": 2.0, "disable_after": mixer_sweeps}
    return options


def default_product_state(n_atoms):
    """Alternating occupations with both end atoms excited.

    Odd n: 1010...101.  Even n: 1010... on the left half and ...0101 on the right
    half, i.e. one 00 domain wall near the centre, which the boundary fields of
    the ordered phase favour.  Returned as a list of 0/1 per atom.
    """
    left = (n_atoms + 1) // 2 if n_atoms % 2 else 2 * ((n_atoms // 2 + 1) // 2)
    occ = [1 - i % 2 for i in range(left)]
    occ += [1 - (n_atoms - 1 - i) % 2 for i in range(left, n_atoms)]
    return occ


def _run(psi, model, options):
    """Run DMRG in place; return the energy and convergence diagnostics."""
    info = tenpy_dmrg.run(psi, model, options)
    stats = info["sweep_statistics"]
    dE = abs(stats["Delta_E"][-1]) if len(stats["E"]) > 1 else float("nan")
    dS = abs(stats["Delta_S"][-1]) if len(stats["E"]) > 1 else float("nan")
    return float(info["E"]), {
        "n_sweeps": len(stats["E"]),
        "last_sweep_dE": float(dE),
        "last_sweep_dS": float(dS),
        "converged": bool(dE < options["max_E_err"] and dS < options["max_S_err"]),
        "max_trunc_err": float(max(stats["max_trunc_err"])),
        "final_trunc_err": float(stats["max_trunc_err"][-1]),
        "max_chi": int(max(psi.chi)) if len(psi.chi) else 1,
    }


# ----------------------------------------------------------------- atom chain

def atom_chain_model(x, g, delta, omega=1.0, rmax=None, local_delta=None):
    """TeNPy MPOModel with one BosonSite(Nmax=1) per atom at positions x."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    V = coupling_matrix(x, g, rmax)
    det = np.full(n, float(delta))
    if local_delta is not None:
        det = det + np.asarray(local_delta, dtype=float).reshape(n)
    site = BosonSite(Nmax=1, conserve=None)
    lat = Chain(n, site, bc="open", bc_MPS="finite")
    model = CouplingModel(lat)
    for i in range(n):
        model.add_onsite_term(0.5 * omega, i, "B")
        model.add_onsite_term(0.5 * omega, i, "Bd")
        model.add_onsite_term(-det[i], i, "N")
    for i in range(n):
        for j in range(i + 1, n):
            if V[i, j] != 0.0:
                model.add_coupling_term(V[i, j], i, j, "N", "N")
    return MPOModel(lat, model.calc_H_MPO())


def atom_chain_ground_state(x, g, delta, omega=1.0, chi=64, rmax=None, local_delta=None,
                            init=None, correlations=False, variance=False, **dmrg_kwargs):
    """Finite DMRG ground state of atoms at positions x (one site per atom).

    init: list of 0/1 occupations (default: default_product_state).
    dmrg_kwargs are passed to dmrg_options (max_sweeps, max_E_err, mixer, ...).
    obs keys: 'occupations', 'density', 'pair_order' (n_L - n_R for atoms
    (2k, 2k+1)), 'half_chain_entropy' (bits, left block of 'entropy_cut' = n//2
    atoms), 'bond_entropies' (bits, cuts 1..n-1), and the DMRG diagnostics
    'n_sweeps', 'converged', 'last_sweep_dE', 'last_sweep_dS', 'max_trunc_err'
    (all sweeps, including the mixer phase), 'final_trunc_err' (last sweep),
    'max_chi'.  correlations=True adds 'nn' (<n_i n_j>) and
    'staggered_order_squared'; variance=True adds 'variance' (<H^2> - <H>^2).
    """
    n = len(x)
    M = atom_chain_model(x, g, delta, omega, rmax, local_delta)
    init = default_product_state(n) if init is None else init
    psi = MPS.from_product_state(M.lat.mps_sites(), [str(int(v)) for v in init], bc="finite",
                              unit_cell_width=M.lat.mps_unit_cell_width)
    energy, obs = _run(psi, M, dmrg_options(chi, **dmrg_kwargs))
    occ = np.real(psi.expectation_value("N"))
    S = np.asarray(psi.entanglement_entropy()) / LN2
    m = 2 * (n // 2)
    obs.update(occupations=occ, density=float(occ.mean()), pair_order=occ[0:m:2] - occ[1:m:2],
               bond_entropies=S, entropy_cut=n // 2,
               half_chain_entropy=float(S[n // 2 - 1]) if n > 1 else 0.0)
    if correlations:
        nn = np.real(psi.correlation_function("N", "N"))
        obs.update(nn=nn, staggered_order_squared=_staggered_from_nn(nn, occ))
    if variance:
        obs["variance"] = float(np.real(M.H_MPO.variance(psi)))
    return energy, psi, obs


# ----------------------------------------------------------------- pair chain

def pair_site():
    """One strong pair as a 4-state TeNPy Site with basis 00, 10, 01, 11 (n_L n_R).

    Operators: nL, nR, XL, XR (flip one atom), nLnR, tau = nL - nR and the
    projectors P00, P10, P01, P11.
    """
    nL = np.diag([0.0, 1.0, 0.0, 1.0])
    nR = np.diag([0.0, 0.0, 1.0, 1.0])
    XL = np.zeros((4, 4))
    XR = np.zeros((4, 4))
    for a, b in [(0, 1), (2, 3)]:          # 00<->10 and 01<->11 flip L
        XL[a, b] = XL[b, a] = 1.0
    for a, b in [(0, 2), (1, 3)]:          # 00<->01 and 10<->11 flip R
        XR[a, b] = XR[b, a] = 1.0
    proj = {f"P{label}": np.diag(np.eye(4)[k]) for k, label in enumerate(["00", "10", "01", "11"])}
    leg = npc.LegCharge.from_trivial(4)
    return Site(leg, ["00", "10", "01", "11"], nL=nL, nR=nR, XL=XL, XR=XR,
                nLnR=nL @ nR, tau=nL - nR, **proj)


def pair_chain_model(eps, g, delta, n_pairs, omega=1.0, mmax=None, local_delta=None):
    """TeNPy MPOModel of the dimerized chain of 2*n_pairs atoms, one site per strong pair.

    Pair a holds atoms at 2a (L, offset 0) and 2a + 1 - eps (R, offset 1 - eps).
    Inter-pair couplings g/d^6 with d = 2m + off_b - off_a are kept for pair
    separations m <= mmax (default: all).  local_delta (length 2*n_pairs, atom
    order) is added to delta.
    """
    s = 1.0 - eps
    V_s = g / s ** 6
    mmax = n_pairs - 1 if mmax is None else min(mmax, n_pairs - 1)
    det = np.full(2 * n_pairs, float(delta))
    if local_delta is not None:
        det = det + np.asarray(local_delta, dtype=float).reshape(2 * n_pairs)
    lat = Chain(n_pairs, pair_site(), bc="open", bc_MPS="finite")
    model = CouplingModel(lat)
    offsets = [("nL", 0.0), ("nR", s)]
    for a in range(n_pairs):
        model.add_onsite_term(0.5 * omega, a, "XL")
        model.add_onsite_term(0.5 * omega, a, "XR")
        model.add_onsite_term(-det[2 * a], a, "nL")
        model.add_onsite_term(-det[2 * a + 1], a, "nR")
        model.add_onsite_term(V_s, a, "nLnR")
        for m in range(1, mmax + 1):
            b = a + m
            if b >= n_pairs:
                break
            for op_a, off_a in offsets:
                for op_b, off_b in offsets:
                    model.add_coupling_term(g / (2 * m + off_b - off_a) ** 6, a, b, op_a, op_b)
    return MPOModel(lat, model.calc_H_MPO())


def pair_chain_ground_state(eps, g, delta, n_pairs, chi=64, mmax=None, omega=1.0,
                            local_delta=None, init=None, correlations=False, variance=False,
                            **dmrg_kwargs):
    """Finite DMRG ground state of 2*n_pairs atoms using 4-state pair sites.

    init: list of pair labels ('00', '10', '01', '11'); default is
    default_product_state regrouped into pairs.
    obs keys: 'occupations' (atom order), 'density', 'pair_order' (<n_L - n_R>),
    'pair_probabilities' (n_pairs x 4, order 00, 10, 01, 11), 'bond_entropies'
    (bits, at atom cuts 2, 4, ..., 2*n_pairs - 2), 'half_chain_entropy' (bits,
    at atom cut 'entropy_cut' = 2*(n_pairs//2)), plus the DMRG diagnostics of
    atom_chain_ground_state; correlations=True adds atom-order 'nn' and
    'staggered_order_squared'; variance=True adds 'variance'.
    """
    M = pair_chain_model(eps, g, delta, n_pairs, omega, mmax, local_delta)
    if init is None:
        occ0 = default_product_state(2 * n_pairs)
        init = [f"{occ0[2 * a]}{occ0[2 * a + 1]}" for a in range(n_pairs)]
    psi = MPS.from_product_state(M.lat.mps_sites(), list(init), bc="finite",
                              unit_cell_width=M.lat.mps_unit_cell_width)
    energy, obs = _run(psi, M, dmrg_options(chi, **dmrg_kwargs))
    nL = np.real(psi.expectation_value("nL"))
    nR = np.real(psi.expectation_value("nR"))
    occ = np.empty(2 * n_pairs)
    occ[0::2], occ[1::2] = nL, nR
    probs = np.column_stack([np.real(psi.expectation_value(op)) for op in ("P00", "P10", "P01", "P11")])
    S = np.asarray(psi.entanglement_entropy()) / LN2
    half = n_pairs // 2
    obs.update(occupations=occ, density=float(occ.mean()), pair_order=nL - nR,
               pair_probabilities=probs, bond_entropies=S, entropy_cut=2 * half,
               half_chain_entropy=float(S[half - 1]) if n_pairs > 1 else 0.0)
    if correlations:
        nn = np.empty((2 * n_pairs, 2 * n_pairs))
        for alpha, op1 in enumerate(("nL", "nR")):
            for beta, op2 in enumerate(("nL", "nR")):
                nn[alpha::2, beta::2] = np.real(psi.correlation_function(op1, op2))
        obs.update(nn=nn, staggered_order_squared=_staggered_from_nn(nn, occ))
    if variance:
        obs["variance"] = float(np.real(M.H_MPO.variance(psi)))
    return energy, psi, obs


def _staggered_from_nn(nn, occ):
    """<Q^2>, Q = (1/n) sum_i (-1)^i (2 n_i - 1), from <n_i n_j> and <n_i>."""
    n = len(occ)
    sgn = (-1.0) ** np.arange(n)
    # <s_i s_j> with s = 2n - 1: 4<n_i n_j> - 2<n_i> - 2<n_j> + 1
    ss = 4 * nn - 2 * occ[:, None] - 2 * occ[None, :] + 1
    return float(sgn @ ss @ sgn / n ** 2)
