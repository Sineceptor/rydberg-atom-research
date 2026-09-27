"""Infinite dimerized chain: iDMRG on pair sites, with a finite-DMRG bulk reference.

The chain is regrouped exactly into pair sites (``dmrg.pair_site``, basis 00, 10,
01, 11), so the strong-bond energy V_s is on-site and only inter-pair couplings
g/d^6 with pair separation m <= mmax enter the MPO.  The MPS unit cell holds
``cell_pairs`` pairs (default 2, four atoms).  Energies are per atom.

Why the pilot iDMRG failed (validation/IDMRG_NOTES.md has the experiments):

1. TeNPy's Lanczos stops when (RitzRes/gap)^2 < P_tol with the gap read off the
   Krylov Ritz values.  |11> states at energy ~V_s put the second Ritz value far
   up, so the test passes after two vectors; otherwise N_max = 20 is hit.  The
   two-site problem was left unconverged by 1e-3 (eps = 0.3) to 7 (eps = 0.4).
2. svd_min = 1e-12 keeps Schmidt values at the noise floor.  iDMRG inverts them
   to build the next two-site guess, amplifying noise by up to 1e12.
3. Even with 1 and 2 fixed, a spectral width of ~2 V_s (2.7e3 at eps = 0.55)
   makes Lanczos too slow near the Z2 transition, and the state drifts.
4. TeNPy's transfer-matrix energy (``MPO.expectation_value_TM``) is unreliable for
   these MPOs (errors up to 3e-4); the growth estimate, a term-by-term sum and the
   power-method evaluation agree with each other and with finite DMRG.

What this module does:

* ``idmrg_options``: energy-based Lanczos stop with more iterations, svd_min 1e-8,
  a short mixer phase, at most 120 sweeps.  The energy is <H> of the final state
  (``energy_per_atom``), checked against a term-by-term sum (``energy_from_terms``).
* For large V_s the local space is truncated variationally to the three lowest
  eigenstates of the single-pair Hamiltonian (``dressed_pair_site``), which drops
  only the dressed |11> level near V_s.  The on-site physics stays exact to all
  orders in Omega; the error is second order in the inter-pair coupling of the
  discarded level, about Omega^2 V_w^2 / (8 V_s^3) per atom, and the energy stays
  an upper bound.  Results obtained this way carry ``obs['truncated'] = True``.
* ``infinite_ground_state`` runs both Z2 product starts (optionally a symmetric
  one) and keeps the lowest valid energy; ``pinned_chain_ground_state`` and
  ``finite_bulk_energy`` give the finite-DMRG bulk reference (and fallback).

Example::

    e, psi, obs = infinite_ground_state(0.4, 11.39, 2.5)     # energy per atom
    obs['order'], obs['density'], obs['xi'], obs['entropy_weak'], obs['truncated']
"""

import numpy as np
import tenpy.linalg.np_conserved as npc
from tenpy.algorithms import dmrg as tenpy_dmrg
from tenpy.models.lattice import Chain
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPS
from tenpy.networks.site import Site

from .dmrg import dmrg_options, pair_site

LN2 = np.log(2.0)
PAIR_LABELS = ("00", "10", "01", "11")
PAIR_OPS = ("nL", "nR", "XL", "XR", "nLnR", "tau", "P00", "P10", "P01", "P11")


# ------------------------------------------------------------------ local basis

def pair_hamiltonian(eps, g, delta, omega=1.0):
    """4 x 4 Hamiltonian of one isolated strong pair in the basis 00, 10, 01, 11."""
    site = pair_site()
    op = {k: site.get_op(k).to_ndarray() for k in ("XL", "XR", "nL", "nR", "nLnR")}
    return (0.5 * omega * (op["XL"] + op["XR"]) - delta * (op["nL"] + op["nR"])
            + g / (1.0 - eps) ** 6 * op["nLnR"])


def dressed_pair_basis(eps, g, delta, omega=1.0):
    """Eigenvalues (4,) of the single-pair Hamiltonian and its three lowest eigenvectors W (4, 3)."""
    levels, vectors = np.linalg.eigh(pair_hamiltonian(eps, g, delta, omega))
    return levels, vectors[:, :3]


def dressed_pair_site(eps, g, delta, omega=1.0):
    """Pair site restricted to the three lowest single-pair eigenstates ('d0', 'd1', 'd2').

    Every operator of ``dmrg.pair_site`` is replaced by W^T O W, so an MPS on these
    sites is exactly the 4-state MPS obtained by applying W on every site, and all
    expectation values refer to the full Hilbert space.  ``site.pair_basis`` is W.
    """
    levels, W = dressed_pair_basis(eps, g, delta, omega)
    if not levels[3] - levels[2] > 10.0 * (abs(delta) + abs(omega)):
        raise ValueError("the highest pair level is not separated; truncation is not controlled")
    full = pair_site()
    ops = {name: W.T @ full.get_op(name).to_ndarray() @ W for name in PAIR_OPS}
    site = Site(npc.LegCharge.from_trivial(3), ["d0", "d1", "d2"], **ops)
    site.pair_basis = W
    return site


def full_pair_site():
    """``dmrg.pair_site`` with ``pair_basis`` set to the identity."""
    site = pair_site()
    site.pair_basis = np.eye(4)
    return site


def truncation_error_estimate(eps, g, delta, omega=1.0):
    """Energy error per atom of the dressed truncation, Omega^2 V_w^2 / (8 V_s^3).

    Second-order coupling of the discarded level (at ~V_s) through the
    inter-pair density terms, whose matrix elements are ~(Omega/2V_s) V_w.
    Exact diagonalization of 8 atoms at g = 11.39, delta = 2.5 gives 1.0e-8,
    1.4e-9, 1.8e-10, 1.9e-11 per atom at eps = 0.4, 0.45, 0.5, 0.55 against
    estimates of 2.0e-8, 2.8e-9, 3.2e-10, 3.2e-11.
    """
    V_s, V_w = g / (1.0 - eps) ** 6, g / (1.0 + eps) ** 6
    return float(omega ** 2 * V_w ** 2 / (8.0 * V_s ** 3))


def use_truncation(eps, g, delta, omega=1.0, truncate="auto", max_error=1e-8):
    """True/False, or "auto": truncate when the estimated error is below max_error
    and the dressed |11> level is well separated."""
    if truncate != "auto":
        return bool(truncate)
    levels, _ = dressed_pair_basis(eps, g, delta, omega)
    separated = levels[3] - levels[2] > 10.0 * (abs(delta) + abs(omega))
    return bool(separated and truncation_error_estimate(eps, g, delta, omega) < max_error)


# ------------------------------------------------------------------ models

def pair_model(eps, g, delta, n_pairs, omega=1.0, mmax=8, truncate=False, infinite=True,
               local_delta=None):
    """TeNPy MPOModel of the dimerized chain with one site per strong pair.

    infinite=True: unit cell of n_pairs pairs, couplings from every pair of the
    cell to the next mmax pairs (crossing the cell boundary).  infinite=False: an
    open chain of n_pairs pairs, the same Hamiltonian as
    ``dmrg.pair_chain_model(eps, g, delta, n_pairs, omega, mmax, local_delta)``.
    Pair a holds atoms at 2a (L) and 2a + 1 - eps (R); the coupling between pairs
    a and a + m is g/d^6 with d = 2m + off_b - off_a.  truncate=True uses
    ``dressed_pair_site`` (variational, see module docstring).
    """
    if mmax < 1:
        raise ValueError("mmax must be at least 1")
    if infinite and n_pairs < 2:
        raise ValueError("two-site iDMRG needs a unit cell of at least 2 pairs")
    if infinite and local_delta is not None:
        raise ValueError("local_delta is only defined for finite chains")
    s = 1.0 - eps
    site = dressed_pair_site(eps, g, delta, omega) if truncate else full_pair_site()
    if infinite:
        lat = Chain(n_pairs, site, bc="periodic", bc_MPS="infinite")
    else:
        lat = Chain(n_pairs, site, bc="open", bc_MPS="finite")
    det = np.full(2 * n_pairs, float(delta))
    if local_delta is not None:
        det = det + np.asarray(local_delta, dtype=float).reshape(2 * n_pairs)
    model = CouplingModel(lat)
    offsets = [("nL", 0.0), ("nR", s)]
    for a in range(n_pairs):
        model.add_onsite_term(0.5 * omega, a, "XL")
        model.add_onsite_term(0.5 * omega, a, "XR")
        model.add_onsite_term(-det[2 * a], a, "nL")
        model.add_onsite_term(-det[2 * a + 1], a, "nR")
        model.add_onsite_term(g / s ** 6, a, "nLnR")
        for m in range(1, mmax + 1):
            if not infinite and a + m >= n_pairs:
                break
            for op_a, off_a in offsets:
                for op_b, off_b in offsets:
                    model.add_coupling_term(g / (2 * m + off_b - off_a) ** 6, a, a + m, op_a, op_b)
    return MPOModel(lat, model.calc_H_MPO())


def infinite_pair_model(eps, g, delta, omega=1.0, mmax=8, cell_pairs=2, truncate=False):
    """``pair_model`` with infinite boundary conditions."""
    return pair_model(eps, g, delta, cell_pairs, omega, mmax, truncate, infinite=True)


def tail_energy_per_atom(eps, g, mmax, n_left, n_right):
    """Classical estimate of the couplings beyond mmax pairs, per atom.

    Uses uncorrelated occupations <n_L>, <n_R>; connected correlations beyond
    mmax pairs are negligible away from criticality.  Add it to an mmax-truncated
    energy to estimate the full 1/r^6 model.
    """
    s = 1.0 - eps
    m = np.arange(mmax + 1, mmax + 20001, dtype=float)
    total = 0.0
    for na, off_a in ((n_left, 0.0), (n_right, s)):
        for nb, off_b in ((n_left, 0.0), (n_right, s)):
            total += na * nb * np.sum(g / (2 * m + off_b - off_a) ** 6)
    return float(total / 2.0)


# ------------------------------------------------------------------ states and measurements

def pair_vector(site, label):
    """Local vector of the pair state label (00, 10, 01, 11 or 'sym' = (10 + 01)/sqrt 2).

    On a dressed site the 4-state vector is projected with W^T and normalized.
    """
    v = np.zeros(4)
    if label == "sym":
        v[1] = v[2] = 1.0 / np.sqrt(2.0)
    else:
        v[PAIR_LABELS.index(label)] = 1.0
    u = site.pair_basis.T @ v
    return u / np.linalg.norm(u)


def initial_states(model):
    """Product starts per unit cell: 'L' = |10>^n (tau = +1), 'R' = |01>^n (tau = -1)
    and 'sym' = ((|10> + |01>)/sqrt 2)^n (I_s symmetric)."""
    sites = model.lat.mps_sites()
    return {name: [pair_vector(site, label) for site in sites]
            for name, label in (("L", "10"), ("R", "01"), ("sym", "sym"))}


def energy_per_atom(model, psi):
    """<H> per atom of an infinite MPS, summing every MPO term that starts in one unit cell.

    Uses TeNPy's power-method evaluation (exact for a finite-range MPO).  The
    transfer-matrix variant ``expectation_value_TM`` is not used: for these MPOs it
    returned energies off by up to 3e-4 that jump with its tolerance.
    """
    return float(np.real(model.H_MPO.expectation_value_power(psi, tol=1e-14))) / 2.0


def energy_from_terms(psi, eps, g, delta, omega=1.0, mmax=8):
    """<H> per atom summed term by term from local expectation values and correlations.

    Independent of the MPO; used to check energy_per_atom.
    """
    L = psi.L
    s = 1.0 - eps
    ex = {op: np.real(psi.expectation_value(op)) for op in ("XL", "XR", "nL", "nR", "nLnR")}
    total = float(np.sum(0.5 * omega * (ex["XL"] + ex["XR"]) - delta * (ex["nL"] + ex["nR"])
                         + g / s ** 6 * ex["nLnR"]))
    offsets = {"nL": 0.0, "nR": s}
    for op_a, off_a in offsets.items():
        for op_b, off_b in offsets.items():
            C = np.real(psi.correlation_function(op_a, op_b, sites1=range(L), sites2=range(L + mmax)))
            for a in range(L):
                for m in range(1, mmax + 1):
                    total += g / (2 * m + off_b - off_a) ** 6 * C[a, a + m]
    return total / (2 * L)


def strong_bond_entropy(psi, a):
    """Entanglement entropy (bits) of the cut inside pair a, between its L and R atom."""
    theta = psi.get_theta(a, n=1).to_ndarray()             # vL, p, vR; Schmidt values on the left
    theta = np.einsum("kp,lpr->lkr", psi.sites[a].pair_basis, theta)   # to the 00, 10, 01, 11 basis
    chi_L, _, chi_R = theta.shape
    t = np.zeros((chi_L, 2, 2, chi_R), dtype=theta.dtype)  # vL, nL, nR, vR
    for k, label in enumerate(PAIR_LABELS):
        t[:, int(label[0]), int(label[1]), :] = theta[:, k, :]
    sv = np.linalg.svd(t.reshape(2 * chi_L, 2 * chi_R), compute_uv=False)
    w = sv ** 2
    w = w[w > 1e-16] / w.sum()
    return float(-np.sum(w * np.log2(w)))


def infinite_observables(psi):
    """Observables of an infinite pair-site MPS.

    'pair_order' <n_L - n_R> per pair of the cell, 'order' |cell mean|,
    'order_spread' max - min of pair_order over the cell, 'n_left', 'n_right'
    cell averages, 'density' per atom, 'pair_probabilities' per pair in the order
    00, 10, 01, 11, 'xi' correlation length in units of the mean spacing a (twice
    the value in pairs), 'entropy_weak' (bits, cuts between pairs),
    'entropy_strong' (bits, cuts inside pairs), 'chi'.
    """
    tau = np.real(psi.expectation_value("tau"))
    nL = np.real(psi.expectation_value("nL"))
    nR = np.real(psi.expectation_value("nR"))
    probs = np.column_stack([np.real(psi.expectation_value(f"P{lab}")) for lab in PAIR_LABELS])
    return {
        "pair_order": tau.tolist(), "order": float(abs(tau.mean())),
        "order_spread": float(tau.max() - tau.min()),
        "n_left": float(nL.mean()), "n_right": float(nR.mean()),
        "density": float((nL.mean() + nR.mean()) / 2.0),
        "pair_probabilities": probs.tolist(),
        "xi": 2.0 * float(psi.correlation_length2()),
        "entropy_weak": (np.asarray(psi.entanglement_entropy()) / LN2).tolist(),
        "entropy_strong": [strong_bond_entropy(psi, a) for a in range(psi.L)],
        "chi": [int(c) for c in psi.chi],
    }


# ------------------------------------------------------------------ solver

def lanczos_params(N_max=100, E_tol=1e-12):
    """Lanczos settings that do not stop early when V_s makes the spectrum wide.

    The residual/gap test is disabled (P_tol = inf); Lanczos stops when the
    ground-state Ritz value changes by less than E_tol, or at N_max, with full
    reorthogonalization.
    """
    return {"N_min": 2, "N_max": N_max, "E_tol": E_tol, "P_tol": np.inf, "reortho": True}


def idmrg_options(chi=64, max_sweeps=120, min_sweeps=30, max_E_err=1e-12, max_S_err=1e-9,
                  mixer_sweeps=10, svd_min=1e-8, lanczos_N_max=100):
    """Options for tenpy.algorithms.dmrg.TwoSiteDMRGEngine on an infinite MPS.

    Differences from TeNPy's defaults and from the pilot runs: lanczos_params()
    with the adaptive P_tol/E_tol updates switched off; svd_min = 1e-8 instead of
    1e-12 (the discarded weight per bond stays below chi * svd_min^2); the mixer
    runs for mixer_sweeps sweeps only (TeNPy's default for infinite MPS is 50 with
    slow decay) and at least 20 sweeps follow without it.  Convergence means
    |dE| < max_E_err and |dS| < max_S_err per sweep, checked every 10 sweeps.

    max_sweeps is kept moderate on purpose.  Near a critical point the entropy
    criterion is not met, and the environments keep growing: after ~150 sweeps
    their energy is ~1e3 in magnitude, the absolute Lanczos tolerance 1e-12 is at
    the rounding level, and the state starts to degrade (one run at eps = 0.55
    collapsed after 310 sweeps).  The energy converges well before that.
    """
    return {
        "trunc_params": {"chi_max": chi, "svd_min": svd_min},
        "max_sweeps": max_sweeps, "min_sweeps": max(min_sweeps, mixer_sweeps + 20),
        "max_E_err": max_E_err, "max_S_err": max_S_err,
        "mixer": mixer_sweeps > 0,
        "mixer_params": {"amplitude": 1e-5, "decay": 2.0, "disable_after": mixer_sweeps},
        "lanczos_params": lanczos_params(lanczos_N_max),
        "P_tol_to_trunc": None, "E_tol_to_trunc": None,
        "N_sweeps_check": 10, "update_env": 5,
        "combine": True,
    }


def run_idmrg(model, init, options):
    """One iDMRG run from a product state (list of local vectors); returns (psi, diagnostics)."""
    psi = MPS.from_product_state(model.lat.mps_sites(), list(init), bc="infinite",
                                 unit_cell_width=model.lat.mps_unit_cell_width)
    engine = tenpy_dmrg.TwoSiteDMRGEngine(psi, model, options)
    E_growth, _ = engine.run()
    st = engine.sweep_stats
    n_lanczos = np.asarray(engine.update_stats["N_lanczos"])
    n_lanczos = n_lanczos[n_lanczos > 0]
    diag = {
        "n_sweeps": int(st["sweep"][-1]),
        "last_dE": float(abs(st["Delta_E"][-1])) if len(st["E"]) > 1 else float("nan"),
        "last_dS": float(abs(st["Delta_S"][-1])) if len(st["E"]) > 1 else float("nan"),
        "growth_energy_per_atom": float(E_growth) / 2.0,
        "max_trunc_err_last": float(st["max_trunc_err"][-1]),
        "norm_err": float(np.linalg.norm(psi.norm_test())),
        "lanczos_at_N_max": float(np.mean(n_lanczos >= options["lanczos_params"]["N_max"]))
        if len(n_lanczos) else 0.0,
    }
    diag["energy_converged"] = bool(diag["last_dE"] < options["max_E_err"])
    diag["converged"] = bool(diag["energy_converged"] and diag["last_dS"] < options["max_S_err"])
    return psi, diag


def infinite_ground_state(eps, g, delta, omega=1.0, chi=64, mmax=8, cell_pairs=2,
                          starts=("L", "R"), truncate="auto", max_truncation_error=1e-8,
                          consistency_tol=1e-9, **option_kwargs):
    """Ground state of the infinite dimerized chain: (energy per atom, psi, obs).

    Runs iDMRG from each product start (see initial_states).  The default starts
    are the two Z2 crystals.  'sym' can be added as a check for metastable states,
    but in the ordered phase it tends to a cat of both crystals that converges
    slowly or not at all (at g = 4.83, delta = 1.5, eps = 0.2 it ended at +5.8 per
    atom) and is never the lowest.  A candidate is
    'valid' when energy_per_atom and the term-by-term energy (energy_from_terms)
    agree to consistency_tol.  The kept state is the valid candidate with the
    lowest energy (variational choice).  'energy_converged' means |dE| < max_E_err
    per sweep; 'converged' additionally requires the entropy criterion, which near
    a critical point is not reached within max_sweeps.
    truncate: True, False or "auto" (see use_truncation).  option_kwargs go to
    idmrg_options.  obs holds the infinite_observables of the kept state and

      'energy'            <H> per atom (energy_per_atom),
      'start'             name of the kept start,
      'candidates'        per start: energy, energy_terms, valid, converged,
                          growth estimate, sweeps, order, pair_order, order_spread, xi,
      'start_spread'      max - min energy over the valid candidates,
      'mirror_check'      {'dE': |E_L - E_R|, 'dtau': max |tau_L + tau_R|} (I_s maps
                          the L start onto the R start),
      'truncated'         True if the dressed 3-state basis was used, with
                          'truncation_error' its estimate per atom (0 otherwise),
      'converged', 'energy_converged', 'valid', 'mmax', 'chi_max', 'cell_pairs'.
    """
    truncated = use_truncation(eps, g, delta, omega, truncate, max_truncation_error)
    model = infinite_pair_model(eps, g, delta, omega, mmax, cell_pairs, truncated)
    options = idmrg_options(chi, **option_kwargs)
    inits = initial_states(model)
    candidates, states = {}, {}
    for name in starts:
        psi, diag = run_idmrg(model, inits[name], options)
        e = energy_per_atom(model, psi)
        e_terms = energy_from_terms(psi, eps, g, delta, omega, mmax)
        obs = infinite_observables(psi)
        candidates[name] = dict(energy=e, energy_terms=e_terms,
                                valid=bool(abs(e - e_terms) < consistency_tol), **diag,
                                order=obs["order"], pair_order=obs["pair_order"],
                                order_spread=obs["order_spread"], density=obs["density"], xi=obs["xi"])
        states[name] = (psi, obs)

    def rank(name):
        c = candidates[name]
        return (not c["valid"], c["energy"])

    name = min(candidates, key=rank)
    psi, obs = states[name]
    e = candidates[name]["energy"]
    obs = dict(obs)
    valid_energies = [c["energy"] for c in candidates.values() if c["valid"]] or [e]
    obs.update(energy=e, start=name, candidates=candidates,
               start_spread=float(max(valid_energies) - min(valid_energies)),
               converged=bool(candidates[name]["converged"]), valid=bool(candidates[name]["valid"]),
               energy_converged=bool(candidates[name]["energy_converged"]),
               truncated=truncated,
               truncation_error=truncation_error_estimate(eps, g, delta, omega) if truncated else 0.0,
               mmax=int(mmax), chi_max=int(chi), cell_pairs=int(cell_pairs))
    if "L" in candidates and "R" in candidates:
        tL = np.asarray(candidates["L"]["pair_order"])
        tR = np.asarray(candidates["R"]["pair_order"])
        obs["mirror_check"] = {"dE": float(abs(candidates["L"]["energy"] - candidates["R"]["energy"])),
                               "dtau": float(np.max(np.abs(tL + tR)))}
    return e, psi, obs


# ------------------------------------------------------------------ finite reference

def pinned_local_delta(n_pairs, delta, pin=None):
    """local_delta for an open chain that ends in |...10> without a domain wall.

    Both ends of an even open chain favour an excited edge atom, which forces a
    delocalized 00 domain wall into the Z2 crystal and makes long-chain DMRG converge
    very slowly.  Lowering the detuning of the last atom by pin (default 2*delta,
    so Delta_last = -delta) removes the wall.  It is a boundary term and cancels
    in the bulk estimate E(N2) - E(N1).
    """
    loc = np.zeros(2 * n_pairs)
    loc[-1] = -(2.0 * delta if pin is None else pin)
    return loc


def pinned_chain_ground_state(eps, g, delta, n_pairs, omega=1.0, chi=64, mmax=8, truncate=False,
                              max_sweeps=80, min_sweeps=6, max_E_err=1e-12, max_S_err=1e-8,
                              mixer_sweeps=6, svd_min=1e-10, lanczos_N_max=100):
    """Finite pair-site DMRG of a pinned open chain of 2*n_pairs atoms: (energy, psi, info).

    Same Hamiltonian as dmrg.pair_chain_model (mmax cutoff) plus
    pinned_local_delta, started from |10>^n, with lanczos_params() (TeNPy's
    default stopping rule is unreliable for finite chains too when V_s is large).
    svd_min = 1e-10 drops noise-level Schmidt values that would otherwise fill chi.
    truncate=True uses the dressed 3-state basis (for testing that approximation).
    """
    model = pair_model(eps, g, delta, n_pairs, omega, mmax, truncate, infinite=False,
                       local_delta=pinned_local_delta(n_pairs, delta))
    init = [pair_vector(site, "10") for site in model.lat.mps_sites()]
    psi = MPS.from_product_state(model.lat.mps_sites(), init, bc="finite",
                                 unit_cell_width=model.lat.mps_unit_cell_width)
    options = dmrg_options(chi, max_sweeps=max_sweeps, min_sweeps=min_sweeps, max_E_err=max_E_err,
                           max_S_err=max_S_err, svd_min=svd_min, mixer=mixer_sweeps > 0,
                           mixer_sweeps=mixer_sweeps)
    options.update(lanczos_params=lanczos_params(lanczos_N_max), P_tol_to_trunc=None,
                   E_tol_to_trunc=None)
    engine = tenpy_dmrg.TwoSiteDMRGEngine(psi, model, options)
    energy, _ = engine.run()
    st = engine.sweep_stats
    dE = abs(st["Delta_E"][-1]) if len(st["E"]) > 1 else float("nan")
    dS = abs(st["Delta_S"][-1]) if len(st["E"]) > 1 else float("nan")
    tau = np.real(psi.expectation_value("tau"))
    occ = np.real(psi.expectation_value("nL")) + np.real(psi.expectation_value("nR"))
    mid = slice(n_pairs // 3, 2 * n_pairs // 3)
    info = {"energy": float(energy), "n_sweeps": len(st["E"]), "last_sweep_dE": float(dE),
            "last_sweep_dS": float(dS), "converged": bool(dE < max_E_err and dS < max_S_err),
            "max_chi": int(max(psi.chi)), "final_trunc_err": float(st["max_trunc_err"][-1]),
            "bulk_order": float(np.mean(tau[mid])), "bulk_density": float(np.mean(occ[mid]) / 2.0),
            "truncated": bool(truncate)}
    return float(energy), psi, info


def finite_bulk_energy(eps, g, delta, n_pairs=(50, 100), omega=1.0, chi=64, mmax=8, **kwargs):
    """Bulk energy per atom [E(N2) - E(N1)]/(N2 - N1) from two pinned open chains.

    n_pairs are the chain lengths in pairs (N = 2 n_pairs atoms); kwargs go to
    pinned_chain_ground_state.  Returns (e_bulk, info) with info keyed by the atom
    number N.  The difference removes the boundary energy; the remaining
    finite-size error is exponentially small in N/xi in gapped phases and of
    order v/(N1 N2) at a critical point.
    """
    info = {}
    for n in n_pairs:
        _, _, info[2 * n] = pinned_chain_ground_state(eps, g, delta, n, omega, chi, mmax, **kwargs)
    N1, N2 = 2 * n_pairs[0], 2 * n_pairs[1]
    return float((info[N2]["energy"] - info[N1]["energy"]) / (N2 - N1)), info


def bulk_energy_fit(energies):
    """Fit E(N) = e N + b + c/N [+ d/N^2] to energies {N: E}; returns (e, b, c, d).

    Two sizes: c = d = 0, i.e. e = [E(N2) - E(N1)]/(N2 - N1).  Three sizes add the
    1/N (Casimir) term of a critical open chain, which removes the O(c/(N1 N2)) bias
    of the two-size difference; four or more also add 1/N^2.
    """
    N = np.array(sorted(energies), dtype=float)
    E = np.array([energies[k] for k in sorted(energies)], dtype=float)
    cols = [N, np.ones_like(N), 1.0 / N, 1.0 / N ** 2][:min(len(N), 4)]
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), E, rcond=None)
    coef = np.concatenate([coef, np.zeros(4 - len(coef))])
    return tuple(float(c) for c in coef)
