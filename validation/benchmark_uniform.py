"""External benchmarks of the finite-DMRG pipeline on two known critical points.

Part 1, exactly solvable: the open transverse-field Ising chain
    H = -J sum_i Sx_i Sx_{i+1} - h sum_i Sz_i   (Pauli matrices, J = h = 1)
at criticality.  DMRG energies and entanglement entropies are compared with the
exact free-fermion values (2N x 2N Majorana matrix), and the central charge is
fitted with the Calabrese-Cardy open-chain formula.

Part 2, literature value: the uniform Rydberg chain with nearest-neighbour
blockade, i.e. the hard-boson model of Fendley, Sengupta and Sachdev
(PRB 69, 075106 (2004)) at V = 0.  In the conventions of rydberg_chain,
    H = (Omega/2) sum_i X_i - Delta sum_i n_i + V1 sum_i n_i n_{i+1},  V1 -> infinity,
the published Ising point U_c/w = -1.31 maps to Delta_c/Omega = 0.655
(w = Omega/2, U = -Delta; see BENCHMARK_NOTES.md).  The transition is located
in odd open chains, whose free ends pin the same Z2 pattern at both ends, from
the crossing of N^(beta/nu) m_c with beta/nu = 1/8, where m_c is the staggered
density at the central atom.  The exact constrained model (P X P plus a
penalty that lifts the blockade-violating sectors) is compared with finite V1
runs of rydberg_chain.dmrg.atom_chain_model (g = V1, rmax = 1), and c is fitted
at the transition.

Part 3, internal cross-check (optional part 'vdw'): the same crossing analysis
for the van der Waals chain at R_b/a = 1.5 (g = 11.39, couplings to fourth
neighbours), compared with the preliminary ring-ED value of the pilot study.

Run from the repository root (each worker is single-threaded):
    python validation/benchmark_uniform.py [--quick] [--workers 2] [--parts tfim,blockade,vdw]
    python validation/benchmark_uniform.py --reanalyse    # refit from the stored raw data
Results go to validation/results/benchmark_uniform.json; see BENCHMARK_NOTES.md.
"""

import argparse
import json
import multiprocessing
import os
import platform
import sys
import time

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")      # one thread per worker; see --workers

import numpy as np                       # noqa: E402
from scipy.linalg import eigh            # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results",
                       "benchmark_uniform.json")
LITERATURE = {
    "model": "hard-boson chain, H = sum_j [-w (d_j + d_j^dag) + U n_j], n_j n_{j+1} = 0",
    "U_c_over_w": -1.31,
    "Delta_c_over_Omega": 0.655,
    "digits": "three significant figures as published",
    "source": "P. Fendley, K. Sengupta, S. Sachdev, Phys. Rev. B 69, 075106 (2004), V = 0 axis "
              "of the phase diagram; same point as the tilted-lattice Ising transition of "
              "S. Sachdev, K. Sengupta, S. M. Girvin, Phys. Rev. B 66, 075128 (2002) "
              "[(E-U)_c = 1.85 J with hopping sqrt(2) J]",
    "status": "quoted from memory of the papers; the text could not be re-read in this "
              "environment (search quota exhausted, publisher and preprint hosts blocked): verify",
}


# --------------------------------------------------------------- helpers

def _tenpy():
    import tenpy
    tenpy.tools.misc.setup_logging(to_stdout=None)
    return tenpy


def _dmrg_run(psi, model, options):
    """Run TeNPy finite DMRG in place; return energy and a small diagnostics dict."""
    from tenpy.algorithms import dmrg as tenpy_dmrg
    info = tenpy_dmrg.run(psi, model, options)
    stats = info["sweep_statistics"]
    dE = abs(stats["Delta_E"][-1]) if len(stats["E"]) > 1 else float("nan")
    dS = abs(stats["Delta_S"][-1]) if len(stats["E"]) > 1 else float("nan")
    return float(info["E"]), {
        "n_sweeps": len(stats["E"]),
        "converged": bool(dE < options["max_E_err"] and dS < options["max_S_err"]),
        "last_sweep_dE": float(dE),
        "final_trunc_err": float(stats["max_trunc_err"][-1]),
        "max_chi": int(max(psi.chi)),
    }


def calabrese_cardy_fit(S, n, lmin=None, parity=None, lmax=None):
    """Fit S(l) = (c/6) ln[(2n/pi) sin(pi l/n)] + s0 on cuts lmin <= l <= lmax.

    lmax defaults to n - lmin.
    S[l-1] is the entropy (nats) of the first l sites, l = 1..n-1.  parity = 0 or 1
    keeps only even or odd l; parity = 'smooth' fits the three-point average
    (S(l-1) + 2 S(l) + S(l+1))/4, which cancels a bond-alternating term to leading
    order.  Returns (c, s0, rms residual).
    """
    S = np.asarray(S, dtype=float)
    lmin = max(2, n // 8) if lmin is None else lmin
    ell = np.arange(1, n)
    y = S.copy()
    if parity == "smooth":
        y[1:-1] = 0.25 * (S[:-2] + 2 * S[1:-1] + S[2:])
    lmax = n - lmin if lmax is None else lmax
    keep = (ell >= lmin) & (ell <= lmax)
    if parity in (0, 1):
        keep &= ell % 2 == parity
    x = np.log(2 * n / np.pi * np.sin(np.pi * ell[keep] / n)) / 6.0
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y[keep], rcond=None)
    rms = float(np.sqrt(np.mean((A @ coef - y[keep]) ** 2)))
    return float(coef[0]), float(coef[1]), rms


# --------------------------------------------------------------- part 1: TFIM

def tfim_majorana_matrix(n, J=1.0, h=1.0):
    """Real antisymmetric M with H = (i/4) sum_ab M_ab g_a g_b for the open TFIM."""
    M = np.zeros((2 * n, 2 * n))
    M[np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)] = 2 * h
    M[np.arange(1, 2 * n - 1, 2), np.arange(2, 2 * n, 2)] = 2 * J
    return M - M.T


def tfim_exact(n, J=1.0, h=1.0, cuts=None):
    """Exact ground energy and entropies S(l) (nats) of the open TFIM.

    S is returned for l = 1..n-1, or only for the cuts given (cost ~ l^3 each).

    The single-particle energies are the positive eigenvalues eps_k of iM and
    E0 = -sum_k eps_k / 2.  The Majorana correlations <g_a g_b> = delta_ab + i G_ab
    with G = -i sign(iM) give S(l) from the eigenvalues +-nu of i G restricted to
    the first 2l Majoranas.
    """
    w, U = eigh(1j * tfim_majorana_matrix(n, J, h))
    energy = -0.5 * float(np.sum(w[w > 0]))
    G = np.real(-1j * (U * np.sign(w)) @ U.conj().T)
    cuts = range(1, n) if cuts is None else cuts
    S = np.empty(len(cuts))
    for k, ell in enumerate(cuts):
        nu = np.linalg.eigvalsh(1j * G[:2 * ell, :2 * ell])
        p = np.clip((1 + nu[nu > 0]) / 2, 1e-300, 1.0)
        q = np.clip(1 - p, 1e-300, 1.0)
        S[k] = -np.sum(p * np.log(p) + q * np.log(q))
    return energy, S


def tfim_exact_central_charge(n, n_cuts=24):
    """Calabrese-Cardy c from exact entropies at n_cuts even cuts in [n/8, 7n/8]."""
    cuts = np.unique(2 * np.round(np.linspace(n // 8, n - n // 8, n_cuts) / 2).astype(int))
    S = tfim_exact(n, cuts=cuts)[1]
    x = np.log(2 * n / np.pi * np.sin(np.pi * cuts / n)) / 6.0
    return float(np.polyfit(x, S, 1)[0])


def tfim_dmrg(n, chi=64, J=1.0, h=1.0):
    """Finite DMRG of the open TFIM with TeNPy's TFIChain (parity conserved)."""
    _tenpy()
    from tenpy.models.tf_ising import TFIChain
    from tenpy.networks.mps import MPS
    from rydberg_chain.dmrg import dmrg_options
    model = TFIChain({"L": n, "J": J, "g": h, "bc_MPS": "finite", "conserve": "parity"})
    psi = MPS.from_product_state(model.lat.mps_sites(), ["up"] * n, bc="finite",
                                 unit_cell_width=model.lat.mps_unit_cell_width)
    energy, diag = _dmrg_run(psi, model, dmrg_options(chi, max_sweeps=40, mixer_sweeps=8))
    return energy, np.asarray(psi.entanglement_entropy()), diag


# --------------------------------------------------------------- part 2: blockaded chain

def two_level_site():
    """Two-level TeNPy site with operators X, N = |1><1| and P = |0><0|."""
    import tenpy.linalg.np_conserved as npc
    from tenpy.networks.site import Site
    X = np.array([[0.0, 1.0], [1.0, 0.0]])
    N = np.diag([0.0, 1.0])
    return Site(npc.LegCharge.from_trivial(2), ["0", "1"], X=X, N=N, P=np.eye(2) - N)


def blockade_chain_model(n, delta, omega=1.0, penalty=None):
    """Exact nearest-neighbour blockade (hard-boson) chain as a TeNPy MPOModel.

    H = (Omega/2) sum_i P_{i-1} X_i P_{i+1} - Delta sum_i n_i + penalty sum_i n_i n_{i+1}.
    The P X P drive conserves the number of blockade violations, so the
    blockaded sector is exactly the hard-boson model; the penalty (default
    10 + 2|Delta|) only lifts the other sectors above it.
    """
    _tenpy()
    from tenpy.models.lattice import Chain
    from tenpy.models.model import CouplingModel, MPOModel
    penalty = 10.0 + 2.0 * abs(delta) if penalty is None else penalty
    lat = Chain(n, two_level_site(), bc="open", bc_MPS="finite")
    model = CouplingModel(lat)
    for i in range(n):
        model.add_onsite_term(-delta, i, "N")
    model.add_coupling_term(0.5 * omega, 0, 1, "X", "P")
    model.add_coupling_term(0.5 * omega, n - 2, n - 1, "P", "X")
    for i in range(1, n - 1):
        model.add_multi_coupling_term(0.5 * omega, [i - 1, i, i + 1], ["P", "X", "P"], ["Id", "Id"])
    for i in range(n - 1):
        model.add_coupling_term(penalty, i, i + 1, "N", "N")
    return MPOModel(lat, model.calc_H_MPO())


def blockade_chain_exact(n, delta, omega=1.0):
    """Ground energy of the blockaded chain by ED in the Fibonacci basis (small n)."""
    states = [s for s in range(2 ** n) if not s & (s >> 1)]
    index = {s: k for k, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s, k in index.items():
        H[k, k] = -delta * bin(s).count("1")
        for i in range(n):
            t = s ^ (1 << i)
            if t in index:
                H[index[t], k] += 0.5 * omega
    return float(np.linalg.eigvalsh(H)[0])


def central_order(occ):
    """Staggered density at the central atom of an odd chain, sign fixed by the end pattern."""
    n = len(occ)
    c = (n - 1) // 2
    return float((-1) ** c * (occ[c] - 0.5 * (occ[c - 1] + occ[c + 1])))


def blockade_point(n, delta, V1=None, chi=64, max_sweeps=30, rmax=1):
    """Ground state of the uniform chain at one Delta, from the Z2 product state.

    V1 = None uses the exact blockade model; a number uses rydberg_chain's atom
    chain with g = V1 and couplings V1/r^6 up to index distance rmax (default 1:
    nearest neighbours only; rmax = None keeps all pairs).  Each point
    starts afresh with the mixer: warm starts from a neighbouring Delta without
    the mixer relax the order-parameter profile of long critical chains very
    slowly.  Returns E, m_c, bond entropies (nats) and DMRG diagnostics.
    """
    _tenpy()
    from tenpy.networks.mps import MPS
    from rydberg_chain import dmrg, geometry
    if V1 is None:
        model = blockade_chain_model(n, delta)
    else:
        model = dmrg.atom_chain_model(geometry.chain_positions(n), V1, delta, rmax=rmax)
    psi = MPS.from_product_state(model.lat.mps_sites(),
                                 [str(v) for v in dmrg.default_product_state(n)], bc="finite",
                                 unit_cell_width=model.lat.mps_unit_cell_width)
    t0 = time.time()
    energy, diag = _dmrg_run(psi, model, dmrg.dmrg_options(chi, max_sweeps=max_sweeps))
    occ = np.real(psi.expectation_value("N"))
    return {"Delta": float(delta), "E": energy, "m_c": central_order(occ),
            "density": float(occ.mean()),
            "S_bonds": [float(x) for x in psi.entanglement_entropy()],
            "seconds": time.time() - t0, **diag}


def blockade_scan(n, deltas, V1=None, chi=64):
    """blockade_point for each Delta in deltas."""
    return [blockade_point(n, d, V1, chi) for d in deltas]


def crossing(deltas, y1, y2, deg=3):
    """Delta where two smooth curves y1(Delta), y2(Delta) cross (polynomial fit)."""
    d = np.asarray(deltas)
    d0 = d.mean()
    p = np.polyfit(d - d0, np.asarray(y1) - np.asarray(y2), min(deg, len(d) - 1))
    roots = [r.real + d0 for r in np.roots(p) if abs(r.imag) < 1e-12]
    roots = [r for r in roots if d.min() <= r <= d.max()]
    if len(roots) != 1:
        return float("nan")
    return float(roots[0])


def extrapolate_crossings(pairs, powers=(1.0, 1.5, 2.0, 2.5)):
    """Fit Delta*(N1, N2) = Delta_c + a Nbar^-p, Nbar = sqrt(N1 N2).

    One two-parameter fit per fixed p, and with four or more pairs a fit with
    p free ('free').  'estimate' is the free fit (or p = 2) with an uncertainty
    equal to the largest deviation of the fixed-p fits with |p - p_free| <= 0.5.
    """
    Nbar = np.array([np.sqrt(a * b) for a, b, _ in pairs])
    y = np.array([c for _, _, c in pairs])
    if len(pairs) < 2:
        return {"estimate": {"Delta_c": float(y[-1]), "uncertainty": float("nan")}}
    fits = {}
    for p in powers:
        A = np.column_stack([np.ones_like(Nbar), Nbar ** -p])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        fits[f"p={p:g}"] = {"Delta_c": float(coef[0]), "a": float(coef[1]), "p": float(p)}
    best = fits["p=2"]
    if len(pairs) >= 4:
        from scipy.optimize import curve_fit
        popt, _ = curve_fit(lambda N, dc, a, p: dc + a * N ** -p, Nbar, y,
                            p0=[y[-1], (y[0] - y[-1]) * Nbar[0] ** 2, 2.0])
        best = fits["free"] = {"Delta_c": float(popt[0]), "a": float(popt[1]), "p": float(popt[2])}
    near = [f["Delta_c"] for f in fits.values() if abs(f["p"] - best["p"]) <= 0.5 + 1e-12]
    fits["estimate"] = {"Delta_c": best["Delta_c"],
                        "uncertainty": float(max(abs(np.asarray(near) - best["Delta_c"])))}
    return fits


def analyse_crossings(scans, sizes, deltas):
    """Scaled order N^(1/8) m_c per size and the crossings of consecutive sizes."""
    y = {n: [n ** 0.125 * r["m_c"] for r in scans[n]] for n in sizes}
    pairs = [(a, b, crossing(deltas, y[a], y[b])) for a, b in zip(sizes[:-1], sizes[1:])]
    return y, pairs


def blockade_ed_checks():
    """DMRG of the exact blockade model against Fibonacci-basis ED."""
    checks = []
    for n, delta in [(9, 0.655), (13, 0.655), (13, 1.2)]:
        E_ex = blockade_chain_exact(n, delta)
        E = blockade_scan(n, [delta])[0]["E"]
        checks.append({"N": n, "Delta": delta, "E_dmrg": E, "E_exact": E_ex, "dE": E - E_ex})
    return checks


def entropy_fits(S, n):
    """Calabrese-Cardy fits of bond entropies S (nats) of an odd, reflection-symmetric chain.

    Reflection maps cut l to n - l, which has the other parity, so even-only and
    odd-only fits over a symmetric window coincide; they are taken on the left
    half (n/8 <= l <= n/2) instead.  'c_smooth' fits the three-point average on
    n/8 <= l <= 7n/8.
    """
    S = np.asarray(S)
    return {"c_even_left": calabrese_cardy_fit(S, n, parity=0, lmax=n // 2)[0],
            "c_odd_left": calabrese_cardy_fit(S, n, parity=1, lmax=n // 2)[0],
            "c_smooth": calabrese_cardy_fit(S, n, parity="smooth")[0],
            "c_all": calabrese_cardy_fit(S, n)[0]}


def central_entropy_c(bond_entropies):
    """c = 6 [S_mid(N2) - S_mid(N1)] / ln(N2/N1) for consecutive odd sizes.

    S_mid is the entropy of the central cut, where the chord length is 2N/pi and
    the bond-alternating term vanishes by reflection symmetry.
    """
    sizes = sorted(bond_entropies)
    mid = {n: bond_entropies[n][(n - 1) // 2 - 1] for n in sizes}
    return [[a, b, float(6 * (mid[b] - mid[a]) / np.log(b / a))]
            for a, b in zip(sizes[:-1], sizes[1:])]


def _task(spec):
    """Worker entry point: ('tfim', n), ('blockade', n, V1, Delta) or ('vdw', n, g, Delta)."""
    t0 = time.time()
    if spec[0] == "tfim":
        n = spec[1]
        E, S, diag = tfim_dmrg(n)
        E_ex, S_ex = tfim_exact(n)
        out = {"N": n, "chi": 64, "E_dmrg": E, "E_exact": E_ex, "dE": E - E_ex,
               "max_abs_dS": float(np.max(np.abs(S - S_ex))),
               "c_dmrg": calabrese_cardy_fit(S, n)[0],
               "c_dmrg_even": calabrese_cardy_fit(S, n, parity=0)[0],
               "c_dmrg_odd": calabrese_cardy_fit(S, n, parity=1)[0],
               "c_exact": calabrese_cardy_fit(S_ex, n)[0], **diag}
    elif spec[0] == "vdw":
        _, n, g, delta = spec
        out = blockade_point(n, delta, g, rmax=VDW_RMAX)
    else:
        _, n, V1, delta = spec
        out = blockade_point(n, delta, V1)
    return spec, out, time.time() - t0


# --------------------------------------------------------------- driver

VDW_G = 11.39            # (R_b/a)^6 at R_b/a = 1.5, as in research/pilots/2026-09-27
VDW_RMAX = 4             # V_5/V_1 = 1/15625 is dropped
VDW_PILOT_DELTA_C = 1.0785   # preliminary ring ED (8/10-pair crossing), not a literature value


def _settings(quick):
    if quick:
        return dict(tfim_sizes=[16, 24], sizes=[17, 25], deltas=[0.63, 0.645, 0.66, 0.675],
                    v1_values=[50.0], v1_sizes=[17, 25], v1_deltas=[0.63, 0.645, 0.66, 0.675],
                    vdw_sizes=[17, 25], vdw_deltas=[1.05, 1.07, 1.09, 1.11])
    return dict(tfim_sizes=[32, 48, 64, 96, 128], sizes=[33, 49, 65, 97, 129],
                deltas=[0.645, 0.65, 0.655, 0.66, 0.665],
                v1_values=[50.0, 100.0, 200.0], v1_sizes=[33, 49],
                v1_deltas=[0.64, 0.655, 0.67, 0.685, 0.7],
                vdw_sizes=[33, 49, 65], vdw_deltas=[1.06, 1.07, 1.08, 1.09, 1.1])


def _tasks(cfg, parts):
    tasks = []
    if "tfim" in parts:
        tasks += [("tfim", n) for n in cfg["tfim_sizes"]]
    if "blockade" in parts:
        tasks += [("blockade", n, None, d) for n in cfg["sizes"] for d in cfg["deltas"]]
        tasks += [("blockade", n, V1, d) for V1 in cfg["v1_values"] for n in cfg["v1_sizes"]
                  for d in cfg["v1_deltas"]]
    if "vdw" in parts:
        tasks += [("vdw", n, VDW_G, d) for n in cfg["vdw_sizes"] for d in cfg["vdw_deltas"]]
    return sorted(tasks, key=lambda t: -t[1])


def _extrapolate_c(sizes, cs):
    """Linear fit c(N) = c_inf + a/N (the open-chain fits approach c roughly as 1/N)."""
    a, c_inf = np.polyfit(1.0 / np.asarray(sizes, dtype=float), np.asarray(cs), 1)
    return {"c_inf": float(c_inf), "a": float(a)}


def _analyse_tfim(cfg, done):
    rows = [done[("tfim", n)] for n in cfg["tfim_sizes"]]
    exact_large = {str(n): tfim_exact_central_charge(n) for n in (256, 512, 1024)}
    return {"rows": rows, "c_exact_large_N": exact_large,
            "c_dmrg_extrapolated": _extrapolate_c(cfg["tfim_sizes"], [r["c_dmrg"] for r in rows]),
            "c_exact_extrapolated": _extrapolate_c([int(n) for n in exact_large],
                                                   list(exact_large.values()))}


def _analyse_blockade(cfg, done):
    sizes, deltas = cfg["sizes"], cfg["deltas"]
    scans = {n: [done[("blockade", n, None, d)] for d in deltas] for n in sizes}
    y, pairs = analyse_crossings(scans, sizes, deltas)
    fits = extrapolate_crossings(pairs)
    k = int(np.argmin(np.abs(np.asarray(deltas) - fits["estimate"]["Delta_c"])))
    blockade = {"ed_checks": blockade_ed_checks(), "deltas": deltas,
                "scans": {str(n): scans[n] for n in sizes},
                "scaled_order": {str(n): y[n] for n in sizes},
                "pairs": [list(p) for p in pairs], "extrapolations": fits,
                "central_charge": {str(n): {"Delta": deltas[k], "max_chi": scans[n][k]["max_chi"],
                                            **entropy_fits(scans[n][k]["S_bonds"], n)}
                                   for n in sizes}}
    finite_v1 = {}
    for V1 in cfg["v1_values"]:
        v_scans = {n: [done[("blockade", n, V1, d)] for d in cfg["v1_deltas"]]
                   for n in cfg["v1_sizes"]}
        _, v_pairs = analyse_crossings(v_scans, cfg["v1_sizes"], cfg["v1_deltas"])
        finite_v1[str(V1)] = {"pairs": [list(p) for p in v_pairs],
                              "scans": {str(n): v_scans[n] for n in cfg["v1_sizes"]}}
    if all(n in sizes for n in cfg["v1_sizes"]):      # same size pair, V1 = infinity
        _, ref_pairs = analyse_crossings(scans, cfg["v1_sizes"], deltas)
        finite_v1["inf"] = {"pairs": [list(p) for p in ref_pairs]}
    blockade["finite_V1"] = finite_v1
    cc = blockade["central_charge"]
    blockade["c_extrapolated"] = {key: _extrapolate_c(sizes, [cc[str(n)][key] for n in sizes])
                                  for key in ("c_even_left", "c_odd_left", "c_smooth")}
    blockade["c_from_central_entropy"] = central_entropy_c(
        {n: scans[n][k]["S_bonds"] for n in sizes})
    return blockade


def _analyse_vdw(cfg, done):
    sizes, deltas = cfg["vdw_sizes"], cfg["vdw_deltas"]
    scans = {n: [done[("vdw", n, VDW_G, d)] for d in deltas] for n in sizes}
    y, pairs = analyse_crossings(scans, sizes, deltas)
    k = int(np.argmin(np.abs(np.asarray(deltas) - pairs[-1][2])))
    return {"g": VDW_G, "R_b_over_a": VDW_G ** (1 / 6), "rmax": VDW_RMAX,
            "pilot_ring_ed_Delta_c": VDW_PILOT_DELTA_C, "deltas": deltas,
            "scans": {str(n): scans[n] for n in sizes},
            "scaled_order": {str(n): y[n] for n in sizes},
            "pairs": [list(p) for p in pairs], "extrapolations": extrapolate_crossings(pairs),
            "c_from_central_entropy": central_entropy_c(
                {n: scans[n][k]["S_bonds"] for n in sizes}),
            "central_charge": {"N": sizes[-1], "Delta": deltas[k],
                               **entropy_fits(scans[sizes[-1]][k]["S_bonds"], sizes[-1])}}


def _done_from_json(result, cfg, parts):
    """Rebuild the per-task results of a previous run from its JSON file."""
    done = {}
    if "tfim" in parts:
        done.update({("tfim", r["N"]): r for r in result["tfim"]["rows"]})
    if "blockade" in parts:
        blk = result["blockade"]
        for n, scan in blk["scans"].items():
            done.update({("blockade", int(n), None, r["Delta"]): r for r in scan})
        for V1 in cfg["v1_values"]:
            for n, scan in blk["finite_V1"][str(V1)]["scans"].items():
                done.update({("blockade", int(n), V1, r["Delta"]): r for r in scan})
    if "vdw" in parts:
        for n, scan in result["vdw_check"]["scans"].items():
            done.update({("vdw", int(n), VDW_G, r["Delta"]): r for r in scan})
    return done


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true", help="small sizes, for a smoke run")
    parser.add_argument("--parts", default="tfim,blockade,vdw",
                        help="comma list of tfim, blockade, vdw; results of other parts "
                             "already in --out are kept")
    parser.add_argument("--workers", type=int, default=2,
                        help="parallel DMRG processes (each single-threaded by default)")
    parser.add_argument("--out", default=RESULTS)
    parser.add_argument("--reanalyse", action="store_true",
                        help="recompute the fits of --parts from the raw data stored in --out")
    args = parser.parse_args(argv)
    cfg = _settings(args.quick)
    parts = [p.strip() for p in args.parts.split(",") if p.strip()]
    t_start = time.time()

    result = {}
    if os.path.exists(args.out) and (args.reanalyse or set(parts) != {"tfim", "blockade", "vdw"}):
        with open(args.out) as f:
            result = json.load(f)
    if args.reanalyse:
        done = _done_from_json(result, cfg, parts)
    else:
        done = {}
        with multiprocessing.get_context("fork").Pool(args.workers) as pool:
            for spec, out, seconds in pool.imap_unordered(_task, _tasks(cfg, parts), chunksize=1):
                done[spec] = out
                out["seconds"] = seconds
                tag = "" if spec[0] == "tfim" else f" sweeps={out['n_sweeps']} conv={out['converged']}"
                print(f"{spec} {seconds:.0f} s{tag}", flush=True)
    result.update({"script": "validation/benchmark_uniform.py", "python": platform.python_version(),
                   "literature": LITERATURE})
    result.setdefault("runs", []).append({"argv": list(argv if argv is not None else sys.argv[1:]),
                                          "parts": parts, "settings": cfg,
                                          "workers": args.workers, "reanalyse": args.reanalyse,
                                          "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                                          "seconds": time.time() - t_start})
    result["quick"] = bool(result.get("quick")) or args.quick
    analyse = {"tfim": _analyse_tfim, "blockade": _analyse_blockade, "vdw": _analyse_vdw}
    names = {"vdw": "vdw_check"}
    for part in parts:
        result[names.get(part, part)] = analyse[part](cfg, done)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=1)

    if "tfim" in parts:
        for r in result["tfim"]["rows"]:
            print(f"TFIM N={r['N']}: dE={r['dE']:.1e} max|dS|={r['max_abs_dS']:.1e} "
                  f"c={r['c_dmrg']:.4f}")
    if "blockade" in parts:
        blk = result["blockade"]
        for a, b, c in blk["pairs"]:
            print(f"blockade crossing ({a},{b}): {c:.5f}")
        for V1, d in blk["finite_V1"].items():
            print(f"V1={V1}: crossings {[round(p[2], 5) for p in d['pairs']]}")
        est = blk["extrapolations"]["estimate"]
        cc = blk["c_extrapolated"]["c_smooth"]["c_inf"]
        print(f"Delta_c = {est['Delta_c']:.5f} +- {est['uncertainty']:.5f}; "
              f"c (1/N extrapolation) = {cc:.4f}; central-entropy c = "
              f"{[round(x[2], 4) for x in blk['c_from_central_entropy']]}")
    if "vdw" in parts:
        v = result["vdw_check"]
        print(f"vdW R_b/a=1.5 crossings {[(a, b, round(c, 5)) for a, b, c in v['pairs']]}; "
              f"estimate {v['extrapolations']['estimate']}; c(N={v['central_charge']['N']}) = "
              f"{v['central_charge']['c_smooth']:.4f}")
    print(f"{time.time() - t_start:.0f} s -> {args.out}", flush=True)
    return result


if __name__ == "__main__":
    main()
