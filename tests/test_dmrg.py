"""Finite DMRG (atom sites and 4-state pair sites) against exact diagonalization."""

import numpy as np
import pytest

from rydberg_chain import ed, geometry, observables as ob

dmrg = pytest.importorskip("rydberg_chain.dmrg")

E0_SAVED = -18.399305357538882     # N = 12, g = 8, Delta = 3, uniform, all pairs
G, DELTA = 11.39, 2.5              # R_b/a = 1.5


def _exact(x, g, delta, V=None, local_delta=None):
    V = geometry.coupling_matrix(x, g) if V is None else V
    E, v, res = ed.lowest_states(ed.build_hamiltonian(V, 1.0, delta, local_delta), k=2)
    assert res.max() < 1e-9
    return E, v[:, 0]


def test_default_product_state():
    assert dmrg.default_product_state(5) == [1, 0, 1, 0, 1]
    assert dmrg.default_product_state(8) == [1, 0, 1, 0, 0, 1, 0, 1]
    assert dmrg.default_product_state(6) == [1, 0, 1, 0, 0, 1]      # wall on a weak bond
    assert dmrg.default_product_state(2) == [1, 0]


def test_pair_site_operators():
    site = dmrg.pair_site()
    op = {name: site.get_op(name).to_ndarray() for name in ("nL", "nR", "XL", "XR", "nLnR", "tau")}
    basis = {"00": (0, 0), "10": (1, 0), "01": (0, 1), "11": (1, 1)}
    index = {label: site.state_labels[label] for label in basis}
    for label, (nl, nr) in basis.items():
        k = index[label]
        assert op["nL"][k, k] == nl and op["nR"][k, k] == nr and op["nLnR"][k, k] == nl * nr
        flipped_L = f"{1 - nl}{nr}"
        flipped_R = f"{nl}{1 - nr}"
        assert op["XL"][index[flipped_L], k] == 1 and op["XL"][:, k].sum() == 1
        assert op["XR"][index[flipped_R], k] == 1 and op["XR"][:, k].sum() == 1
    assert np.allclose(op["tau"], op["nL"] - op["nR"])
    assert np.allclose(sum(site.get_op(f"P{s}").to_ndarray() for s in basis), np.eye(4))


@pytest.mark.parametrize("n,eps", [(10, 0.0), (11, 0.3), (12, 0.3)])
def test_atom_dmrg_matches_ed(n, eps):
    x = geometry.chain_positions(n, eps)
    E, psi = _exact(x, G, DELTA)
    energy, mps, obs = dmrg.atom_chain_ground_state(x, G, DELTA, chi=64, correlations=True)
    assert obs["converged"] and obs["final_trunc_err"] < 1e-12
    assert abs(energy - E[0]) < 1e-9
    assert np.allclose(obs["occupations"], ob.site_occupations(psi, n), atol=1e-5)
    assert np.allclose(obs["pair_order"], ob.pair_order(psi, n), atol=1e-5)
    assert obs["half_chain_entropy"] == pytest.approx(ob.entanglement_entropy(psi, n, n // 2), abs=1e-5)
    assert np.allclose(obs["nn"], ob.density_correlations(psi, n), atol=1e-5)
    assert obs["staggered_order_squared"] == pytest.approx(ob.staggered_order_squared(psi, n), abs=1e-5)


def test_atom_dmrg_reproduces_saved_uniform_energy():
    energy, _, obs = dmrg.atom_chain_ground_state(geometry.chain_positions(12), 8.0, 3.0, chi=64,
                                                  variance=True)
    assert abs(energy - E0_SAVED) < 1e-9
    assert abs(obs["variance"]) < 1e-9


def test_atom_dmrg_local_detuning_and_cutoff():
    n = 10
    x = geometry.chain_positions(n, 0.2)
    local = np.random.default_rng(5).normal(scale=0.5, size=n)
    E, _ = _exact(x, G, DELTA, V=geometry.coupling_matrix(x, G, rmax=3), local_delta=local)
    energy, _, _ = dmrg.atom_chain_ground_state(x, G, DELTA, rmax=3, local_delta=local)
    assert abs(energy - E[0]) < 1e-9


@pytest.mark.parametrize("n_pairs,eps,tight", [(6, 0.0, False), (6, 0.3, True), (7, 0.6, False)])
def test_pair_dmrg_matches_ed(n_pairs, eps, tight):
    """Default tolerances give energies to 1e-8; max_S_err=1e-10 gives states to ~1e-7."""
    n = 2 * n_pairs
    E, psi = _exact(geometry.chain_positions(n, eps), G, DELTA)
    kwargs = {"max_S_err": 1e-10} if tight else {}
    energy, mps, obs = dmrg.pair_chain_ground_state(eps, G, DELTA, n_pairs, chi=64,
                                                    correlations=tight, **kwargs)
    atol = 1e-7 if tight else 1e-5
    assert obs["converged"]
    assert abs(energy - E[0]) < 1e-8
    assert np.allclose(obs["occupations"], ob.site_occupations(psi, n), atol=atol)
    assert np.allclose(obs["pair_order"], ob.pair_order(psi, n), atol=atol)
    assert np.allclose(obs["pair_probabilities"], ob.pair_state_probabilities(psi, n), atol=atol)
    assert obs["entropy_cut"] == 2 * (n_pairs // 2)
    assert obs["half_chain_entropy"] == pytest.approx(
        ob.entanglement_entropy(psi, n, obs["entropy_cut"]), abs=atol)
    if tight:
        assert np.allclose(obs["nn"], ob.density_correlations(psi, n), atol=atol)
        assert obs["staggered_order_squared"] == pytest.approx(ob.staggered_order_squared(psi, n), abs=atol)


def test_pair_dmrg_range_cutoff():
    """mmax keeps couplings between pairs at most mmax apart (intra-pair V_s always kept)."""
    n_pairs, eps, mmax = 5, 0.4, 1
    n = 2 * n_pairs
    x = geometry.chain_positions(n, eps)
    V = geometry.coupling_matrix(x, G)
    pair = np.arange(n) // 2
    V[abs(pair[:, None] - pair[None, :]) > mmax] = 0.0
    E, _ = _exact(x, G, DELTA, V=V)
    energy, _, _ = dmrg.pair_chain_ground_state(eps, G, DELTA, n_pairs, mmax=mmax)
    assert abs(energy - E[0]) < 1e-9
