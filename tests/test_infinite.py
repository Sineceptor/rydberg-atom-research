"""Infinite pair-site iDMRG: MPO, energy estimators, dressed truncation, finite-bulk agreement."""

import numpy as np
import pytest

from rydberg_chain import ed, geometry

inf = pytest.importorskip("rydberg_chain.infinite")
from tenpy.networks.mps import MPS  # noqa: E402  (TeNPy is present once the import above succeeded)

G, DELTA = 11.39, 2.5


def _classical_energy_per_atom(labels, eps, g, delta, mmax):
    """-delta sum n + V_s n_L n_R + inter-pair couplings for a periodic product pattern."""
    s = 1.0 - eps
    L = len(labels)
    occ = [(int(lab[0]), int(lab[1])) for lab in labels]
    E = 0.0
    for a in range(L):
        nL, nR = occ[a]
        E += -delta * (nL + nR) + g / s ** 6 * nL * nR
        for m in range(1, mmax + 1):
            mL, mR = occ[(a + m) % L]
            for na, fa in ((nL, 0.0), (nR, s)):
                for nb, fb in ((mL, 0.0), (mR, s)):
                    E += na * nb * g / (2 * m + fb - fa) ** 6
    return E / (2 * L)


@pytest.mark.parametrize("cell_pairs,eps,mmax", [(2, 0.3, 6), (3, 0.55, 5), (4, 0.85, 3)])
def test_infinite_mpo_matches_classical_energies(cell_pairs, eps, mmax):
    """Diagonal part of the MPO, including couplings that cross the unit cell."""
    M = inf.infinite_pair_model(eps, G, DELTA, mmax=mmax, cell_pairs=cell_pairs)
    rng = np.random.default_rng(cell_pairs)
    for _ in range(4):
        labels = [inf.PAIR_LABELS[k] for k in rng.integers(0, 4, cell_pairs)]
        psi = MPS.from_product_state(M.lat.mps_sites(), labels, bc="infinite", unit_cell_width=cell_pairs)
        assert inf.energy_per_atom(M, psi) == pytest.approx(
            _classical_energy_per_atom(labels, eps, G, DELTA, mmax), abs=1e-9)


@pytest.mark.parametrize("truncate", [False, True])
def test_transfer_matrix_energy_matches_term_sum(truncate):
    """Off-diagonal terms too: random product superpositions, MPO vs term-by-term sum."""
    eps, mmax = 0.5, 5
    M = inf.infinite_pair_model(eps, G, DELTA, mmax=mmax, cell_pairs=3, truncate=truncate)
    d = 3 if truncate else 4
    rng = np.random.default_rng(3)
    vecs = [v / np.linalg.norm(v) for v in rng.normal(size=(3, d))]
    psi = MPS.from_product_state(M.lat.mps_sites(), vecs, bc="infinite", unit_cell_width=3)
    assert inf.energy_per_atom(M, psi) == pytest.approx(
        inf.energy_from_terms(psi, eps, G, DELTA, mmax=mmax), abs=1e-11)


def _dressed_ed_energy(n, eps):
    """Lowest energy of the full n-atom Hamiltonian, and of it restricted to dressed pair states."""
    H = ed.build_hamiltonian(geometry.coupling_matrix(geometry.chain_positions(n, eps), G), 1.0, DELTA)
    _, W = inf.dressed_pair_basis(eps, G, DELTA)
    U = np.ones((1, 1))
    for _ in range(n // 2):        # pair k contributes (n_L + 2 n_R) 4^k to the ED index
        U = np.kron(U, W)
    return np.linalg.eigvalsh(H.toarray())[0], np.linalg.eigvalsh(U.T @ (H @ U))[0]


def test_dressed_truncation_is_exact_for_one_pair_and_variational():
    E_full, E_tr = _dressed_ed_energy(2, 0.6)
    assert E_tr == pytest.approx(E_full, abs=1e-12)
    for eps in (0.4, 0.5):
        E_full, E_tr = _dressed_ed_energy(8, eps)
        assert -1e-12 < (E_tr - E_full) / 8 < inf.truncation_error_estimate(eps, G, DELTA)


def test_dressed_site_operators_are_projected_pair_operators():
    site = inf.dressed_pair_site(0.6, G, DELTA)
    full = inf.full_pair_site()
    W = site.pair_basis
    for name in ("nL", "nR", "XL", "tau", "P11"):
        assert np.allclose(site.get_op(name).to_ndarray(), W.T @ full.get_op(name).to_ndarray() @ W)
    H1 = inf.pair_hamiltonian(0.6, G, DELTA)
    assert np.allclose(W.T @ H1 @ W, np.diag(np.linalg.eigvalsh(H1)[:3]))


def test_use_truncation_auto():
    assert not inf.use_truncation(0.3, G, DELTA)
    assert inf.use_truncation(0.6, G, DELTA)
    assert inf.use_truncation(0.3, G, DELTA, truncate=True)
    assert not inf.use_truncation(0.8, G, DELTA, truncate=False)


def test_pinned_open_chain_has_no_domain_wall():
    """With the pinned last atom the Z2 crystal keeps one sign of <n_L - n_R>."""
    energy, psi, info = inf.pinned_chain_ground_state(0.2, G, DELTA, 10, chi=32, mmax=4)
    assert info["converged"]
    assert np.all(np.real(psi.expectation_value("tau")) > 0.8)


@pytest.fixture(scope="module")
def ordered_point():
    eps, mmax, chi = 0.2, 4, 32
    e, psi, obs = inf.infinite_ground_state(eps, G, DELTA, chi=chi, mmax=mmax, starts=("L", "R"))
    return eps, mmax, chi, e, psi, obs


def test_idmrg_matches_finite_bulk_energy(ordered_point):
    eps, mmax, chi, e, psi, obs = ordered_point
    assert obs["converged"] and obs["valid"] and not obs["truncated"]
    e_fin, info = inf.finite_bulk_energy(eps, G, DELTA, n_pairs=(10, 20), chi=chi, mmax=mmax)
    assert all(v["converged"] for v in info.values())
    assert e == pytest.approx(e_fin, abs=1e-8)
    assert e == pytest.approx(inf.energy_from_terms(psi, eps, G, DELTA, mmax=mmax), abs=1e-10)
    assert obs["candidates"]["L"]["growth_energy_per_atom"] == pytest.approx(e, abs=1e-8)


def test_idmrg_order_is_translation_and_mirror_consistent(ordered_point):
    eps, mmax, chi, e, psi, obs = ordered_point
    assert obs["order_spread"] < 1e-8
    assert obs["order"] > 0.9
    assert obs["mirror_check"]["dE"] < 1e-10 and obs["mirror_check"]["dtau"] < 1e-6
    assert 0.4 < obs["density"] < 0.5
    assert obs["xi"] > 0 and min(obs["entropy_weak"]) >= 0


def test_truncated_idmrg_matches_finite_chains():
    """Large V_s (eps = 0.7, disordered): dressed iDMRG against pinned finite chains."""
    eps, mmax, chi = 0.7, 3, 24
    e, psi, obs = inf.infinite_ground_state(eps, G, DELTA, chi=chi, mmax=mmax, starts=("L",))
    assert obs["truncated"] and obs["converged"]
    e_tr, _ = inf.finite_bulk_energy(eps, G, DELTA, n_pairs=(16, 32), chi=chi, mmax=mmax, truncate=True)
    assert e == pytest.approx(e_tr, abs=1e-8)
    E_tr, _, _ = inf.pinned_chain_ground_state(eps, G, DELTA, 8, chi=chi, mmax=mmax, truncate=True)
    E_full, _, info = inf.pinned_chain_ground_state(eps, G, DELTA, 8, chi=chi, mmax=mmax)
    assert info["converged"] and abs(E_tr - E_full) / 16 < 1e-10
    assert obs["order"] < 1e-6
    assert sum(obs["pair_probabilities"][0]) == pytest.approx(1.0, abs=1e-10)


def test_bulk_energy_fit():
    f = {N: -1.2 * N + 3.0 + 0.5 / N + 7.0 / N ** 2 for N in (100, 200, 400, 800)}
    e, b, c, d = inf.bulk_energy_fit(f)
    assert (e, b, c, d) == pytest.approx((-1.2, 3.0, 0.5, 7.0), rel=1e-6)
    two = {N: f[N] for N in (100, 200)}
    assert inf.bulk_energy_fit(two)[0] == pytest.approx((f[200] - f[100]) / 100, abs=1e-13)
