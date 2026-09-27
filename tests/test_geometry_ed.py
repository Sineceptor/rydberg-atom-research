"""Geometry, exact-diagonalization builder and state observables."""

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

from rydberg_chain import ed, geometry, observables as ob

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "study"))
sys.dont_write_bytecode = True     # study/ is a checksummed record: never write caches there
import rydberg_model as study  # noqa: E402  (historical reference implementation)
sys.dont_write_bytecode = False

E0_SAVED = -18.399305357538882     # N = 12, g = 8, Delta = 3, uniform, all pairs


def _saved_row(n, g, delta):
    with open(ROOT / "study" / "results" / "size_comparison_V8.csv") as f:
        for row in csv.DictReader(f):
            if int(row["N"]) == n and float(row["V_over_omega"]) == g and float(row["delta_over_omega"]) == delta:
                return {k: float(v) for k, v in row.items()}
    raise LookupError("saved row not found")


# ------------------------------------------------------------------ geometry

def test_chain_positions():
    assert np.array_equal(geometry.chain_positions(5), np.arange(5.0))
    x = geometry.chain_positions(7, 0.3)
    assert np.allclose(np.diff(x), [0.7, 1.3, 0.7, 1.3, 0.7, 1.3])
    assert x[0] == 0.0 and np.isclose(x[-1], 6.0)          # mean spacing 1
    assert geometry.chain_positions(1, 0.4).tolist() == [0.0]


def test_coupling_matrix():
    x = geometry.chain_positions(6, 0.2)
    V = geometry.coupling_matrix(x, 5.0)
    assert np.array_equal(V, V.T) and np.all(np.diag(V) == 0)
    assert np.isclose(V[0, 1], 5.0 / 0.8 ** 6) and np.isclose(V[1, 2], 5.0 / 1.2 ** 6)
    assert np.isclose(V[0, 3], 5.0 / (x[3] - x[0]) ** 6)
    Vc = geometry.coupling_matrix(x, 5.0, rmax=2)
    i, j = np.indices(V.shape)
    assert np.array_equal(Vc, np.where(abs(i - j) <= 2, V, 0.0))
    points = np.column_stack([x, np.zeros_like(x)])
    assert np.allclose(geometry.coupling_matrix(points, 5.0), V, rtol=1e-14)
    assert geometry.strong_weak_couplings(0.2, 5.0) == pytest.approx((V[0, 1], V[1, 2]))
    with pytest.raises(ValueError):
        geometry.coupling_matrix([0.0, 1.0, 1.0], 1.0)


# ----------------------------------------------------------- ED Hamiltonian

@pytest.mark.parametrize("seed", range(6))
def test_builder_equals_study_hamiltonian(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 10))
    x = np.cumsum(rng.uniform(0.5, 1.5, n))
    g, omega, delta = rng.uniform(0.5, 20), rng.uniform(0.1, 3), rng.uniform(-3, 5)
    rmax = None if seed % 2 else int(rng.integers(1, max(2, n)))
    V = geometry.coupling_matrix(x, g, rmax)
    V_study = study.interactions(np.column_stack([x, np.zeros(n)]), g)
    if rmax is None:
        assert np.allclose(V, V_study, rtol=1e-14, atol=0)     # equal up to pow() rounding
    new = ed.build_hamiltonian(V, omega, delta)
    old = study.hamiltonian(V, omega, delta)
    assert new.shape == old.shape and new.has_sorted_indices
    assert np.array_equal(new.toarray(), old.toarray())      # bit-for-bit


def test_two_atom_matrix():
    H = ed.build_hamiltonian(geometry.coupling_matrix([0.0, 1.0], 8.0), omega=0.6, delta=1.3).toarray()
    w, d, V = 0.3, 1.3, 8.0
    expected = np.array([[0, w, w, 0],          # |00>
                         [w, -d, 0, w],         # atom 0 excited (bit 0)
                         [w, 0, -d, w],         # atom 1 excited (bit 1)
                         [0, w, w, V - 2 * d]])
    assert np.array_equal(H, expected)


def test_local_detuning():
    rng = np.random.default_rng(3)
    n = 6
    V = geometry.coupling_matrix(geometry.chain_positions(n, 0.3), 4.0)
    local = rng.normal(size=n)
    H = ed.build_hamiltonian(V, 1.0, 2.0, local_delta=local)
    bits = ed.basis_occupations(n).astype(float)
    classical = -bits @ (2.0 + local) + np.einsum("si,ij,sj->s", bits, np.triu(V, 1), bits)
    assert np.allclose(H.diagonal(), classical, atol=1e-12)
    difference = (H - ed.build_hamiltonian(V, 1.0, 2.0)).toarray()
    assert np.allclose(difference, np.diag(-bits @ local), atol=1e-12)


@pytest.mark.parametrize("delta,omega", [(-1.5, 1.0), (0.0, 1.0), (0.7, 1.0), (2.0, 2.5)])
def test_noninteracting_density(delta, omega):
    n = 6
    H = ed.build_hamiltonian(np.zeros((n, n)), omega, delta)
    E, v, res = ed.lowest_states(H, k=1)
    expected = (1 + delta / np.sqrt(delta ** 2 + omega ** 2)) / 2
    assert np.allclose(ob.site_occupations(v[:, 0], n), expected, atol=1e-12)
    assert E[0] == pytest.approx(n * (-delta - np.sqrt(delta ** 2 + omega ** 2)) / 2, abs=1e-12)
    assert res[0] < 1e-12


def test_eigensolvers_agree():
    H = ed.chain_hamiltonian(10, 0.3, 11.39, 2.5)
    E_d, v_d, r_d = ed.lowest_states(H, k=3, method="dense")
    assert np.all(np.diff(E_d) >= 0) and r_d.max() < 1e-12
    for method in ("eigsh", "lobpcg"):
        E, v, r = ed.lowest_states(H, k=3, method=method)
        assert np.allclose(E, E_d, atol=1e-11) and np.all(np.diff(E) >= 0)
        assert r.max() < 1e-10
        assert abs(abs(v_d[:, 0] @ v[:, 0]) - 1) < 1e-10
        assert np.array_equal(E, ed.lowest_states(H, k=3, method=method)[0])   # seeded
    assert np.array_equal(ed.lowest_states(H, k=3)[0], E_d)                  # auto -> dense at 1024


def test_saved_uniform_reference():
    saved = _saved_row(12, 8.0, 3.0)
    H = ed.chain_hamiltonian(12, 0.0, 8.0, 3.0)
    assert abs(ed.lowest_states(H, k=1, method="lobpcg")[0][0] - E0_SAVED) < 1e-9
    E, v, res = ed.lowest_states(H, k=2)
    psi = v[:, 0]
    assert abs(E[0] - E0_SAVED) < 1e-9
    assert abs(E[0] - saved["energy_over_omega"]) < 1e-9
    assert abs((E[1] - E[0]) - saved["gap_over_omega"]) < 1e-8
    assert abs(ob.mean_density(psi, 12) - saved["density"]) < 1e-9
    assert abs(ob.staggered_order_squared(psi, 12) - saved["staggered_spin_squared"]) < 1e-9
    assert abs(ob.entanglement_entropy(psi, 12, 6) - saved["entropy_bits"]) < 1e-8
    assert res.max() < 1e-10


@pytest.mark.parametrize("n", [7, 9])
def test_reflection_invariance_odd_dimerized(n):
    x = geometry.chain_positions(n, 0.3)        # starts with a strong bond, ends with a weak one
    xr = geometry.reflect(x)                    # starts weak, ends strong: a different chain
    assert not np.allclose(np.diff(x), np.diff(xr))
    spectra = []
    states = []
    for pos in (x, xr, -x):
        H = ed.build_hamiltonian(geometry.coupling_matrix(pos, 11.39), 1.0, 2.5)
        spectra.append(np.linalg.eigvalsh(H.toarray()))
        states.append(ed.lowest_states(H, k=1)[1][:, 0])
    assert np.allclose(spectra[0], spectra[1], atol=1e-10)
    assert np.allclose(spectra[0], spectra[2], atol=1e-10)
    occ, occ_r = ob.site_occupations(states[0], n), ob.site_occupations(states[1], n)
    assert np.allclose(occ, occ_r[::-1], atol=1e-10)


# -------------------------------------------------------------- observables

def test_observables_on_crystal():
    n = 6
    psi = np.zeros(2 ** n)
    psi[int("010101", 2)] = 1.0                 # atoms 0, 2, 4 excited (bit i = atom i)
    assert np.array_equal(ob.site_occupations(psi, n), [1, 0, 1, 0, 1, 0])
    assert ob.staggered_order_squared(psi, n) == pytest.approx(1.0)
    assert np.allclose(ob.pair_order(psi, n), 1.0)
    bonds = ob.bond_probabilities(psi, n)
    assert np.allclose(bonds["p10"], [1, 0, 1, 0, 1]) and np.allclose(bonds["p01"], [0, 1, 0, 1, 0])
    assert np.allclose(bonds["p00"] + bonds["p11"], 0)
    assert np.allclose(ob.pair_state_probabilities(psi, n), [[0, 1, 0, 0]] * 3)
    assert ob.entanglement_entropy(psi, n, 3) == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(ob.connected_correlations(psi, n), 0)


def test_bell_entropy_and_correlations():
    n = 4
    psi = np.zeros(2 ** n)
    psi[0b0001] = psi[0b0010] = 1 / np.sqrt(2)  # pair (0, 1) shares one excitation
    assert ob.entanglement_entropy(psi, n, 1) == pytest.approx(1.0)
    assert ob.entanglement_entropy(psi, n, 2) == pytest.approx(0.0, abs=1e-12)
    assert ob.pair_order(psi, n) == pytest.approx([0.0, 0.0])
    C = ob.connected_correlations(psi, n)
    assert C[0, 1] == pytest.approx(-0.25) and C[0, 0] == pytest.approx(0.25)


def test_observables_against_brute_force():
    n = 7
    rng = np.random.default_rng(11)
    psi = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    psi /= np.linalg.norm(psi)
    p = abs(psi) ** 2
    bits = ed.basis_occupations(n).astype(float)
    assert np.allclose(ob.site_occupations(psi, n), p @ bits)
    assert np.allclose(ob.density_correlations(psi, n), bits.T @ (p[:, None] * bits))
    Q = ((2 * bits - 1) @ (-1.0) ** np.arange(n)) / n
    assert ob.staggered_order_squared(psi, n) == pytest.approx(p @ Q ** 2)
    bonds = ob.bond_probabilities(psi, n)
    for i in range(n - 1):
        for a in (0, 1):
            for b in (0, 1):
                mask = (bits[:, i] == a) & (bits[:, i + 1] == b)
                assert bonds[f"p{a}{b}"][i] == pytest.approx(p[mask].sum())
    probs = ob.pair_state_probabilities(psi, n)
    assert np.allclose(probs.sum(axis=1), 1)
    assert np.allclose(probs[:, 1] + probs[:, 3], (p @ bits)[0:6:2])
    for cut in range(n + 1):
        keep = list(range(cut))
        expected = study.entropy(psi, keep) if 0 < cut < n else 0.0
        assert ob.entanglement_entropy(psi, n, cut) == pytest.approx(expected, abs=1e-10)
