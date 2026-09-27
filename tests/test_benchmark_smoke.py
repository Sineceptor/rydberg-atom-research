"""Fast checks of validation/benchmark_uniform.py (exact references and fit helpers)."""

import importlib.util
import json
import os

import numpy as np
import pytest

from rydberg_chain import ed, geometry

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    "benchmark_uniform", os.path.join(ROOT, "validation", "benchmark_uniform.py"))
bench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bench)


def _tfim_dense(n, J, h):
    """Dense ED of -J sum Z_i Z_{i+1} - h sum X_i: ground energy and S(l) in nats."""
    X, Z, I = np.array([[0.0, 1.0], [1.0, 0.0]]), np.diag([1.0, -1.0]), np.eye(2)

    def op(o, i):
        m = np.ones((1, 1))
        for k in range(n):
            m = np.kron(m, o if k == i else I)
        return m

    H = -J * sum(op(Z, i) @ op(Z, i + 1) for i in range(n - 1)) - h * sum(op(X, i) for i in range(n))
    w, v = np.linalg.eigh(H)
    S = []
    for ell in range(1, n):
        p = np.linalg.svd(v[:, 0].reshape(2 ** ell, -1), compute_uv=False) ** 2
        p = p[p > 1e-15]
        S.append(-np.sum(p * np.log(p)))
    return w[0], np.array(S)


@pytest.mark.parametrize("n,h", [(6, 1.0), (8, 1.0), (8, 0.6)])
def test_tfim_free_fermions_match_dense_ed(n, h):
    E, S = bench.tfim_exact(n, 1.0, h)
    E_ed, S_ed = _tfim_dense(n, 1.0, h)
    assert E == pytest.approx(E_ed, abs=1e-12)
    assert np.allclose(S, S_ed, atol=1e-12)


def test_tfim_dmrg_matches_free_fermions():
    pytest.importorskip("tenpy")
    E, S, diag = bench.tfim_dmrg(16, chi=32)
    E_ex, S_ex = bench.tfim_exact(16)
    assert diag["converged"]
    assert E == pytest.approx(E_ex, abs=1e-10)
    assert np.allclose(S, S_ex, atol=1e-7)


def test_calabrese_cardy_fit_recovers_synthetic_c():
    n = 101
    ell = np.arange(1, n)
    chord = np.log(2 * n / np.pi * np.sin(np.pi * ell / n))
    S = (0.5 / 6) * chord + 0.3 + 0.05 * (-1.0) ** ell
    for parity in (0, 1, "smooth"):
        assert bench.calabrese_cardy_fit(S, n, parity=parity)[0] == pytest.approx(0.5, abs=2e-3)


def test_crossing_and_extrapolation():
    d = np.linspace(0.6, 0.7, 6)
    assert bench.crossing(d, 1 + 2 * (d - 0.65), 1 + 5 * (d - 0.65)) == pytest.approx(0.65, abs=1e-12)
    pairs = [(n, 2 * n, 0.655 + 3.0 / (2 ** 0.5 * n) ** 2) for n in (20, 40, 80)]
    assert bench.extrapolate_crossings(pairs)["p=2"]["Delta_c"] == pytest.approx(0.655, abs=1e-12)


@pytest.mark.parametrize("n,delta", [(9, 0.655), (10, 1.5)])
def test_blockade_dmrg_matches_constrained_ed(n, delta):
    pytest.importorskip("tenpy")
    r = bench.blockade_scan(n, [delta])[0]
    assert r["converged"]
    assert r["E"] == pytest.approx(bench.blockade_chain_exact(n, delta), abs=1e-9)


def test_large_v1_approaches_blockade_limit():
    n, delta = 8, 0.655
    E_c = bench.blockade_chain_exact(n, delta)
    x = geometry.chain_positions(n)
    gaps = []
    for V1 in (50.0, 100.0, 200.0):
        H = ed.build_hamiltonian(geometry.coupling_matrix(x, V1, rmax=1), 1.0, delta)
        gaps.append(E_c - ed.lowest_states(H, k=1)[0][0])
    # virtual double occupation lowers the energy by O(Omega^2/V1)
    assert all(g > 0 for g in gaps)
    assert gaps[0] / gaps[1] == pytest.approx(2.0, rel=0.05)
    assert gaps[1] / gaps[2] == pytest.approx(2.0, rel=0.05)


def test_saved_results_agree_with_literature():
    path = bench.RESULTS
    if not os.path.exists(path):
        pytest.skip("run validation/benchmark_uniform.py first")
    with open(path) as f:
        res = json.load(f)
    if res.get("quick") or "--quick" in res.get("argv", []):
        pytest.skip("saved results come from a --quick run")
    rows = res["tfim"]["rows"]
    for row in rows:
        assert abs(row["dE"]) < 1e-8 and row["max_abs_dS"] < 1e-6
    assert rows[-1]["c_dmrg"] == pytest.approx(0.5, abs=0.02)
    blk = res["blockade"]
    assert all(abs(c["dE"]) < 1e-9 for c in blk["ed_checks"])
    assert all(r["converged"] for scan in blk["scans"].values() for r in scan)
    assert blk["extrapolations"]["estimate"]["Delta_c"] == pytest.approx(0.655, abs=0.005)
    largest = max(blk["central_charge"], key=int)
    assert blk["central_charge"][largest]["c_smooth"] == pytest.approx(0.5, abs=0.05)
    assert blk["c_extrapolated"]["c_smooth"]["c_inf"] == pytest.approx(0.5, abs=0.02)
