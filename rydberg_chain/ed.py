"""Exact diagonalization in the full 2^N occupation basis.

Bit i of a basis index is the occupation of atom i (atom 0 is the least
significant bit), the same convention as ``study/rydberg_model.py``.
"""

import warnings

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.linalg import LinearOperator, eigsh, lobpcg

from .geometry import chain_positions, coupling_matrix

DENSE_LIMIT = 1024          # dense eigh up to this dimension (eigsh is ~100x faster at 4096)
_BLOCK = 1 << 18            # basis states processed per vectorized block


def basis_occupations(n, start=0, stop=None, dtype=np.uint8):
    """Occupations of basis states start..stop-1: array (stop-start, n), column i = atom i."""
    stop = 2 ** n if stop is None else stop
    idx = np.arange(start, stop, dtype=np.int64)
    return ((idx[:, None] >> np.arange(n, dtype=np.int64)) & 1).astype(dtype)


def _site_detunings(n, delta, local_delta):
    det = np.full(n, float(delta))
    if local_delta is not None:
        local_delta = np.asarray(local_delta, dtype=float)
        if local_delta.shape != (n,):
            raise ValueError(f"local_delta must have shape ({n},)")
        det = det + local_delta
    return det


def diagonal_energies(V, delta=0.0, local_delta=None):
    """Classical energies -sum_i Delta_i n_i + sum_{i<j} V_ij n_i n_j of every basis state.

    Delta_i = delta + local_delta[i].  Without local_delta the terms are added in
    the same order as ``study/rydberg_model.hamiltonian``, so the diagonal agrees
    with it bit for bit.
    """
    V = np.asarray(V, dtype=float)
    n = len(V)
    dim = 2 ** n
    pairs = [(i, j, V[i, j]) for i in range(n) for j in range(i + 1, n) if V[i, j] != 0.0]
    det = None if local_delta is None else _site_detunings(n, delta, local_delta)
    energy = np.empty(dim)
    for start in range(0, dim, _BLOCK):
        stop = min(dim, start + _BLOCK)
        bits = basis_occupations(n, start, stop)
        if det is None:
            e = -delta * bits.sum(axis=1).astype(float)
        else:
            e = -(bits @ det)
        both = np.empty(stop - start, dtype=np.uint8)
        for i, j, v in pairs:
            np.bitwise_and(bits[:, i], bits[:, j], out=both)
            e += v * both
        energy[start:stop] = e
    return energy


def build_hamiltonian(V, omega=1.0, delta=0.0, local_delta=None):
    """Sparse Hamiltonian (Omega/2) sum X_i - sum Delta_i n_i + sum_{i<j} V_ij n_i n_j.

    V is the symmetric coupling matrix (only i < j is read); local_delta is an
    optional per-site detuning added to delta, Delta_i = delta + local_delta[i].
    Returns a real symmetric CSR matrix with sorted indices and N + 1 stored
    entries per row (the diagonal is stored even where it is zero).
    """
    V = np.asarray(V, dtype=float)
    n = len(V)
    dim = 2 ** n
    diag = diagonal_energies(V, delta, local_delta)
    per_row = n + 1
    index_type = np.int32 if dim * per_row < 2 ** 31 else np.int64
    indices = np.empty(dim * per_row, dtype=index_type)
    data = np.empty(dim * per_row)
    flips = np.left_shift(1, np.arange(n, dtype=np.int64))
    for start in range(0, dim, _BLOCK):
        stop = min(dim, start + _BLOCK)
        rows = np.arange(start, stop, dtype=np.int64)
        cols = np.empty((stop - start, per_row), dtype=np.int64)
        cols[:, 0] = rows
        cols[:, 1:] = rows[:, None] ^ flips
        cols.sort(axis=1)
        vals = np.where(cols == rows[:, None], diag[start:stop, None], 0.5 * omega)
        indices[start * per_row:stop * per_row] = cols.ravel()
        data[start * per_row:stop * per_row] = vals.ravel()
    indptr = np.arange(0, dim * per_row + 1, per_row, dtype=index_type)
    H = csr_matrix((data, indices, indptr), shape=(dim, dim))
    H.has_sorted_indices = True
    return H


def chain_hamiltonian(n_atoms, eps, g, delta, omega=1.0, rmax=None, local_delta=None):
    """Hamiltonian of an open dimerized chain; see geometry.chain_positions."""
    V = coupling_matrix(chain_positions(n_atoms, eps), g, rmax)
    return build_hamiltonian(V, omega, delta, local_delta)


def lowest_states(H, k=2, seed=2026, tol=1e-12, dense_limit=DENSE_LIMIT, method="auto"):
    """Lowest k eigenpairs of a real symmetric H.

    method
      'auto'    dense ``scipy.linalg.eigh`` when dim <= dense_limit, else 'eigsh';
      'dense'   dense eigh;
      'eigsh'   ARPACK ``eigsh(which='SA', tol=tol)`` from a seeded Gaussian vector;
      'lobpcg'  block LOBPCG (block size k + 2, seeded) with the Jacobi preconditioner
                (diag(H) - sigma)^-1; tol is the residual-norm target.  Typically 4-5x
                faster than eigsh for Rydberg chains, whose spectrum is wide but whose
                high-energy states are nearly diagonal.  Falls back to eigsh if the
                residual exceeds 1e-8 max(1, |E0|).

    Returns (energies, vectors, residuals): energies ascending (shape (k,)),
    normalized eigenvectors as columns (dim, k), and residual norms
    ||H v - E v|| for each pair (shape (k,)).
    """
    dim = H.shape[0]
    if not 1 <= k <= dim:
        raise ValueError("need 1 <= k <= dim")
    if method == "auto":
        method = "dense" if dim <= dense_limit else "eigsh"
    if method not in ("dense", "eigsh", "lobpcg"):
        raise ValueError(f"unknown method {method!r}")
    if method == "lobpcg" and 5 * (k + 2) >= dim:
        method = "dense"
    if method == "eigsh" and k >= dim - 1:
        method = "dense"
    if method == "dense":
        dense = H.toarray() if issparse(H) else np.asarray(H)
        energies, vectors = eigh(dense, subset_by_index=[0, k - 1])
    elif method == "lobpcg":
        energies, vectors = _lobpcg(H, k, seed, tol)
        residuals = np.linalg.norm(H @ vectors - vectors * energies, axis=0)
        if not np.all(residuals <= 1e-8 * max(1.0, abs(energies[0]))):
            return lowest_states(H, k, seed, tol, dense_limit, method="eigsh")
        return energies, vectors, residuals
    else:
        v0 = np.random.default_rng(seed).normal(size=dim)
        energies, vectors = eigsh(H, k=k, which="SA", tol=tol, v0=v0)
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    residuals = np.linalg.norm(H @ vectors - vectors * energies, axis=0)
    return energies, vectors, residuals


def _lobpcg(H, k, seed, tol, maxiter=1000):
    diag = H.diagonal()
    inverse = 1.0 / (diag - (diag.min() - 1.0))

    def precondition(R):
        return R * (inverse if R.ndim == 1 else inverse[:, None])

    M = LinearOperator(H.shape, matvec=precondition, matmat=precondition, dtype=float)
    X = np.random.default_rng(seed).normal(size=(H.shape[0], k + 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)      # non-convergence is checked by the caller
        energies, vectors = lobpcg(H, X, M=M, largest=False, tol=tol, maxiter=maxiter)
    order = np.argsort(energies)[:k]
    vectors = vectors[:, order]
    return energies[order], vectors / np.linalg.norm(vectors, axis=0)
