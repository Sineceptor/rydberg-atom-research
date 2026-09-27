"""Ground-state observables of full state vectors (bit i of the index = atom i).

All functions take a normalized state of length 2^n and the atom number n.
Entropies are in bits.
"""

import numpy as np


def _probabilities(state, n):
    state = np.asarray(state)
    if state.shape != (2 ** n,):
        raise ValueError(f"state must have shape ({2 ** n},)")
    return np.abs(state) ** 2


def _diagonal_sum(site_values):
    """Array over basis states of sum_i site_values[i, n_i] (site_values has shape (n, 2))."""
    total = np.zeros(1)
    for values in site_values[::-1]:          # most significant bit first
        total = (total[:, None] + np.asarray(values, dtype=float)[None, :]).ravel()
    return total


def site_occupations(state, n):
    """<n_i> for i = 0..n-1."""
    p = _probabilities(state, n)
    return np.array([p.reshape(2 ** (n - 1 - i), 2, 2 ** i)[:, 1, :].sum() for i in range(n)])


def mean_density(state, n):
    """(1/n) sum_i <n_i>."""
    return float(site_occupations(state, n).mean())


def bond_probabilities(state, n):
    """Joint probabilities of neighbouring atoms (i, i+1), i = 0..n-2.

    Returns a dict with arrays 'p00', 'p01', 'p10', 'p11' of length n-1, where
    'pab' is the probability that n_i = a and n_{i+1} = b.
    """
    p = _probabilities(state, n)
    out = {key: np.empty(n - 1) for key in ("p00", "p01", "p10", "p11")}
    for i in range(n - 1):
        block = p.reshape(2 ** (n - 2 - i), 2, 2, 2 ** i)   # axes: high, n_{i+1}, n_i, low
        for a in (0, 1):
            for b in (0, 1):
                out[f"p{a}{b}"][i] = block[:, b, a, :].sum()
    return out


def staggered_order_squared(state, n):
    """<Q^2> with Q = (1/n) sum_i (-1)^i (2 n_i - 1); both ideal crystals give 1."""
    p = _probabilities(state, n)
    signs = (-1.0) ** np.arange(n)
    Q = _diagonal_sum(np.column_stack([-signs, signs])) / n
    return float(p @ Q ** 2)


def pair_order(state, n):
    """<n_L - n_R> for each strong pair (atoms 2k, 2k+1), k = 0..n//2 - 1."""
    occ = site_occupations(state, n)
    m = 2 * (n // 2)
    return occ[0:m:2] - occ[1:m:2]


def pair_state_probabilities(state, n):
    """Probabilities of |n_L n_R> = 00, 10, 01, 11 for each strong pair; shape (n//2, 4)."""
    p = _probabilities(state, n)
    out = np.empty((n // 2, 4))
    for k in range(n // 2):
        i = 2 * k
        block = p.reshape(2 ** (n - 2 - i), 2, 2, 2 ** i)   # axes: high, n_R, n_L, low
        out[k] = [block[:, 0, 0, :].sum(), block[:, 0, 1, :].sum(),
                  block[:, 1, 0, :].sum(), block[:, 1, 1, :].sum()]
    return out


def entanglement_entropy(state, n, cut):
    """Von Neumann entropy (bits) of the left block of atoms 0..cut-1."""
    if not 0 <= cut <= n:
        raise ValueError("cut must lie in 0..n")
    state = np.asarray(state)
    if state.shape != (2 ** n,):
        raise ValueError(f"state must have shape ({2 ** n},)")
    # Row index = atoms cut..n-1 (high bits), column index = atoms 0..cut-1.
    s = np.linalg.svd(state.reshape(2 ** (n - cut), 2 ** cut), compute_uv=False)
    w = s ** 2
    w = w[w > 1e-15]
    w = w / w.sum()
    return float(-np.sum(w * np.log2(w)))


def density_correlations(state, n):
    """Matrix <n_i n_j> (diagonal = <n_i>)."""
    p = _probabilities(state, n)
    nn = np.diag(site_occupations(state, n))
    for i in range(n):
        for j in range(i + 1, n):
            block = p.reshape(2 ** (n - 1 - j), 2, 2 ** (j - 1 - i), 2, 2 ** i)
            nn[i, j] = nn[j, i] = block[:, 1, :, 1, :].sum()
    return nn


def connected_correlations(state, n):
    """Matrix <n_i n_j> - <n_i><n_j>."""
    nn = density_correlations(state, n)
    occ = np.diag(nn)
    return nn - np.outer(occ, occ)
