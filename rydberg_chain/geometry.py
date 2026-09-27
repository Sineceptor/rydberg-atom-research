"""Atom positions and van der Waals interaction matrices."""

import numpy as np


def chain_positions(n_atoms, eps=0.0):
    """Positions of an open dimerized chain with mean spacing 1.

    The gap after atom i is 1 - eps for even i (strong bond inside a pair) and
    1 + eps for odd i (weak bond between pairs).  eps = 0 is the uniform chain.
    """
    if n_atoms < 1:
        raise ValueError("n_atoms must be positive")
    if not -1.0 < eps < 1.0:
        raise ValueError("eps must lie in (-1, 1) so that all gaps are positive")
    gaps = np.where(np.arange(n_atoms - 1) % 2 == 0, 1.0 - eps, 1.0 + eps)
    return np.concatenate([[0.0], np.cumsum(gaps)])


def strong_weak_couplings(eps, g):
    """Nearest-neighbour energies (V_s, V_w) = (g/(1-eps)^6, g/(1+eps)^6)."""
    return g / (1.0 - eps) ** 6, g / (1.0 + eps) ** 6


def coupling_matrix(x, g, rmax=None):
    """Symmetric V_ij = g/|x_i - x_j|^6 with zero diagonal.

    x may be a 1D array of positions on a line or an (n, d) array of points.
    If rmax is given, couplings with |i - j| > rmax (index distance) are set
    to zero; this is the interaction cutoff used by truncated models.
    """
    x = np.asarray(x, dtype=float)
    points = x[:, None] if x.ndim == 1 else x
    n = len(points)
    diff = points[:, None, :] - points[None, :, :]
    dist = np.sqrt(np.sum(diff ** 2, axis=-1))
    off = ~np.eye(n, dtype=bool)
    if np.any(dist[off] == 0.0):
        raise ValueError("two atoms occupy the same position")
    V = np.zeros((n, n))
    V[off] = g / dist[off] ** 6
    if rmax is not None:
        idx = np.arange(n)
        V[np.abs(idx[:, None] - idx[None, :]) > rmax] = 0.0
    return V


def reflect(x):
    """Mirror a 1D configuration and relabel atoms left to right.

    Atom i of the result sits where atom n-1-i was, measured from the other end,
    so the returned positions are again ascending and start at zero.
    """
    x = np.asarray(x, dtype=float)
    return x.max() - x[::-1]
