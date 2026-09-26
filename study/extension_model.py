"""Ground states of small open Rydberg chains, with controllable interaction range.

Units: hbar = Omega = nearest-neighbour spacing = 1.
Atom 0 is the rightmost bit, as in the unchanged original rydberg_model.py.
No periodic boundary and no projection onto a blockade-only basis are used.
"""

from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from rydberg_model import chain, interactions, occupations, hamiltonian


def chain_couplings(number, strength, cutoff="all"):
    """Keep separation <= cutoff lattice spacings; 'all' retains every pair.

    'nearest1' keeps g, 'nearest2' also keeps g/64; remaining pairs are
    g/729, g/4096, ... . The laser and the number of atoms never change.
    """
    if cutoff not in ("nearest1", "nearest2", "all"):
        raise ValueError("cutoff must be nearest1, nearest2, or all")
    coupling = interactions(chain(number), strength)
    if cutoff != "all":
        distance_limit = 1 if cutoff == "nearest1" else 2
        for i in range(number):
            for j in range(number):
                if abs(i-j) > distance_limit:
                    coupling[i, j] = 0.0
    return coupling


def discarded_tail_bound(number, strength, cutoff):
    """Exact operator norm of discarded positive interactions in this chain.

    The all-excited bit string contains every discarded pair. There are N-r
    pairs separated by r spacings, each with interaction g/r**6.
    This bounds each energy shift and the gap change, not the density error.
    """
    if cutoff == "all":
        return 0.0
    if cutoff not in ("nearest1", "nearest2") or strength < 0:
        raise ValueError("This bound requires a recognised cutoff and nonnegative g")
    keep = 1 if cutoff == "nearest1" else 2
    return float(strength * sum((number-r)/r**6 for r in range(keep+1, number)))


@lru_cache(maxsize=None)
def basis_observables(number):
    """Counts and two staggered observables, keeping the earlier one for audit.

    M = (2/N) sum (-1)^i n_i. For odd N, this has an offset because
    there is one more even-index site. Its maximum magnitude exceeds 1.
    S = (1/N) sum (-1)^i (2 n_i - 1) removes that offset.
    Both perfect alternating strings have S**2 = 1 for either parity.
    Neither diagnostic on a small chain by itself establishes a bulk phase.
    """
    bits = occupations(number)
    count = bits.sum(axis=1)
    imbalance = (2.0 / number) * (bits @ (-1.0)**np.arange(number))
    centered = imbalance - np.sum((-1.0)**np.arange(number))/number
    return bits, count, imbalance**2, centered


def lowest_two(matrix, method="auto", seed=230923):
    """Find the two lowest energies, then explicitly check both eigenvectors.

    A deterministic random starting vector has no deliberately imposed
    reflection symmetry. A symmetric starting vector could miss one sector.
    """
    if method == "dense" or (method == "auto" and matrix.shape[0] <= 64):
        energies, vectors = eigh(matrix.toarray(), subset_by_index=[0, 1])
    else:
        start = np.random.default_rng(seed).normal(size=matrix.shape[0])
        energies, vectors = eigsh(matrix, k=2, which="SA", tol=1e-11,
                                 ncv=min(40, matrix.shape[0]-1),
                                 maxiter=50000, v0=start)
        order = np.argsort(energies)
        energies, vectors = energies[order], vectors[:, order]
    return energies, vectors


def ground_point(number, strength, delta, cutoff="all", method="auto", seed=230923):
    """Return observables and per-site occupations of the ground state.

    Density = mean occupation over all atoms.
    Centre density excludes just the first and last atom, not all edge effects.
    Earlier alternating_order_squared is retained for reproducibility and
    has an odd-N offset. Use staggered_spin_squared for comparisons
    across odd/even N: both ideal alternating strings then give 1.
    """
    coupling = chain_couplings(number, strength, cutoff)
    matrix = hamiltonian(coupling, omega=1.0, delta=delta)
    energies, vectors = lowest_two(matrix, method=method, seed=seed)
    state = vectors[:, 0]
    weights = abs(state)**2
    bits, count, alternating_squared, centered = basis_observables(number)
    site_density = weights @ bits
    residuals = [np.linalg.norm(matrix @ vectors[:, i] - energies[i]*vectors[:, i])
                 for i in range(2)]
    checks = {
        "residual": float(max(residuals)),
        "relative_residual": float(max(residuals)/(1+max(abs(energies)))),
        "orthogonality_error": float(abs(np.vdot(vectors[:, 0], vectors[:, 1]))),
        "normalization_error": float(max(abs(np.linalg.norm(vectors[:, i])-1)
                                         for i in range(2))),
    }
    row = {
        "N": number, "g": float(strength), "delta": float(delta), "cutoff": cutoff,
        "density": float(weights @ count / number),
        "alternating_order_squared": float(weights @ alternating_squared),
        "staggered_spin_squared": float(weights @ centered**2),
        "staggered_spin_mean": float(weights @ centered),
        "gap": float(energies[1]-energies[0]), "ground_energy": float(energies[0]),
        "centre_density": float(site_density[1:-1].mean()),
        "edge_density": float(site_density[[0, -1]].mean()),
        **checks,
    }
    return row, site_density, matrix
