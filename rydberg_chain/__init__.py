"""Exact diagonalization, DMRG and effective models for dimerized Rydberg chains.

Conventions (shared by every module):

    hbar = 1, |0> = ground, |1> = Rydberg,
    H = (Omega/2) sum_i X_i - Delta sum_i n_i + sum_{i<j} C6/|x_i - x_j|^6 n_i n_j.

Lengths are in units of the mean spacing a and energies in units of Omega, so
g = C6/(Omega a^6) = (R_b/a)^6.  In a dimerized chain the gap after atom i is
1 - eps for even i (strong bond, inside a pair) and 1 + eps for odd i (weak
bond, between pairs).  In the exact-diagonalization basis, bit i of the basis
index is the occupation of atom i.

Modules
-------
geometry     positions and interaction matrices
ed           vectorized sparse Hamiltonian and lowest eigenpairs
observables  occupations, bond probabilities, order parameters, entropy
effective    one-excitation-per-pair transverse-field Ising estimates
dmrg         finite DMRG (TeNPy) on atoms or on 4-state pair sites
"""

from . import geometry, ed, observables, effective

# ``dmrg`` needs TeNPy and is imported on request: ``from rydberg_chain import dmrg``.
__all__ = ["geometry", "ed", "observables", "effective"]
__version__ = "0.1.0"
