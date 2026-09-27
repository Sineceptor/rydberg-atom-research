# Notation and conventions

One set of symbols for the code, the data files, the figures and the manuscript. Where older files use a different name, the mapping is given below.

## Units and Hamiltonian

ħ = 1. Energies are in units of the Rabi frequency Ω, lengths in units of the mean spacing a, and times as Ωt.

$$H=\frac{\Omega}{2}\sum_i X_i-\Delta\sum_i n_i+\sum_{i<j}V(|x_i-x_j|)\,n_in_j,\qquad V(r)=\frac{C_6}{r^6},\qquad n_i=|1_i\rangle\langle1_i|.$$

|0⟩ is the ground state and |1⟩ the Rydberg state. Each pair i<j is counted once. The sign of the drive term is a gauge choice (∏ᵢZᵢ maps X → −X), so only |Ω| matters.

## Symbols

| Symbol | Meaning | Code / CSV name |
|---|---|---|
| Ω | Rabi frequency (unit of energy) | `omega` (default 1) |
| Δ | laser detuning; Δ > 0 favours excitation | `delta`, `delta_over_omega` |
| δ = Δ/Ω | dimensionless detuning (equal to Δ in code units) | `delta_over_omega` |
| C₆ | van der Waals coefficient | — |
| R_b | blockade radius, defined by C₆/R_b⁶ = Ω | — |
| g = C₆/(Ωa⁶) = (R_b/a)⁶ | coupling constant; for a uniform chain g = V₀/Ω | `g`, `V_over_omega` (chains) |
| V₀ = C₆/a⁶ | nearest-neighbour interaction of the uniform chain | `strength` in `study/` |
| V/Ω | interaction of an isolated pair of atoms | `V_over_omega` (two-atom files) |
| V_ij | pair interaction matrix | `coupling_matrix`, `interactions` |
| N | number of atoms | `N`, `n_atoms` |
| ε | dimerization, −1 < ε < 1 | `eps` |
| s = a(1−ε), w = a(1+ε) | strong (intra-pair) and weak (inter-pair) gaps | `s` |
| V_s = g/(1−ε)⁶, V_w = g/(1+ε)⁶ | strong- and weak-bond nearest-neighbour energies | `V_s`, `V_w` |
| n̄ = N⁻¹Σᵢ⟨nᵢ⟩ | mean excitation density | `density` |
| Q = N⁻¹Σᵢ(−1)ⁱ(2nᵢ−1) | staggered order; Q² = 1 for either perfect alternating pattern at any N | `staggered_spin_squared` (= ⟨Q²⟩), `staggered_order_squared()` |
| M = (2/N)Σᵢ(−1)ⁱnᵢ | legacy imbalance, offset for odd N (M = Q + (N mod 2)/N) | `alternating_order_squared` (= ⟨M²⟩) |
| τ_k = n_{2k} − n_{2k+1} | pair pseudo-spin of strong pair k (L = atom 2k, R = atom 2k+1) | `tau`, `pair_order()` |
| J_m, B_m, h | effective pseudo-spin couplings at pair separation m and transverse field | `effective.tfim_couplings` |
| T₂, I_s, I_w | translation by one pair; inversion about a strong-bond centre; I_w = T₂I_s | — |
| E₀, E₁, G = E₁ − E₀ | two lowest energies and the gap | `ground_energy`, `gap`, `gap_over_omega` |
| S_A | von Neumann entropy of the reduced state of atoms A, in bits | `entropy_bits`, `S_*_bits` |

The letter S denotes entropy only. Some `study/` code and CSV columns call the staggered order "S" or "staggered spin"; that quantity is Q here.

**Sign of J_m.** `docs/DIMERIZED_THEORY.md` writes the pair interaction as U_m = A_m + B_m(τ_k − τ_{k+m}) − J_m τ_kτ_{k+m} with J_m = ¼[V(2m−s) + V(2m+s) − 2V(2m)] > 0. `rydberg_chain.effective.tfim_couplings` returns the same coefficient with the opposite sign (U_m contains +J_m τ_kτ_{k+m}, so J_m < 0 is ferromagnetic). Both describe ferromagnetic pseudo-spin coupling; compare magnitudes, and state which convention a table uses.

## Interaction ranges

| Name | Meaning |
|---|---|
| `nearest1` | only \|i − j\| = 1 |
| `nearest2` | \|i − j\| ≤ 2 |
| `all` | every pair (the reference model) |
| `rmax` (ED, atom DMRG) | keep \|i − j\| ≤ rmax; default: all pairs |
| `mmax` (pair-site DMRG) | keep inter-pair couplings up to pair separation mmax; intra-pair V_s always kept; default: all |

## Geometry of the dimerized chain

The gap after atom i is 1 − ε for even i and 1 + ε for odd i, so atoms (2k, 2k+1) form strong pair k, and an open chain always starts with a strong bond. Pair k sits at positions 2k and 2k + s. ε = 0 is the uniform chain. For even N an open chain holds N/2 complete pairs; for odd N the last atom is unpaired.

## Basis ordering

- **Package and `study/` code:** bit i of the basis index is the occupation of atom i, so index k = Σᵢ nᵢ 2ⁱ and atom 0 is the least significant bit. Example for N = 3: index 1 is atom 0 excited, index 4 is atom 2 excited.
- **Qiskit labels** print qubit 0 on the right: the label `001` means atom 0 excited. The Qiskit operators in `study/rydberg_model.py` use the same bit order as the package, so vectors agree element by element.
- **Written configurations** in the manuscript, figures and `study/results/parity_top_bitstrings.csv` (column `bitstring_atom0_left`) list atoms left to right starting from atom 0, as they sit in the chain.
- **Pair sites** (pair-site DMRG) use the local basis |n_L n_R⟩ = 00, 10, 01, 11 in that order.

## Other conventions

- Entropies are in bits (log₂). A Bell pair gives 1 bit.
- "Open" chains have free ends; "rings" are periodic, with the minimum-image distance for V.
- Energies per atom are E/N. Densities, order parameters and gaps are ground-state expectation values unless a figure states otherwise.
- Interaction switch-off controls (`isolated`, `split` in `study/rydberg_model.interactions`) set couplings to zero while keeping every atom and its drive. They are not vacancies.
