# Results and limitations

The corrected study uses |0⟩ for the ground state, |1⟩ for a Rydberg excitation and ℏ=Ω=a=1. Open-chain interactions are g/|i−j|⁶. Three models retain nearest neighbours, neighbours through index separation two, or every pair. Comparisons use matching sizes, drive, detuning and state definitions.

## Evidence available

| Study | Saved evidence | Interpretation |
|---|---|---|
| Blockade dynamics | Two atoms, a common initial state and time interval; `blockade_scan.csv` | Increasing repulsion suppresses simultaneous excitation in the resonant examples. A finite interaction does not make its probability exactly zero. |
| Defects and geometry | Disconnected interaction controls and four-atom shapes | Entanglement depends on which couplings and partitions are chosen. A removed interaction is an ideal model intervention. |
| Interaction cutoffs | 6,561 coarse ground states and 312 local-refinement states, N=4 through 12, g=2,8,16 | Discarding distant pairs changes predictions. Accuracy depends on size, parameters, observable and cutoff. |
| Parity and low states | Site profiles, bond probabilities, configuration weights and several low energies | Odd and even open chains can behave differently. A low gap needs a state-level explanation. |
| Optional kagome pilot | Twelve-site patch plus a same-size strip control | Mean density can hide local redistribution. The shape comparison changes coordination and triangle/bond counts as well as boundaries. |

## Selected quantitative examples

At g=16, the largest sampled nearest-neighbour-only density discrepancy is approximately **8.0916 percentage points for N=11**, compared with **4.8217 points for N=12**. These maxima occur at their respective sampled detunings. This demonstrates why an even-only size sequence cannot justify a general monotonic trend. The maximum sampled second-neighbour density discrepancy across the 27 size/strength combinations is about **0.1863 points**, at N=4, g=16; a smaller selected value is not a global bound.

At N=8, g=8, δ=3, the two perfectly alternating configurations together carry only about **6.35%** of the ground-state probability. The group with exactly one 00 bond and no 11 bond carries about **76.68%**. Consequently the small gap at that setting should not automatically be explained as tunnelling between only two perfect strings. Spatial profiles and additional low states support a boundary-sensitive interpretation for this example; they do not identify the mechanism of every small gap.

The underlying numbers are in `maximum_errors.csv`, `parity_low_states.csv`, `parity_top_bitstrings.csv` and `parity_summary.json`. Atom 0 is the rightmost bit in conventional Qiskit labels, while the exported `parity_top_bitstrings.csv` explicitly prints atom 0 on the left.

## A corrected order diagnostic

The historical occupation-imbalance quantity M=(2/N)Σᵢ(−1)ⁱnᵢ has unequal ideal-pattern magnitudes at odd N. The corrected quantity is

$$Q=\frac1N\sum_i(-1)^i(2n_i-1).$$

Either perfect alternating string has Q²=1 for both parities. The saved column `staggered_spin_squared` represents ⟨Q²⟩; some code calls Q “S”. Entanglement entropy is a different quantity.

For independent sites with uniform excitation probability p,

$$\langle Q^2\rangle=\frac{4p(1-p)}N+\frac{(2p-1)^2(N\bmod2)}{N^2}.$$

A nonzero value in a short system therefore does not identify a phase. Use configurations, correlations and size dependence alongside it. Cropping the two endpoints from an average does not remove every boundary effect.

## A bound with a precise scope

For repulsive interactions in a uniform open chain, discarding pairs beyond separation R removes a positive diagonal operator of norm

$$B=g\sum_{r=R+1}^{N-1}\frac{N-r}{r^6}.$$

Each ordered eigenvalue of the full model can shift upward relative to the truncated model by between zero and B; hence the absolute change in the lowest gap is at most B. This does **not** give a bound on density or guarantee that wavefunctions remain close, particularly when levels compete.

## What is not established

These calculations omit experimental decoherence, state-preparation and detection errors. They run on a classical computer. Small open chains and a twelve-site patch do not establish bulk phase boundaries, critical exponents, topological protection or a quantum spin liquid. A numerical residual checks an eigenpair, not every physical assumption or a claim of novelty.

Published work already addresses interaction range and boundaries. Relevant context includes [Bernien et al.](https://arxiv.org/html/1707.04344v2), [Rader and Läuchli](https://arxiv.org/abs/1908.02068), [Slagle et al.](https://arxiv.org/html/2108.09309v2) and [Wang et al.](https://arxiv.org/abs/2502.07886). Model conventions and approximations differ across these papers; quantitative comparisons require matching them.
