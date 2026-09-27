# Results and limitations

This page summarizes what the saved evidence in `study/results/` supports, and what it does not. Conventions: |0⟩ ground, |1⟩ Rydberg, ħ = Ω = a = 1, open chains with couplings g/|i − j|⁶ unless stated; symbols as in [NOTATION.md](NOTATION.md). Publication versions of the figures, with exact parameters and captions, are listed in [`paper/figures/FIGURES.md`](../paper/figures/FIGURES.md).

## Evidence available

| Study | Saved evidence | What it shows |
|---|---|---|
| Two-atom blockade | `qiskit_baseline.csv`, `blockade_scan.csv` | Repulsion suppresses double excitation, but never completely: the time-averaged P₁₁ over 0 ≤ Ωt ≤ 2π falls from 3/8 at V = 0 to 1.5 × 10⁻⁴ at V/Ω = 40 and approaches 0.2356(Ω/V)². |
| Interaction switch-off controls (idealised) | `defect_*.csv` | Setting chosen couplings to zero, while every atom and its drive are kept, removes the corresponding entanglement. These controls are not vacancies, which would remove an atom's drive and Hilbert space as well. |
| Four-atom geometries | `geometry_*.csv` | Entanglement follows the short, strongly blockaded bonds and depends on the chosen bipartition. |
| Interaction-range cutoffs | 6,561 coarse and 312 refinement ground states, N = 4–12, g = 2, 8, 16 (`ground_states_*.csv`, `cutoff_errors_coarse.csv`, `maximum_errors.csv`) | Truncating the 1/r⁶ tail changes densities, order and gaps by amounts that grow with g and, for nearest-neighbour truncation, with N. This is the reason later calculations keep the tail. |
| Odd/even chains and low states | `parity_*.csv`, `parity_summary.json` | Odd and even open chains differ qualitatively in the ℤ₂ regime; even chains hold a boundary-forced domain wall. |
| Kagome patch (exploratory appendix) | `kagome_*.csv`, `kagome_*summary.json` | A twelve-site patch and a same-size strip control. Mean density can hide local redistribution, but the shape comparison also changes coordination and bond counts. Not part of the chain question. |

## Interaction-range cutoffs: why the full tail is kept

Each truncated model keeps nearest neighbours (`nearest1`) or neighbours up to index separation two (`nearest2`), and is compared with the all-pairs ground state at the same N, g and Δ.

- With nearest neighbours only, the largest sampled density error at g = 16 is **8.0916 percentage points for N = 11** and **4.8217 points for N = 12**, each at its own sampled detuning (Δ/Ω = 0.75 and 0.6875). Every odd chain exceeds both neighbouring even chains, and at g = 8 and 16 the error is still growing with N within each parity at N = 12. An even-only size sequence therefore understates it.
- Keeping second neighbours reduces the largest error over all 27 size/strength combinations to **0.1863 points** (N = 4, g = 16).
- The errors peak near the detuning where the density rises most steeply. Dividing the error by that slope gives an equivalent detuning shift of 1.05–1.47 times the omitted next-nearest repulsion V(2a) = g/64 for N = 11 and 12 at g = 8 and 16.

These maxima come from a grid of 81 detunings (−3 ≤ Δ/Ω ≤ 7) plus local refinement; they are sampled values, not bounds. The dimerized-chain calculations therefore keep all pairs in ED. The pair-site DMRG benchmarks keep couplings up to eight pair separations, which shifts the N = 20 ground-state energy by 2.4 × 10⁻⁷ relative to all pairs.

**A bound on energies.** For repulsive interactions in a uniform open chain, discarding pairs beyond separation R removes a positive diagonal operator of norm

$$B=g\sum_{r=R+1}^{N-1}\frac{N-r}{r^6}.$$

Each ordered eigenvalue of the full model lies between that of the truncated model and that value plus B, so the lowest gap changes by at most B. This does not bound densities or guarantee that eigenvectors stay close, particularly when levels compete. The values are in `discarded_tail_bounds.csv` and the `tail_norm_bound` column of `maximum_errors.csv`.

## Interaction switch-off controls

In a six-atom chain (V₀/Ω = 5, Δ = 0, initial state |000000⟩, 0 ≤ Ωt ≤ 4):

- with every interaction of atom 2 set to zero, atom 2 Rabi-oscillates as sin²(Ωt/2) to 1.2 × 10⁻¹⁵ and stays unentangled;
- with every coupling across the central cut set to zero, the half-chain entropy stays at rounding level (below 10⁻¹⁴ bits), against 1.25 bits at Ωt = 4 in the intact chain;
- with only atom 2 decoupled, the half-chain entropy still reaches 0.017 bits by Ωt = 4, because atoms 1 and 3 remain coupled through the tail V₀/2⁶ = 0.078Ω.

These are model interventions that isolate the role of particular couplings. They do not simulate atom loss.

## Four-atom geometries

With V₀/Ω = 5 at the nearest spacing, the rectangle's partition {0,2}|{1,3}, which cuts both short bonds, reaches 1.97 bits at Ωt = 2.23, while {0,1}|{2,3}, which cuts neither, stays below 0.031 bits. In the chain the interleaved partition carries 1.52 bits at Ωt = 2, twice the contiguous value of 0.75 bits. The square's two partitions are identical by symmetry.

## Odd and even chains: the boundary-forced domain wall

At g = 8 and Δ/Ω = 3 (all pairs), odd chains N = 5–11 keep ⟨Q²⟩ = 0.938–0.942 and a gap of 2.27–2.56Ω. In even chains ⟨Q²⟩ falls from 0.32 to 0.17 and the gap from 0.26Ω to 0.079Ω between N = 4 and 12 (`parity_fixed_setting_scan.csv`).

The cause is the boundary. Both end atoms favour excitation, which a perfect alternating pattern of an even chain cannot provide, so the ground state carries one 00 bond. At N = 8 the two perfect alternating strings together carry only **6.35%** of the ground-state probability, while configurations with exactly one 00 bond and no 11 bond carry **76.68%**. The two lowest states place the 00 bond on the three bonds between pairs with probabilities 0.29, 0.42, 0.29 and 0.38, 0.03, 0.38 (`parity_bond_probabilities.csv`), the pattern of the two lowest standing waves of one wall. The small gap at this setting is therefore the splitting of a delocalized wall rather than tunnelling between two perfect patterns. The same wall slows finite DMRG on long even chains (see [METHODS.md](METHODS.md)).

The numbers are in `parity_low_states.csv`, `parity_top_bitstrings.csv` (configurations printed with atom 0 on the left) and `parity_summary.json`.

**Order diagnostic.** Q = N⁻¹Σᵢ(−1)ⁱ(2nᵢ − 1) gives Q² = 1 for either perfect pattern at any N; the saved column `staggered_spin_squared` is ⟨Q²⟩. The older M = (2/N)Σᵢ(−1)ⁱnᵢ (column `alternating_order_squared`) is offset for odd N and is kept only for comparison. For independent sites with uniform excitation probability p, ⟨Q²⟩ = 4p(1−p)/N + (2p−1)²(N mod 2)/N², so ⟨Q²⟩ alone does not identify a phase in short chains.

**N = 8 map.** Over −3 ≤ Δ/Ω ≤ 7 and 0.5 ≤ g ≤ 12, ⟨Q²⟩ of the eight-atom chain never exceeds 0.237, for the same boundary reason. The map describes one finite size and is not a phase diagram.

## Dimerized chains

The current research question concerns dimerized chains. Its pilot results (critical detunings, gap ratios, the ℤ₄ splitting and the ℤ₃ control) are preliminary and are tracked in [RESEARCH_PLAN.md](RESEARCH_PLAN.md) and [DIMERIZED_THEORY.md](DIMERIZED_THEORY.md). No dimerized-chain number is final until it has been recomputed with `rydberg_chain/` and, where a restricted basis was used, checked against the full Hilbert space.

## Scope

- The calculations model ideal atoms: no decoherence, state-preparation or detection errors, and no position disorder. They ran on a classical computer.
- The completed study covers N ≤ 12 with open ends, and the kagome appendix one twelve-site patch. Neither establishes bulk phase boundaries, critical exponents or topological order.
- Numerical agreement checks the implementation, not the physical assumptions.

Interaction range, boundaries and ordering in Rydberg chains have been studied before; relevant context includes [Bernien et al.](https://arxiv.org/abs/1707.04344), [Rader and Läuchli](https://arxiv.org/abs/1908.02068), [Slagle et al.](https://arxiv.org/abs/2108.09309), [Wang et al.](https://arxiv.org/abs/2502.07886) and, for dimerized chains with exchange interactions, [Zeybek, Mukherjee and Schmelcher](https://arxiv.org/abs/2303.17470). Conventions differ between these works; quantitative comparisons require matching them.
