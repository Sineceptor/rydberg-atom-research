# Dimerized chain: effective Ising theory (working notes)

These notes derive the analytic picture behind the proposed "alternating coupling" results section. They are meant to be re-derived by hand; every step is second-year quantum mechanics (perturbation theory) plus the transverse-field Ising model.

## Geometry and Hamiltonian

Atoms sit on a line with alternating gaps: a short gap s = a(1 − ε) inside each pair and a long gap w = a(1 + ε) between pairs. The mean spacing a is fixed, so s + w = 2a. Pair ("dimer") k holds atoms L = 2k and R = 2k + 1 at positions 2ka and 2ka + s.

$$H=\frac{\Omega}{2}\sum_i X_i-\Delta\sum_i n_i+\sum_{i<j}\frac{C_6}{|x_i-x_j|^6}n_in_j .$$

Units: ℏ = Ω = a = 1 and g = C₆/(Ωa⁶) = (R_b/a)⁶. The strong and weak nearest-neighbour energies are V_s = g/(1 − ε)⁶ and V_w = g/(1 + ε)⁶.

**Exact regrouping.** Each pair can be treated as one site with four states |00⟩, |10⟩, |01⟩, |11⟩. The large V_s then becomes an on-site energy, and only the weaker inter-pair terms connect sites. This is the same Hamiltonian, not an approximation. It is what makes DMRG stable at strong dimerization and halves the chain length (validated against exact diagonalization to 10⁻⁹ or better).

## Low-energy manifold: one excitation per pair

For V_s ≫ Δ > 0 and moderate inter-pair interactions, the lowest states put exactly one excitation in every pair. Each pair then has two states, |L⟩ (left atom excited) and |R⟩ (right atom excited). Call this a pseudo-spin: τᶻ = +1 for L and −1 for R.

**Zeroth order (Ω = 0).** Excitations in pairs k and k + m are separated by d = 2m + s(τ_k − τ_{k+m})/2. Write u₀ = g/(2m)⁶ and u_± = g/(2m ± s)⁶. The pair energy then decomposes exactly as

$$U_m=A_m+B_m(\tau_k-\tau_{k+m})+J_m\tau_k\tau_{k+m},\qquad J_m=\tfrac12\Big[u_0-\tfrac12(u_++u_-)\Big],\quad B_m=\tfrac14(u_+-u_-).$$

- J_m < 0 because 1/r⁶ is convex, so the pseudo-spins interact ferromagnetically. Both perfect crystals (…1010… and …0101…) are the two ferromagnetic states.
- For m = 1, u₋ = V_w. So |J₁| = [V_w + g/(2+s)⁶ − 2g/64]/4. The next-nearest repulsion g/64 **reduces** the coupling, because a "00" domain wall moves two excitations farther apart. This correction is not small near melting; at g = 11.39 and ε = 0.6 it reduces the coupling by about 44% relative to V_w alone.
- The B terms telescope in the bulk and leave only **boundary fields**. These favour exciting the outermost atom at both ends. In an even chain, that forces one "00" domain wall into the interior. This explains the earlier parity result that at N = 8, g = 8, δ = 3 about 77% of the ground-state weight has exactly one 00 bond.

**Second order in Ω.** Flipping one atom takes |L⟩ or |R⟩ to |00⟩ (energy cost Δ) or to |11⟩ (cost V_s − Δ). Both routes connect L and R equally, giving

$$H_\text{eff}=\sum_{k,m}J_m\tau^z_k\tau^z_{k+m}-h\sum_k\tau^x_k+\text{(boundary fields)},\qquad h=\frac{\Omega^2}{4}\Big[\frac1\Delta+\frac1{V_s-\Delta}\Big].$$

The transverse field h is the pair's laser-driven "superatom" delocalization: the excitation tunnels between the two atoms of a pair through a virtual empty or doubly excited state.

## Prediction

The nearest-neighbour transverse-field Ising chain orders for |J₁| > h and is critical at |J₁| = h. Hence the Z2 crystal melts at approximately

$$V_w+\frac{g}{(2+s)^6}-\frac{2g}{64}\approx\frac{\Omega^2}{\Delta}\qquad(V_s\gg\Delta).$$

Consequences to test numerically:

1. Increasing ε melts the crystal into a **pair paramagnet**. That state has density ≈ 1/2, but each excitation is shared symmetrically within its pair, so ⟨n_L − n_R⟩ = 0.
2. The melting dimerization ε_c increases with Δ, because h ∝ 1/Δ.
3. The transition should be in the 2D Ising universality class (central charge c = 1/2), where this effective model holds.
4. Corrections: second-order energy denominators shift with neighbour interactions, and longer-range J_m and higher orders in Ω/Δ contribute. Numerics decide how accurate the leading estimate is.

At small Δ, empty pairs |00⟩ become cheap. The relevant effective model then has three states per pair (empty, L, R), like a quantum Blume–Capel chain. Such chains can host tricritical points, where a continuous Ising transition becomes first order. Whether that happens here is an open question for the numerics; no claim is made yet.

## Pilot evidence (26 Sept 2026)

- Exact diagonalization, N = 16, g = 11.39 (R_b/a = 1.5), Δ = 2.5: interior |⟨n_L − n_R⟩| = 0.21, 0.19, 0.17, 0.10, 0.02, 0.00 at ε = 0.25, 0.4, 0.5, 0.6, 0.7, 0.85. The density stays near 0.46–0.47. The dominant configurations contain one central 00 wall, as predicted.
- Leading-order estimate at the same settings: |J₁|/h ≈ 0.96 at ε = 0.6, i.e. melting near ε ≈ 0.6.
- Infinite-chain DMRG (pair sites, bond dimension 64): ordered with ⟨n_L − n_R⟩ = 0.95, 0.95, 0.94 at ε = 0, 0.2, 0.3. **Beyond ε ≈ 0.3 the infinite-chain results are not converged:** the energy per atom jumps nonphysically (to −0.80 at ε = 0.55). No phase boundary may be quoted from them until the solver passes the validation gate in [RESEARCH_PLAN.md](RESEARCH_PLAN.md).
- Finite-chain checks: atom-level DMRG matches the existing exact code to about 10⁻¹² (including the saved N = 12 uniform-chain energy); pair-site DMRG matches exact diagonalization at N = 16 to 2 × 10⁻⁹ or better.
- Exact-diagonalization cost on a 4-core, 15 GB machine: 4 s (N = 16), 25 s (N = 18), 196 s and 2.8 GB (N = 20); the practical ceiling is N ≈ 22.
