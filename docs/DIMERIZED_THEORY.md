# Dimerized chain: pair pseudo-spin theory (working notes, revised 27 September)

These notes give the analytic picture behind the alternating-coupling results. They are written to be re-derived by hand; the tools are second-order perturbation theory and the transverse-field Ising model. Section 6 lists corrections to the 26 September version.

## 1. Geometry, Hamiltonian and symmetries

Atoms sit on a line with alternating gaps: a short gap s = a(1 − ε) inside each pair and a long gap w = a(1 + ε) between pairs. Because s + w = 2a, pair k holds atoms L = 2k and R = 2k + 1 at positions 2ka and 2ka + s.

$$H=\frac{\Omega}{2}\sum_i X_i-\Delta\sum_i n_i+\sum_{i<j}V(|x_i-x_j|)\,n_in_j,\qquad V(r)=\frac{C_6}{r^6}.$$

Units: ℏ = Ω = a = 1 and g = C₆/(Ωa⁶) = (R_b/a)⁶. Strong and weak nearest-neighbour energies are V_s = V(s) and V_w = V(w).

Symmetries of the infinite chain:
- translation by one pair, T₂;
- inversion about the centre of a strong bond, I_s, which swaps L and R inside that pair;
- inversion about a weak-bond centre, I_w = T₂I_s.

The uniform chain (ε = 0) also has translation by one atom, which dimerization removes.

**Exact regrouping.** Treat each pair as one site with four states |00⟩, |10⟩, |01⟩, |11⟩. V_s becomes an on-site energy and only inter-pair terms connect sites. This is the same Hamiltonian, not an approximation. Pair-site DMRG matches exact diagonalization to 10⁻⁸ or better in the package tests (N = 12–14), and to 4 × 10⁻¹³ for the 20-atom benchmark (validation/results/scaling_benchmark.json).

## 2. Pseudo-spin couplings

An occupied pair has two classical states: |L⟩ (left atom excited) and |R⟩. Define τᶻ = +1 for L and −1 for R, so I_s flips τ.

**Zeroth order in Ω.** Take excitations in pairs k and k + m. With u₀ = V(2m), u₊ = V(2m + s) and u₋ = V(2m − s), their interaction decomposes exactly as

$$U_m=A_m+B_m(\tau_k-\tau_{k+m})-J_m\tau_k\tau_{k+m},\qquad J_m=\tfrac14\left[u_-+u_+-2u_0\right],\qquad B_m=\tfrac14\left(u_+-u_-\right).$$

- **J_m > 0 (ferromagnetic):** it is a second difference of the convex function 1/r⁶.
- **J_m → 0 as ε → 1:** u₊ and u₋ approach u₀.
- **m = 1:** J₁ = [V_w + V(2a+s) − 2V(2a)]/4. The next-nearest repulsion V(2a) reduces the coupling (by about 44% relative to V_w/4 at g = 11.39, ε = 0.6).
- **m = 2:** J₂ = [V(3+ε) + V(5−ε) − 2V(4)]/4 couples every-other-pair occupations. This is the relevant coupling in the Z₄ regime.
- **Boundary fields:** the B_m terms telescope in the bulk and leave only fields at the ends. These favour exciting the outermost atom. An even open chain in the Z₂ regime is therefore forced to hold one "00" domain wall. This is standard kink-in-a-box physics, and it explains the earlier finding that about 77% of the N = 8 (g = 8, δ = 3) ground-state weight has exactly one 00 bond.

**Second order in Ω.** Flipping one atom takes |L⟩ or |R⟩ to |00⟩ or |11⟩, and both routes connect L with R. With bare energy denominators this gives

$$H_\text{eff}=-\sum_{k,m}J_m\tau^z_k\tau^z_{k+m}-h\sum_k\tau^x_k+\ldots,\qquad h_0=\frac{\Omega^2}{4}\Big[\frac1\Delta+\frac1{V_s-\Delta}\Big].$$

The denominators should include the interaction energy of the removed excitation with its neighbours. In the Z₂ crystal that energy is I_L = 2V(2a) for an L excitation and I_R = V_w + V(2a+s) for an R excitation, giving

$$h\approx\frac{\Omega^2}{8}\Big[\frac1{\Delta-I_L}+\frac1{\Delta-I_R}\Big]+\frac{\Omega^2}{4(V_s-\Delta)}.$$

This dressed form is empirical so far. A controlled Schrieffer–Wolff derivation is planned, including the renormalization of J and the symmetry-allowed odd term τˣ(τᶻ − τᶻ).

## 3. Z₂ regime: one Ising line

The ordered state (all τ = +1 or all τ = −1) breaks I_s and keeps T₂. The nearest-neighbour Ising chain is critical at J₁ = h. The pilot ring exact diagonalization at g = 11.39 gives the following (size pairs agree to 1e-4; gap ratios R = 7.9–8.0, Ising). Dressed values are recomputed from the formula in section 2 with nearest-pair J₁; other reasonable dressing prescriptions change the ε ≥ 0.6 errors within 0.2–1.1%:

| ε | Δ_c (ED) | leading-order h₀ | dressed h |
|---|---|---|---|
| 0.3 | 1.326 | 0.494 (−63%) | 2.67 (+101%, fails) |
| 0.45 | 1.805 | 1.101 (−39%) | 2.08 (+15%) |
| 0.6 | 3.149 | 2.616 (−17%) | 3.18 (+0.9%) |
| 0.7 | 5.611 | 5.186 (−7.6%) | 5.64 (+0.5%) |

- For ε ≲ 0.45 the transition is driven by empty pairs, not τ flips, and no pseudo-spin model describes it.
- At Δ = 2.5 the bulk crystal melts at ε_c ≈ 0.53–0.55.

The regime beyond the line, where each pair shares one excitation, is **not** a separate phase. There is no U(1) symmetry, and each pair's ground state (a|00⟩ + b(|10⟩ + |01⟩)/√2) is I_s-even. It is adiabatically connected to the ordinary disordered phase: a crossover. This differs from the U(1)-protected bond-order phase of Zeybek et al.

## 4. Even-period crystals split, odd-period crystals are protected (headline hypothesis)

**Z₄ regime (R_b/a ≈ 3–3.5).** At ε = 0 there are four period-4 patterns, permuted by D₄. Dimerization leaves the Klein group Z₂×Z₂, generated by T₂ and I_s. All four patterns remain degenerate, but the order parameter separates into a pair-density-wave part u (which pairs are occupied) and a pseudo-spin part τ. The single Z₄ melting can therefore split into two transitions:
- Z₄ crystal → **pair density wave (PDW)**: u ordered, τ disordered; breaks T₂, preserves I_s;
- PDW → disordered.

The τ transition is governed by J₂ against h. At large ε the PDW is simply the Z₂ crystal of a chain of two-atom superatoms. It resembles the rung-symmetric Z_p⁺ phases of two-leg ladders. Here τ is reversed by a spatial inversion rather than by an internal leg exchange, and a derivative (Huse–Fisher chiral) coupling between u and the I_s-odd components of the order parameter is symmetry-allowed near ε = 0. Chiral transitions are also reported in ladders, so this is a difference in mechanism, not a claim that ladders forbid chirality.

**Ashkin–Teller view (added 28 September; re-derive by hand).** Label the four patterns by k = 0…3, the position of the excitation modulo 4, and write them as two Ising variables: (σ₁, σ₂) = (+,+), (−,+), (−,−), (+,−) for k = 0, 1, 2, 3. Then u = σ₂ and τ = σ₁σ₂.
- Translation by one atom, T₁ (a symmetry only at ε = 0), maps (σ₁, σ₂) → (−σ₂, σ₁). It exchanges the two Ising copies, which is why they order together in the uniform chain.
- T₂ flips both variables, and I_s (k → 1 − k) flips σ₁ alone. After dimerization these two flips are separate symmetries, and no remaining operation exchanges σ₁ with σ₂.
- T₁ also maps ε → −ε. An anisotropy λ(ε)(ε₁ − ε₂) between the two Ising energy operators is therefore allowed, with λ odd in ε: dimerization enters at first order. A term in τ = σ₁σ₂ itself is forbidden, because I_s flips it.
- An anisotropic Ashkin–Teller model generically melts through two Ising transitions, with an intermediate phase in which only one variable orders. The PDW is σ₂ ordered and σ₁ disordered, so ⟨τ⟩ = 0. The alternative, σ₁ ordered alone, would be a PDW centred on the weak bonds. The pilot spectra (both states I_s-even, T₂ = ±1) point to the strong-bond PDW for ε > 0.
- Other models follow the same two-Ising route out of a Z₄ phase, with the two lines merging at an Ashkin–Teller point: multicomponent chains (arXiv:2308.12838), quantum loop models on a zig-zag ladder (arXiv:2406.20093) and two-leg ladders (arXiv:2411.05494, 2412.11201). What is specific here is that the anisotropy comes from geometry, is linear in ε, and vanishes at the uniform chain.
- A standard reference for the anisotropic Ashkin–Teller model must be found and cited [verify].

**Comparison with the ladder crystals (full texts read 28 September).**
- *Square ladder (arXiv:2411.05494).* Every crystal except the checkerboard appears twice. In one version each occupied rung shares its excitation between the legs; in the other the leg symmetry is broken. An Ising transition separates the two, and the two Ising lines fuse into an extended Ashkin–Teller interval, followed by a Z₄-chiral transition and a floating phase.
- *Rectangular ladder, rungs half the leg spacing (arXiv:2412.11201).* The leg-symmetric Z_p⁺ crystals give way through Ising transitions to leg-broken crystals that are staggered, which doubles the period to 2p. They are staggered because atoms on opposite legs of neighbouring occupied rungs are farther apart than atoms on the same leg. Here J_m > 0, so τ orders uniformly, and the τ-ordered state is the uniform chain's Z₄ crystal.
- *In the ladders the leg swap is an internal symmetry.* It forbids the odd couplings B_m of section 2. Here they are allowed, and they produce the boundary fields and the forced domain wall in even open chains.

**Z₃ regime (R_b/a ≈ 2.7).** T₂ still cycles the three period-3 patterns and I_s exchanges two of them, so D₃ survives. Dimerization should add no new transition, only non-universal shifts. In the square ladder a period-3 crystal of rungs does split (arXiv:2411.05494), because every ladder crystal is built from whole rungs. The period-3 crystal of the dimerized chain straddles the pairs, so it cannot. A crystal with one excitation on every third pair has period 6 and falls under the even case.

**Pilot evidence (constrained ring ED; must be validated against full ED, gate P1):**
- **R_b/a = 3.3, Δ = 3.** The inversion-odd (τ-channel) gap at N = 16/24/32 is:
  - ε = 0: 0.011 / 0.0020 / 0.0004 (closing);
  - ε = 0.3: 0.046 / 0.029 / 0.022 (N·gap ≈ 0.72, critical);
  - ε = 0.7: 0.166 / 0.164 / 0.164 (gapped).

  For ε ≥ 0.4 the two lowest states have T₂ = ±1 and are both I_s-even.
- **R_b/a = 2.7, Δ = 5, N = 24.** The Z₃ threefold degeneracy is exact at ε = 0, 0.3 and 0.6.

## 5. Why dimerization does not produce a tricritical point (classical-limit argument)

At small Δ each pair is a three-state superatom (empty, L, R), which maps onto a spin-1 Blume–Emery–Griffiths chain.
- **Condition for a first-order jump.** A direct jump from empty to ordered, which is the Blume–Capel route to tricriticality, requires net attraction between aligned occupied pairs.
- **Why the condition fails.** Aligned occupied pairs are always 2a apart, so their interaction V(2a) > 0 does not depend on ε.
- **Trend with ε.** Dimerization lowers J₁ while V(2a) stays fixed, so the biquadratic-to-bilinear ratio K/J = −1 − V(2a)/J₁ moves away from the tricritical value.
- **Numerics.** Ring ED gap ratios stay at R ≈ 8 along the line for R_b/a from 1.3 to 1.85 and ε ≤ 0.7.

This is a classical-limit criterion backed by numerics, not a theorem. It is consistent with the need for an inverted sweep (effective V₂ < 0) in the uniform-chain tricritical experiment (arXiv:2601.16275), and for attractive interspecies interactions in ladder proposals.

## 6. Corrections to the 26 September notes

- The interior open-chain values of |⟨n_L − n_R⟩| (0.21 … 0.00) are contaminated by the boundary-forced domain wall. They do not measure bulk order.
- The bulk melting point at Δ = 2.5, g = 11.39 is ε_c ≈ 0.53–0.55, not 0.6.
- The leading-order effective model is quantitative only at strong dimerization (see the table in section 3).
- "Pair paramagnet" was a misnomer: the regime is a crossover within the disordered phase.
- The early infinite-DMRG scan was unconverged beyond ε ≈ 0.3 (energies up to −0.80 per atom where ring ED gives −1.2284). No number from it is used.
