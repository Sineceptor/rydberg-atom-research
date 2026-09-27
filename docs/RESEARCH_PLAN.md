# Research plan to 15 October 2026 (revised 27 September)

This revision replaces the 26 September draft. An adversarial literature review and new pilot calculations changed the headline. The direction must be confirmed with Prof. Rubenstein before production runs.

## 1. Why the plan changed

- **Zeybek, Mukherjee & Schmelcher** (PRL 131, 203003 (2023); arXiv:2303.17470), read in full. Their model uses two Rydberg levels, dipolar hopping, fixed particle number and the laser off. In that dimerized chain they already find a density-wave → bond-order transition at half filling. "Dimerization melts the Z2 crystal into shared pairs" is therefore not a new idea, although our laser-driven, single-level model is a different system.
- **Tricritical Ising (TCI)** has already been measured in a uniform single-species vdW chain (Sun et al., arXiv:2601.16275, Nature 2026, using an inverted sweep). Pilot exact diagonalization also finds only Ising behaviour along the dimerized Z2 line: gap ratios 7.87–8.15 against the Ising value 8. There is also a classical-limit reason. Tricriticality needs a net attraction V₂ < 0 between aligned occupied pairs. Aligned pairs always sit 2a apart, so dimerization cannot supply it. TCI is therefore a pre-registered negative test, not the headline.
- **The single boundary domain wall** in even open chains is textbook kink-in-a-box physics. It becomes a methods/system-size subsection.
- **Degenerate-packing selection** is trivial here: the unit-disk graph is disjoint edges, and the vdW tails lift the degeneracy at first order. At most one motivating sentence.

## 2. Direction

**Working title:** *Pair pseudo-spins in dimerized Rydberg chains: even-period crystals split, odd-period crystals are protected.*

**Model:** single Rydberg level, laser-driven, van der Waals only. Gaps alternate s = a(1−ε) and w = a(1+ε); units Ω = a = 1 and g = (R_b/a)⁶. Symmetries: translation by one pair T₂, inversion about a strong-bond centre I_s, and I_w = T₂I_s.

**Unifying idea.** An occupied strong pair carries a pseudo-spin τ = n_L − n_R, which I_s flips. Occupied pairs m pair-spacings apart couple through the second difference

$$J_m(\varepsilon)=\tfrac14\left[V(2ma-s)+V(2ma+s)-2V(2ma)\right],$$

- **Sign:** positive by convexity of 1/r⁶, so the coupling is ferromagnetic in τ.
- **Limit:** it vanishes as ε → 1, while the laser-induced transverse field h stays finite. Every crystal whose pairs carry a free L/R choice must therefore melt its τ order at some ε.

**The even/odd rule (headline hypothesis):**

- **Even-period crystals split.** In the Z₄ regime (R_b/a ≈ 3–3.5), dimerization reduces the symmetry acting on the four period-4 patterns from D₄ to the Klein group Z₂×Z₂. The single Z₄ melting then separates into two transitions, with an intermediate **pair-density-wave (PDW) phase** between them. In the PDW, every other strong pair holds one excitation shared symmetrically between its atoms. It breaks T₂ and preserves I_s, so it is a distinct phase.
- **Odd-period crystals are protected.** In the Z₃ regime, T₂ still cycles the three patterns and I_s swaps two of them. D₃ survives, so there is no new transition, only shifts.
- **Z₂ backbone.** The Z₂ crystal melts through a single Ising line Δ_c(ε, g) into the ordinary disordered phase.
  - The regime where each pair shares one excitation is a crossover, not a phase: there is no U(1), and the inversion index is trivial.
  - An effective transverse-field Ising model with a dressed energy denominator predicts this line quantitatively at strong dimerization.

**Pilot evidence** (exploratory scripts: [research/pilots/2026-09-27](../research/pilots/2026-09-27/README.md)):

| Question | Result so far |
|---|---|
| Z₂ line, g = 11.39 (R_b/a = 1.5), ring ED | Δ_c = 1.0785, 1.326, 1.805, 3.149, 5.611 at ε = 0, 0.3, 0.45, 0.6, 0.7 (size pairs agree to 1e-4). R = 7.9–8.0. |
| Effective theory | Leading order misses Δ_c by 63, 39, 17, 7.6% at ε = 0.3, 0.45, 0.6, 0.7. Dressed denominators give +0.9% and +0.5% at ε = 0.6 and 0.7, but +15% at 0.45 and fail at 0.3. |
| Z₄ → PDW, R_b/a = 3.3, Δ = 3, constrained ring ED | The inversion-odd gap closes for ε ≤ 0.2, is critical near 0.3 (N·gap ≈ 0.72), and saturates for ε ≥ 0.4. The PDW ground states have T₂ = ±1 and are both I_s-even. |
| Z₃ control, R_b/a = 2.7, Δ = 5, N = 24 | Threefold degeneracy exact at ε = 0, 0.3 and 0.6 |
| Bulk melting at Δ = 2.5, g = 11.39 | ε_c ≈ 0.53–0.55 (the earlier estimate of 0.6 was wrong) |

**Novelty status.** Two independent searches found no study of a spacing-dimerized, single-level vdW chain, and none of its Z₃/Z₄ regimes. The web search was incomplete (arXiv full text was unreachable).

**Must read before committing (gate P0, by 29 Sept):**
- Chepiga, arXiv:2308.12838 (PRL 132, 076505)
- arXiv:2408.15965 (alternating-species chain)
- arXiv:2412.11201 and 2411.05494 (ladders)
- arXiv:2406.20093
- arXiv:2403.06649
- arXiv:2601.16275 and 2506.21963

**Not to be claimed:** "first", "simplest route to TCI", "emergent supersymmetry", a new phase for the shared-pair disordered regime, "exact for hundreds of atoms".

## 3. Decisive computations (4 cores; about 75–95 CPU-hours in total)

| Priority | Computation | Decides | Due |
|---|---|---|---|
| P0 | Literature gate (list above), read in full | Whether the Z₄/PDW split is already published | 29 Sept |
| P1 | Validate the constrained basis against full ED (N = 16, 20; periodic and open; R_b/a = 3.3; ε ∈ {0, 0.3, 0.5}; Δ ∈ {2, 3, 5}) | Go/no-go for the headline | 29 Sept |
| P2 | Z₄/PDW gap map: constrained periodic ED, N = 16–36 on a 12×8 (Δ, ε) grid, N = 40 on three cuts. Count degeneracies using T₂ and I_s quantum numbers; locate N·gap crossings | Extent of the PDW and its phase boundaries | 2 Oct |
| P3 | Z₂ line Δ_c(ε), gap ratio R, and c from E₀(N) at g ∈ {4.83, 11.39, 40} | Backbone boundary and universality | 2 Oct |
| P4 | Effective theory for both lines: leading order, dressed, and second-order Schrieffer–Wolff; classical Ω → 0 skeleton; BEG K/J table | Quantitative theory and validity map; classical TCI criterion | 3 Oct |
| P5 | Z₃ control: N = 24, 30, 36; ε from 0 to 0.7 | Odd-period protection | 3 Oct |
| P6 | ~~Infinite/finite DMRG repair~~ **Done 27 Sept.** Infinite DMRG passes its gate at 18 of 18 points (≤ 3.4 × 10⁻¹⁰ per atom against pinned 100/200-atom chains; ε = 0.55 critical, four-size fit). Bulk order curve at R_b/a = 1.5, Δ = 2.5 is in the draft. Fixes: Lanczos stopping rule, svd_min = 10⁻⁸, variational three-state pair basis at large V_s. See [IDMRG_NOTES.md](../validation/IDMRG_NOTES.md) | DMRG usable for large systems | done |
| P7 | Large-system confirmation: pinned open chains (last atom detuned by −2Δ; converge in 7–15 sweeps at 100–200 atoms) or odd chains; Calabrese–Cardy central charge on the Z₂ line and two Z₄ cuts. An odd number of pairs removes only the PDW wall, not the τ wall | c = 1/2 at hundreds of atoms (the professor's "larger systems" comment) | 6 Oct |
| P8 | Boundary-wall subsection (reuse the N = 12/16 data; add N = 20) | Methods warning; explains the earlier N = 8 finding | 5 Oct |
| P9 | Hardware table in lab units (QuEra Aquila, 4 µm minimum spacing) | Which predictions current hardware can test | 5 Oct |

**Gates:**
- 29 Sept: P0 and P1 pass.
- 2 Oct: P2 shows the predicted degeneracy counts beyond 12 pairs.
- 6 Oct: numerics frozen.
- 6–13 Oct: figures and manuscript.
- 13–15 Oct: self-audit, with every quoted number recomputed, then send to the professor.

**Fallback.** If P0 finds the Z₄ split published, or P1/P2 show the PDW is an artefact of the restricted basis, the paper becomes the Z₂ backbone plus the negative TCI result plus the boundary wall. Working title: *Dimerized Rydberg chains as quantum order–disorder ferroelectrics: pseudo-spin theory, phase boundary, and why geometry alone does not make them tricritical.* That is still an unstudied model with a validated quantitative theory, and most of it is already in hand.

## 4. Mapping onto the paper outline

| Heading | Content |
|---|---|
| Methods: algorithms for solving the Hamiltonian | Dense and Lanczos ED, symmetry-resolved ring ED, the exact pair-site regrouping, DMRG/iDMRG, observables, gap ratios and entanglement scaling |
| Methods: computational details in Qiskit | Pauli-operator construction and split-step circuits as independent small-N checks |
| Results: coupling constant | R_b/a selects the Z₂, Z₃ or Z₄ regime; existing blockade, defect-control and geometry results, condensed |
| Results: alternating coupling constants | The even/odd rule, the PDW phase, the Z₂ line and the effective theory, and the negative TCI test |
| Results: system size | ED cost versus DMRG; symmetry-resolved ED to N ≈ 40; DMRG to hundreds of atoms; even-chain boundary wall; interaction-cutoff errors |

## 5. Existing-material fixes (status, 27 September)

The historical `study/` package is frozen and checksum-verified. New code lives in `rydberg_chain/`, with `tests/` (66 tests) and `validation/`.

| Item | Status |
|---|---|
| 1. Infinite DMRG | Fixed and validated (P6 above) |
| 2. External benchmark | Blockaded-chain Z₂ point Δ_c/Ω = 0.65515 ± 0.00016 against the published 0.655, with c = 0.50 ± 0.01. Critical Ising chain exact to 10⁻¹³. Full vdW chain at R_b/a = 1.5: 1.0794 (open DMRG) against 1.0785 (ring ED) ([BENCHMARK_NOTES.md](../validation/BENCHMARK_NOTES.md)); the literature value still needs a full-text check |
| 3. ED speed | Vectorized builder, identical to the study Hamiltonian; N = 20 in 55 s with LOBPCG |
| 4. Consolidation and tests | `rydberg_chain/` package with pytest suite |
| 5. Requirements | Minimal root `requirements.txt`; study pins unchanged |
| 6. Notation | [NOTATION.md](NOTATION.md) |
| 7. Figures | Seven publication figures in `paper/figures/`, regenerated from saved data |
| 8. Docs tone | README and docs rewritten |
| 9. Cutoff study | Reframed as the reason all tails are kept |
| 10. Defect label | Labelled as an idealised interaction switch-off control |
| 11. Kagome | Marked as an exploratory appendix |

**Manuscript draft:** `paper/DRAFT.md` and `paper/DRAFT.docx` (see `paper/README.md`).

## 6. Ownership and disclosure

The student derives J_m, h, the dressing, the symmetry table and the gap-ratio method himself, and reproduces the key numbers from the maintained package. AI assistance in code, pilots and drafting is disclosed, in line with the professor's and the journal's policy.
