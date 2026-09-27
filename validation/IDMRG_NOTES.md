# Infinite DMRG for the dimerized chain: diagnosis, fix and validation

Code: `rydberg_chain/infinite.py`. Tests: `tests/test_infinite.py` (12 tests, about 30 s). Gate script: `validation/idmrg_validation.py`. Results: `validation/results/idmrg_validation.json`.

```
OMP_NUM_THREADS=1 python validation/idmrg_validation.py --workers 2   # all 18 points
python validation/idmrg_validation.py --table                         # Markdown table of saved results
```

Conventions are those of the package: Ω = 1, mean spacing a = 1, g = (R_b/a)⁶, energies per atom, entropies in bits, ξ in units of a.

## 1. The failure

The pilot iDMRG (`research/pilots/2026-09-27/scripts/dimer_model.py`) used pair sites, a two-pair unit cell, χ = 64, TeNPy's default Lanczos, the default infinite-MPS mixer (50 sweeps) and svd_min = 1e-12. At g = 11.39, Δ = 2.5 it gave −1.219414 per atom at ε = 0.3 and then −1.201444, −1.181529 and −0.798889 at ε = 0.4, 0.5 and 0.55. The two pairs of the unit cell had different ⟨n_L − n_R⟩: −0.756 and −0.924 at ε = 0.4. An earlier atom-site iDMRG reached +2396 per atom at ε = 0.85.

## 2. Controlled experiments

All runs use g = 11.39 and Δ = 2.5. Scripts are in the session scratchpad; the numbers are quoted here.

| # | Question | Result |
|---|---|---|
| 1 | Is the infinite MPO wrong, e.g. couplings that cross the unit cell? | No. For product states in unit cells of 1–4 pairs at ε = 0.3, 0.55 and 0.85 (m_max = 6, 10, 3), ⟨H⟩ per site equals the directly summed classical energy to ≤ 6e-11, which is 1e-16 relative to V_s = 1e6. For random product superpositions, which test the Ω terms, the agreement is ≤ 1e-12. |
| 2 | Is the energy estimator wrong? | Not in the pilot. Reproducing it at ε = 0.4 gives a growth estimate of −1.2014442 and ⟨H⟩ of the final state of −1.2014389. The *state* is wrong: the correct value is −1.2208317. The run was unconverged after 110 sweeps, with an energy change of −0.036 per sweep. |
| 3 | Is the local eigensolver converged? | No. See below. |
| 4 | Do tiny Schmidt values matter? | Yes. See below. |
| 5 | What remains once 3 and 4 are fixed? | The spectral width at large V_s. See below. |
| 6 | Is TeNPy's transfer-matrix energy reliable here? | No. See below. |

**Experiment 3: the local eigensolver.**
- Pilot state after 20 sweeps, two-site effective Hamiltonian on bond (0, 1).
- ε = 0.3 (dimension 21904):
  - TeNPy's default Lanczos gives −901.221791; ARPACK gives −901.222737, an error of 9.5e-4.
  - Over the run, 53% of the updates stopped after 2 Krylov vectors and 47% hit N_max = 20.
- ε = 0.4 (dimension 65536):
  - The environment already holds E = +7832; the |11⟩ weight has leaked in.
  - The default Lanczos is off by 7.0.
  - N_max = 100 leaves 5e-3; N_max = 400 is exact.
- Cause: TeNPy stops when (RitzRes/gap)² < P_tol, with the gap taken from the Krylov Ritz values. |11⟩ states near V_s put the second Ritz value far up, so the test passes almost at once. P_tol is also raised automatically to as much as 1e-4 from the truncation error.

**Experiment 4: tiny Schmidt values.**
- ε = 0.4, fixed Lanczos, mixer off after sweep 20.
- svd_min = 1e-12:
  - The state is converged at sweep 20 (χ = 12, S_min = 4e-8).
  - One sweep after the mixer stops, noise-level Schmidt values (1e-12) enter and χ grows 12 → 34 → 48 → 64.
  - The two-site guess error 1 − |⟨guess|θ⟩| grows 5e-8 → 3e-6 → 2e-2, and Lanczos hits N_max.
  - A 150-sweep run shows growth-energy jumps of ±21 and ±251 per sweep.
- Cause: iDMRG builds the next guess with the inverse Schmidt values, so noise is amplified by up to 1e12.
- svd_min = 1e-8: the state is stable for 70 sweeps (χ = 12, two Lanczos vectors per update, energy constant at −1.2208316815135). This holds even with TeNPy's default Lanczos, so at ε = 0.4 this was the decisive fix.

**Experiment 5: spectral width.**
- ε = 0.55, V_s = 1372, svd_min = 1e-8, near the Z₂ transition.
- Default Lanczos: diverges. The growth energies range from −0.6 to −2043 and the final growth value is +275 per atom.
- Lanczos with N_max = 100: every update from sweep 2 on hits N_max; the guess error is O(1); the energy wanders between −1.12 and −1.26.
- The same physics with the |11⟩ level removed (the dressed basis of Section 4) converges monotonically: −1.2276829 → −1.2283418 within 70 sweeps, with 2–40 Lanczos vectors per update. The cause is therefore the spectral width ~2V_s, not the physics.
- Finite DMRG with 4-state sites is affected too. At ε = 0.7 (V_s = 1.6e4), a 64-atom chain does not converge in 80 sweeps and stays 1.1e-9 above the dressed-basis energy.

**Experiment 6: the transfer-matrix energy.**
- ε = 0.6, χ = 32. The growth estimate, a term-by-term sum and TeNPy's power-method ⟨H⟩ agree to 1e-13 at −1.2345684269.
- `MPO.expectation_value_TM` gives −1.2346110 in the dressed basis. In the full basis it gives −1.2345684267 at tol = 1e-10 but −1.2348276 at tol = 1e-12.
- In the first validation run it also produced the only large "disagreements", of 3e-7 to 1e-3.

## 3. Root cause

The pilot failure has three layers. The first two are numerical settings; the third is intrinsic to the exact pair-site formulation at large V_s:

1. **Unconverged Lanczos.** The on-site energy V_s = g/(1−ε)⁶ makes the local spectrum ~2V_s wide, and TeNPy's residual/gap stopping test is fooled by it.
2. **Inverted noise.** svd_min = 1e-12 keeps Schmidt values at the noise floor, and iDMRG inverts them.
3. **Slow Lanczos near the transition.** Even when Lanczos is not stopped early, a width of 10³–10⁶ against gaps of 10⁻² makes it too slow once the state is critical or nearly so.

A fourth problem is independent of these: TeNPy's transfer-matrix energy `expectation_value_TM` is unreliable for these MPOs. It would have hidden or faked problems regardless of 1–3.

The following were ruled out as causes:
- The mixer, at ε ≤ 0.4. The TeNPy default for infinite MPS (50 sweeps with slow decay) only delays convergence.
- The initial state.
- The unit-cell size.
- The MPO.

## 4. The fix (`rydberg_chain.infinite`)

- **Lanczos (`lanczos_params`).**
  - The residual test is disabled (P_tol = ∞); Lanczos stops when the ground-state Ritz value changes by less than 1e-12, or at N_max = 100.
  - Full reorthogonalization.
  - No automatic P_tol/E_tol adaptation.
- **Truncation, mixer and sweeps (`idmrg_options`).**
  - svd_min = 1e-8; the discarded weight per bond is below χ·1e-16.
  - The mixer runs for 10 sweeps, followed by at least 20 sweeps without it.
  - Convergence is judged every 10 sweeps: |ΔE| < 1e-12 per sweep (`energy_converged`), plus |ΔS| < 1e-9 (`converged`).
  - At most 120 sweeps. Near a critical point the entropy criterion is never met. After ~150 sweeps the environment energy (~10³) puts the absolute Lanczos tolerance at the rounding level, and one run at ε = 0.55 collapsed after 310 sweeps. The energy has converged long before that point.
- **Dressed basis for large V_s (`dressed_pair_site`).** Each pair is restricted to the three lowest eigenstates of the single-pair Hamiltonian (Ω/2)(X_L + X_R) − Δ(n_L + n_R) + V_s n_L n_R. Only the dressed |11⟩ level near V_s is dropped. Every operator is replaced by Wᵀ O W, so:
  - the MPS is exactly a state of the full Hilbert space;
  - all observables, including P(|11⟩), refer to it;
  - the energy is a variational upper bound;
  - the on-site physics is exact to all orders in Ω.

  The error comes from the inter-pair coupling of the dropped level, about Ω² V_w² / (8 V_s³) per atom (`truncation_error_estimate`). Exact diagonalization of 8 atoms gives 1.0e-8, 1.4e-9, 1.8e-10, 1.9e-11, 1.8e-12 and 3e-13 per atom at ε = 0.4, 0.45, 0.5, 0.55, 0.6 and 0.7, against estimates of 2.0e-8, 2.8e-9, 3.2e-10, 3.2e-11, 2.7e-12 and 1.5e-14. At ε ≥ 0.7 the comparison is limited by the roundoff of the full ED.

  `truncate="auto"` (the default) uses the dressed basis when the estimate is below 1e-8 per atom and the dropped level lies more than 10(|Δ| + Ω) above the others. At g = 11.39, Δ = 2.5 that means ε ≥ 0.45. Such results carry `obs['truncated'] = True` and `obs['truncation_error']`.

  This replaces the suggested "project out |11⟩ and add the second-order correction". That option neglects terms of order Ω²V_w/V_s², which are larger, and it is not variational.
- **Energy (`energy_per_atom`).** ⟨H⟩ is computed with TeNPy's power method, which sums every MPO term that starts in the unit cell. It is cross-checked against `energy_from_terms`, a sum over local expectation values and correlation functions that is independent of the MPO.
- **Starts and selection (`infinite_ground_state`).**
  - Default starts: |10⟩ⁿ (τ = +1) and |01⟩ⁿ (τ = −1). A third start, 'sym' = ((|10⟩ + |01⟩)/√2)ⁿ, is opt-in.
  - A candidate is *valid* when `energy_per_atom` and `energy_from_terms` agree to 1e-9.
  - The kept state is the valid candidate with the lowest energy (variational choice).
  - The I_s mirror maps the τ = +1 start onto the τ = −1 start. `mirror_check` reports |E_L − E_R| and max|τ_L + τ_R|; `candidates` reports convergence and the within-cell order spread of every start.
  - The 'sym' start is unreliable in the ordered phase, where it tends to a cat of both crystals. At g = 4.83, Δ = 1.5, ε = 0.2 it ended at +5.83 per atom, with Lanczos at N_max in 98% of the updates. It was never the lowest.
- **Finite reference (`pinned_chain_ground_state`, `finite_bulk_energy`, `bulk_energy_fit`).**
  - Open chains with the last atom's detuning lowered by 2Δ. This removes the forced 00 domain wall that kept N = 100–200 pair DMRG from converging.
  - Same Lanczos settings; svd_min = 1e-10; the same local basis as the iDMRG run.
  - The shift is a boundary term and cancels in [E(N₂) − E(N₁)]/(N₂ − N₁).
  - Pinned N = 100 and 200 chains converge in 7–15 sweeps. At ε = 0.3 this took 238 s, against 812 s with the default Lanczos and svd_min = 1e-12; the two bulk estimates agree to 1e-13.
  - `bulk_energy_fit` fits E(N) = eN + b + c/N (+ d/N²) when three (four) sizes exist, for critical points.

## 5. Validation gate

Gate (`validation/idmrg_validation.py`) at 14 values of ε ∈ [0, 0.8] at g = 11.39, Δ = 2.5, plus 4 other (g, Δ, ε). A point passes when all of the following hold:
- **Energy.** |e_inf − e_fin| < 1e-6 per atom, where e_fin = [E(200) − E(100)]/100 from pinned chains with the same m_max = 8, χ = 64 and local basis. Near a critical point (iDMRG ξ > 20a) e_fin is replaced by the fit over N = 100, 200, 400, 800; this is flagged.
- **Validity.** The kept state is valid.
- **Translation consistency.** ⟨n_L − n_R⟩ is the same on both pairs of the cell (spread < 1e-6).
- **Mirror consistency.** The τ = ±1 starts are I_s images: |ΔE| < 1e-6 and |τ_L + τ_R| < 1e-3.

The three starts L, R and 'sym' were all run. Entropies are in bits and ξ is in units of a. S_weak and S_strong are the cuts between pairs and inside a pair. e_fin − e_inf > 0 means the finite estimate lies above.

VALIDATION_TABLE

VALIDATION_DETAILS

## 6. Open problems and fallback

- **Critical point ε = 0.55 (Δ = 2.5).** This is the only point that fails the gate as literally specified, with two chain lengths. [E(200) − E(100)]/100 is 6.6e-6 above e_inf, because a critical open chain carries a Casimir term that the two-size difference does not remove. [E(400) − E(200)]/200 is still 2.3e-6 above. CRITICAL_FIT

  On the infinite side the energy is converged in χ: −1.2283418318551, −1.2283418318451 and −1.2283418292830 at χ = 32, 64 and 96. The χ = 96 run had not converged after 310 sweeps (dE = 2e-8 per sweep). At criticality ⟨n_L − n_R⟩ and ξ are set by χ: |τ| = 6e-7, 2e-3 and 1.3e-2 and ξ = 118, 123 and 143 a at χ = 32, 64 and 96. They are not physical values; critical properties need finite-size or finite-entanglement scaling.
- **Fallback.** Use the finite reference when an iDMRG result fails its own checks: `valid` false, `energy_converged` false, a large `order_spread`, or a mirror mismatch.
  - In gapped phases: `finite_bulk_energy` (pinned chains, N = 100 and 200). It agrees with iDMRG to ≤ 3.4e-10 at every gapped point here, and to ≤ 1.6e-12 at all but ε = 0.6.
  - Near criticality: add N = 400 and 800 and use `bulk_energy_fit`.
- **The dressed basis is an approximation.** It is variational and flagged, and is used here for ε ≥ 0.45 at g = 11.39, Δ = 2.5, and for (11.39, 4, 0.6) and (24, 3.5, 0.45). The gate compares with finite chains in the same basis, so it validates the solver. The truncation error is validated separately by 10-atom ED at every dressed point: ≤ 1.6e-9 per atom (ε = 0.45), ≤ 2.1e-10 for ε ≥ 0.5 at Δ = 2.5, and 7.9e-10 at (24, 3.5, 0.45).
- **4-state finite DMRG at large V_s** (`dmrg.pair_chain_ground_state`) has the same Lanczos problem. At ε = 0.7 a 64-atom chain did not converge in 80 sweeps and stayed 1.1e-9 above the dressed-basis energy. `dmrg.py` is outside this task and is unchanged. For ε ≳ 0.5 at g = 11.39, use `infinite.pinned_chain_ground_state(..., truncate=True)`, or pass `lanczos_params()` to TeNPy.
- **`MPO.expectation_value_TM` should not be used** for these MPOs anywhere in the project.
- **Interaction cutoff.** All energies use m_max = 8 pairs. The classical tail beyond it is +1.3e-7 to +7.8e-7 per atom (`tail_correction` in the JSON); add it before comparing with full-range calculations.
- **ED references.** The N = 16 open-chain ED values quoted in the task (−1.25057 … −1.27203) are not bulk references. They include an O(1/N) edge energy: at ε = 0.3, 16 × (−1.25022 + 1.21941) ≈ −0.49 per chain. The ring-ED value −1.2284 at ε = 0.55 (research plan) agrees with e_inf + tail = −1.2283415 to the quoted digits.
- **Unit cell.** Two pairs suffice for the Z₂ crystal and the disordered regime; the within-cell spread of ⟨n_L − n_R⟩ is ≤ 1.1e-8 at every point. Phases that double the pair cell, such as the Z₄ or pair-density-wave regime at R_b/a ≈ 3.3, need `cell_pairs ≥ 4` and were not validated.
- **'sym' start.** It is kept opt-in because it can diverge in ordered phases (Section 4). A general safeguard against a Z₂ start diverging has not been built. Diverged runs so far show up as large energies, `energy_converged = False` and within-cell order spreads of 0.04–0.46, so check these flags.
- **Cost** (one core of a shared machine):
  - iDMRG with three starts takes 5–240 s per point; the diverging 'sym' run at g = 4.83 took 1500 s.
  - Pinned N = 100 + 200 chains take 80–1100 s.
  - At ε = 0.55 the N = 400 chain took ~30 min and the N = 800 chain N800_TIME.
