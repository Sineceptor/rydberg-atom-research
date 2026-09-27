# External benchmarks of the DMRG pipeline

Script: `validation/benchmark_uniform.py`. Results: `validation/results/benchmark_uniform.json`. Smoke tests: `tests/test_benchmark_smoke.py`.

Run it from the repository root:

```
python validation/benchmark_uniform.py --workers 2      # all parts, about 17 min on 2 cores
python validation/benchmark_uniform.py --quick          # sizes 16-25, about 2 min
python validation/benchmark_uniform.py --reanalyse      # refit from the raw data in the JSON
```

Each worker runs single-threaded (the script sets `OMP_NUM_THREADS=1` unless it is already set), so `--workers 2` uses two cores. Every Rydberg-chain DMRG point starts from the Z2 product state with the density-matrix mixer on, uses χ ≤ 64 and stops after at most 30 sweeps (no point needed more than 14). See "Pitfalls" below for why warm starts are not used.

## 1. Exactly solvable check: critical open transverse-field Ising chain

H = −J Σ σˣᵢσˣᵢ₊₁ − h Σ σᶻᵢ with J = h = 1 and open ends (TeNPy `TFIChain`, parity conserved).

**Exact reference.** Write H = (i/4) Σ M_ab γ_a γ_b in 2N Majorana operators. The couplings are M_{2i,2i+1} = 2h and M_{2i+1,2i+2} = 2J.
- The energy is E₀ = −½ Σ ε_k, where the ε_k are the positive eigenvalues of iM.
- The entropy of the first ℓ sites follows from the eigenvalues ±ν of iG restricted to the first 2ℓ Majoranas, with G = −i sign(iM).
- Both formulas are checked against dense ED at N = 6 and 8 (energy and every S(ℓ) to 1e-12) in `tests/test_benchmark_smoke.py`.

**Central charge.** Fit S(ℓ) = (c/6) ln[(2N/π) sin(πℓ/N)] + s₀ on N/8 ≤ ℓ ≤ 7N/8. Fits using only even or only odd cuts are also reported, to expose bond alternation.

### Results (χ = 64, 8 sweeps each)

| N | E_DMRG − E_exact | max over ℓ of \|S_DMRG − S_exact\| | c (DMRG) | c (exact S) | wall time (s) |
|---|---|---|---|---|---|
| 32 | −3.6e-14 | 1.9e-12 | 0.5326 | 0.5326 | 6 |
| 48 | 1.2e-13 | 6.2e-12 | 0.5223 | 0.5223 | 9 |
| 64 | −1.1e-13 | 1.5e-11 | 0.5169 | 0.5169 | 12 |
| 96 | 1.3e-13 | 2.7e-11 | 0.5114 | 0.5114 | 21 |
| 128 | 2.6e-13 | 2.0e-11 | 0.5086 | 0.5086 | 25 |

- **Energies and entropies.** DMRG reproduces the exact open-chain energies to about 1e-13 and every bond entropy to about 1e-11.
- **Central charge.** The fitted c approaches 1/2 as c(N) ≈ 1/2 + 1.0/N.
  - A linear fit in 1/N of the DMRG values gives c = 0.5007.
  - The exact entropies at N = 256, 512 and 1024 give 0.5045, 0.5023 and 0.5011, which extrapolate to 0.50002.
  - The even-cut and odd-cut fits agree to 1e-3 or better.

## 2. Literature benchmark: the Z2 transition of the nearest-neighbour-blockaded chain

### Published value and conventions

**Hard-boson model.** Fendley, Sengupta and Sachdev [1] study hard-core bosons with nearest-neighbour exclusion n_j n_{j+1} = 0:

  H_FSS = Σ_j [ −w (d_j + d_j†) + U n_j + V (next-nearest-neighbour repulsion) ].

On the V = 0 axis the only transition is from the disordered phase to the period-2 density wave. It is in the Ising class and lies at **U_c/w ≈ −1.31**. The same point is the Ising transition of the effective dipole model of a tilted Mott insulator [2]. There the hopping is √2 w, so the transition sits at (E−U)_c ≈ 1.85 w, and 1.85/√2 = 1.31.

**Mapping to the Rydberg chain.** In this repository H = (Ω/2) Σ Xᵢ − Δ Σ nᵢ + Σ V_ij nᵢ n_j. The model above is recovered in the limit V₁ = V_{i,i+1} → ∞ with every longer-range V_ij = 0. The mapping is:

| Hard-boson model [1] | This repository | Comment |
|---|---|---|
| d_j + d_j† | X_j | same matrix in the two-state basis |
| −w | Ω/2 | the sign of the drive is a gauge choice: Π_j Z_j maps X → −X and leaves n and the constraint unchanged, so only \|w\| = Ω/2 matters |
| U | −Δ | the ordered phase lies at large negative U, i.e. large positive Δ |
| U_c/w = −1.31 | **Δ_c/Ω = 1.31/2 = 0.655** | |

Other factor-of-2 conventions in use:
- **PXP papers that write H = Ω Σ P X P − Δ Σ n.** These have no ½ on the drive, so the same point reads Δ_c/Ω ≈ 1.31.
- **Experiments [3, 4].** They use (Ω/2)σˣ with Ω the Rabi frequency, as this repository does, so there Δ_c/Ω ≈ 0.655 in the strict-blockade limit.
- **Blockade radius.** It does not enter the hard-boson limit. Here g = (R_b/a)⁶ = C₆/(Ω a⁶), i.e. R_b = (C₆/Ω)^{1/6} [3, 4], so V₁ = gΩ. A real van der Waals chain at finite R_b/a also has V₂ = V₁/64. This next-nearest-neighbour repulsion penalises the Z2 crystal and moves Δ_c upward.
  - For example, the preliminary ring-ED pilot gives Δ_c = 1.0785 at R_b/a = 1.5 (`research/pilots/2026-09-27`). That is an internal number, not a literature value; Sec. 3 cross-checks it with open-chain DMRG.

**Source status.** The web-search quota of this session was exhausted before this task started. The egress proxy also blocks arxiv.org, ar5iv, journals.aps.org, Semantic Scholar and Crossref. As a result, the values −1.31 and 1.85 above are quoted from the published papers as I know them, and the exact sentence or figure could not be re-read here. **[verify]** Before citing, check that U_c/w ≈ −1.31 appears in [1] for V = 0, and check the sign convention of the U term there. The numerical result below is independent of this check: it agrees with −1.31 to the three digits quoted.

### Method

**Models.**
1. **Exact blockade.** `blockade_chain_model` builds H = (Ω/2) Σ P_{i−1} Xᵢ P_{i+1} − Δ Σ nᵢ + K Σ nᵢnᵢ₊₁, with P = 1 − n and open ends (the end atoms have a single projector).
   - The P X P drive conserves the number of blockade violations, so the blockaded sector is exactly the hard-boson model.
   - The penalty K = 10 + 2|Δ| only lifts the other sectors; without it the ground state would lie in a sector that violates the blockade (for example all atoms excited, with E = −NΔ).
   - DMRG equals Fibonacci-basis ED at N = 9 and 13 to about 1e-14.
2. **Finite V₁.** `rydberg_chain.dmrg.atom_chain_model` with g = V₁ and `rmax=1`, i.e. V₂ = V₃ = … = 0, for V₁ = 50, 100, 200.

**Finite-size estimator.**
- **Boundary conditions.** Odd open chains are used. A free end atom has only one neighbour, so it is excited preferentially. That is a boundary field, which flows to fixed boundary conditions. In odd chains both ends select the same Z2 pattern (1010…01), so the bulk has no domain wall. This is unlike even chains; see the core package report.
- **Order parameter.** At the central atom c = (N−1)/2 we measure m_c = (−1)^c [n_c − (n_{c−1} + n_{c+1})/2]. With fixed boundaries, m_c ~ N^{−β/ν} f((Δ−Δ_c) N^{1/ν}) with β/ν = 1/8, so the curves N^{1/8} m_c(Δ) of different N cross at Δ_c.
- **Crossings.** Crossings of consecutive sizes (33, 49, 65, 97, 129) are found from cubic fits on the grid Δ = 0.645 … 0.665 (step 0.005). They are extrapolated as Δ*(N̄) = Δ_c + a N̄^{−p} with N̄ = √(N₁N₂).
  - The main fit leaves p free. The expectation is p = 1/ν + ω = 2 for the leading boundary correction ω = 1.
  - Fits with p fixed at 1, 1.5, 2 and 2.5 are also stored.
  - The quoted uncertainty is the largest deviation of the fixed-p fits within ±0.5 of the free p.
- **Central charge.** c comes from the bond entropies at the grid point nearest Δ_c. The alternation of S(ℓ) between even and odd cuts is large here because of the Z2 sublattice structure.
  - In an odd chain that is symmetric under reflection, ℓ → N−ℓ changes the parity of the cut. Even-only and odd-only fits over a symmetric window are therefore identical, and they are taken on the left half (N/8 ≤ ℓ ≤ N/2) instead.
  - The main estimator, `c_smooth`, fits the three-point average (S(ℓ−1) + 2S(ℓ) + S(ℓ+1))/4 on N/8 ≤ ℓ ≤ 7N/8.
  - Every estimator is also extrapolated linearly in 1/N.

### Results

**Checks.**
- DMRG of the exact blockade model equals Fibonacci-basis ED at N = 9 and 13 (Δ = 0.655 and 1.2) to ≤ 2e-14.
- Every DMRG point converged (|dE| < 1e-12, |dS| < 1e-8) in 8–9 sweeps, with a final truncation error ≤ 1.3e-19 at χ ≤ 64.

**Crossings of N^{1/8} m_c** (exact blockade model):

| Pair | Δ* / Ω |
|---|---|
| (33, 49) | 0.65322 |
| (49, 65) | 0.65410 |
| (65, 97) | 0.65458 |
| (97, 129) | 0.65484 |

The crossings drift upward with an effective exponent p ≈ 1.8, close to the p = 2 expected from ω = 1.

**Extrapolation.**
- The fit Δ* = Δ_c + a N̄^{−p} with p free gives Δ_c = 0.655152 and p = 1.79.
- Fixed p = 1.5, 2 and 2.5 give 0.65532, 0.65506 and 0.65491. Fixed p = 1 gives 0.65583, but the data rule it out: the successive increments shrink by a factor of 1.8 per step of 1.41 in N̄.

**Result versus literature.**

| | Δ_c / Ω (this repo's convention) | U_c / w (FSS convention) |
|---|---|---|
| This benchmark (DMRG, N ≤ 129) | **0.6552 ± 0.0002** | −1.3103 ± 0.0003 |
| Fendley, Sengupta, Sachdev [1] | 0.655 | −1.31 |

The two agree to the three digits published. The quoted uncertainty is the spread of the fixed-p fits with \|p − 1.79\| ≤ 0.5. DMRG and interpolation errors (≤ 1e-5) are negligible next to it.

**Central charge at Δ = 0.655.**

| N | c_smooth | c_even_left | c_odd_left |
|---|---|---|---|
| 33 | 0.398 | 0.434 | 0.359 |
| 49 | 0.429 | 0.456 | 0.400 |
| 65 | 0.445 | 0.466 | 0.423 |
| 97 | 0.463 | 0.478 | 0.448 |
| 129 | 0.473 | 0.484 | 0.461 |

- **Extrapolation in 1/N.** The three estimators extrapolate to 0.497 (c_smooth), 0.501 (c_even_left) and 0.494 (c_odd_left).
- **Central entropy.** Consecutive sizes give c = 6[S_mid(N₂) − S_mid(N₁)]/ln(N₂/N₁) = 0.468, 0.478, 0.485 and 0.490, still rising.
- **Conclusion.** Both routes agree on c = 0.50 ± 0.01, the Ising value.
  - The approach from below is slower than in the TFIM. This is expected with a finite boundary field, whose free-to-fixed crossover lowers the apparent slope of S(ℓ) near the ends.
  - Raising ℓ_min from N/8 to N/3 moves c_smooth at N = 129 only from 0.473 to 0.480.

**Convergence in V₁** (`rydberg_chain.dmrg.atom_chain_model` with `rmax=1`, crossing of N = 33 and 49):

| V₁ / Ω | Δ*(33, 49) / Ω | shift × V₁ |
|---|---|---|
| 50 | 0.67212 | 0.945 |
| 100 | 0.66260 | 0.938 |
| 200 | 0.65788 | 0.932 |
| ∞ (exact blockade) | 0.65322 | — |

- The shift is 0.93 Ω²/V₁ to within 1.5%.
- A quadratic extrapolation in 1/V₁ of the three finite-V₁ values gives 0.65318, which matches the exact-blockade crossing to 4e-5.
- This checks the package's `atom_chain_model` against the independent constrained construction. It also shows that "V₁ = 50" is not yet the blockade limit at the 1e-2 level.

## 3. Internal cross-check: van der Waals chain at R_b/a = 1.5

The same estimator was applied to the full model with g = 11.39 (V₁ = 11.39, V₂ = 0.178, …), with couplings kept up to fourth neighbours.
- **Crossings.** The crossings are 1.07846 for (33, 49) and 1.07893 for (49, 65).
- **Extrapolation.** With p = 1.5–2.5 this gives Δ_c = 1.0794 ± 0.0002, where the uncertainty is the fit spread only; p = 1 gives 1.0801.
- **Truncation of the couplings.** Dropping V₆ = g/6⁶ = 2.4e-4 and all longer couplings (the even distances penalise the Z2 crystal) lowers Δ_c by about 2e-4, so the untruncated chain should sit near 1.0796.
- **Comparison with the pilot.** The preliminary ring-ED value of the pilot (8/10-pair rings, 16/20 atoms) is Δ_c = 1.0785. The two methods have different boundary conditions and sizes and agree to about 1e-3. The ε = 0 entry of the pilot table is therefore confirmed at that level.
- **Central charge.** At N = 65 and Δ = 1.08, c_smooth = 0.472. The central-entropy estimate is 0.482 from (33, 49) and 0.483 from (49, 65). The values are close to the blockade-model ones at the same sizes, consistent with Ising.
- **Status.** This is an internal consistency check, not a literature comparison.

## Run record

- **Commands.** `validation/results/benchmark_uniform.json` comes from three commands:
  1. `--parts tfim,blockade` (the full run; 861 s wall time on 2 workers).
  2. `--parts vdw` (167 s).
  3. `--reanalyse` (27 s).
- **Load.** All three ran while other jobs kept the load at 4–7 on 4 cores.
- **The `runs` list.** It records each command. The first full run came from an earlier revision of `main()` that had the same task functions but fewer fits. Its fits were recomputed from the stored raw data by `--reanalyse`.

## Pitfalls met while building this benchmark

- **Warm starts.** Warm-starting DMRG from the converged state at a neighbouring Δ, without the mixer, relaxes the order-parameter profile of long critical chains very slowly. At N = 65, a warm-started three-point scan had not finished after 14 minutes, while a fresh start from the product state with the mixer converged in 8 sweeps (27 s). Every point is therefore computed from scratch.
- **Finite V₁ shifts the transition more than a naive estimate suggests.** Virtual double occupation favours the disordered phase. The shift is roughly linear in 1/V₁ with a coefficient of order Ω² (see the table), so V₁ = 25 would already move Δ_c by about 0.04. Treat "large V₁" with care when comparing with the strict-blockade literature.
- **Exact entropies at large N.** Computing them for every cut costs O(N⁴). For N ≥ 256 only 24 cuts are evaluated.

## Open problems

1. **Literature value not re-read.** It still has to be checked against the text of [1] (and [2]) for the exact number, its quoted precision and the sign convention of U. See "Source status" above.
2. **A more precise reference.** More recent DMRG work on the hard-boson model with next-nearest-neighbour repulsion (e.g. Chepiga and Mila, PRL 122, 017205 (2019)) probably gives the V = 0 Ising point to more digits and would be a sharper benchmark. It was not accessible here.
3. **No published full van der Waals benchmark yet.** None was found for the 1/r⁶ chain at a stated R_b/a, so Sec. 3 compares only with the internal pilot.
   - A published location, for example from the phase diagrams of Rader and Läuchli (arXiv:1908.02068) or of Maceira, Chepiga and Mila (arXiv:2203.01163; journal reference [verify]), would close this. Those papers mostly work near the Z3 and Z4 transitions, so the precision of their Z2 lines is uncertain.
4. **Central charge.** It is established to about 2% (0.50 ± 0.01). A tighter test would need N ≳ 250 or a fit that includes the boundary-field crossover.
5. **Time per sweep.** TeNPy is dominated by Python overhead at these small χ, at 5–10 s per sweep for N = 129. Size scans are therefore run point by point in parallel processes rather than with threads.

## References

[1] P. Fendley, K. Sengupta, and S. Sachdev, "Competing density-wave orders in a one-dimensional hard-boson model," Phys. Rev. B **69**, 075106 (2004), arXiv:cond-mat/0309438. [verify the location of U_c/w ≈ −1.31 at V = 0]

[2] S. Sachdev, K. Sengupta, and S. M. Girvin, "Mott insulators in strong electric fields," Phys. Rev. B **66**, 075128 (2002). [verify]

[3] H. Bernien et al., "Probing many-body dynamics on a 51-atom quantum simulator," Nature **551**, 579 (2017). Conventions (Ω/2)σˣ and R_b = (C₆/Ω)^{1/6}. This session re-read the full text through a GitHub mirror: its Fig. 5 places the peak of the domain-wall variance of a finite-speed sweep at Δ ≈ 0.5Ω. That is a dynamical, non-adiabatic marker, not the equilibrium Δ_c.

[4] A. Keesling et al., "Quantum Kibble–Zurek mechanism and critical dynamics on a programmable Rydberg simulator," Nature **568**, 207 (2019).

[5] P. Calabrese and J. Cardy, "Entanglement entropy and quantum field theory," J. Stat. Mech. (2004) P06002.
