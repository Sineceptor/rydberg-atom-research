# Validation

What has been checked, how, and with what observed error. Status as of 27 September 2026.

## Maintained package (`rydberg_chain/`)

Run `python -m pytest` from the repository root with `OMP_NUM_THREADS=2` (see [README](../README.md)). The 43 tests of the core modules (`test_geometry_ed.py`, `test_effective.py`, `test_dmrg.py`) passed in 43 s on 27 September. `test_infinite.py` (13 tests) and `test_benchmark_smoke.py` (10 tests) cover the infinite-DMRG module and the benchmark references; the full suite of 66 tests passes.

| Check | Observed agreement |
|---|---|
| ED Hamiltonian against `study/rydberg_model.hamiltonian`, six random cases including cutoffs | identical matrices; coupling matrices to 10⁻¹⁴ relative |
| Hand-written two-atom matrix, `local_delta`, non-interacting density formula | exact |
| Dense `eigh`, `eigsh` and LOBPCG on the same matrix | 10⁻¹¹ |
| Saved uniform reference N = 12, g = 8, Δ/Ω = 3 (E₀ = −18.399305357538882) | E₀ to 4 × 10⁻¹⁴; density, ⟨Q²⟩, entropy and gap to ≤ 10⁻⁹ |
| Reflection of odd dimerized chains (N = 7, 9) | spectrum unchanged; occupations mirrored |
| Atom-site DMRG against ED, N = 10–12, ε = 0 and 0.3, with `local_delta` and `rmax` | energies ≤ 10⁻⁹ |
| Pair-site DMRG against ED, N = 12–14, ε = 0, 0.3, 0.6, with `mmax` | energies ≤ 10⁻⁸; observables 10⁻⁵ (default settings), 10⁻⁷ with `max_S_err = 1e-10` |
| Effective model: J_m, B_m against a fit to direct 1/r⁶ pair energies; h against the exact single-pair doublet splitting; Pfeuty limits of the Ising chain | J_m, B_m consistent; h to 10⁻³ |

## Benchmarks (`validation/`)

**Solver cost and accuracy** (`scaling_benchmark.py` → `results/scaling_benchmark.json`; ε = 0.3, g = 11.39, Δ/Ω = 2.5, even open chains; figure `paper/figures/fig07_solver_scaling`):

- ED with `eigsh` at N = 16, 18, 20: 6.6, 40 and 362 s; peak memory 0.11–0.56 GB. LOBPCG: 1.8, 8.7 and 55 s, with the two lowest energies equal to the `eigsh` values within 1.1 × 10⁻¹².
- Pair-site DMRG (χ = 64) at N = 20 reproduces the ED ground-state energy to 4 × 10⁻¹³ with all couplings; keeping couplings up to eight pair separations shifts it by −2.4 × 10⁻⁷.
- At N = 50 DMRG is close to converged but meets neither criterion within 30 sweeps (last-sweep energy change 1.6 × 10⁻⁹, variance 3.7 × 10⁻⁹). The N = 100 and 200 runs are not converged ground states (variance 1.9 × 10⁻⁴ and 4.3 × 10⁻⁴) and give only the cost per sweep (24 and 33 s). The slow mode is the delocalized domain wall of the even open chain.
- The machine was shared with other jobs (load average up to 9); times are upper bounds and vary by about a factor of two.

**External benchmarks** (`benchmark_uniform.py`, [notes](../validation/BENCHMARK_NOTES.md)): the critical open transverse-field Ising chain against its exact free-fermion solution, and the ℤ₂ transition of the nearest-neighbour-blockaded chain against the published value.

**Infinite DMRG gate** (`idmrg_validation.py`, [notes](../validation/IDMRG_NOTES.md)): energy per atom within 10⁻⁶ of the bulk estimate from pinned open chains of 100 and 200 atoms, with consistent order parameters. All 18 points pass (`results/idmrg_validation.json`): 17 on the two-size estimate, with differences ≤ 3.4 × 10⁻¹⁰ per atom, and the near-critical point ε = 0.55 (g = 11.39, Δ/Ω = 2.5, ξ ≈ 127a) only on a four-size fit (N = 100–800), at −4.1 × 10⁻⁷. Points with ε ≥ 0.45 at this g and Δ use the variational three-state pair basis (error ≤ 1.6 × 10⁻⁹ per atom against 10-atom ED).

## Publication figures (`paper/make_figures.py`)

The figures are drawn from the saved CSV/JSON files without recomputing any physics. The manifest `paper/figures/FIGURES.md` records the SHA-256 prefix of every source file. Numbers quoted in the draft captions are computed from the same data when the script runs, and the script asserts the qualitative statements it makes (for example, that every odd chain's cutoff error exceeds both even neighbours). Regenerating the figures reproduces the PNG and PDF files byte for byte.

## Historical study (`study/`)

The study is a frozen record. `python study/verify_delivery.py` checks all 104 delivered files against `study/SHA256SUMS.txt`; it verifies stored bytes, not conclusions.

- **Regression.** The bundled reference recomputes all 62 historical N = 8 rows; the maximum differences are 1.8 × 10⁻¹³ (density), 4.7 × 10⁻¹² (order) and 1.2 × 10⁻¹³ (gap), below the required 10⁻⁸ ([reference_validation.json](reference_validation.json)).
- **Model check.** `tools/check_current_model.py` tests a hand-built two-atom Hamiltonian, the non-interacting density, dense against sparse energies and gaps, and split-step circuits at both signs of Δ, writing its outputs to a temporary folder.
- **Saved checks.** `extension_summary.json` (6,561 coarse and 312 refinement states, eigenpair residuals ≤ 5.5 × 10⁻¹⁰), `parity_summary.json` (dense and independent-seed repeats, all below 10⁻⁸), `circuit_controls_summary.json` and `summary.json` (Qiskit against SciPy states to 10⁻¹⁴; norm and energy conserved to 5 × 10⁻¹⁵).
- **Notebooks.** Both notebooks ran from fresh kernels without errors on 26 September (Foundations: 16 code cells, 7 figures; Finite Chain Study: 15 code cells, 6 figures). Hashes of the 88 retained data, figure and reference files matched before and after execution.

The larger historical scans were not all rerun when the study was frozen; their data were compared byte for byte with the corrected baseline.
