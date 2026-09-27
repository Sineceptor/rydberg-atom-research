# Methods and code map

This page describes the model, the observables, the numerical methods and which code implements each. Symbols follow [NOTATION.md](NOTATION.md).

## Repository map

| Location | Role | Change policy |
|---|---|---|
| `study/` | The completed finite-chain study (N ≤ 12): model, scripts, notebooks, `results/`, `figures/`, `reference_data/` | Frozen historical record. `python study/verify_delivery.py` checks every file against `study/SHA256SUMS.txt`. Do not edit; reproduce or extend in the maintained code instead. |
| `rydberg_chain/` | Maintained package used for all new calculations | Changes need passing tests |
| `tests/` | pytest suite for the package | Extend with every new feature |
| `validation/` | Benchmark scripts with saved JSON results (`validation/results/`) and [notes](../validation/BENCHMARK_NOTES.md) | Scripts record parameters, versions and machine load |
| `paper/` | Section drafts; `make_figures.py` writes `paper/figures/` and its manifest from saved data | Figures are regenerated, never edited by hand |
| `research/pilots/` | Exploratory pilot scripts and raw outputs | Numbers must be recomputed with the package before publication |
| `tools/check_current_model.py` | Short standalone check of the historical model | |

### Package modules

| Module | Main functions | Purpose |
|---|---|---|
| `geometry` | `chain_positions(n, eps)`, `coupling_matrix(x, g, rmax)`, `strong_weak_couplings`, `reflect` | Positions of uniform or dimerized chains and V_ij = g/\|x_i − x_j\|⁶, optionally truncated |
| `ed` | `build_hamiltonian`, `chain_hamiltonian`, `lowest_states`, `basis_occupations` | Vectorized sparse Hamiltonian in the 2ᴺ basis; lowest eigenpairs by dense `eigh` (dimension ≤ 1024), `eigsh` or LOBPCG, with residuals |
| `observables` | `site_occupations`, `mean_density`, `bond_probabilities`, `staggered_order_squared`, `pair_order`, `pair_state_probabilities`, `entanglement_entropy`, `connected_correlations` | Expectation values of full state vectors |
| `effective` | `tfim_couplings`, `transverse_field`, `critical_eps` | Pseudo-spin (transverse-field Ising) couplings for one excitation per strong pair |
| `dmrg` | `atom_chain_ground_state`, `pair_chain_ground_state`, `dmrg_options` | Finite two-site DMRG (TeNPy) with one atom per site, or one strong pair per four-state site |
| `infinite` | `infinite_ground_state`, `pinned_chain_ground_state`, `finite_bulk_energy` | Infinite DMRG on pair sites and domain-wall-free open chains with a pinned end atom |
| `plotstyle` | `use`, `subplots`, `panel_labels`, `save` | Shared publication figure style |

`study/` modules keep their original roles: `rydberg_model.py` (model, dynamics, Qiskit equivalents), `extension_model.py` and `run_extension.py` (interaction-cutoff scan), `run_august_experiments.py` (blockade, switch-off controls, geometry, scaling, circuits), `phase_diagram_bridge.py` (N = 8 map), `parity_diagnostics.py` (odd/even analysis), `circuit_controls.py`, `kagome_pilot.py` and `kagome_shape_control.py` (exploratory appendix), `run_project.py` and `run_support.py` (managed reruns into `study/runs/`).

## Model

$$H=\frac{\Omega}{2}\sum_iX_i-\Delta\sum_i n_i+\sum_{i<j}\frac{C_6}{|x_i-x_j|^6}n_i n_j,\qquad n_i=|1_i\rangle\langle1_i|=\tfrac12(I-Z_i).$$

With ħ = Ω = a = 1 the only couplings are δ = Δ/Ω and g = C₆/(Ωa⁶) = (R_b/a)⁶, plus the dimerization ε for the alternating chain. Report all three, and τ = Ωt for dynamics.

**From geometry to matrix.** `geometry.chain_positions` places the atoms (gap 1 − ε after even i, 1 + ε after odd i). `geometry.coupling_matrix` gives V_ij for every pair, or for |i − j| ≤ rmax. `ed.build_hamiltonian` stores the diagonal (−Δ times the number of excitations plus the occupied-pair energies) and the drive, which flips one bit with matrix element Ω/2. The builder is vectorized; at N = 20 it takes about 2 s for 2.2 × 10⁷ stored entries. It reproduces `study/rydberg_model.hamiltonian` exactly (see [VALIDATION.md](VALIDATION.md)).

**Exact regrouping into pair sites.** Atoms (2k, 2k+1) form a four-state site |n_L n_R⟩. The strong-bond energy V_s n_L n_R becomes an on-site term and only inter-pair couplings enter the matrix product operator. This is the same Hamiltonian, not an approximation, unless `mmax` truncates couplings beyond a pair separation.

## Observables

For a state ψ with basis weights |a_s|²:

$$\langle n_i\rangle=\sum_s|a_s|^2 n_i(s),\qquad \bar n=\frac1N\sum_i\langle n_i\rangle,\qquad C_{ij}=\langle n_in_j\rangle-\langle n_i\rangle\langle n_j\rangle,$$

$$Q=\frac1N\sum_i(-1)^i(2n_i-1),\qquad \tau_k=n_{2k}-n_{2k+1},\qquad G=E_1-E_0.$$

Either perfect alternating string has Q² = 1 for both parities of N. For independent sites with uniform excitation probability p, ⟨Q²⟩ = 4p(1−p)/N + (2p−1)²(N mod 2)/N², so a nonzero ⟨Q²⟩ in a short chain does not by itself indicate order; use configurations, correlations and size dependence alongside it.

For a pure state and a set of atoms A, ρ_A = Tr_B|ψ⟩⟨ψ| and S_A = −Tr ρ_A log₂ρ_A measures bipartite entanglement in bits (Bell pair: 1; product state: 0).

For dynamics under constant H from a stated initial state, |ψ(t)⟩ = e^{−iHt}|ψ(0)⟩. A ground state and the state reached after a sudden quench are different objects.

## Numerical methods

**Exact diagonalization.** `ed.lowest_states(H, k)` returns (energies, vectors, residuals). It uses dense `eigh` up to dimension 1024, and otherwise SciPy `eigsh` (implicitly restarted Lanczos, `which='SA'`, tol 10⁻¹², fixed seed) or, with `method='lobpcg'`, Jacobi-preconditioned block LOBPCG with a fallback to `eigsh` if the residual is too large. In the N = 16–20 benchmark LOBPCG solved 4–7 times faster than `eigsh` and agreed with it to 10⁻¹²; `eigsh` at N = 20 took about 6 minutes on a shared machine. Close levels need several eigenpairs and more than one start vector.

**Finite DMRG.** Two-site DMRG (TeNPy 1.1) with a density-matrix mixer for the first sweeps, bond dimension χ (default 64) and convergence when both the energy change (10⁻¹²) and the mean bond-entropy change (10⁻⁸) between sweeps are small. Report χ, sweeps, truncation error, the convergence flag and the energy variance ⟨H²⟩ − ⟨H⟩². The pair-site form is preferred for dimerized chains; `mmax = 8` shifts E₀ at N = 20 by 2.4 × 10⁻⁷.

**Even open chains.** In the ℤ₂ regime the boundary fields favour exciting both end atoms, which forces one 00 domain wall into an even chain. The wall is delocalized, the gap closes like a particle in a box, and DMRG moves it slowly: at ε = 0.3, g = 11.39, Δ = 2.5 the N ≥ 50 runs did not meet the entropy criterion in 30 sweeps. Bulk quantities should come from odd chains, a pinned end (`infinite.pinned_chain_ground_state`, via `local_delta`) or infinite DMRG.

**Infinite DMRG.** `infinite.infinite_ground_state` runs iDMRG on a two-pair unit cell of pair sites. Its acceptance gate compares the energy per atom with the bulk estimate [E(N₂) − E(N₁)]/(N₂ − N₁) from pinned open chains (N₁ = 100, N₂ = 200), requiring agreement to 10⁻⁶ ([`validation/idmrg_validation.py`](../validation/idmrg_validation.py)).

**Dynamics and Qiskit.** The historical study evolves states with SciPy `expm_multiply`, and independently with Qiskit's Pauli-operator representation; split-step circuits approximate the same evolution and were checked for convergence with step size at both signs of Δ. All are classical calculations; nothing ran on quantum hardware.

**What to record.** Atom count and parity, positions or ε, boundary conditions, g, Δ, interaction range, the eigenstate or initial state, time grid, solver and tolerances, χ and sweeps for DMRG, library versions and the source file of each figure. Cite the [SciPy `eigsh`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html), [SciPy `lobpcg`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html), [SciPy `expm_multiply`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html), [TeNPy](https://tenpy.readthedocs.io/) and [Qiskit quantum information](https://quantum.cloud.ibm.com/docs/en/api/qiskit/quantum_info) documentation for the methods actually used.

## Checks to run

Hand-built two-atom matrices, the g = 0 limit, agreement between independent constructions and solvers, conservation of norm and energy in dynamics, reflection symmetry, eigenpair residuals and, for DMRG, the energy variance and comparison with ED where both are possible. [VALIDATION.md](VALIDATION.md) lists what has been checked and the observed errors; [RESULTS_AND_LIMITS.md](RESULTS_AND_LIMITS.md) states what the results do and do not support.
