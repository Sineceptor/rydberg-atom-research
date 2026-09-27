<!--
paper/sections/02_methods.md: Methods ("Algorithms for Solving the Hamiltonian", "Computational Details of Modeling in Qiskit", and the software package).
Status: draft of 27 September 2026, for Ethan to edit before moving it into the shared Google Doc.
Addresses the professor's comments 2 (what is in the package and how it works), 3 (what solves the Hamiltonian, with key equations) and 5 (how easy it is to go to larger systems).

Conventions used in this file
- Equations are numbered (2.1), (2.2), ... by section, matching the Results files; Eq. (1) is the Hamiltonian in the Introduction.
- Figures and tables are labelled M1, M2, ... because 03_results_coupling.md already uses Figures 2-7 on the assumption that the Methods have no figures. At merge, Figures M1-M3 become Figures 2-4 (or move to an appendix) and the Results figures shift accordingly.
- Plain integers in brackets, e.g. [3], are citations to the list at the end of this file. Renumber on merge.
- Image paths are relative to the repository root. Run pandoc from the root, or pass --resource-path=.
- "Preliminary" marks anything computed with the exploratory scripts in research/pilots/2026-09-27/; each such item names the computation (P1-P9, docs/RESEARCH_PLAN.md) that must validate it.
- Every quoted number has a "source" comment directly after it. "Scratch check, 27 Sept" means a short calculation run while drafting this section; its script and output are in the drafting session's scratchpad, not in the repository, and must be promoted into validation/ before submission (see the open items at the end).
- Section cross-references to the Results (Sec. 3.1-3.3) follow the outline order: coupling constant, alternating couplings, system size.
-->

## 2. Methods: Algorithms for Solving the Hamiltonian

All results in this paper are classical numerical calculations for the ideal, fully coherent Hamiltonian of Eq. (1); none was run on quantum hardware. We represent the same Hamiltonian in three ways: as a sparse matrix in the occupation basis for exact diagonalization (ED), as a matrix product operator for the density-matrix renormalization group (DMRG), and as a sum of Pauli strings for independent checks in Qiskit. Figure M1 shows how these representations, the solvers and the observables fit together, and Sec. 2.8 describes the code. Calculations that rely on exploratory pilot scripts rather than on the maintained, tested package are labelled preliminary, together with the validation step that must precede their use.

[FIGURE PENDING: workflow schematic of the calculation. Left to right: (i) inputs N, ε, g, Δ/Ω, boundary condition and interaction cutoff; (ii) `geometry.chain_positions` and `geometry.coupling_matrix` giving $x_i$ and $V_{ij}$; (iii) three representations of the same H: sparse CSR matrix (`ed.build_hamiltonian`), matrix product operator (`dmrg.pair_chain_model` / `dmrg.atom_chain_model`), Pauli sum (`study/rydberg_model.qiskit_hamiltonian`); (iv) solvers: dense `eigh`, ARPACK Lanczos or LOBPCG; two-site DMRG (TeNPy); `expm_multiply`, `HamiltonianGate` and split-step circuits; (v) observables (`observables.py`, DMRG `obs` dictionary); (vi) CSV/JSON files in `study/results/` or `validation/results/` with parameters, solver settings, versions and residuals; (vii) figure scripts that read only saved files. Mark the pilot ring-ED branch (`research/pilots/2026-09-27/`) as exploratory. Data source: `rydberg_chain/*.py`, `study/rydberg_model.py`. Computation ID: none (schematic).]

*Figure M1. How a calculation proceeds from atom positions to a figure. The coupling matrix $V_{ij}=g/\lvert x_i-x_j\rvert^6$ is built once and then converted into a sparse matrix for exact diagonalization, a matrix product operator for DMRG, or a Pauli-string operator for Qiskit. Each solver returns states together with their numerical diagnostics (eigenpair residuals, DMRG truncation error and energy variance, circuit step-size convergence), which are saved with the observables. Figures are drawn only from saved files. Dashed boxes mark the exploratory pilot scripts used for the preliminary ring calculations.*

### 2.1 Model, units and geometry

We set $\hbar=1$ and measure energies in units of the Rabi frequency $\Omega$, lengths in units of the mean interatomic spacing $a$, and time as $\Omega t$. Equation (1) then becomes

$$\frac{H}{\Omega}=\frac12\sum_{i=0}^{N-1}X_i-\delta\sum_{i=0}^{N-1}n_i+\sum_{i<j}V_{ij}\,n_in_j,\qquad V_{ij}=\frac{g}{\lvert x_i-x_j\rvert^{6}},\qquad (2.1)$$

$$\delta=\frac{\Delta}{\Omega},\qquad g=\frac{C_6}{\Omega a^6}=\left(\frac{R_b}{a}\right)^{6},\qquad (2.2)$$

where $X_i=\sigma^x_i$, $n_i=|1\rangle\langle1|_i$ projects atom $i$ onto the Rydberg state, the positions $x_i$ are in units of $a$, and $R_b=(C_6/\Omega)^{1/6}$ is the blockade radius, the distance at which the van der Waals energy equals $\Omega$ [1, 2]. Each pair is counted once. The dimerized-chain calculations use mainly $g=11.39$, 387.4 and 1291.5, that is, $R_b/a=1.50$, 2.70 and 3.30, which select the $\mathbb Z_2$, $\mathbb Z_3$ and $\mathbb Z_4$ regimes (Sec. 3.1); a few use $g=4.83$ ($R_b/a=1.30$), and the historical uniform-chain scans use $g=2$, 8 and 16<!-- source: 1.5^6 = 11.39, 2.7^6 = 387.42, 3.3^6 = 1291.47, 1.3^6 = 4.83; values as used in research/pilots/2026-09-27/outputs/crit_*.json, z3_24_*.txt and out24_*.txt; g = 2, 8, 16 from study/results/extension_summary.json grid.g. Checker, 27 Sept: the earlier sentence listed only the first three values, but 04_results_alternating.md also uses g = 4.83 and 05_results_size.md uses g = 2, 8, 16. -->.

**Dimerized geometry.** The chain has alternating gaps,

$$x_0=0,\qquad x_{i+1}-x_i=\begin{cases}1-\varepsilon, & i\ \text{even (strong bond)},\\ 1+\varepsilon, & i\ \text{odd (weak bond)},\end{cases}\qquad 0\le\varepsilon<1.\qquad (2.3)$$

Atoms $2k$ and $2k+1$ form strong pair $k$, located at $x=2k$ and $2k+s$ with $s=1-\varepsilon$. Because a strong and a weak gap always sum to $2$, the mean spacing, and with it the meaning of $g$, does not depend on $\varepsilon$; $\varepsilon=0$ is the uniform chain. The nearest-neighbour energies are $V_s=g/(1-\varepsilon)^6$ inside a pair and $V_w=g/(1+\varepsilon)^6$ between pairs. An open chain with even $N$ consists of $N/2$ complete pairs and starts and ends with a strong bond.

All pairs of atoms are retained unless a cutoff is stated. Truncated models keep couplings with index separation $\lvert i-j\rvert\le r_{\max}$ (atom sites) or pair separation $\le m_{\max}$ (pair sites, Sec. 2.3). Periodic rings (Sec. 2.2) use the minimum-image distance in the pair-basis scripts and a sum over the images $x_j+mN$, $\lvert m\rvert\le3$, in the constrained-basis script.

**A bound on interaction cutoffs.** Because every $V_{ij}$ is positive, discarding a set $\mathcal D$ of pairs removes the diagonal, positive-semidefinite operator $D=\sum_{(i,j)\in\mathcal D}V_{ij}n_in_j$, whose largest eigenvalue is $\lVert D\rVert=\sum_{(i,j)\in\mathcal D}V_{ij}$ (all atoms excited). Weyl's inequality [3] then bounds every ordered eigenvalue,

$$E_k^{\text{trunc}}\le E_k\le E_k^{\text{trunc}}+\lVert D\rVert,\qquad (2.4)$$

so no energy or gap changes by more than $\lVert D\rVert$. The bound says nothing about densities or correlations, which can change more when two levels are close. For the pair-site cutoff $m_{\max}=8$ in a 20-atom chain with $\varepsilon=0.3$, $g=11.39$, the bound is $\lVert D\rVert=1.4\times10^{-6}$<!-- source: scratch check, 27 Sept: sum of V_ij over atom pairs in pairs more than 8 apart, from rydberg_chain.geometry.coupling_matrix(chain_positions(20, 0.3), 11.39) = 1.361e-6 -->, and the ground-state energy actually shifts by $2.4\times10^{-7}$<!-- source: validation/results/scaling_benchmark.json, dmrg_pair entry N = 20, mmax = 8: E0_minus_ED = -2.44e-7 -->.

**Basis and bit convention.** We use the occupation basis $|b\rangle$ labelled by the integer $b=\sum_ib_i2^i$, where $b_i\in\{0,1\}$ is the occupation of atom $i$, so atom 0 is the least significant bit. This is Qiskit's little-endian qubit order: in a printed Qiskit label atom 0 is the rightmost character. With $|0\rangle$ the atomic ground state and $Z|0\rangle=+|0\rangle$,

$$n_i=\tfrac12\,(I-Z_i).\qquad (2.5)$$

### 2.2 Exact diagonalization

**Matrix elements.** The Hilbert space of $N$ atoms has dimension $2^N$. In the occupation basis the detuning and interaction terms are diagonal and the drive flips one atom:

$$\langle b|H|b\rangle=-\delta\sum_ib_i+\sum_{i<j}V_{ij}\,b_ib_j,\qquad \langle b\oplus2^i|H|b\rangle=\tfrac12,\qquad (2.6)$$

where $b\oplus2^i$ (bitwise exclusive OR) is $b$ with atom $i$ flipped. Each row thus holds one diagonal element and $N$ off-diagonal elements equal to $\Omega/2$. The package builds this matrix directly in compressed sparse row (CSR) format, with $N+1$ stored entries per row, and vectorizes the construction over blocks of $2^{18}$ basis states. The result is identical, element by element, to the historical builder in `study/rydberg_model.py`<!-- source: tests/test_geometry_ed.py::test_builder_equals_study_hamiltonian (np.array_equal on six random chains, N up to 9) -->. At $N=16$ it is built in 0.05 s instead of 0.74 s, and the whole calculation, eigensolver included, peaks at 0.11 GB of memory, whereas the historical builder alone peaks at 0.20 GB<!-- source: validation/results/scaling_benchmark.json: ed[N=16].build_s = 0.0516, peak_rss_GB = 0.110 (build plus eigensolver); study_builder[N=16].build_s = 0.738, peak_rss_GB = 0.201 (build only) -->.

The matrix is real and symmetric. Conjugating it by $\prod_iZ_i$ reverses the sign of every off-diagonal element and leaves the diagonal unchanged, so it is unitarily equivalent to a matrix with non-positive off-diagonal elements in which every pair of basis states is connected by a sequence of single flips. By the Perron–Frobenius theorem [3] the ground state of any finite array with $\Omega\neq0$ is therefore non-degenerate. Ground-state degeneracies of the infinite chain, such as the two $\mathbb{Z}_2$ crystals, appear at finite $N$ as small splittings, which we treat as gaps and follow as functions of $N$.

**Eigensolvers.** For dimensions up to 1024 ($N\le10$) we use dense LAPACK diagonalization (`scipy.linalg.eigh`). For larger $N$ we compute only the few lowest eigenpairs with Krylov-space methods. The Lanczos method [4] builds the space

$$\mathcal K_m(H,v_0)=\operatorname{span}\{v_0,Hv_0,H^2v_0,\dots,H^{m-1}v_0\},\qquad (2.7)$$

in which $H$ is tridiagonal, and approximates the extreme eigenpairs of $H$ by those of this small projected matrix. Each step costs one sparse matrix-vector product, $O(N2^N)$ operations. We use ARPACK's implicitly restarted Lanczos method [5] through `scipy.sparse.linalg.eigsh` [6], requesting the smallest algebraic eigenvalues with relative tolerance $10^{-12}$, and start from a seeded Gaussian random vector. A symmetric start vector would confine the Krylov space to one symmetry sector and could miss the lowest state of another. The package also offers the locally optimal block preconditioned conjugate gradient method (LOBPCG) [7] with the diagonal preconditioner $(\operatorname{diag}H-\sigma)^{-1}$, $\sigma=\min\operatorname{diag}H-1$, and a block of $k+2$ vectors for $k$ requested states; it falls back to ARPACK if any residual exceeds $10^{-8}\max(1,\lvert E_0\rvert)$. Every eigenpair is returned with its residual $r_k=\lVert Hv_k-E_kv_k\rVert$. The test suite requires the dense, ARPACK and LOBPCG spectra of a 10-atom dimerized chain to agree to $10^{-11}$<!-- source: tests/test_geometry_ed.py::test_eigensolvers_agree (N = 10, eps = 0.3, g = 11.39, Delta = 2.5, k = 3) -->.

The historical finite-chain scans in `study/` used the same approach: dense diagonalization up to dimension 64 and ARPACK (tolerance $10^{-11}$, seeded start vector) above it. Over 6,561 coarse-grid and 312 refinement ground states ($N=4$–12) the largest eigenpair residual is $5.5\times10^{-10}$, dense and sparse energies differ by at most $4.3\times10^{-14}$, and a bundled 62-row regression reference is reproduced to $4.7\times10^{-12}$<!-- source: study/results/extension_summary.json: grid.coarse_ground_states = 6561, grid.extra_refinement_states = 312; checks.maximum_absolute_eigenpair_residual = 5.49e-10; checks.dense_sparse_energy_error = 4.26e-14; checks.prior_csv_matching.maximum_errors (largest: alternating_order_squared = 4.70e-12, 62 points) -->.

**Time evolution.** For dynamics from a stated initial state (in this work the all-ground state $|0\cdots0\rangle$) we compute

$$|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle\qquad (2.8)$$

with `scipy.sparse.linalg.expm_multiply`, which evaluates the action of the matrix exponential on a vector by a truncated Taylor series with scaling [8] and therefore has no time-step (Trotter) error. The norm is conserved to $4.4\times10^{-15}$ or better up to $N=14$<!-- source: study/results/scaling.csv, column norm_error, maximum 4.44e-15 (N = 10) over N = 4-14 -->. Time evolution is implemented in the historical module `study/rydberg_model.py` (function `evolve`); the maintained package currently treats stationary states only.

**Symmetry-resolved ring ED (preliminary).** Bulk quantities such as critical points and gap ratios are computed on periodic rings, which have no boundary fields. In the four-state pair basis of Sec. 2.3 a ring of $N_p=N/2$ pairs is symmetric under translation by one pair, $T_2$, and under inversion about a strong-bond centre, $I_s$, which maps pair $j$ to pair $-j \bmod N_p$ and exchanges its two atoms. We work at zero momentum of $T_2$ [9]. For each orbit of pair configurations under $T_2$ we keep one representative $|a\rangle$ (the smallest integer label), with orbit length $P_a$, and form

$$|\bar a\rangle=\frac{1}{\sqrt{P_a}}\sum_{t=0}^{P_a-1}T_2^{\,t}|a\rangle,\qquad \langle\bar b|H|\bar a\rangle=\frac{\Omega}{2}\sum_{a\to b'}\sqrt{\frac{P_a}{P_b}}\quad(\bar b\neq\bar a),\qquad (2.9)$$

where the sum runs over the single-atom flips that take $|a\rangle$ to a configuration $|b'\rangle$ in the orbit of $|b\rangle$, and the diagonal elements are the classical energies of Eq. (2.6). By Burnside's lemma the number of orbits is $N_p^{-1}\sum_{t=0}^{N_p-1}d^{\gcd(t,N_p)}\approx d^{N_p}/N_p$ for local dimension $d$. For 10 pairs with $d=4$ this gives 104,968 states instead of $4^{10}=1{,}048{,}576$<!-- source: scratch check, 27 Sept: Burnside count 104968; equals nrep printed by research/pilots/2026-09-27/scripts/pairk0.py build(10, 4, ...) in the d = 3 vs d = 4 check below -->. Since $I_sT_2I_s=T_2^{-1}$, inversion maps the zero-momentum space to itself. The inversion parity $p_n=\langle\psi_n|I_s|\psi_n\rangle=\pm1$ of each of the eight lowest eigenvectors (ARPACK, tolerance $10^{-10}$) is measured after $I_s$ has been diagonalized within near-degenerate multiplets. These calculations use the pilot scripts `pairk0.py` and `crit.py` in `research/pilots/2026-09-27/scripts/`, which will be reimplemented and tested in the package (P1, P3). Two shortcuts in these scripts turned out to be harmless at the precision quoted in this paper: they accumulate interaction energies in single precision, with errors up to $9\times10^{-5}\,\Omega$ in individual high-energy diagonal elements, and they start Lanczos from an inversion-symmetric vector. Repeating three test points in double precision from random start vectors changed the lowest inversion-odd gap by less than $10^{-8}\,\Omega$ and the gap ratio of Sec. 2.5 by less than $3\times10^{-6}$<!-- source: scratch check, 27 Sept (methods_pairk0_check.py): Np = 8, d = 4, eps = 0.3, Delta = 1.3257; Np = 10, d = 3, eps = 0.3, Delta = 1.3257; Np = 8, d = 4, eps = 0, Delta = 1.0785; g = 11.39. Max |Eint(float32) - Eint(float64)| = 9.2e-5; |d gap_sigma| <= 7.5e-9; |d R| <= 2.6e-6 -->. The production code will nevertheless use double precision and random start vectors.

**Constrained bases (preliminary).** Two restrictions remove configurations of high energy and extend ED to larger rings. (i) The three-state pair basis ($d=3$) excludes doubly excited strong pairs $|11\rangle$, which cost at least $V_s-\Delta$. (ii) At $R_b/a=2.7$ and 3.3 the constrained-ring script keeps only configurations with no two excitations within $\ell_c=1$ and $\ell_c=2$ sites, respectively, a hard-boson constraint [10]. The dimension then grows as $\varphi^N$, with $\varphi=1.618$ (the golden ratio) for $\ell_c=1$ and $\varphi=1.466$ for $\ell_c=2$, instead of $2^N$: 103,682 states for $N=24$ with $\ell_c=1$, and 946,781 states for $N=36$ with $\ell_c=2$, compared with $2^{36}\approx6.9\times10^{10}$ in the full basis<!-- source: research/pilots/2026-09-27/outputs/z3_24_0.0.txt (N = 24, dim = 103682) and out36_5.0.txt (N = 36, dim = 946781). Growth rates: roots of x^2 = x + 1 and x^3 = x^2 + 1 -->. The eigenvalues of $T_2$ and $I_s$ are measured as expectation values in the lowest states. These restrictions are approximations, not symmetries: virtual transitions into excluded configurations shift energies at second order in $\Omega$. For $\ell_c=2$ at $R_b/a=3.3$ and $\Delta/\Omega=5$, the cheapest excluded configurations cost about $V(2a)-\Delta\approx15\,\Omega$<!-- source: g/2^6 = 1291.5/64 = 20.18 (scratch check, 27 Sept) minus Delta = 5 -->, so each blocked neighbour contributes shifts of order $(\Omega/2)^2/15\Omega\approx0.02\,\Omega$, comparable to the $10^{-2}\,\Omega$ gaps that distinguish the candidate phases (Sec. 3.2). The $d=3$ restriction also matters. At $\varepsilon=0.3$ and $g=11.39$, restoring $|11\rangle$ moves the finite-size estimate of the $\mathbb{Z}_2$ critical point from 16- and 20-atom rings from $\Delta_c/\Omega=1.3257$ to 1.3312, a shift of 0.4%, while the gap ratio stays at 8.00<!-- source: scratch check, 27 Sept (methods_d3_vs_d4.py, pilot pairk0.build/gaps, scaled-gap crossing of Np = 8 and 10, brentq xtol 1e-6): d = 3: Dc = 1.325677, R(Np=10) = 7.9985; d = 4: Dc = 1.331185, R(Np=10) = 7.9972; difference +0.00551. The d = 3 value agrees with research/pilots/2026-09-27/outputs/crit_eps0.3_g11.39_d3.json (Dc = 1.325676 for sizes [8, 10]) -->. This systematic shift is about seventy times larger than the $8\times10^{-5}$ difference between the estimates from two size pairs (Sec. 2.5), and it should decrease roughly as $1/V_s$ at stronger dimerization (an estimate, not yet checked). Every constrained-basis result is therefore preliminary until it is reproduced in the full basis (computation P1: $N=16$ and 20, periodic and open, $R_b/a=3.3$, $\varepsilon\in\{0,0.3,0.5\}$, $\Delta/\Omega\in\{2,3,5\}$; the $\mathbb{Z}_2$ line is to be recomputed with $d=4$ in P3).

[TABLE PENDING: P1 validation of the constrained bases. For each (N, boundary, ε, Δ/Ω) of the P1 grid: lowest four energies, their T2 and I_s quantum numbers, and gaps in the constrained and the full basis, with the differences. Data source: new validation script to be written in validation/ using rydberg_chain.ed; the constrained results from research/pilots/2026-09-27/scripts/ced.py. Computation ID: P1.]

**Cost and practical limits.** Table M1 lists measured costs of full-basis ED for dimerized open chains with all pairs retained.

*Table M1. Cost of full-basis ED for the two lowest eigenpairs of an open dimerized chain ($\varepsilon=0.3$, $g=11.39$, $\Delta/\Omega=2.5$, all pairs). Times include building the matrix; each case ran in its own process with two threads on a shared 4-core Intel Xeon (2.1 GHz, 16 GB), so times are upper bounds that varied by about a factor of two between runs. Peak memory is the resident-set maximum of the process.*

| $N$ | Dimension | Stored matrix entries | ARPACK time (s) | LOBPCG time (s) | Peak memory, ARPACK / LOBPCG (GB) | Largest ARPACK residual |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 65,536 | 1,114,112 | 6.6 | 1.8 | 0.11 / 0.12 | $4.5\times10^{-12}$ |
| 18 | 262,144 | 4,980,736 | 39.7 | 8.7 | 0.21 / 0.26 | $5.2\times10^{-12}$ |
| 20 | 1,048,576 | 22,020,096 | 362 | 55 | 0.56 / 0.85 | $3.5\times10^{-12}$ |

<!-- source: validation/results/scaling_benchmark.json: ed[*].dim, nnz, total_s (6.64, 39.72, 361.99), peak_rss_GB (0.110, 0.210, 0.561), residual_max (4.47e-12, 5.19e-12, 3.46e-12); ed_lobpcg[*].total_s (1.83, 8.70, 54.96), peak_rss_GB (0.119, 0.263, 0.853), residual_max <= 9.4e-13, E0_minus_eigsh and E1_minus_eigsh <= 1.2e-12 in magnitude. Machine block: 4 cores, Intel Xeon 2.10 GHz, MemTotal 16481980 kB, OMP/OPENBLAS/MKL threads = 2, load average 4-9 (notes field). -->

ARPACK and LOBPCG agree to $1.2\times10^{-12}$ in the two lowest energies, and LOBPCG was 3.6 to 6.6 times faster; the ARPACK runs had the heavier machine load, so part of this ratio reflects load<!-- source: validation/results/scaling_benchmark.json, ed_lobpcg[*].E0_minus_eigsh, E1_minus_eigsh (max magnitude 1.12e-12); time ratios 6.64/1.83 = 3.6, 39.72/8.70 = 4.6, 361.99/54.96 = 6.6; machine.load_average_at_start = 4.2 (ARPACK pass) vs load_average_at_lobpcg_pass = 1.8; summary field: "the latter under heavier load" -->. The solve time grew by a factor of five to nine for every two added atoms. The stored matrix alone takes $12(N+1)2^N$ bytes (8-byte values and 4-byte column indices): 0.26 GB at $N=20$, 1.2 GB at $N=22$, 5.0 GB at $N=24$ and 22 GB at $N=26$. Each Krylov vector adds $8\cdot2^N$ bytes (0.13 GB at $N=24$), and ARPACK keeps about twenty of them<!-- source: formula checked against scaling_benchmark.json nnz = (N+1) 2^N and rydberg_chain/ed.py (float64 data, int32 indices while (N+1) 2^N < 2^31, true up to N = 26); byte counts recomputed by the checker, 27 Sept; ARPACK's default ncv = max(2k+1, 20) (SciPy documentation) -->. On the benchmark machine, full-basis ED is therefore routine up to $N=20$ and practical at $N=22$. At $N=24$ it is possible but untested: scaling the measured peak memory by four and the measured times by six to nine per two atoms gives about 9 GB and several hours with ARPACK, or about 14 GB and under an hour with LOBPCG, close to the 16 GB of the machine<!-- source: extrapolation by the checker, 27 Sept, from scaling_benchmark.json: ARPACK 0.56 GB x 16 = 9.0 GB, 362 s x 9^2 = 8 h; LOBPCG 0.85 GB x 16 = 13.6 GB, 55 s x 6.3^2 = 36 min. An estimate, not a measurement. -->. $N=26$ would need a matrix-free implementation that stores only the diagonal (0.54 GB) and applies the spin flips on the fly, with a solver that keeps few vectors of 0.54 GB each. Symmetry resolution reduces the 20-atom ring to 104,968 states, and the restricted bases above reach rings of 24 atoms ($d=3$) and 36 atoms ($\ell_c=2$; 40 atoms are planned in P2). DMRG removes the exponential cost altogether for weakly entangled states (Sec. 2.4 and Table M2). Figure M2 compares the two approaches.

[FIGURE PENDING: cost of ground-state calculations versus N. Left: wall time; right: peak memory (both logarithmic). Series: full-basis ED with ARPACK and with LOBPCG (N = 16, 18, 20) and extrapolated memory 12(N+1)2^N bytes to N = 26 (dashed); pair-site DMRG at χ = 64, m_max = 8 (N = 20, 50, 100, 200), with open symbols for runs that stopped at the sweep limit. Optional inset: the historical time-evolution benchmark, study/results/scaling.csv (also plotted in study/figures/08_scaling.png). Data source: validation/results/scaling_benchmark.json (all cases present as of 27 Sept; add converged large-N runs from P7). Computation IDs: P6, P7.]

*Figure M2. Cost of computing the ground state of an open dimerized chain ($\varepsilon=0.3$, $g=11.39$, $\Delta/\Omega=2.5$) by full-basis exact diagonalization and by pair-site DMRG. The ED memory grows as $2^N$ and the solve time by a factor of five to nine per two atoms, which limits full-basis ED to about 22–24 atoms on a 16 GB machine. DMRG at fixed bond dimension $\chi=64$ needs at most 0.26 GB up to 200 atoms, and its time per sweep grows roughly linearly with $N$ once $\chi$ saturates; its limitation is convergence in the ordered phase of even open chains (Sec. 2.4), not memory. Open symbols mark runs that reached the sweep limit before meeting the entropy criterion.*

### 2.3 Exact pair-site regrouping

Grouping the two atoms of each strong pair, $L=2k$ and $R=2k+1$, into one site with the four states $|n_Ln_R\rangle\in\{|00\rangle,|10\rangle,|01\rangle,|11\rangle\}$ rewrites Eq. (2.1) as

$$\frac{H}{\Omega}=\sum_{k}h_k+\sum_{k<l}\ \sum_{\alpha,\beta\in\{L,R\}}\frac{g}{\big[2(l-k)+o_\beta-o_\alpha\big]^6}\,n_{k\alpha}n_{l\beta},\qquad h_k=\tfrac12\big(X_{kL}+X_{kR}\big)-\delta\big(n_{kL}+n_{kR}\big)+V_s\,n_{kL}n_{kR},\qquad (2.10)$$

with offsets $o_L=0$ and $o_R=s=1-\varepsilon$. This is a relabelling of the same Hilbert space, not an approximation. It has four practical advantages. (i) The largest coupling, $V_s$ ($96.8\,\Omega$ at $\varepsilon=0.3$, $g=11.39$<!-- source: rydberg_chain.geometry.strong_weak_couplings(0.3, 11.39) = (96.81, 2.36), scratch check, 27 Sept -->), becomes an on-site energy. (ii) The number of sites, and the interaction range measured in sites, are halved. (iii) The pseudo-spin $\tau_k=n_{kL}-n_{kR}$ and the pair-state projectors are on-site operators. (iv) $T_2$ becomes a one-site translation, and $I_s$ a reflection of the site order combined with the on-site exchange $|10\rangle\leftrightarrow|01\rangle$. Pair-site DMRG reproduces atom-site ED to within $10^{-11}\,\Omega$ (Sec. 2.4).

### 2.4 Density-matrix renormalization group

**Matrix product states.** A matrix product state (MPS) on $L$ sites with local dimension $d$ ($d=2$ for atoms, 4 for pairs) is

$$|\psi\rangle=\sum_{s_1,\dots,s_L}A_1^{s_1}A_2^{s_2}\cdots A_L^{s_L}\,|s_1s_2\cdots s_L\rangle,\qquad (2.11)$$

where each $A_\ell^{s}$ is a matrix of dimension at most $\chi\times\chi$, the bond dimension [11–13]. In canonical form, the MPS across the bond after site $\ell$ is a Schmidt decomposition $|\psi\rangle=\sum_{\alpha=1}^{\chi}\lambda_\alpha|\alpha\rangle_{\text{left}}|\alpha\rangle_{\text{right}}$, so the entanglement entropy of the left block is $S_\ell=-\sum_\alpha\lambda_\alpha^2\log_2\lambda_\alpha^2\le\log_2\chi$. Ground states of gapped one-dimensional Hamiltonians obey an area law [14] and are represented accurately at fixed $\chi$ for any length. At a critical point $S_\ell$ grows logarithmically with $\ell$ (Sec. 2.5), so the required $\chi$ grows as a power of the system size. This is why DMRG reaches hundreds of atoms at modest cost in gapped phases, and why critical points require convergence checks in $\chi$.

**Variational optimization.** DMRG minimizes $E[\psi]=\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ over MPS of bond dimension $\chi$. The two-site algorithm sweeps back and forth along the chain. At each bond it solves the eigenvalue problem $H_{\text{eff}}\Theta=E\Theta$ for the two-site tensor $\Theta$ in the basis of the current left and right blocks (by Lanczos), splits $\Theta$ by a singular value decomposition, and keeps the $\chi$ largest singular values. The discarded weight

$$\epsilon_{\text{trunc}}=\sum_{\alpha>\chi}\lambda_\alpha^2\qquad (2.12)$$

measures the truncation error. Because every step is variational, the DMRG energy is an upper bound, $E_{\text{DMRG}}\ge E_0$. We use the TeNPy library [15, 16] with two-site updates, singular values below $10^{-12}$ discarded, and a density-matrix perturbation ("mixer") [17] of amplitude $10^{-5}$ that decays by a factor of two per sweep and is switched off after ten sweeps; the mixer lets the optimization escape poor local minima. A run is declared converged when the energy and the mean bond entropy change between sweeps by less than $10^{-12}$ and $10^{-8}$, respectively, with between 6 and 60 sweeps<!-- source: rydberg_chain/dmrg.py, dmrg_options defaults (chi = 64, max_sweeps = 60, min_sweeps = 6, max_E_err = 1e-12, max_S_err = 1e-8, svd_min = 1e-12, mixer amplitude 1e-5, decay 2, disable_after 10). TeNPy 1.1.1 TwoSiteDMRGEngine.DefaultMixer = DensityMatrixMixer (checked 27 Sept) -->. As an independent check we evaluate the energy variance

$$\sigma_E^2=\langle H^2\rangle-\langle H\rangle^2,\qquad (2.13)$$

which vanishes for an exact eigenstate. By Temple's inequality [18], $E-E_0\le\sigma_E^2/(E_1-E)$ when $E<E_1$, so a small variance certifies the energy whenever the gap to the next level is not also small.

**Initial state.** DMRG starts from an alternating product state with both end atoms excited. For even $N$ the two halves are joined by a single 00 domain wall at the centre, which the boundary fields of the ordered phase favour (Sec. 3.3).

**Long-range interactions as a matrix product operator.** The Hamiltonian is encoded as a matrix product operator (MPO) that contains every two-body term up to the cutoff explicitly; we do not fit $r^{-6}$ by a sum of exponentials [19, 20]. With explicit terms, each MPO bond carries one channel for every operator that has been placed to its left and still awaits its partner. The MPO bond dimension is therefore $D_W=2+(\text{number of distinct coupling operators})\times(\text{range})$: $2+r_{\max}$ on atom sites (one operator, $n$) and $2+2m_{\max}$ on pair sites ($n_L$ and $n_R$). The benchmark MPOs have $D_W=18$ for $m_{\max}=8$ and $D_W=20$ for all pairs in a 10-pair chain<!-- source: validation/results/scaling_benchmark.json, dmrg_pair[*].mpo_bond_dim = 18 (mmax = 8) and 20 (mmax = null, 10 pairs) -->. Because $r^{-6}$ decays quickly, the couplings discarded by $m_{\max}=8$ (atom separations above $17a$) are small: the bound of Eq. (2.4) is $1.4\times10^{-6}\,\Omega$ for 20 atoms and grows in proportion to $N$, to $2.9\times10^{-4}\,\Omega$ ($1.4\times10^{-6}\,\Omega$ per atom) for 200 atoms<!-- source: checker, 27 Sept: sum of V_ij over atoms in pairs more than 8 apart, rydberg_chain.geometry.coupling_matrix(chain_positions(N, 0.3), 11.39): 1.36e-6 (N = 20), 4.58e-5 (50), 1.26e-4 (100), 2.85e-4 (200). The earlier text gave "about 1e-6" without the N dependence. -->. At fixed $\chi$ and $D_W$ the cost of a sweep grows linearly with the number of sites.

**Finite and infinite DMRG.** All DMRG results in this paper use finite DMRG on open chains, with one site per atom or one site per pair. Infinite DMRG (iDMRG) [21] would optimize a repeating unit cell directly in the thermodynamic limit and remove boundary effects. An exploratory iDMRG scan with a two-atom unit cell ($\chi=48$, interaction range 12 atoms, $g=11.39$, $\Delta/\Omega=2.5$) failed from $\varepsilon=0.3$ upward: its energy per atom varies smoothly from $-1.2219$ to $-1.2193\,\Omega$ for $\varepsilon\le0.25$, then jumps to $-0.725\,\Omega$ at 0.3 and changes non-monotonically beyond. Several values are positive, which violates the variational bound $E_0\le0$ set by the empty state, whose energy is exactly zero<!-- source: research/pilots/2026-09-27/outputs/iscan_g11.39_D2.5.json: chi = 48, cell = 2; E_per_site = -1.2219 at eps = 0 ... -1.2193 at eps = 0.25 (smooth), -0.725 at 0.3, -0.631 at 0.35, +9.48 at 0.4, +26.2 at 0.5, +2.65 at 0.55, +1.36 at 0.6, ... +2396 at 0.85; interaction range rmax = 12 from scripts/dimer_model.py build_infinite default. docs/DIMERIZED_THEORY.md section 6 records the failure -->. No iDMRG number is used in this paper. A repaired implementation on pair sites with a two-pair unit cell (`rydberg_chain/infinite.py`) is in development and does not yet pass its own validation against finite DMRG; before it is used it must also reproduce the preliminary ring-ED energy of $-1.2284\,\Omega$ per atom at $\varepsilon=0.55$, $\Delta/\Omega=2.5$, $g=11.39$ (P6)<!-- source: rydberg_chain/infinite.py (created 27 Sept, after the first draft of this section); tests/test_infinite.py: test_idmrg_matches_finite_bulk_energy fails (-1.23261 vs finite bulk -1.21962) and test_idmrg_order_is_translation_and_mirror_consistent fails (order spread 1.8e-8 > 1e-8) in a checker run on 27 Sept; validation/idmrg_validation.py has no results file yet. -1.2284: docs/RESEARCH_PLAN.md P6 row and docs/DIMERIZED_THEORY.md Sec. 6. -->.

**Validation against ED.** With the package defaults ($\chi=64$), atom-site DMRG reproduces the ED ground-state energies of open chains with $N=10$–12 ($g=11.39$, $\Delta/\Omega=2.5$, $\varepsilon=0$ and 0.3) to $5\times10^{-14}\,\Omega$ and the site occupations to $2\times10^{-13}$. Pair-site DMRG reproduces them for $N=12$ and 14 ($\varepsilon=0$, 0.3 and 0.6; for $\varepsilon=0.3$ with the tighter entropy criterion $10^{-10}$) to $8\times10^{-12}\,\Omega$ in energy and $1.2\times10^{-11}$ in occupations<!-- source: scratch check, 27 Sept (methods_dmrg_check.py; same cases and settings as tests/test_dmrg.py): atom |dE| = 5.3e-15, 3.4e-14, 5.3e-14, max |d n_i| = 1.7e-13; pair |dE| = 2.1e-14, 6.4e-14, 7.5e-12, max |d n_i| = 1.2e-11. The tests assert |dE| < 1e-9 (atom) and < 1e-8 (pair). docs/DIMERIZED_THEORY.md quotes "2e-9 or better" for pair sites, from an earlier setting -->, and at $N=20$ ($\varepsilon=0.3$, all pairs) to $4\times10^{-13}\,\Omega$ with an energy variance of order $10^{-12}$<!-- source: validation/results/scaling_benchmark.json, dmrg_pair[N = 20, mmax = null]: E0_minus_ED = -4.16e-13, variance = -1.0e-12 (rounding) -->. Atom-site DMRG also reproduces the saved uniform-chain energy $E_0/\Omega=-18.399305357538882$ ($N=12$, $g=8$, $\Delta/\Omega=3$) to $5\times10^{-14}$<!-- source: E0 constant in tests/test_dmrg.py and tests/test_geometry_ed.py; saved row in study/results/size_comparison_V8.csv (N = 12, V = 8, delta = 3: energy_over_omega = -18.3993053575389). DMRG - saved = -5.3e-14 (scratch check, 27 Sept) -->. The test suite checks all of these cases except the 20-atom benchmark, with energy tolerances of $10^{-9}$ (atom sites) and $10^{-8}$ (pair sites).

**Reach at large N.** Table M2 shows how far pair-site DMRG reaches.

*Table M2. Pair-site DMRG for open dimerized chains ($\varepsilon=0.3$, $g=11.39$, $\Delta/\Omega=2.5$, $\chi=64$, $m_{\max}=8$, at most 30 sweeps; same machine as Table M1). $\Delta E_{\text{sweep}}$ is the energy change in the last sweep and $\sigma_E^2$ the energy variance, Eq. (2.13).*

| $N$ | $E_0/(N\Omega)$ | Sweeps | Time per sweep (s) | Total time (s) | Peak memory (GB) | $\sigma_E^2/\Omega^2$ | $\lvert\Delta E_{\text{sweep}}\rvert/\Omega$ | Converged |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 20 | $-1.24436$ | 16 | 1.4 | 22 | 0.16 | $3\times10^{-13}$ | $7\times10^{-15}$ | yes |
| 50 | $-1.22960$ | limit | 11.3 | 351 | 0.23 | $3.7\times10^{-9}$ | $1.6\times10^{-9}$ | no (entropy) |
| 100 | $-1.22453$ | limit | 24.2 | 751 | 0.24 | $1.9\times10^{-4}$ | $2.1\times10^{-5}$ | no |
| 200 | $-1.22197$ | limit | 32.6 | 1012 | 0.26 | $4.3\times10^{-4}$ | $9.3\times10^{-5}$ | no |

<!-- source: validation/results/scaling_benchmark.json, dmrg_pair entries with mmax = 8: E_per_atom (-1.2443602, -1.2296048, -1.2245256, -1.2219652), n_sweeps (16, 31, 31, 31; max_sweeps = 30, so "limit"), seconds_per_sweep (1.39, 11.31, 24.23, 32.64), seconds (22.3, 350.6, 751.2, 1012.0), peak_rss_GB (0.157, 0.234, 0.243, 0.259), variance (3.4e-13, 3.69e-9, 1.88e-4, 4.27e-4), last_sweep_dE (7.1e-15, 1.58e-9, 2.14e-5, 9.35e-5), last_sweep_dS (1.4e-9, 7.2e-6, 4.3e-3, 2.4e-3), converged (true, false, false, false). NOTE: the file's "notes" field says "the energy variance shows the energy is converged regardless"; that holds for N = 50 but not for N = 100 and 200. Correct the note. -->

Memory stays at or below 0.26 GB up to 200 atoms, and the time per sweep grows roughly linearly with $N$ once $\chi$ saturates at 64 (11.3, 24.2 and 32.6 s at $N=50$, 100 and 200 on the shared machine, whose timings vary by about a factor of two). The limitation is convergence. In even open chains in the ordered phase, the boundary fields force a single domain wall whose position is a soft mode: the ED gap of this chain falls from $0.051\,\Omega$ at $N=16$ to $0.036\,\Omega$ at $N=20$<!-- source: validation/results/scaling_benchmark.json, ed[*].gap = 0.0512 (N = 16), 0.0427 (N = 18), 0.0362 (N = 20) -->, and local two-site updates move the wall only slowly. At $N=50$ the energy has converged (variance $3.7\times10^{-9}$) but the entropy criterion was not met within 30 sweeps; at $N=100$ and 200 neither has converged, and the energy variance is of order $10^{-4}\,\Omega^2$. Before any DMRG result with $N\ge50$ is used, P7 will pin the chain ends with a local detuning (the `local_delta` argument), use odd atom numbers with atom sites, and raise $\chi$ and the sweep limit.

### 2.5 Observables and diagnostics

For a state with amplitudes $\psi_b$ in the occupation basis, or an MPS, we compute the following.

*Occupations and density,*

$$\langle n_i\rangle=\sum_b\lvert\psi_b\rvert^2\,b_i,\qquad \bar n=\frac1N\sum_i\langle n_i\rangle.\qquad (2.14)$$

*Staggered order.* With

$$Q=\frac1N\sum_i(-1)^i\,(2n_i-1),\qquad (2.15)$$

both perfectly alternating configurations have $Q^2=1$ for either parity of $N$. This corrects the occupation imbalance $M=(2/N)\sum_i(-1)^in_i$ used in earlier versions of the study, whose two ideal values differ in magnitude for odd $N$. For independent sites with excitation probability $p$, $\langle Q^2\rangle=4p(1-p)/N+(2p-1)^2(N\bmod2)/N^2$, so a nonzero $\langle Q^2\rangle$ in a short chain does not by itself indicate order<!-- source: docs/RESULTS_AND_LIMITS.md, "A corrected order diagnostic"; implemented in rydberg_chain/observables.py::staggered_order_squared -->.

*Pair observables.* The pair order $\langle\tau_k\rangle=\langle n_{2k}-n_{2k+1}\rangle$ and the probabilities of the four pair states $|n_Ln_R\rangle$ of each strong pair.

*Structure factors and ring order parameters.* With the pair occupation $\mathcal N_k=n_{kL}+n_{kR}$,

$$S_{\mathcal O}(q)=\frac1{N_p}\Big\langle\Big\lvert\sum_{k=0}^{N_p-1}e^{-iqk}\mathcal O_k\Big\rvert^2\Big\rangle,\qquad \mathcal O\in\{\mathcal N,\tau\}.\qquad (2.16)$$

Uniform pseudo-spin order ($\mathbb Z_2$ crystal of the dimerized chain, which breaks $I_s$) makes $S_\tau(0)/N_p$ finite in the limit of large $N_p$; a period-two pair density wave (the candidate phase of Sec. 3.2, which would break $T_2$) makes $S_{\mathcal N}(\pi)/N_p$ finite. The constrained-ring script uses the equivalent order parameters $u=N_p^{-1}\sum_k(-1)^k\mathcal N_k$, $\bar\tau=N_p^{-1}\sum_k\tau_k$ and $v=N_p^{-1}\sum_k(-1)^k\tau_k$ and reports $\langle u^2\rangle=S_{\mathcal N}(\pi)/N_p$, $\langle\bar\tau^2\rangle=S_\tau(0)/N_p$ and $\langle v^2\rangle=S_\tau(\pi)/N_p$. The four $\mathbb Z_4$ patterns have $(u,v)=(\pm\tfrac12,\pm\tfrac12)$, whereas a pair density wave has $u\neq0$ with $v$ and $\bar\tau$ disordered<!-- source: research/pilots/2026-09-27/scripts/ced.py (sgn = (-1)^(i//2), lr = +1 for L and -1 for R; u = B sgn/Np, tau = B lr/Np, v = B (sgn lr)/Np) and sfac.py (S_N, S_tau from FFT over pairs). Pattern values worked out by hand -->.

*Connected correlations,* $C_{ij}=\langle n_in_j\rangle-\langle n_i\rangle\langle n_j\rangle$.

*Entanglement entropy.* For a bipartition of the chain into a left block $A$ and the rest $B$,

$$S_A=-\operatorname{Tr}\rho_A\log_2\rho_A,\qquad \rho_A=\operatorname{Tr}_B|\psi\rangle\langle\psi|,\qquad (2.17)$$

measured in bits. For ED we reshape $\psi$ into a $2^{N-\ell}\times2^{\ell}$ matrix, whose singular values are the Schmidt coefficients $\lambda_\alpha$, and use $S_A=-\sum_\alpha\lambda_\alpha^2\log_2\lambda_\alpha^2$; DMRG provides the $\lambda_\alpha$ directly. A two-atom Bell pair gives one bit and a product state zero, which the tests check.

*Central charge from entanglement.* At a critical point described by a conformal field theory (CFT) with central charge $c$, the entropy of the first $\ell$ sites of an open chain of length $L$ is [22]

$$S(\ell)=\frac c6\log_2\!\Big[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}\Big]+s_1+\ldots,\qquad (2.18)$$

and $S(\ell)=(c/3)\log_2[(L/\pi)\sin(\pi\ell/L)]+s_1'$ for a ring. The omitted terms include oscillations that alternate from bond to bond near open boundaries [23]; in the dimerized chain the planned fits (P7) use cuts of one type only, between pairs, which pair-site DMRG provides directly. The Ising and tricritical Ising CFTs have $c=1/2$ and $7/10$ [24, 25].

*Finite-size spectra, critical points and the gap ratio.* On a critical ring of length $L$, the low-lying levels obey [26–28]

$$E_0(L)=e_\infty L-\frac{\pi vc}{6L}+\ldots,\qquad E_n(L)-E_0(L)=\frac{2\pi v}{L}\,x_n+\ldots,\qquad (2.19)$$

where $v$ is the velocity of excitations and $x_n$ the scaling dimension of the corresponding operator. With dynamic exponent $z=1$, the scaled gap $L\,[E_\sigma(L)-E_0(L)]$ is independent of $L$ at the critical point, so we locate $\Delta_c$ as the crossing of this quantity for two ring sizes, found by Brent's method (requested tolerance $10^{-4}$ in $\Delta/\Omega$; a repeat at $10^{-6}$ changed the $\varepsilon=0.3$ value below by $10^{-6}$). Here $E_\sigma$ is the lowest zero-momentum level that is odd under $I_s$. The velocity cancels in the gap ratio

$$R=\frac{E_\epsilon-E_0}{E_\sigma-E_0}\ \xrightarrow{\,L\to\infty\,}\ \frac{x_\epsilon}{x_\sigma},\qquad (2.20)$$

where $E_\epsilon$ is the first excited zero-momentum level that is even under $I_s$. For the Ising CFT $x_\sigma=1/8$ and $x_\epsilon=1$, so $R=8$; for the tricritical Ising CFT $x_\sigma=3/40$ and $x_\epsilon=1/5$, so $R=8/3$ [24, 25]. Because the two values differ by a factor of three and require no fitted velocity, $R$ discriminates the two universality classes at modest ring sizes. Finite-size spectra of this kind have been used to identify the Ising CFT in numerical studies of the uniform Rydberg chain [29], and universal ratios of CFT energy levels, at both Ising and tricritical Ising transitions, have been measured in a Rydberg-chain experiment [30]. As an illustration (preliminary, pilot $d=3$ basis, to be recomputed in P3), at $\varepsilon=0.3$ and $g=11.39$ the crossings of 16/20-atom and 20/24-atom rings give $\Delta_c/\Omega=1.32568$ and 1.32575 with $R=8.00$<!-- source: research/pilots/2026-09-27/outputs/crit_eps0.3_g11.39_d3.json: pairs [8, 10] Dc = 1.325676, R_b = 7.998; [10, 12] Dc = 1.325755, R_b = 7.998. crit.py uses brentq xtol = 1e-4; the scratch repeat methods_d3_vs_d4.json (xtol = 1e-6, d = 3, sizes [8, 10]) gives 1.325677, i.e. 1.0e-6 from the pilot value. See the d = 3 systematic shift of +0.0055 in Sec. 2.2 -->; the systematic shift from the $d=3$ restriction (Sec. 2.2) is larger than the difference between these two estimates.

### 2.6 Effective pseudo-spin theory

The analytic picture used in Sec. 3.2 is derived in the theory notes and summarized here. When each strong pair holds at most one excitation, an occupied pair is a two-level system with pseudo-spin $\tau^z=+1$ (left atom excited) or $-1$ (right atom excited). The interaction of excitations in pairs $k$ and $k+m$ decomposes exactly as

$$U_m=A_m+B_m\big(\tau^z_k-\tau^z_{k+m}\big)-J_m\,\tau^z_k\tau^z_{k+m},\qquad J_m=\tfrac14\big[V(2m-s)+V(2m+s)-2V(2m)\big],\qquad B_m=\tfrac14\big[V(2m+s)-V(2m-s)\big],\qquad (2.21)$$

with $V(r)=g/r^6$ and lengths in units of $a$. Because $r^{-6}$ is convex, $J_m>0$, a ferromagnetic coupling of the pseudo-spins; $J_m\to0$ as $\varepsilon\to1$ ($s\to0$), and the $B_m$ terms cancel in the bulk and leave fields at the chain ends. The drive connects $\tau^z=+1$ and $-1$ at second order through the empty ($|00\rangle$) and doubly excited ($|11\rangle$) pair states, giving the effective transverse-field Ising model

$$H_{\text{eff}}=-\sum_{k,m}J_m\,\tau^z_k\tau^z_{k+m}-h\sum_k\tau^x_k+\ldots,\qquad h\simeq h_0=\frac{\Omega^2}{4}\left[\frac1\Delta+\frac1{V_s-\Delta}\right]\quad(0<\Delta<V_s).\qquad (2.22)$$

The nearest-neighbour model is critical at $J_1=h$ [31]. The leading-order field $h_0$ neglects the interaction of the flipped excitation with its neighbours in the energy denominators; the dressed form and its accuracy are discussed in Sec. 3.2. The module `rydberg_chain/effective.py` evaluates $J_m$, $B_m$, $h_0$ and the dimerization at which $J_1=h_0$. Its tests check Eq. (2.21) against the $r^{-6}$ energies of explicit configurations, and $h_0$ against half the exact splitting of the one-excitation doublet of a single driven pair (to a relative accuracy of $10^{-3}$ at $\Omega=0.02$)<!-- source: tests/test_effective.py::test_couplings_match_direct_pair_energies, test_transverse_field_matches_single_pair. Sign convention: effective.tfim_couplings returns J_m <= 0 as the coefficient of +tau tau; Eq. (2.21) uses J_m >= 0 with an explicit minus sign, as in docs/DIMERIZED_THEORY.md -->.

### 2.7 Computational details of modeling in Qiskit

<!-- Outline heading: "Computational Details of Modeling in Qiskit". It can be promoted to its own section heading at merge. -->

Qiskit [32] provides an independent construction of Eq. (2.1) and of its time evolution. Substituting Eq. (2.5) into Eq. (2.1) gives the Pauli-string form

$$H=c_0\,I+\frac{\Omega}{2}\sum_iX_i+\sum_ih_i^zZ_i+\sum_{i<j}\frac{V_{ij}}{4}Z_iZ_j,\qquad c_0=-\frac{N\Delta}{2}+\frac14\sum_{i<j}V_{ij},\qquad h_i^z=\frac{\Delta}{2}-\frac14\sum_{j\neq i}V_{ij}.\qquad (2.23)$$

The function `qiskit_hamiltonian` in `study/rydberg_model.py` writes these terms as a `SparsePauliOp`, placing atom $i$ at character $N-1-i$ of each label, and combines duplicate strings. Its matrix equals the occupation-basis matrix of Eq. (2.6) to $2.2\times10^{-16}$ for a three-atom test chain, and to $3.6\times10^{-15}$ for four-atom chains at three couplings and three interaction cutoffs<!-- source: study/results/summary.json, checks.Pauli_vs_bit_basis_H_error = 2.2e-16 (3 atoms, V_nn = 3, Omega = 0.7, Delta = 0.4 in the code's energy units); study/results/extension_summary.json, checks.qiskit_operator_max_entry_error = 3.6e-15 (study/run_extension.py: N = 4, (g, delta) = (2, -1.5), (8, 1), (16, 5), cutoffs nearest1, nearest2, all) -->.

**Exact evolution.** `HamiltonianGate(H, t)` appends the unitary $e^{-iHt}$, obtained by dense matrix exponentiation, to a circuit, which we simulate as a state vector (`Statevector`). For two atoms at $V/\Omega=0$ and 20, sampled at 161 times in $0\le\Omega t\le2\pi$, these states agree with the SciPy evolution of Eq. (2.8) to $1.0\times10^{-14}$<!-- source: study/results/summary.json, checks.qiskit_scipy_state_error_V0 = 1.02e-14, checks.qiskit_scipy_state_error_V20 = 5.2e-15 -->. The entanglement entropy of a random five-atom state computed by our code agrees with Qiskit's `partial_trace` and `entropy` to $2.2\times10^{-16}$ bits, and a two-qubit Bell circuit gives 1.0000 bit<!-- source: study/results/summary.json, checks.entropy_vs_Qiskit_error = 2.2e-16; checks.bell_entropy_bits = 1.0000000000000002 -->.

**Split-step circuits.** To represent the evolution with standard gates we split $H=H_X+H_D$ into the drive $H_X=(\Omega/2)\sum_iX_i$ and the diagonal part $H_D=-\Delta\sum_in_i+\sum_{i<j}V_{ij}n_in_j$, and use the symmetric second-order (Strang) product [33–35]

$$U_n=\Big[e^{-iH_X\,dt/2}\,e^{-iH_D\,dt}\,e^{-iH_X\,dt/2}\Big]^n,\qquad dt=\frac tn,\qquad \big\lVert U_n-e^{-iHt}\big\rVert=O\big(t\,dt^2\big).\qquad (2.24)$$

The prefactor of the error is set by nested commutators of $H_X$ and $H_D$ [36]. Each factor is exact in gates:

$$e^{-iH_X\,dt/2}=\prod_i R_X^{(i)}\!\Big(\frac{\Omega\,dt}{2}\Big),\qquad e^{-iH_D\,dt}=\prod_iP^{(i)}(\Delta\,dt)\prod_{i<j}CP^{(ij)}(-V_{ij}\,dt),\qquad (2.25)$$

with $R_X(\theta)=e^{-i\theta X/2}$, $P(\lambda)=\operatorname{diag}(1,e^{i\lambda})$ and $CP(\lambda)=\operatorname{diag}(1,1,1,e^{i\lambda})$. The diagonal factor is exact because all terms of $H_D$ commute. One step uses $2N$ $R_X$ gates, $N$ phase gates and one controlled-phase gate per retained pair: 18 gates for four atoms with all pairs<!-- source: study/results/circuit_convergence.csv, gates = 180 for 10 steps (N = 4, all 6 pairs); stepped_circuit in study/rydberg_model.py -->.

**Convergence tests.** For four atoms with all pair interactions (nearest-neighbour $V/\Omega=5$), $\Delta=0$ and $\Omega t=1$, the infidelity with respect to the exact state falls from $1.3\times10^{-5}$ at 10 steps to $3.1\times10^{-9}$ at 80 steps, by a factor of 16 per halving of $dt$, as expected for a state error proportional to $dt^2$<!-- source: study/results/circuit_convergence.csv: infidelity 1.301e-5 (10 steps), 8.07e-7 (20), 5.04e-8 (40), 3.15e-9 (80); successive ratios 16.1, 16.0, 16.0 -->. Figure M3 repeats the test at both signs of the detuning against a reference built independently from Kronecker products, which agrees with the sparse matrix to $8.9\times10^{-16}$<!-- source: study/results/circuit_controls_summary.json, maximum_independent_tensor_vs_sparse_matrix_error = 8.88e-16 -->. The fitted order is 2.002 for $\Delta/\Omega=-1.3$ and 2.001 for $+1.3$. A separate test fixes the sign convention: with the laser off, an atom prepared in $(|0\rangle+|1\rangle)/\sqrt2$ must acquire $\langle Y\rangle=\sin(\Delta t)$; the circuit reproduces this to $2\times10^{-16}$, while the opposite sign would give a state error of 0.88<!-- source: study/results/circuit_controls_summary.json, convergence[*].observed_state_error_order_last_three = 2.0017 (-1.3), 2.0006 (+1.3); phase_sign: expected_y = -0.7895037396899504, circuit_y = -0.7895037396899502 (Delta = -1.3, t = 0.7), deliberately_reversed_sign_state_error = 0.879 -->.

![](study/figures/circuit_detuning_convergence.png)

*Figure M3. Convergence of the second-order split-step Qiskit circuit, Eqs. (2.24)–(2.25), to exact evolution. Three atoms in an open chain with all pair interactions (nearest-neighbour $V/\Omega=5$), initial state $|000\rangle$, $\Omega t=1.7$, and $\Delta/\Omega=-1.3$ or $+1.3$, for 8 to 128 steps. Left: state error $\lVert\psi_{\text{circuit}}-e^{i\phi}\psi_{\text{exact}}\rVert$ after removing the global phase, with an $n^{-2}$ reference line (dashed); right: infidelity $1-\lvert\langle\psi_{\text{exact}}|\psi_{\text{circuit}}\rangle\rvert^2$. The reference state is the dense matrix exponential of a Hamiltonian built independently from Kronecker products; the finest-step state errors are $8.4\times10^{-5}$ and $4.4\times10^{-5}$ for the two signs.*<!-- source: study/results/circuit_detuning_convergence.csv (step_counts 8-128; phase_aligned_state_error at 128 steps: 8.39e-5 for -1.3, 4.43e-5 for +1.3); figure generated by study/circuit_controls.py -->

**Scope.** Qiskit serves as an independent check of the operator, the bit ordering, the sign of the detuning and the time evolution for two to five atoms. It is not a route to large systems: `HamiltonianGate` exponentiates a dense $2^N\times2^N$ matrix, whose memory grows as $4^N$, and state-vector simulation of the split-step circuits stores $2^N$ complex amplitudes, the same exponential cost as ED. All circuits were simulated classically and without noise; none was executed on a quantum processor.

### 2.8 The software package

**Contents.** The repository separates a frozen record, a maintained package and exploratory work (Table M3).

*Table M3. Repository layout and the role of each part.*

| Directory | Contents | Status |
|---|---|---|
| `study/` | Historical model (`rydberg_model.py`: couplings, sparse Hamiltonian, time evolution, observables, Qiskit operator and circuits), scan scripts, two teaching notebooks, 61 result files, 25 figures and two regression references | Frozen; verified against a SHA-256 manifest by `study/verify_delivery.py` |
| `rydberg_chain/` | Maintained package (version 0.1.0): `geometry`, `ed`, `observables`, `effective`, `dmrg`; `infinite` (iDMRG on pair sites) | Core modules maintained and tested; `infinite` in development, not yet validated |
| `tests/` | 43 `pytest` tests of the core modules: builder identical to the historical one, hand-built two-atom matrix, non-interacting limit, eigensolver agreement, saved references, reflection symmetry, observables against brute force, DMRG against ED, effective couplings. Tests for `infinite` and for the external benchmarks are being added | Core tests all pass |
| `validation/` | Cost benchmark of ED and DMRG (`scaling_benchmark.py`, results in `results/scaling_benchmark.json`); DMRG benchmarks against the critical transverse-field Ising chain and the hard-boson $\mathbb Z_2$ transition (`benchmark_uniform.py`); iDMRG validation (`idmrg_validation.py`) | In progress; only the cost benchmark has results so far |
| `research/pilots/2026-09-27/` | Exploratory ring ED (pair basis, momentum and inversion resolved), constrained ED, the failed iDMRG scan, boundary-wall checks, raw outputs | Exploratory; numbers preliminary |
| `docs/` | Methods notes, validation record, results and limitations, theory notes, research plan | Documentation |

<!-- source: study/ counts from docs/VALIDATION.md ("61 result files, 25 figures and two regression-reference files"); package version from rydberg_chain/__init__.py; test count from a pytest run on 27 Sept (43 passed in 23.8 s, OMP_NUM_THREADS=1; rerun by the checker: 43 passed in 25.5 s). Checker, 27 Sept, later: the whole tests/ directory now holds 63 tests; 60 pass, 1 is skipped (needs validation/results/benchmark_uniform.json) and 2 in tests/test_infinite.py fail (see Sec. 2.4). Update this table when infinite.py and the benchmarks are finished. Module lists from the files themselves. -->

The maintained modules are small and documented. `geometry` returns positions (Eq. 2.3) and the coupling matrix $V_{ij}$, with optional cutoff, for points on a line or in the plane. `ed` builds the diagonal energies and the CSR Hamiltonian of Eq. (2.6), optionally with a site-dependent detuning (`local_delta`, used for pinning or disorder), and returns the lowest eigenpairs with their residuals (`lowest_states`, methods `dense`, `eigsh` and `lobpcg`). `observables` evaluates Eqs. (2.14), (2.15) and (2.17), connected correlations, bond probabilities and pair-state probabilities on full state vectors, using the bit convention of Sec. 2.1; the structure factors of Eq. (2.16) are so far computed only by the pilot ring scripts. `effective` implements Eqs. (2.21)–(2.22) and Pfeuty's order parameter of the transverse-field Ising chain. `dmrg` builds the atom-site and pair-site MPOs with TeNPy, runs two-site DMRG and returns the energy, the MPS and a dictionary of observables and convergence diagnostics. `infinite`, still in development, will provide iDMRG on pair sites (Sec. 2.4).

**From positions to a figure.** A calculation proceeds in six steps (Fig. M1). (1) Positions are generated from $N$ and $\varepsilon$, and $V_{ij}$ from $g$ and the cutoff. (2) The Hamiltonian is represented as a sparse matrix, an MPO or a Pauli operator. (3) A solver returns eigenpairs with residuals, a DMRG state with its truncation error, variance and convergence flags, or a time-evolved state. (4) Observables are evaluated. (5) Results are written to CSV or JSON together with the parameters, solver settings, library versions, machine and diagnostics; the benchmark file also stores the source of the script that produced it. (6) Figures are drawn from the saved files only. A minimal ED calculation for a 12-atom dimerized chain reads:

```python
from rydberg_chain import geometry, ed, observables as ob

x = geometry.chain_positions(12, eps=0.3)           # open dimerized chain, a = 1
V = geometry.coupling_matrix(x, g=11.39)             # all pairs, g / r^6
H = ed.build_hamiltonian(V, omega=1.0, delta=2.5)    # sparse CSR, dimension 2^12
E, v, res = ed.lowest_states(H, k=2)                 # energies, vectors, residuals
n_i = ob.site_occupations(v[:, 0], 12)               # <n_i>
S = ob.entanglement_entropy(v[:, 0], 12, cut=6)      # half-chain entropy in bits
```

<!-- Checked 27 Sept (scratch check): runs and gives E = (-15.114421, -15.037077), residuals <= 1.2e-12. -->

The corresponding DMRG call for the 50-atom benchmark of Table M2 is `dmrg.pair_chain_ground_state(eps=0.3, g=11.39, delta=2.5, n_pairs=25, chi=64, mmax=8, max_sweeps=30, variance=True)`<!-- source: validation/scaling_benchmark.py, lines 77-78 (max_sweeps=MAX_SWEEPS = 30, variance=True); the package default is max_sweeps = 60 -->.

**Running the code.** With Python 3.11 or later, `pip install -r requirements.txt` installs the pinned dependencies of the package (NumPy, SciPy, Matplotlib, TeNPy and pytest), and `python -m pytest tests/test_geometry_ed.py tests/test_dmrg.py tests/test_effective.py` runs the 43 core tests from the repository root in about 25 s on one thread<!-- source: requirements.txt; pytest run 27 Sept: 43 passed in 23.76 s with OMP_NUM_THREADS=1; checker rerun 25.5 s. A bare `python -m pytest` now also runs the in-development infinite-DMRG and benchmark tests (about 105 s, 2 failures as of 27 Sept) -->. The historical study uses its own pinned environment (`study/requirements.txt`, Python 3.12, including Qiskit 2.5.2), is checked with `python study/verify_delivery.py`, and is reproduced with `python study/run_project.py --section chains` (other sections: `foundations`, `optional`, `all`); each managed run writes to a new directory under `study/runs/` and never overwrites the delivered data. Random start vectors are seeded, so repeated ED runs give identical results. The maintained package was tested with Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1 and TeNPy 1.1.1 on the 4-core Intel Xeon machine of Table M1; the historical study was run with Python 3.12.14, NumPy 2.5.3, SciPy 1.18.1 and Qiskit 2.5.2 on an Apple arm64 computer<!-- source: validation/results/scaling_benchmark.json, machine block; study/results/summary.json, versions and benchmark_machine; study/results/extension_summary.json, versions -->.

**Scope and limitations of the methods.** All calculations treat the ideal, coherent Hamiltonian of Eq. (1): they omit spontaneous decay from the Rydberg state, laser phase noise, position disorder of the atoms and state-preparation and detection errors, which a comparison with experiment would have to add. Full-basis ED is limited to about 22–24 atoms. Open chains carry boundary fields, which in even chains force a domain wall. The ring calculations that enter Sec. 3.2 currently come from pilot scripts, some in restricted bases whose systematic errors (Sec. 2.2) are not yet controlled (P1, P3). DMRG runs of 50 atoms and more have not yet met the convergence criteria in the ordered phase of even chains, and no validated infinite-DMRG result is available (P6, P7). The pseudo-spin theory is a second-order perturbative expansion whose range of validity is established only by comparison with the numerics.

<!--
Open items for this section (for Ethan and for the computations in docs/RESEARCH_PLAN.md)
1. The d = 3 pair basis shifts Delta_c at eps = 0.3, g = 11.39 by +0.0055 (1.3257 -> 1.3312, 16/20-atom rings; scratch check methods_d3_vs_d4.py). The pilot Delta_c values at eps = 0.3, 0.45, 0.6, 0.7 quoted in docs/, 00_front.md and the Results used d = 3; recompute them with d = 4 (or with the package) in P3, and quantify the shift versus eps (expected to fall roughly as 1/V_s; not yet checked).
2. Promote the scratch checks cited above (DMRG vs ED differences, pairk0 precision check, d = 3 vs d = 4 crossing, tail bound) into validation/ with saved outputs. The benchmark script is now saved as validation/scaling_benchmark.py (relative ROOT path); the copy embedded in scaling_benchmark.json as script_source still has the old hard-coded path and differs only in that line, so regenerate the JSON from the saved script at the next benchmark run.
3. Make Figures M1 and M2 (all Table M1/M2 data are now in scaling_benchmark.json, including N = 200). Correct the "notes" field of that JSON (and the same text in validation/scaling_benchmark.py, lines 27-35), which claims the energy is converged regardless of the entropy criterion; the N = 100 and 200 variances (1.9e-4, 4.3e-4) contradict it. The bulk |pair order| also drifts with N (0.42, 0.53, 0.79 at N = 50, 100, 200), another sign that these runs are unconverged.
4. rydberg_chain/infinite.py now exists (created after this section was first drafted) but does not pass its validation: in a checker run on 27 Sept, tests/test_infinite.py::test_idmrg_matches_finite_bulk_energy (-1.23261 vs -1.21962 per atom) and test_idmrg_order_is_translation_and_mirror_consistent failed, and validation/idmrg_validation.py has no results yet (P6). Update Sec. 2.4, Table M3 and the test counts once it passes. Symmetry-resolved and constrained ring ED are still pilot code (pairk0.py, crit.py, ced.py); reimplement in rydberg_chain with double-precision diagonals and random start vectors, and add tests against full ED (P1).
5. The maintained package has no time-evolution function; dynamics use study/rydberg_model.evolve. Either add one or keep citing the frozen module.
6. docs/DIMERIZED_THEORY.md says pair-site DMRG matches ED "to 2e-9 or better"; with the current defaults the difference is <= 8e-12 (scratch check) and 4e-13 at N = 20 (benchmark). Update the doc.
7. Check the [verify] references below and harmonize citation numbers with the other sections.
-->

### References for this section

<!-- Entries marked [verify] need a bibliographic check before submission. arXiv identifiers 2405.08810, 2408.02010 and 2601.16275 were confirmed against the local arXiv abstract archive (titles only). -->

[1] H. Bernien, S. Schwartz, A. Keesling, H. Levine, A. Omran, H. Pichler, S. Choi, A. S. Zibrov, M. Endres, M. Greiner, V. Vuletić, and M. D. Lukin, "Probing many-body dynamics on a 51-atom quantum simulator," Nature **551**, 579 (2017), arXiv:1707.04344.

[2] A. Browaeys and T. Lahaye, "Many-body physics with individually controlled Rydberg atoms," Nat. Phys. **16**, 132 (2020), arXiv:2002.07413.

[3] R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed. (Cambridge University Press, Cambridge, 2013): Weyl's inequalities (Sec. 4.3) and the Perron–Frobenius theorem (Ch. 8).

[4] C. Lanczos, "An iteration method for the solution of the eigenvalue problem of linear differential and integral operators," J. Res. Natl. Bur. Stand. **45**, 255 (1950).

[5] R. B. Lehoucq, D. C. Sorensen, and C. Yang, *ARPACK Users' Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods* (SIAM, Philadelphia, 1998).

[6] P. Virtanen et al., "SciPy 1.0: fundamental algorithms for scientific computing in Python," Nat. Methods **17**, 261 (2020), arXiv:1907.10121. Documentation of `eigsh`, `lobpcg` and `expm_multiply`: https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html.

[7] A. V. Knyazev, "Toward the optimal preconditioned eigensolver: locally optimal block preconditioned conjugate gradient method," SIAM J. Sci. Comput. **23**, 517 (2001).

[8] A. H. Al-Mohy and N. J. Higham, "Computing the action of the matrix exponential, with an application to exponential integrators," SIAM J. Sci. Comput. **33**, 488 (2011).

[9] A. W. Sandvik, "Computational studies of quantum spin systems," AIP Conf. Proc. **1297**, 135 (2010), arXiv:1101.3281.

[10] P. Fendley, K. Sengupta, and S. Sachdev, "Competing density-wave orders in a one-dimensional hard-boson model," Phys. Rev. B **69**, 075106 (2004), arXiv:cond-mat/0309438 [verify arXiv id].

[11] S. R. White, "Density matrix formulation for quantum renormalization groups," Phys. Rev. Lett. **69**, 2863 (1992).

[12] S. R. White, "Density-matrix algorithms for quantum renormalization groups," Phys. Rev. B **48**, 10345 (1993).

[13] U. Schollwöck, "The density-matrix renormalization group in the age of matrix product states," Ann. Phys. **326**, 96 (2011), arXiv:1008.3477.

[14] M. B. Hastings, "An area law for one-dimensional quantum systems," J. Stat. Mech. (2007) P08024, arXiv:0705.2024.

[15] J. Hauschild and F. Pollmann, "Efficient numerical simulations with Tensor Networks: Tensor Network Python (TeNPy)," SciPost Phys. Lect. Notes **5** (2018), arXiv:1805.00055.

[16] J. Hauschild et al., "Tensor Network Python (TeNPy) version 1," arXiv:2408.02010 (2024) [verify journal version; believed to be SciPost Phys. Codebases **41** (2024)].

[17] S. R. White, "Density matrix renormalization group algorithms with a single center site," Phys. Rev. B **72**, 180403(R) (2005), arXiv:cond-mat/0508709 [verify arXiv id].

[18] G. Temple, "The theory of Rayleigh's principle as applied to continuous systems," Proc. R. Soc. Lond. A **119**, 276 (1928) [verify].

[19] G. M. Crosswhite, A. C. Doherty, and G. Vidal, "Applying matrix product operators to model systems with long-range interactions," Phys. Rev. B **78**, 035116 (2008), arXiv:0804.2504 [verify arXiv id].

[20] B. Pirvu, V. Murg, J. I. Cirac, and F. Verstraete, "Matrix product operator representations," New J. Phys. **12**, 025012 (2010), arXiv:0804.3976 [verify arXiv id].

[21] I. P. McCulloch, "Infinite size density matrix renormalization group, revisited," arXiv:0804.2509 (2008).

[22] P. Calabrese and J. Cardy, "Entanglement entropy and quantum field theory," J. Stat. Mech. (2004) P06002, arXiv:hep-th/0405152.

[23] N. Laflorencie, E. S. Sørensen, M.-S. Chang, and I. Affleck, "Boundary effects in the critical scaling of entanglement entropy in 1D systems," Phys. Rev. Lett. **96**, 100603 (2006).

[24] D. Friedan, Z. Qiu, and S. Shenker, "Conformal invariance, unitarity, and critical exponents in two dimensions," Phys. Rev. Lett. **52**, 1575 (1984).

[25] P. Di Francesco, P. Mathieu, and D. Sénéchal, *Conformal Field Theory* (Springer, New York, 1997).

[26] J. L. Cardy, "Conformal invariance and universality in finite-size scaling," J. Phys. A **17**, L385 (1984).

[27] H. W. J. Blöte, J. L. Cardy, and M. P. Nightingale, "Conformal invariance, the central charge, and universal finite-size amplitudes at criticality," Phys. Rev. Lett. **56**, 742 (1986).

[28] I. Affleck, "Universal term in the free energy at a critical point and the conformal anomaly," Phys. Rev. Lett. **56**, 746 (1986).

[29] K. Slagle, D. Aasen, H. Pichler, R. S. K. Mong, P. Fendley, X. Chen, M. Endres, and J. Alicea, "Microscopic characterization of Ising conformal field theory in Rydberg chains," Phys. Rev. B **104**, 235109 (2021), arXiv:2108.09309 [verify authors, volume and page].

[30] X. Sun, Y. Le, S. Naus, R. B.-S. Tsai, L. R. B. Picard, S. Murciano, M. Knap, J. Alicea, and M. Endres, "Experimental observation of conformal field theory spectra," arXiv:2601.16275 (2026) [verify authors and journal version; project notes list Nature (2026)].

[31] P. Pfeuty, "The one-dimensional Ising model with a transverse field," Ann. Phys. **57**, 79 (1970).

[32] A. Javadi-Abhari et al., "Quantum computing with Qiskit," arXiv:2405.08810 (2024). API documentation for `SparsePauliOp`, `HamiltonianGate`, `RXGate`, `PhaseGate` and `CPhaseGate`: https://quantum.cloud.ibm.com/docs/en/api/qiskit.

[33] H. F. Trotter, "On the product of semi-groups of operators," Proc. Am. Math. Soc. **10**, 545 (1959).

[34] G. Strang, "On the construction and comparison of difference schemes," SIAM J. Numer. Anal. **5**, 506 (1968).

[35] M. Suzuki, "Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems," Commun. Math. Phys. **51**, 183 (1976).

[36] A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, and S. Zhu, "Theory of Trotter error with commutator scaling," Phys. Rev. X **11**, 011020 (2021), arXiv:1912.08854.
