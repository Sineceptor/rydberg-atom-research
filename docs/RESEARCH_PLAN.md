# Research plan to 15 October 2026

Status: proposal, drafted 26 September 2026. Nothing below has been run beyond the pilot checks listed in section 1. The direction must be confirmed with Prof. Rubenstein before work starts.

## 1. Where the project stands

**Existing work is careful but not new.** The package validates exact diagonalization for N ≤ 12 against Qiskit and studies blockade, defects and entanglement, small geometries, interaction cutoffs, odd/even chains and a 12-site kagome patch. All of this is established physics (Bernien et al. 2017; Rader & Läuchli 2019). It supplies a foundations section and methods, not a novelty claim.

**The paper document** has an introduction and headings. The abstract, project-direction paragraph, Methods and all three Results sections are empty. Open comments from Prof. Rubenstein:

1. Insert plots, with captions.
2. Explain what is in the package and how it works.
3. Explain what solves the Hamiltonian, with key equations.
4. Expand the introduction: phases and quantum computing (Lukin group papers).
5. Check how easy it is to go to larger systems.

**Pilot checks (26 Sept, not committed as code):**

- New DMRG reproduces the saved N = 12 energy (−18.399305357539) and matches exact diagonalization to about 10⁻¹², including dimerized chains.
- Exact-diagonalization cost: 4 s (N = 16), 25 s (N = 18), 196 s and 2.8 GB (N = 20). The practical ceiling is N ≈ 22; DMRG reaches hundreds of atoms. This answers comment 5.
- Exact diagonalization at N = 16, R_b/a = 1.5, Δ/Ω = 2.5: the crystal melts as dimerization grows. Interior |⟨n_L − n_R⟩| is 0.21 → 0.10 → 0.02 at ε = 0.25 → 0.6 → 0.7. The effective-model estimate gives melting near ε ≈ 0.6.
- **Blocking issue:** infinite-chain DMRG is unconverged beyond ε ≈ 0.3.

## 2. Recommended direction

**Working title:** *Melting a Rydberg crystal by dimerization: effective Ising theory and large-scale simulations of alternating-spacing atom chains.*

This fills the outline's existing heading "How the Physics Varies with Alternating Coupling Constants".

| Option | Verdict |
|---|---|
| Kagome / spin liquid | Crowded field; needs large 2D simulations; not feasible by 15 Oct |
| Interaction-cutoff study | Not new; reuse it to justify keeping all tails |
| End-detuning boundary study | Incremental; kept as fallback |
| **Dimerized van der Waals chain** | Matches the outline, feasible, apparently open |

**Novelty case.** Published dimerized Rydberg studies add dipolar exchange between two or three Rydberg levels (Zeybek, Mukherjee & Schmelcher, PRL 131, 203003 (2023); Mögerle et al., PRX Quantum 2025; de Léséleuc et al., Science 2019) or use ladders. No phase diagram or effective theory was found for the single-level, van der Waals–only chain. This is based on a search-engine scan only; arXiv was unreachable from the analysis environment. Prof. Rubenstein must confirm it, and the Zeybek et al. paper must be read in full.

**Scientific content:**

1. Exact pair-site formulation (each strong pair becomes one 4-state site).
2. Effective transverse-field Ising model ([DIMERIZED_THEORY.md](DIMERIZED_THEORY.md)): |J| = [V_w + V(2a+s) − 2V(2a)]/4 and h = (Ω²/4)[1/Δ + 1/(V_s − Δ)].
   - Test the exact Ising order parameter m = (1 − h²/J²)^{1/8} against ⟨n_L − n_R⟩ from the simulations.
3. The same model's boundary fields explain the earlier finding that ~77% of the N = 8 weight has exactly one 00 bond.
4. Phase diagram in (Δ/Ω, ε) compared with the prediction.
5. Nature of the transition: Ising (central charge c = 1/2) is expected. At small Δ the problem resembles a Blume–Capel chain, so a tricritical point is possible. This is a stretch goal, not a promise.
6. Hardware window. A 4 µm minimum spacing limits ε to roughly 0.3–0.45 for realistic R_b, so scans must include regimes (smaller Δ or R_b/a) where melting occurs at accessible ε.

**Not to be claimed:** a new phase of matter, topological order or a spin liquid. The pair paramagnet is symmetric, so no symmetry sets it apart from the ordinary disordered phase.

## 3. Fixes to existing material

**Paper document:**

- Write the project-direction paragraph.
- Expand the introduction (phases; blockade gates and logical qubits; key papers: Jaksch 2000, Lukin 2001, Bernien 2017, Keesling 2019, Levine 2019, Ebadi 2021 and 2022, Semeghini 2021, Evered 2023, Bluvstein 2022 and 2024).
- Write Methods: Hamiltonian, exact diagonalization, DMRG/infinite DMRG, observables, the Calabrese–Cardy entanglement formula, validation, package map.
- Captioned figures; one numbered citation style; abstract last.

**Repository** (🔴 = blocks new results):

| # | Problem | Fix |
|---|---|---|
| 🔴1 | Infinite DMRG unconverged past ε ≈ 0.3 | Project out the doubly excited pair state with a perturbative correction; use several initial states and keep the lowest energy; tune the mixer schedule; compare with finite DMRG at N = 100–200. Gate: agreement to 10⁻⁶ at 10 points spanning the transition |
| 🔴2 | No external benchmark of the DMRG pipeline | Reproduce the published uniform-chain Z2 transition before any new claim |
| 3 | Slow exact-diagonalization builder (Python lists) | Vectorize; optionally use reflection symmetry |
| 4 | Duplicated ground-state helpers; legacy M diagnostic | One module plus automated tests; keep the reference regression |
| 5 | Requirements pin Mac-only and unrelated packages | Minimal core plus notebook extras; add TeNPy |
| 6 | Notation drift (Q, S, M; bit order) | One symbol, one stated convention, R_b/a alongside g |
| 7 | Figure titles carry caveats | Publication figures with panel labels and a consistent colorblind-safe style |
| 8 | Docs dominated by disclaimers | One limitations section; README states question and results |
| 9 | Cutoff study framed as a headline result | Present it as the reason all tails are kept; the V(2a) cancellation shows why |
| 10 | "Defect" is an ideal interaction switch-off | Label it as an idealized control, or redo it as a vacancy or local detuning |
| 11 | Kagome pilot off-topic | Keep in the repository as an appendix, not in the paper |

## 4. Work plan with checkpoints

- **A. Fix and benchmark** (items 1–2). If infinite DMRG is not validated by 1 Oct, use finite DMRG with N up to 200 and finite-size scaling.
- **B. Theory.** Ethan re-derives the effective model by hand, adds next-order corrections and plots the predicted melting line.
- **C. Phase diagram.** R_b/a ≈ 1.2, 1.35, 1.5; Δ/Ω from 0 to 5; ε from 0 to 0.7; adaptive grid. Record order, density, correlation length and entanglement. About 2,000 points, run overnight.
- **D. Transition character.** Central-charge fits at about 4 boundary points (N up to 200); first-order tests at small Δ.
- **E. System size and boundaries.** Exact diagonalization to N = 20 plus DMRG; odd/even chains; quantitative test of the boundary-field prediction.
- **F. Hardware mapping.** Accessible region on current analog hardware and the snapshot signature of melting. Optional: simulated adiabatic ramp.

**Planned figures:**

1. Setup and blockade
2. Uniform-chain benchmark
3. Pair-to-pseudo-spin mapping
4. Order versus ε, compared with the Ising formula
5. Phase diagram with the predicted line
6. Entanglement scaling
7. Finite-chain boundaries and parity
8. Hardware window

## 5. Schedule

| Dates | Work | Checkpoint |
|---|---|---|
| 27–28 Sept | Proposal email to Prof. Rubenstein; item A; read Zeybek 2023, Bernien 2017, Browaeys & Lahaye 2020 | Go-ahead |
| 29 Sept–1 Oct | Derivation (Ethan); code consolidation and tests; Methods and Introduction drafts | Infinite DMRG validated, or fallback |
| 2–5 Oct | Phase-diagram runs; figures 4–5 | Professor check-in, 5 Oct |
| 6–8 Oct | Items D–E; all captioned figures in the paper document | Figures final |
| 9–11 Oct | Results, discussion, abstract | Full draft to professor, 11 Oct |
| 12–15 Oct | Revisions, repository cleanup, reproducibility; arXiv only with professor approval | Final |

**Cut order if time slips:** ramp simulation, then the tricritical search, then the third blockade radius.

**If the direction is already published:** move the same machinery to the end-detuning study in [FUTURE_WORK.md](FUTURE_WORK.md).

## 6. Ownership and disclosure

Ethan should personally derive the effective model, make the Ising comparison and write the introduction and discussion. The paper and any application materials disclose AI assistance, as the repository already does. Peer review cannot finish by 15 October; a complete manuscript (and an arXiv preprint if the professor approves) is the realistic outcome.
