# Next steps as ready-to-use prompts (state as of 28 September 2026)

Each prompt below is self-contained. Paste one into a new Claude session (Cowork or Claude Code) that has this repository open. Run them in order, because each one depends on the results of the one before.

**Where things stand.**
- **Done:**
  - Package, tests and benchmarks.
  - Infinite-DMRG validation (P6).
  - Literature gate (P0).
  - Full manuscript draft in `paper/DRAFT.md`.
- **Open:** P1–P5 and P7–P9 in `docs/RESEARCH_PLAN.md` §3.
- **Headline decision:** P1 decides whether the headline, the Z₄ → PDW split, stands.
- **Where the work lives:** until PR #1 is merged, it is on branch `claude/vigilant-mccarthy-xvgdqs`. After a squash-merge it is on `main`.

**Which tool.**
- Heavy numerics (P1, P2, P5, P7) need a machine that can run 2²⁰-state ED and long DMRG jobs. Claude Code in a terminal is the safest choice for those.
- Writing, meeting preparation and Google-Doc text suit either tool.

**Rules that apply to every prompt** (also in `CLAUDE.md`):
- Never modify `study/`.
- Give every new number in `paper/DRAFT.md` a `<!-- source: ... -->` comment.
- Mark anything unchecked `[verify]`.
- Do not write into the Google Doc.
- Keep personal material out of this public repository.
- Commit under Ethan's GitHub identity, keeping the AI co-author line.

---

## 0. Orientation (use this first in any new chat)

```text
You are helping me (Ethan) with a research project on dimerized Rydberg-atom chains, supervised by a professor. The repository is Sineceptor/rydberg-atom-research. If PR #1 is not merged yet, work on branch claude/vigilant-mccarthy-xvgdqs; otherwise work on main.

Before doing anything else, read these in order: CLAUDE.md, docs/RESEARCH_PLAN.md, docs/DIMERIZED_THEORY.md, the "Notes for Ethan" preface at the top of paper/DRAFT.md, and docs/NEXT_STEPS_PROMPTS.md. Then tell me in plain language:
(1) what the project claims and what evidence exists so far;
(2) which computations (P1–P9) are still open and their due dates;
(3) the three biggest risks to finishing by 15 October.

Do not change any files in this step. Follow every rule in CLAUDE.md.
```

## 1. P1: check the restricted basis against full ED (go/no-go, due 29–30 Sept)

```text
Task P1 from docs/RESEARCH_PLAN.md: decide whether the Z4 -> pair-density-wave (PDW) evidence survives without the restricted basis.

Background: research/pilots/2026-09-27/scripts/ced.py diagonalizes the periodic dimerized ring in a basis that forbids two excitations within two sites (excl=2). Its ring couplings sum periodic images m = -3..3 (Vmat). It labels states by T2 (translation by one pair) and I_s (inversion about a strong-bond centre, site i -> 1 - i mod N). The pilot numbers are in paper/DRAFT.md contribution 3 and Sec. 3.2.3, and in research/pilots/2026-09-27/drafting_checks/ced_check_D3.txt.

Do this:
1. Write validation/p1_constrained_vs_full.py. For each point, build the SAME ring Hamiltonian in the full 2^N basis and in the excl=2 basis, reusing ced.py's positions/Vmat convention. Use the sparse builder in rydberg_chain/ed.py where possible. Resolve T2 and I_s, either by projecting into symmetry sectors or by measuring them on the eigenvectors.
2. Grid: N = 16 and 20; R_b/a = 3.3 (g = 3.3**6); eps in {0, 0.3, 0.5}; Delta/Omega in {2, 3, 5}. Add open chains (N = 16, 20) at the same points. Add the Z3 control at R_b/a = 2.7, Delta = 5, eps in {0, 0.3, 0.6}, on an N = 18 ring (compatible with both pairs and period 3). Compare it against the adjacent-only exclusion basis (excl=1) that the Z3 pilot runs used.
3. For each point, record the lowest 6 energies with their (T2, I_s) labels in both bases, the inversion-odd gap, and the splitting of the two lowest I_s-even states with T2 = +1 and -1. Save to validation/results/p1_constrained_vs_full.json. Run long jobs in the background, and use LOBPCG for N = 20 (an earlier N = 16 attempt timed out at 600 s with the default solver).
4. Pass criteria: the same level ordering and symmetry labels in both bases; the inversion-odd gap agreeing to within 10% or 0.005 (whichever is larger); the PDW signature at eps = 0.5 (two lowest states I_s-even with T2 = +/-1) present in the full basis; Z3 threefold degeneracy present at N = 18.
5. Write validation/P1_NOTES.md with a table and a clear verdict: PASS, PARTIAL (say which parts), or FAIL.
6. Update the P1 row in docs/RESEARCH_PLAN.md and the preface table in paper/DRAFT.md. Replace the [TABLE PENDING] in Sec. 2.2 with the comparison table, with source comments. If the verdict is FAIL, do not edit the claims yet: list what would change under the fallback title in RESEARCH_PLAN.md.
7. Run the test suite, then commit and push.

Report the verdict and the three most important numbers in plain language.
```

## 2. P2: Z₄/PDW gap map (due 2 Oct; only if P1 passes)

```text
Task P2 from docs/RESEARCH_PLAN.md: map the extent of the PDW phase at R_b/a = 3.3.

Using the basis validated in P1 (see validation/P1_NOTES.md), compute on periodic rings of N = 16, 20, 24, 28, 32 and 36 atoms, on a 12 x 8 grid of (Delta/Omega, eps). Choose Delta between 1.5 and 7 and eps between 0 and 0.7, and refine the grid near the boundaries you find. Add N = 40 on three cuts: eps = 0.3, eps = 0.5, and Delta = 3.
- Resolve T2 and I_s quantum numbers for the lowest states. Count degeneracies to classify each point: disordered; PDW (two quasi-degenerate I_s-even states, T2 = +/-1); or Z4 crystal (fourfold, including an I_s-odd state).
- Locate phase boundaries from crossings of N * gap in the relevant symmetry channel. For the tau-ordering line, use the inversion-odd gap. For the PDW -> disordered line, use the T2-odd gap.
- Save validation/results/p2_gap_map.json and a script that makes Figure 10b and Figure 12 of paper/DRAFT.md (see the [FIGURE PENDING] blocks for what each panel should show). Use rydberg_chain/plotstyle.py.
- Check whether the two Ising lines approach each other as eps -> 0, as the Ashkin-Teller argument in docs/DIMERIZED_THEORY.md section 4 predicts (anisotropy linear in eps).
- Update Table 5, Sec. 3.2.3, the bracketed [R_b/a = 3.3] in the abstract, and the P2 row, all with source comments. Rebuild paper/DRAFT.docx with the pandoc command in CLAUDE.md.
- Commit and push. Summarize the phase diagram in five sentences.
```

## 3. P3: the Z₂ line with the maintained package (due 2 Oct)

```text
Task P3 from docs/RESEARCH_PLAN.md: recompute the Z2 melting line with the maintained package instead of the pilot scripts.

For g in {4.83, 11.39, 40} (R_b/a = 1.30, 1.50, 1.85) and eps in {0, 0.15, 0.3, 0.45, 0.5, 0.55, 0.6, 0.7}:
- Locate Delta_c from N * gap crossings on rings of 16, 20 and 24 atoms in the full pair basis (4 states per pair, no truncation), resolved by pair momentum as in research/pilots/2026-09-27/scripts/pairk0.py. Report the gap ratio R at the crossing (Ising value 8), and the central charge from the ground-state energy scaling E0(N).
- Compare with the pilot values in paper/DRAFT.md Table 4. Pilot Delta_c at g = 11.39: 1.0785, 1.326, 1.805, 3.149, 5.611 at eps = 0, 0.3, 0.45, 0.6, 0.7. Quantify the shift caused by the three-state pair basis used in the pilots (about 0.4% at eps = 0.3).
- Save validation/p3_z2_line.py and validation/results/p3_z2_line.json. Make Figures 9a, 10a and 11.
- Update Table 4, contribution 2, the abstract brackets, Secs. 3.2.2, 3.2.4 and 3.2.5, and the P3 row, with source comments. Rebuild the docx, commit and push.
```

## 4. P4: effective theory with a controlled derivation (due 3 Oct)

```text
Task P4 from docs/RESEARCH_PLAN.md: turn the pseudo-spin effective model into a controlled theory.

docs/DIMERIZED_THEORY.md section 2 derives J_m exactly at zeroth order in Omega. It gives a leading-order transverse field h0 and an empirical "dressed" field that fits the Z2 line within 1% at eps = 0.6-0.7, but is 15% off at eps = 0.45 and fails at 0.3.
1. Do a second-order Schrieffer-Wolff transformation onto the sector with one excitation per occupied pair. Keep the interaction energy of the flipped excitation with its neighbours in the denominators, and derive the renormalized J and h and the symmetry-allowed odd term. Write each step so a high-school student can follow it and re-derive it by hand.
2. Implement the result in rydberg_chain/effective.py (keep the sign convention documented in CLAUDE.md), with tests in tests/test_effective.py.
3. Compare the predicted Delta_c(eps) with the P3 results. Produce the validity map (where the theory is within 1%, 5% and 20%) and the classical Omega -> 0 skeleton. Also produce the Blume-Emery-Griffiths K/J table that supports the "no tricritical point" argument.
4. Update docs/DIMERIZED_THEORY.md, Table 4's theory columns, contribution 2, Eq. (28) and its text (no longer "empirical" if the derivation works), and Figure 11. Commit and push.
Also list which equations I (Ethan) must re-derive myself before the paper goes out.
```

## 5. P5: Z₃ control (due 3 Oct)

```text
Task P5 from docs/RESEARCH_PLAN.md: confirm that the Z3 crystal is protected (no extra transition) under dimerization.

At R_b/a = 2.7 (g = 2.7**6) and Delta/Omega = 5, and on one extra cut in Delta through the Z3 melting transition:
- compute rings of N = 24, 30 and 36 atoms (multiples of 6, so both pairs and period 3 fit), for eps from 0 to 0.7 in steps of 0.1. Use the basis validated in P1.
- Check that the threefold degeneracy persists: report the splitting of the three lowest states and the gap to the fourth, with T2 and I_s labels. Test whether any additional level crossing appears as eps grows.
- Save the script and JSON under validation/, and make Figure 13 and panel (b) of Figure 9.
- Update Sec. 3.2.3 (last paragraph), contribution 3's [for eps up to 0.6], and the P5 row. Commit and push.
```

## 6. P7: large systems and central charge (due 6 Oct)

```text
Task P7 from docs/RESEARCH_PLAN.md: answer the professor's "larger systems" question with DMRG on hundreds of atoms.

Even open chains in the Z2 regime carry a boundary-forced domain wall, so use pinned chains (rydberg_chain.infinite.pinned_chain_ground_state; the last atom is detuned by -2 Delta) or odd numbers of atoms. Note that an odd number of pairs removes only the PDW wall, not the tau wall.
1. Z2 line at g = 11.39: at the P3 critical points for eps = 0.3 and 0.6, run chains of 100, 200 and 400 atoms. Extract the central charge from the Calabrese-Cardy formula after removing Friedel oscillations (as in arXiv:2411.05494, Sec. III C). Expected c = 1/2.
2. Two Z4 cuts at R_b/a = 3.3 through the transitions found in P2: extract c at each transition, and the order-parameter scaling dimension (Ising d = 1/8). Use Ornstein-Zernike fits for xi and q, and the product |q - pi/2| * xi to test for chiral behaviour near eps = 0 (see arXiv:2406.20093, Sec. 2.3).
3. Record convergence (energy change, variance, bond dimension) for every run. Report failures honestly.
4. Save scripts and JSON under validation/. Make Figures 14 and 19, and update Table 2, Table 6, Sec. 3.3.4, "Reach at large N" in Sec. 2.4, and outlook (i) of the Discussion. Commit and push.
Long runs: start them in the background and check progress; do not block for hours.
```

## 7. P8 and P9: boundary wall and hardware table (due 5 Oct)

```text
Two smaller tasks from docs/RESEARCH_PLAN.md.

P8, boundary-forced domain wall: research/pilots/2026-09-27/scripts/dwbox_check.py prints its results but saves nothing. Promote it to validation/p8_boundary_wall.py. Run ED for even open chains up to N = 20 in the Z2 regime (g = 11.39, a few eps and Delta values). Save the wall-position probabilities and the fitted hopping t to JSON, and compare the hard-wall box model with the projected model with end states (eff_model). Update Sec. 3.3.3 (Eqs. 30-31) and panel (a) of Figure 19.

P9, hardware table: for QuEra Aquila, check the parameters used in Sec. 3.2.5 against Ref. [77] and the current public device specification: minimum spacing 4 um, C6/h = 862,690 MHz um^6, Omega/2pi up to 2.5 MHz. Mark anything you cannot confirm [verify]. Tabulate, for R_b/a = 1.5, 2.7 and 3.3, the lattice constant a, the short gap s and the largest reachable eps at each Rabi frequency. State plainly which predictions current hardware can test (the weakly dimerized Z2 line) and which it cannot (the Z4/PDW regime needs gaps of about 1.3-2.1 um). Update Sec. 3.2.5 and outlook (ii).

Commit and push after each task.
```

## 8. Literature and citation clean-up (any time before 13 Oct)

```text
Clean up the references in paper/DRAFT.md.
1. List every reference marked [verify] and every inline "[verify ...]" note. For each, find the published journal version (volume, pages, year) and check the specific claim it supports. Arxiv full texts may need to be uploaded by me if the web is blocked; ask me for the PDFs you need.
2. Find a standard reference for the anisotropic Ashkin-Teller model (two Ising transitions when the symmetry exchanging the two Ising copies is broken). Cite it where docs/DIMERIZED_THEORY.md section 4 and the Discussion ("Ladders and tricritical points") use that argument.
3. Check the Fendley-Sengupta-Sachdev benchmark value U_c/w = -1.31 and its sign convention (Ref. [11], validation/BENCHMARK_NOTES.md Sec. 2).
4. Run one more search for papers from 2025-2026 on dimerized or superlattice Rydberg chains, "pair density wave" in Rydberg arrays, and anisotropic Ashkin-Teller physics in Rydberg systems. Report anything that overlaps with the claims in docs/RESEARCH_PLAN.md ("Novelty status").
Update the reference list and the P0 notes. Commit and push.
```

## 9. After the meeting with the professor

```text
I just met my professor about this project. Here are my notes from the meeting: [PASTE NOTES].

1. Turn each piece of feedback into a concrete change: to the plan (docs/RESEARCH_PLAN.md), to the draft (paper/DRAFT.md), or to a computation (P1-P9). Put them in a table with owner (me or Claude) and due date.
2. If she changed the direction, headline or deadline, update the plan and the "Title, and what to delete if a gate fails" section of the draft preface.
3. Do not put private details of the meeting into the public repository. Keep them in meeting_notes/ (gitignored).
Show me the table before editing any files.
```

## 10. Final assembly and self-audit (13–15 Oct)

```text
Prepare paper/DRAFT.md for sending to my professor.
1. Replace every bracketed number and [FIGURE PENDING]/[TABLE PENDING] with final values or figures from the saved results. Each number needs a source comment naming the JSON or script that produced it.
2. Self-audit: for every number in the abstract, the contributions and the tables, recompute or re-read it from its source file and list any mismatch. Run the test suite and python study/verify_delivery.py.
3. Check the claims against the "Not to be claimed" list in docs/RESEARCH_PLAN.md and the rules in CLAUDE.md. The PDW must be presented as the chain counterpart of the rung-resonant crystals of two-leg ladders.
4. Check that the Acknowledgements disclose the AI assistance accurately.
5. Produce a clean copy without the "Notes for Ethan" preface and without HTML comments, as paper/DRAFT_clean.md and paper/DRAFT_clean.docx, for me to paste into the Google Doc.
Commit and push, and list anything I must still check by hand.
```
