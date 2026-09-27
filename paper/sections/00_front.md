<!--
paper/sections/00_front.md: title, author line, abstract (draft) and project direction.
Status: draft of 27 September 2026, for Ethan to edit before moving it into the shared Google Doc.
Every number was re-checked against its source file on 27 September.

Conventions used in this file
- A phrase in square brackets that contains a number, e.g. "[at ε = 0.6 and 0.7]", is a pilot
  value (research/pilots/2026-09-27). Replace it with the production value from the maintained package
  (computations P1–P7 in docs/RESEARCH_PLAN.md) before submission, then remove the brackets.
  Plain integers in brackets, e.g. [3], are citations to the list at the end of this file.
  Renumber the citations when the sections are merged.
- Every quoted number has a "source" comment directly after it.
- Units: ħ = Ω = a = 1 in the source files; g = C6/(Ω a^6) = (R_b/a)^6, so g = 4.83, 11.39, 40, 387.4, 1291.5
  are R_b/a = 1.30, 1.50, 1.85, 2.70, 3.30. Pair counts in crit_*.json are pairs (atoms = 2 × pairs).
- Numbers from study/ (historical, checksum-protected) and validation/results/ (maintained package) are not
  bracketed, but still carry a source comment.

Title options
  (a) The provisional title below. It calls the pair-density wave a "candidate" because that result is preliminary
      (constrained basis, rings of at most 36 atoms). Drop "candidate" only if P1 and P2 pass.
  (b) The plan's working title, to be used only if gates P0, P1 and P2 pass:
      "Pair pseudo-spins in dimerized Rydberg chains: even-period crystals split, odd-period crystals are protected"
  (c) The fallback title, if P0 finds the Z4 split published or P1/P2 show the PDW is a basis artefact:
      "Dimerized Rydberg chains as quantum order-disorder ferroelectrics: pseudo-spin theory, phase boundary,
      and why geometry alone does not make them tricritical"
      With (c), delete contribution 3, the even/odd sentence of the abstract, and the second prediction
      in the project direction.

The acknowledgments and the AI-use disclosure (docs/RESEARCH_PLAN.md, section 6) go at the end of the manuscript.
-->

# Geometric dimerization of laser-driven van der Waals Rydberg chains: pseudo-spin theory and a candidate pair-density-wave phase

Ethan Zhang$^{1}$ and Brenda M. Rubenstein$^{2}$

$^{1}$[Affiliation to be added]
$^{2}$Department of Chemistry, Brown University, Providence, Rhode Island 02912, USA [confirm affiliations, author order and corresponding author]

## Abstract

Rydberg atoms in optical tweezers form chains of laser-driven two-level systems with repulsive van der Waals interactions. As the blockade radius grows, a uniform chain at large detuning orders into $\mathbb{Z}_2$, $\mathbb{Z}_3$ or $\mathbb{Z}_4$ density-wave crystals. We study how geometric dimerization, with gaps alternating between $a(1-\varepsilon)$ and $a(1+\varepsilon)$, reshapes this phase diagram, using exact diagonalization (ED) and density-matrix renormalization group (DMRG) calculations. The pseudo-spins $\tau = n_L - n_R$ of occupied strong pairs couple ferromagnetically through a second difference of $r^{-6}$, which vanishes as $\varepsilon \to 1$, while the laser acts as a transverse field that stays finite. In the $\mathbb{Z}_2$ regime this gives a single Ising melting line $\Delta_c(\varepsilon)$, which a transverse-field Ising model with dressed energy denominators reproduces [to within 1% at $\varepsilon = 0.6$ and 0.7]<!-- source: ED Δ_c from research/pilots/2026-09-27/outputs/crit_eps0.6_g11.39_d3.json (3.1494) and crit_eps0.7_g11.39_d3.json (5.6107), largest size pair (10/12 pairs = 20/24 atoms). Dressed prediction recomputed on 27 Sept (twice, independently) from the formula in docs/DIMERIZED_THEORY.md §2 with nearest-pair J_1: 3.178 (+0.9%) and 5.640 (+0.5%). The doc quotes "about 0.5%"; the error depends on the dressing prescription (0.2–1.1%). Only ε = 0.6 and 0.7 were computed, so the claim is not stated as "ε ≥ 0.6". Finalize in P4. --> but not [for $\varepsilon \le 0.45$]<!-- source: recomputed on 27 Sept: dressed error +15% at ε = 0.45 (crit_eps0.45_g11.39_d3.json, Δ_c = 1.8052; dressed 2.076) and +101% at ε = 0.3 (Δ_c = 1.3258; dressed 2.666). docs/DIMERIZED_THEORY.md §3 says "about +6%" at ε = 0.45 and "fails" at 0.3; the +6% could not be reproduced. Resolve in P4. -->. Preliminary ED in a blockade-constrained basis, still to be validated against full ED and on larger rings, suggests an even/odd rule: at [$R_b/a = 3.3$] an intermediate pair-density-wave phase, which breaks pair translation but preserves inversion about strong-bond centres, appears to separate the $\mathbb{Z}_4$ crystal from the disordered phase, whereas at [$R_b/a = 2.7$] the $\mathbb{Z}_3$ crystal keeps its threefold degeneracy<!-- source: g = 1291.5 = 3.3^6 in research/pilots/2026-09-27/outputs/out24_*.txt and out36_*.txt (ε = 0.5; Δ/Ω = 1.5 disordered, 2.0 PDW, 6.5 Z4-like) and the ced.py rerun at Δ/Ω = 3 (see contribution 3); g = 387.4 = 2.7^6 in z3_24_*.txt. -->. We find no tricritical point: the ratio of the two lowest critical gaps stays near the Ising value 8, far from the tricritical Ising value 8/3, [for $R_b/a = 1.3$–$1.85$ and $\varepsilon \le 0.7$]<!-- source: research/pilots/2026-09-27/outputs/crit_*.json (g = 4.83, 11.39) and results.jsonl (g = 11.39, 40); 4.83^(1/6) = 1.30, 40^(1/6) = 1.85. The results.jsonl rows with g = 100 and 150 (R_b/a = 2.15, 2.31) lie outside this range; they give R = 8.71 and 10.05 at 20 atoms and still drift with size, so they are not used. -->, as a classical-limit argument predicts.

## Project direction

<!-- This fills the outline placeholder "[Your Exact Project Direction]". -->

In the uniform Rydberg chain, the crystal that forms at large detuning changes from $\mathbb{Z}_2$ to $\mathbb{Z}_3$ and $\mathbb{Z}_4$ as the blockade radius $R_b = (C_6/\Omega)^{1/6}$ grows relative to the spacing $a$ [1]. The transitions out of these crystals include Ising, three-state Potts and chiral transitions, as well as intervening floating phases [2–4]. Tweezer arrays place each atom independently, so the spacing itself is a control parameter. This project asks: *how does geometric dimerization, with alternating gaps $s = a(1-\varepsilon)$ and $w = a(1+\varepsilon)$, reshape the ground-state phase diagram of the single-level, laser-driven van der Waals chain, and can one pair pseudo-spin rule predict where each crystal melts?* Dimerization keeps the form of the Hamiltonian and the atomic species. It changes the couplings and lowers the symmetry: translation by one atom is lost, and what remains is translation by one pair, $T_2$, and inversion about the centre of a strong bond, $I_s$.

The pseudo-spin picture makes three predictions, which we test against numerically exact results. First, every crystal whose occupied pairs carry a free left/right choice must lose its left/right order at some $\varepsilon < 1$ for fixed $\Delta$ and $R_b/a$. Second, even-period crystals acquire an additional transition, whereas odd-period crystals do not. Third, dimerization alone should not make the $\mathbb{Z}_2$ melting tricritical.

Related systems show a reorganization of density-wave order: number-conserving dimerized chains, with dipolar exchange between two Rydberg states [5] or with extended Bose–Hubbard hopping [6]; two-species chains [7]; and multicomponent arrays [8]. The closest analogues of the proposed splitting are two-leg ladders, in which a crystal that breaks the symmetry between the legs can be separated from the disordered phase by a leg-symmetric crystal and an additional Ising transition [9, 10], and quantum loop models, in which the transition out of a $\mathbb{Z}_4$ phase can proceed through a pair of Ising transitions [11]. Rydberg chains tuned to tricritical Ising points have been analysed theoretically [12], and tricritical Ising spectra have been observed experimentally [13]<!-- The plan (docs/RESEARCH_PLAN.md §1) states that [13] reached the tricritical point with an inverted sweep (effective V_2 < 0). The abstract in the local arXiv archive does not say this. Check the full text (P0) before stating it in the paper. -->. Our literature search is not yet complete, and so far it has found no study of the spacing-dimerized, single-level chain. The Results are organized in three parts: the coupling constant $R_b/a$, which selects the $\mathbb{Z}_2$, $\mathbb{Z}_3$ or $\mathbb{Z}_4$ regime; the alternating coupling set by $\varepsilon$; and system size. Qiskit is used only for independent checks on small systems. All results are zero-temperature ground-state properties of the ideal Hamiltonian, without decoherence, positional disorder or state-preparation errors.

**Contributions.**

1. **Pair pseudo-spin theory.** Grouping each strong pair into one four-state site rewrites the chain exactly. When no strong pair is doubly excited, two occupied pairs $m$ pair spacings apart couple, to zeroth order in $\Omega$, through
   $$J_m(\varepsilon) = \tfrac14\left[V(2ma-s) + V(2ma+s) - 2V(2ma)\right], \qquad V(r) = \frac{C_6}{r^6}.$$
   This coupling is positive (ferromagnetic in $\tau$) because $r^{-6}$ is convex, and it vanishes as $\varepsilon \to 1$. The laser-induced transverse field stays finite. Because $I_s$ flips $\tau$ while $T_2$ moves the pattern by one pair, the remaining symmetries act differently on crystals of even and odd period (contribution 3).

2. **The $\mathbb{Z}_2$ melting line and where the theory is valid.** Symmetry-resolved ED of rings of [up to 24 atoms] locates $\Delta_c(\varepsilon)$ at $R_b/a = 1.5$. It rises [from $\Delta_c/\Omega = 1.079$ at $\varepsilon = 0$]<!-- source: research/pilots/2026-09-27/outputs/crit_eps0.0_g11.39_d4.json, crossing of 8- and 10-pair rings (16 and 20 atoms), Δ_c = 1.0785 --> [to $5.61$ at $\varepsilon = 0.7$]<!-- source: research/pilots/2026-09-27/outputs/crit_eps0.7_g11.39_d3.json, crossing of 10- and 12-pair rings (20 and 24 atoms), Δ_c = 5.6107. The d3 runs (ε ≥ 0.3) use a pair basis without doubly excited strong pairs; validate in P3. -->, and the gap ratios along it are consistent with Ising universality. At leading order the effective model underestimates $\Delta_c$, [by 63% at $\varepsilon = 0.3$ and by 7.6% at $\varepsilon = 0.7$]<!-- source: docs/DIMERIZED_THEORY.md §3 table; recomputed on 27 Sept from the leading-order formula (J_1 = h_0) against crit_eps0.3_g11.39_d3.json and crit_eps0.7_g11.39_d3.json: −62.7% and −7.6% (also −39.0% at ε = 0.45 and −16.9% at ε = 0.6) -->. With dressed energy denominators it agrees [to within 1% at $\varepsilon = 0.6$ and 0.7]<!-- source: as in the abstract; +0.9% at ε = 0.6 and +0.5% at ε = 0.7 -->. It fails at weaker dimerization, where empty pairs rather than pseudo-spin flips appear to drive the transition<!-- source: docs/DIMERIZED_THEORY.md §3 (interpretation, not yet tested directly) -->. On the disordered side of the line at strong dimerization, each strong pair shares one excitation between its two atoms; this regime is a crossover within the disordered phase, not a separate phase<!-- source: docs/DIMERIZED_THEORY.md §3: no U(1) symmetry, and the pair ground state a|00⟩ + b(|10⟩ + |01⟩)/√2 is I_s-even -->.

3. **An even/odd rule and a candidate pair-density-wave phase (preliminary).** ED in a basis that forbids simultaneous excitation of atoms one or two sites apart was run at $R_b/a = 3.3$ and $\Delta/\Omega = 3$ on rings of 16, 24 and 32 atoms. The inversion-odd gap closes with system size [at $\varepsilon = 0$ and 0.2]<!-- source: docs/DIMERIZED_THEORY.md §4 and research/pilots/2026-09-27/scripts/ced.py (unchanged), rerun on 27 Sept; output in research/pilots/2026-09-27/drafting_checks/ced_check_D3.txt. Inversion-odd gap at N = 16/24/32: 0.0111/0.0020/0.0004 (ε = 0) and 0.0270/0.0111/0.0053 (ε = 0.2). Re-verified line by line in the check pass. -->, scales as $1/N$ [at $\varepsilon = 0.3$], and saturates [at $\varepsilon = 0.4$ and 0.7]<!-- source: same rerun: 0.0460/0.0294/0.0219 (ε = 0.3, N·gap = 0.74/0.71/0.70); 0.0749/0.0641/0.0606 (ε = 0.4) and 0.1660/0.1644/0.1643 (ε = 0.7). -->. Where it saturates, two quasi-degenerate ground states remain. Both are inversion-even and they have $T_2 = \pm 1$, as expected for a pair density wave<!-- source: same rerun: the T2 = −1, Is = +1 state lies 0.0059/0.0007/0.00009 (ε = 0.4) and 0.0036/0.0003/0.00003 (ε = 0.7) above the ground state at N = 16/24/32. -->. At $R_b/a = 2.7$ and $\Delta/\Omega = 5$ the threefold degeneracy of the $\mathbb{Z}_3$ crystal on a 24-atom ring persists [for $\varepsilon$ up to 0.6]<!-- source: research/pilots/2026-09-27/outputs/z3_24_0.0.txt, z3_24_0.3.txt and z3_24_0.6.txt (N = 24, Δ/Ω = 5, g = 387.4). These runs exclude only adjacent excitations (basis dimension 103,682 = Lucas number L_24), unlike the Z4 runs, which also exclude excitations 2a apart. The three lowest levels agree to the printed precision of 1e-5; the next level is 2.79, 2.34 and 0.60 above. Include the Z3 point in P1 anyway (the excluded adjacent pairs cost V_s ≥ 387Ω). -->. These results remain preliminary until the constrained basis is checked against full ED (computation P1) and the degeneracy counts are confirmed on larger rings and with DMRG (P2, P5, P7).

4. **No tricritical point along the $\mathbb{Z}_2$ line.** At the largest rings the gap ratios lie [between 7.9 and 8.15]<!-- source: largest-size values in research/pilots/2026-09-27/outputs/crit_*.json (R_b = 7.915–7.998 at 16–24 atoms) and results.jsonl rows with g ≤ 40 (N10 = 20 atoms: 7.926–8.147). Over all sizes the range is 7.84–8.32. --> [for $R_b/a = 1.3$–$1.85$ and $\varepsilon \le 0.7$]<!-- source: same files -->. This is close to the Ising value 8 and far from the tricritical Ising value 8/3. A classical-limit argument accounts for this. Aligned occupied pairs always sit $2a$ apart, whatever the value of $\varepsilon$, so they repel with $V(2a) > 0$. Dimerization therefore cannot supply the net attraction between them that a first-order segment, and hence a tricritical point, would require. This is a criterion for the classical limit, not a proof.

The system-size section reports two supporting results. The first is the accuracy and cost of each solver as the chain grows. For example, pair-site DMRG reproduces the full-ED ground-state energy of a 20-atom open chain to within $5\times10^{-13}\,\Omega$<!-- source: validation/results/scaling_benchmark.json, dmrg_pair[0].E0_minus_ED = −4.16e-13 (ε = 0.3, g = 11.39, Δ/Ω = 2.5, χ = 64, all-range couplings). The 50-atom case met the energy criterion (variance 3.7e-9) but not the entropy criterion within 30 sweeps. The 100- and 200-atom runs have finished but are not converged ground states (variance 1.9e-4 and 4.3e-4); they give only the cost per sweep. -->. The second is a domain wall that the boundaries force into even open chains in the $\mathbb{Z}_2$ regime. It accounts for an earlier result for a uniform eight-atom chain at $R_b/a = 1.41$ and $\Delta/\Omega = 3$: 77% of the ground-state weight<!-- source: study/results/parity_summary.json, N = 8, g = 8 (R_b/a = 8^(1/6) = 1.414), Δ/Ω = 3, cutoff "all", level 0: probability_exactly_one_00_and_no_11 = 0.7668 --> lies on configurations with exactly one 00 bond and no 11 bond.

### References for this section

[1] H. Bernien, S. Schwartz, A. Keesling, H. Levine, A. Omran, H. Pichler, S. Choi, A. S. Zibrov, M. Endres, M. Greiner, V. Vuletić, and M. D. Lukin, "Probing many-body dynamics on a 51-atom quantum simulator," *Nature* **551**, 579–584 (2017), arXiv:1707.04344.

[2] P. Fendley, K. Sengupta, and S. Sachdev, "Competing density-wave orders in a one-dimensional hard-boson model," *Phys. Rev. B* **69**, 075106 (2004), arXiv:cond-mat/0309438.

[3] A. Keesling, A. Omran, H. Levine, H. Bernien, H. Pichler, S. Choi, R. Samajdar, S. Schwartz, P. Silvi, S. Sachdev, P. Zoller, M. Endres, M. Greiner, V. Vuletić, and M. D. Lukin, "Quantum Kibble–Zurek mechanism and critical dynamics on a programmable Rydberg simulator," *Nature* **568**, 207–211 (2019), arXiv:1809.05540.

[4] M. Rader and A. M. Läuchli, "Floating phases in one-dimensional Rydberg Ising chains," arXiv:1908.02068 (2019).

[5] Z. Zeybek, R. Mukherjee, and P. Schmelcher, "Quantum phases from competing van der Waals and dipole-dipole interactions of Rydberg atoms," *Phys. Rev. Lett.* **131**, 203003 (2023), arXiv:2303.17470. [verify volume and article number]

[6] Z. Zeybek, P. Schmelcher, and R. Mukherjee, "Bond-order density wave phases in dimerized extended Bose–Hubbard models," arXiv:2403.06649 (2024). [verify whether a journal version exists]

[7] L.-Y.-N. Liu, S.-Y. Yu, S.-R. Peng, J. Sheng, S. Yi, P. Xu, S.-S. Gong, T. Shi, and J. Cui, "Novel ground states and emergent quantum many-body scars in a two-species Rydberg atom array," arXiv:2408.15965 (2024).

[8] N. Chepiga, "Tunable quantum criticality in multicomponent Rydberg arrays," *Phys. Rev. Lett.* **132**, 076505 (2024), arXiv:2308.12838. [verify journal reference; not yet read in full (P0)]

[9] J. Soto Garcia and N. Chepiga, "Numerical investigation of quantum phases and phase transitions in a two-leg ladder of Rydberg atoms," arXiv:2411.05494 (2024).

[10] [Authors to be added], "Phase diagram of Rydberg atoms in a two-leg rectangular ladder," arXiv:2412.11201 (2024). [verify authors]

[11] B. M. La Riviere and N. Chepiga, "$\mathbb{Z}_4$ transitions in quantum loop models on a zig-zag ladder," arXiv:2406.20093 (2024).

[12] S. Naus, Y. Liu, S. Murciano, P. Sala, M. Endres, and J. Alicea, "Practical roadmap to measurement-altered criticality in Rydberg arrays," arXiv:2506.21963 (2025).

[13] X. Sun, Y. Le, S. Naus, R. B.-S. Tsai, L. R. B. Picard, S. Murciano, M. Knap, J. Alicea, and M. Endres, "Experimental observation of conformal field theory spectra," arXiv:2601.16275 (2026). [verify: docs/RESEARCH_PLAN.md lists it as published in *Nature* (2026); volume and pages not checked]

<!-- Reference check, 27 Sept: titles of [6], [7], [9]–[13] match the local arXiv abstract archive; authors of [4]–[9], [11]–[13] match the arXiv listing in the drafting scratchpad (kaneko/MyTodaysArxiv_*.txt). [1]–[3] are standard and were checked from memory. Authors of [10] are not in either local source. -->
