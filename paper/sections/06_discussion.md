<!--
paper/sections/06_discussion.md: Discussion, Conclusion, Acknowledgements, Data and code availability.
Status: draft of 27 September 2026, for Ethan to edit before moving it into the shared Google Doc.
Adversarial check, 27 Sept: every number re-read from its source file; changes are listed in the open items at the end.

Conventions used in this file
- Follows "## 3. Results" (03-05), so the Discussion is "## 4" and the Conclusion "## 5". Cross-references:
  Sec. 3.2.2 (Z2 line and shared-pair crossover), 3.2.3 (even/odd rule, Table 2, chiral term), 3.2.4 (tricritical
  test, Eq. 3.3), 3.2.5 (hardware), 3.3.1 (cost and reach), 3.3.3 (boundary-forced wall); Eq. (2.16) structure
  factors, Eq. (2.21) J_m, Eq. (3.2) dressed field. Renumber on merge.
- Plain integers in brackets, e.g. [3], are citations to the list at the end of this file. Renumber on merge.
- Every quoted number has a "source" comment. "Scratch" means a file in the drafting session's scratchpad that the
  Results drafts also cite; move it into research/pilots/2026-09-27/outputs/ or validation/ before submission.
- No new calculation was run for this section. It contains no figures; all figures are in Sections 1-3.
- The Discussion is written for the plan's main direction. If gate P0, P1 or P2 fails (docs/RESEARCH_PLAN.md,
  "Fallback"), delete the second result under "What is new", the PDW sentences of the Zeybek paragraph, of the
  ladder paragraph and of outlook item (ii), and rewrite the Conclusion around the Z2 line, the negative tricritical
  test and the boundary wall.
-->

## 4. Discussion

**What is new.** Subject to the literature check (P0) and to the validation computations (P1–P3), three results go beyond the uniform chain. First, an exact rewriting of the $r^{-6}$ interaction between singly occupied strong pairs (Eq. 2.21) turns dimerization into a competition between two quantities. One is a ferromagnetic coupling $J_m$ between pseudo-spins, which record whether the left or the right atom of a pair is excited; it vanishes as $\varepsilon\to1$. The other is a transverse field induced by the laser, which does not vanish. Second, the symmetries that survive dimerization decide what happens to each crystal. The $\mathbb{Z}_2$ crystal melts through a single Ising line. The $\mathbb{Z}_4$ crystal can melt in two steps, through a pair-density wave (PDW), and the pilot spectra point that way. The $\mathbb{Z}_3$ crystal keeps its threefold degeneracy. This even/odd rule rests on constrained bases and small rings and remains preliminary. Third, the $\mathbb{Z}_2$ line is located quantitatively. At $R_b/a=1.5$ it rises from $\Delta_c=1.0785\,\Omega$ at $\varepsilon=0$ to $5.611\,\Omega$ at $\varepsilon=0.7$<!-- source: research/pilots/2026-09-27/outputs/crit_eps0.0_g11.39_d4.json (sizes [8, 10] pairs = 16/20 atoms, Dc = 1.078508) and crit_eps0.7_g11.39_d3.json (sizes [10, 12] pairs = 20/24 atoms, Dc = 5.610715); g = 11.39 --><!-- The eps = 0 point agrees with an independent open-chain DMRG estimate, 1.0794 +/- 0.0002 (validation/results/benchmark_uniform.json, vdw_check.extrapolations, p = 1.5-2.5; validation/BENCHMARK_NOTES.md Sec. 3), to 0.1%. That check is not yet quoted in Sec. 3.2.2; add it to Table 1 or its caption and then cite it here. -->. A pseudo-spin model with dressed energy denominators (Eq. 3.2) reproduces the line to within 1% at $\varepsilon=0.6$ and 0.7, but is 15% off at $\varepsilon=0.45$<!-- source: scratch eff_check.py (formulas of docs/DIMERIZED_THEORY.md section 2), rerun during the check, against crit_eps0.6_g11.39_d3.json and crit_eps0.7_g11.39_d3.json: +0.9% (3.1776 vs 3.1494) and +0.5% (5.6404 vs 5.6107); +15.0% at eps = 0.45 (2.0761 vs 1.8052, crit_eps0.45_g11.39_d3.json). Same values in Table 1 of 04_results_alternating.md. Only eps = 0.6 and 0.7 were computed above 0.45, so the claim is not stated as "eps >= 0.6" (as in 00_front.md) -->.

**Number-conserving dimerized chains.** Zeybek et al. [1] studied a dimerized chain of atoms with two Rydberg levels, coupled by van der Waals and dipolar exchange interactions. Their microwave drive is switched off after the state is prepared, so the number of excitations is conserved. At half filling they find a density wave and a bond-order phase in which each dimer shares one excitation; at quarter filling they find a bond-order density wave with a shared excitation on every other dimer<!-- source: Z. Zeybek et al., arXiv:2303.17470v2, full text (drafting scratchpad zeybek.txt), re-read in the check: Eq. (3) H_dim ("After reaching a given filling rho, the microwave laser is switched off"); Fig. 3(a): BO and DW at rho = 1/2, BODW1 at rho = 1/4 and BODW2 at rho = 1/3; BODW1 is approximately the product over i of (b+_4i + b+_4i+1)/sqrt(2) acting on the empty chain --> (see also [2]). In pseudo-spin language their half-filled density wave is our $\mathbb{Z}_2$ crystal, and their hopping within a dimer plays the role of our transverse field. Their bond-order phase is a distinct phase because the model has a U(1) symmetry (conserved excitation number) and a fixed filling. Our laser breaks U(1), and the state with one excitation shared in each pair has the same symmetries as the nearly empty disordered state at negative detuning. No symmetry separates the two, so the change between them is a crossover, not a transition (Sec. 3.2.2). The candidate PDW is different: it breaks translation by one pair, a spatial symmetry that does not depend on U(1). It would be the laser-driven analogue of their quarter-filling phase, with the sharing within a pair generated at second order in $\Omega$ rather than by hopping.

**Ladders and tricritical points.** In two-leg ladders, crystals in which excitations are shared across rungs are separated by Ising transitions from crystals that break the symmetry between the legs [3,4]; the candidate PDW resembles the first kind. Two features are specific to the dimerized chain. First, the two atoms of a pair are exchanged by a spatial inversion rather than by an internal leg swap. Near $\varepsilon=0$ this allows a Huse–Fisher chiral coupling [5] between the pair occupation and the pseudo-spin (Sec. 3.2.3); ladders also show chiral transitions [3,4]. Second, only even-period crystals carry a free pseudo-spin: in the period-3 crystal the pattern itself fixes which atom of each occupied pair is excited. The tricritical test that the research plan fixed in advance was negative. The gap ratio stays between 7.91 and 8.15, close to the Ising value 8 and far from the tricritical Ising value 8/3, for $R_b/a=1.30$–1.85 and $\varepsilon\le0.7$<!-- source: largest rings of research/pilots/2026-09-27/outputs/crit_*.json (R_b of the largest size pair: 7.915 at g = 4.83, eps = 0.45, to 7.998 at g = 11.39, eps = 0.3) and results.jsonl (N10 ratio, 20 atoms: 7.926 at g = 11.39, eps = 0.6, to 8.147 at g = 40, eps = 0.3). R_b/a = 1.30, 1.85 are g = 4.83, 40. The results.jsonl rows at g = 100 and 150 (R_b/a = 2.15, 2.31; R = 8.71 and 10.05 at 20 atoms, still drifting with size) lie outside this range and are set aside, as in Sec. 3.2.4; they are not shown to be outside the Z2 regime -->, as the classical criterion of Eq. (3.3) predicts. Tricritical Ising spectra have been measured in a Rydberg chain [6], and Rydberg ladders have been proposed for probing tricritical points [7]. Within the range we tested, dimerizing the spacing of a single-species chain does not provide a further route to such a point.

**Limitations.** All results are classical simulations of the ideal Hamiltonian, without decay, dephasing, position disorder or preparation and detection errors. Bulk quantities come from rings of at most 24 atoms ($\mathbb{Z}_2$) and 36 atoms ($\mathbb{Z}_4$), mostly in restricted bases. Excluding doubly excited pairs shifts $\Delta_c$ by 0.4% at $\varepsilon=0.3$<!-- source: scratch methods_d3_vs_d4.json: Dc = 1.325677 (d = 3) vs 1.331185 (d = 4), 16/20-atom rings, g = 11.39; difference 0.0055 = 0.41%; quoted in 02_methods.md Sec. 2.2 -->, much more than the spread between ring sizes. The $\mathbb{Z}_4$ constraint omits virtual processes that shift levels by about $0.02\,\Omega$<!-- source: 02_methods.md Sec. 2.2 estimate, (Omega/2)^2 / (V(2a) - Delta) with V(2a) = 1291.5/64 = 20.2 and Delta = 5: 0.016 (0.015 at Delta = 3) -->, comparable to the gaps that distinguish the candidate PDW from the $\mathbb{Z}_4$ crystal. Full-basis ED has been run only up to 20 atoms<!-- source: validation/results/scaling_benchmark.json, ed[N = 20]: dim = 1048576, total_s = 362 (ARPACK); ed_lobpcg[N = 20]: total_s = 55. 02_methods.md and Sec. 3.3.1 put the limit at about 22-24 atoms -->. Pair-site DMRG on 100 and 200 atoms did not converge within 30 sweeps at bond dimension 64, because of the boundary-forced domain wall (energy variance up to $4\times10^{-4}\,\Omega^2$; Sec. 3.3.1)<!-- source: validation/results/scaling_benchmark.json, dmrg_pair entries n_atoms = 100 and 200 (eps = 0.3, g = 11.39, Delta = 2.5, chi = 64, mmax = 8): variance = 1.88e-4 and 4.27e-4, converged = false, n_sweeps = 31 > max_sweeps = 30 -->. Finally, the dressed field of Eq. (3.2) is empirical.

**Outlook.** (i) The uniform $\mathbb{Z}_4$ crystal melts through chiral transitions or a floating phase [5,8–10]. Whether the transitions of the dimerized chain inherit this requires DMRG on hundreds of atoms with correlations resolved by wave vector (P7). (ii) The pair occupation and the pseudo-spin are both diagonal in the occupation basis, so projective snapshots [11] give the structure factors of Eq. (2.16) directly. A PDW would show long-range order in the pair occupation with short-range pseudo-spin correlations. With the device parameters of Sec. 3.2.5 [12], a chain at $R_b/a=1.5$ can be dimerized only up to $\varepsilon=0.28$–0.39, which covers the weakly dimerized part of the $\mathbb{Z}_2$ line. The $\mathbb{Z}_4$ regime would need short gaps of 1.3–2.1 µm for $\varepsilon=0.3$–0.5, below the 4 µm minimum spacing<!-- source: research/pilots/2026-09-27/scripts/quick_checks.py, printed output, rerun during the check (C6 = 862690 MHz um^6, Omega/2pi = 1, 1.5, 2.5 MHz, R_b/a = 1.5: eps_max = 0.39, 0.34, 0.28); Z4 gaps: 04_results_alternating.md Sec. 3.2.5 (R_b/a = 3.3: a = 2.96-2.54 um, s = 2.07-1.27 um for eps = 0.3-0.5). Device values to be checked in P9 -->. Open chains also need suitable lengths, pinned end atoms or boundary control [13] to avoid boundary-forced domain walls (Sec. 3.3.3). (iii) For dimerized two-dimensional arrays, the symmetry counting must be redone for each lattice, starting from the known density-wave phases [14,15].

## 5. Conclusion

In a laser-driven van der Waals Rydberg chain with alternating spacing, the pseudo-spin of each strong pair, coupled ferromagnetically through a second difference of $r^{-6}$ and flipped by the laser, organizes the preliminary results. These are a single Ising melting line of the $\mathbb{Z}_2$ crystal, which a dressed transverse-field Ising model reproduces to 1% at $\varepsilon=0.6$ and 0.7; a candidate pair-density wave between the $\mathbb{Z}_4$ crystal and the disordered phase; a $\mathbb{Z}_3$ degeneracy that survives dimerization; and no tricritical point for $R_b/a=1.30$–1.85 and $\varepsilon\le0.7$. The $\mathbb{Z}_2$ line is the firmest result and still needs recomputation in the full pair basis (P3). The even/odd rule becomes a finding only after validation in the full basis and on larger systems (P1, P2, P7).

## Acknowledgements

We thank [names] for discussions. [Funding: agency and grant number, or "This research received no specific funding"; to be completed by B.M.R.]

**Use of AI tools.** AI tools [name each tool and version, e.g. Anthropic's Claude] assisted in writing the code in `study/`, `rydberg_chain/`, `tests/` and `validation/`, in running the pilot calculations in `research/pilots/2026-09-27/`, and in drafting this manuscript. The authors reviewed and tested the code, checked every quoted number against the saved data, and take full responsibility for this article.
<!-- This statement is true only once Ethan has done what docs/RESEARCH_PLAN.md section 6 requires: re-derive J_m, h, the dressing, the symmetry table and the gap-ratio method, and reproduce the key numbers with the maintained package. Do not submit before then. Check the target journal's AI policy for where the disclosure belongs (Acknowledgements or Methods) and whether tools must be named. -->

## Data and code availability

The code, data and figures supporting this study are available at https://github.com/Sineceptor/rydberg-atom-research [release tag or commit, once rydberg_chain/, tests/, validation/ and paper/ are pushed; archived DOI]. `study/` holds the historical calculations, results and figures, checked against a SHA-256 manifest; `rydberg_chain/` and `tests/` hold the maintained package for exact diagonalization and DMRG; `validation/results/` holds the benchmarks; and `research/pilots/2026-09-27/` holds the exploratory scripts and raw outputs behind every result labelled preliminary. No third-party data were used.
<!-- Before submission: (1) rydberg_chain/, tests/, validation/ and paper/ are untracked in git as of 27 Sept and are not yet on GitHub; commit and push them. (2) Move the scratch outputs cited in 02-06 (d3_N32.txt, methods_d3_vs_d4.json, eff_check.py and others) into the repository. (3) Tag a release and archive it (for example on Zenodo) to obtain a DOI. (4) README.md states that no software licence is assigned; choose one, or state the terms of reuse here. (5) Confirm that the repository is public. -->

### References for this section

<!-- arXiv identifiers from 2024 onward were confirmed against the local arXiv abstract archive. Author lists of [2]-[4], [6], [7] and [13] and the journal reference of [7] were confirmed during the check from the arXiv listing metadata (dc:creator, journal_reference) in the drafting scratchpad (rss/rss_data/*/*.xml). The authors and journal of [1] were read from the arXiv v2 full text, which does not give the article number. Entries marked [verify] need a bibliographic check. Cross-references: [1] = Intro [30]; [2] = Intro [31]; [3] = Intro [34]; [4] = Intro [35]; [6] = Intro [36]; [7] = Intro [37]; [8] = Intro [11]; [9] = Intro [14]; [10] = Intro [12]; [11] = Intro [1]; [12] = Sec. 3.2 [13]; [14] = Intro [17]. -->

[1] Z. Zeybek, R. Mukherjee, and P. Schmelcher, "Quantum phases from competing van der Waals and dipole-dipole interactions of Rydberg atoms," *Phys. Rev. Lett.* **131**, 203003 (2023), arXiv:2303.17470. [verify volume and article number]

[2] Z. Zeybek, P. Schmelcher, and R. Mukherjee, "Bond-order density wave phases in dimerized extended Bose–Hubbard models," arXiv:2403.06649 (2024). [verify journal version]

[3] J. Soto Garcia and N. Chepiga, "Numerical investigation of quantum phases and phase transitions in a two-leg ladder of Rydberg atoms," arXiv:2411.05494 (2024). [verify author name form and journal version: arXiv v2 (Feb 2025) lists the first author as "Jose Soto"]

[4] S.-A. Liao, J. Zhang, and L.-P. Yang, "Phase diagram of Rydberg atoms in a two-leg rectangular ladder," arXiv:2412.11201 (2024). [verify journal version]

[5] D. A. Huse and M. E. Fisher, "Commensurate melting, domain walls, and dislocations," *Phys. Rev. B* **29**, 239–270 (1984).

[6] X. Sun, Y. Le, S. Naus, R. B.-S. Tsai, L. R. B. Picard, S. Murciano, M. Knap, J. Alicea, and M. Endres, "Experimental observation of conformal field theory spectra," arXiv:2601.16275 (2026). [verify journal version; project notes list *Nature* (2026)]

[7] H. Wang, X. Li, and C. Li, "Tricritical Kibble–Zurek scaling in Rydberg atom ladders," *Nat. Commun.* **16**, 10584 (2025), arXiv:2505.12979.

[8] P. Fendley, K. Sengupta, and S. Sachdev, "Competing density-wave orders in a one-dimensional hard-boson model," *Phys. Rev. B* **69**, 075106 (2004), arXiv:cond-mat/0309438.

[9] N. Chepiga and F. Mila, "Kibble–Zurek exponent and chiral transition of the period-4 phase of Rydberg chains," *Nat. Commun.* **12**, 414 (2021). [verify]

[10] M. Rader and A. M. Läuchli, "Floating phases in one-dimensional Rydberg Ising chains," arXiv:1908.02068 (2019).

[11] H. Bernien, S. Schwartz, A. Keesling, H. Levine, A. Omran, H. Pichler, S. Choi, A. S. Zibrov, M. Endres, M. Greiner, V. Vuletić, and M. D. Lukin, "Probing many-body dynamics on a 51-atom quantum simulator," *Nature* **551**, 579–584 (2017), arXiv:1707.04344.

[12] J. Wurtz, A. Bylinskii, B. Braverman, J. Amato-Grill, S. H. Cantu, F. Huber, A. Lukin, F. Liu, P. Weinberg, J. Long, S.-T. Wang, N. Gemelke, and A. Keesling, "Aquila: QuEra's 256-qubit neutral-atom quantum computer," arXiv:2306.11727 (2023). [verify authors and the device parameters used in Sec. 3.2.5]

[13] Y. M. Lokare and M. J. Coley-O'Rourke, "Stabilization of bulk quantum orders in finite Rydberg atom arrays," arXiv:2604.18890 (2026). [verify journal version]

[14] S. Ebadi, T. T. Wang, H. Levine, A. Keesling, G. Semeghini, A. Omran, D. Bluvstein, R. Samajdar, H. Pichler, W. W. Ho, S. Choi, S. Sachdev, M. Greiner, V. Vuletić, and M. D. Lukin, "Quantum phases of matter on a 256-atom programmable quantum simulator," *Nature* **595**, 227–232 (2021), arXiv:2012.12281.

[15] R. Samajdar, W. W. Ho, H. Pichler, M. D. Lukin, and S. Sachdev, "Complex density wave orders and quantum phase transitions in a model of square-lattice Rydberg atom arrays," *Phys. Rev. Lett.* **124**, 103601 (2020), arXiv:1910.09548. [verify arXiv id; authors, journal, volume and article number match reference [42] of the Zeybek et al. full text]

<!--
Open items after the adversarial check of 27 September (for Ethan):
1. The chiral-term sentence of "Ladders and tricritical points" now follows Sec. 3.2.3: it no longer implies that
   ladders forbid chiral terms (the ladder abstracts report Z3- and Z4-chiral transitions). docs/DIMERIZED_THEORY.md
   section 4 ("Unlike rectangular ladders ...") still uses the older contrast; align it.
2. The eps = 0 cross-check against open-chain DMRG (1.0794 vs 1.0785) is in validation/ but not in Sec. 3.2.2; add it
   there before citing it here.
3. Scratch files cited here (eff_check.py, methods_d3_vs_d4.json) must be moved into the repository.
4. The authors of [13] (Y. M. Lokare, M. J. Coley-O'Rourke) and of [2]-[4], [6], [7] come from the arXiv listing
   metadata; check for published versions of [2]-[4], [6] and [13].
-->
