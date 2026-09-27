## 1. Introduction

<!-- Draft for review by E. Zhang before transfer to the shared Google Doc. Citations are numbered in order of first appearance within this section; renumber when the sections are merged. Figure numbering also needs harmonizing with the other sections. -->

Rydberg atom arrays are a programmable quantum simulation platform in which individual atoms are trapped in optical tweezers, arranged atom by atom into chosen one- and two-dimensional geometries, and laser-excited to Rydberg states of high principal quantum number [1–3]. Because Rydberg atoms have large electronic orbitals, they interact strongly over micrometer distances through long-range van der Waals interactions, whose coefficient $C_6$ grows approximately as the eleventh power of the principal quantum number [4,5]. By controlling laser parameters and atom geometry, researchers can simulate interacting quantum many-body systems and investigate quantum phase transitions, ordered phases, frustrated magnetism and topological phases of matter. Because the interaction geometry and Hamiltonian parameters are tuned experimentally, theoretical models can be compared directly with quantum simulations [1].

The standard model for a Rydberg atom array is

$$
H=\frac{\hbar\Omega}{2}\sum_i\sigma_i^x-\hbar\Delta\sum_i n_i+\sum_{i<j}\frac{C_6}{r_{ij}^6}\,n_in_j, \qquad (1)
$$

where $\sigma_i^x=|g_i\rangle\langle r_i|+|r_i\rangle\langle g_i|$ drives transitions between the ground state $|g\rangle$ and the Rydberg state $|r\rangle$ of atom $i$ at the Rabi frequency $\Omega$, $n_i=|r_i\rangle\langle r_i|$ counts Rydberg excitations, $\Delta$ is the laser detuning and $r_{ij}$ is the distance between atoms $i$ and $j$ [1]. The detuning acts as a chemical potential: for $\Delta>0$ each excitation lowers the energy by $\hbar\Delta$, while the last term penalizes excitations that are close together. When $C_6/r^6\gg\hbar\Omega$, the doubly excited state of two atoms is shifted out of resonance and only one of them can be excited. This Rydberg blockade between two individually trapped atoms was first observed in 2009 [6,7]. The blockade radius $R_b$, defined by $C_6/R_b^6=\hbar\Omega$, sets the distance within which simultaneous excitation is suppressed [1]. For a chain with spacing $a$, the physics is controlled by $\Delta/\Omega$ and the dimensionless coupling constant $g=C_6/(\hbar\Omega a^6)=(R_b/a)^6$. Mapping the resulting phase diagrams requires solving a quantum many-body problem, usually by exact diagonalization (ED) for small systems or by tensor-network methods such as the density-matrix renormalization group (DMRG) for larger ones [8,9].

**Phases in one and two dimensions.** In a uniform chain at positive detunings of a few $\Omega$, too small to overcome the blockade, the excitations crystallize: every second atom is excited in the $\mathbb{Z}_2$ phase, every third in the $\mathbb{Z}_3$ phase and every fourth in the $\mathbb{Z}_4$ phase, for $R_b/a$ roughly between 1 and 2, 2 and 3, and 3 and 4, respectively [1,10]. Bernien et al. prepared these crystals in chains of up to 51 atoms by sweeping the detuning [1], and Keesling et al. used the quantum Kibble–Zurek mechanism to probe the critical dynamics of the transitions into them [10]. The transitions differ in kind [11–15]. Melting of the $\mathbb{Z}_2$ crystal is in the Ising universality class, whose conformal-field-theory fields have been identified with microscopic operators of Rydberg chains [15]. Three-state Potts criticality occurs only at one commensurate point of the $\mathbb{Z}_3$ boundary; away from such special points the $\mathbb{Z}_3$ and $\mathbb{Z}_4$ crystals melt through chiral transitions of Huse–Fisher type or through an intermediate incommensurate floating phase [11–14], which has since been probed experimentally [16]. In two dimensions, Ebadi et al. observed checkerboard, striated and star-ordered phases on a 256-atom square array [17], and Scholl et al. studied antiferromagnetic ordering on square and triangular arrays of up to hundreds of atoms [18]. Verresen, Lukin and Vishwanath predicted that the blockade constraint on the ruby lattice stabilizes a $\mathbb{Z}_2$ topological spin liquid of the toric-code type [19], and Semeghini et al. detected signatures of this state using topological string operators [20].

**Rydberg arrays for quantum computing.** Besides analog Hamiltonian simulation [17], the blockade makes Rydberg arrays a platform for digital quantum computing. Jaksch et al. proposed fast controlled-phase gates between neutral atoms based on Rydberg interactions [21], and Lukin et al. showed how dipole blockade enables quantum information processing in atomic ensembles [22]. The Lukin group and collaborators then implemented parallel multiqubit gates, including a three-qubit Toffoli gate [23], reached two-qubit gate fidelities of 99.5% in parallel across 60 atoms [24], and used coherent transport of entangled atoms to obtain reconfigurable, non-local connectivity [25]. This architecture underlies a logical quantum processor with up to 48 logical qubits [26]. Because the blockade constraint is the independent-set condition on a unit-disk graph, the arrays also encode combinatorial optimization directly, as in the study of the maximum independent set problem on graphs of up to 289 atoms [27].

**Why alternate the spacing?** Tweezers can place atoms with alternating short and long gaps, down to a minimum spacing set by the hardware. Bond alternation is the defining ingredient of the Su–Schrieffer–Heeger (SSH) model [28], which de Léséleuc et al. realized for interacting bosons using dipolar exchange between two Rydberg states, observing a symmetry-protected topological phase [29]. Zeybek, Mukherjee and Schmelcher found a bond-order density-wave phase in a number-conserving dimerized chain with van der Waals and dipolar exchange interactions [30]; related phases occur in dimerized extended Bose–Hubbard models [31]. Other structured chains alternate the atomic species [32], and multicomponent arrays have been proposed as a way to tune quantum criticality [33].
<!-- [verify] the one-line description of Ref. [33] (Chepiga, multicomponent arrays) against the full text in P0; it is paraphrased from the title only. -->
The closest geometric relatives are two-leg ladders, where crystals with one excitation shared across a rung are separated by an Ising transition from crystals that also break the symmetry between the legs [34,35]. Tricritical Ising points, where a continuous Ising line turns first order, appear in theoretical phase diagrams of constrained chains [11]; Ising and tricritical Ising finite-size spectra have been measured in a Rydberg chain [36], and tricritical Kibble–Zurek scaling has been proposed for Rydberg ladders [37]. The same alternating geometry has been proposed for constrained four-body dynamics, with blockade across short bonds and antiblockade across long ones [38].
<!-- Ref. [15] (Slagle et al.) was removed from the tricritical citation because it could not be confirmed that it treats the tricritical point; Fendley, Sengupta and Sachdev [11] do. Ref. [38] (arXiv:2608.20233) uses the same geometry for antiblockade dynamics, not for the ground-state phase diagram; confirmed from the abstract in the local arXiv archive. -->

To our knowledge, the ground-state phase diagram of the simplest member of this family has not been mapped: a chain with a single Rydberg level, laser driving and van der Waals interactions only, in which the gaps alternate between $s=a(1-\varepsilon)$ and $w=a(1+\varepsilon)$ with $0\le\varepsilon<1$ (Figure 1).
<!-- Novelty statement is contingent on literature gate P0 (docs/RESEARCH_PLAN.md, section 3). Chepiga PRL 132, 076505 and arXiv:2408.15965, 2412.11201, 2411.05494, 2403.06649, 2406.20093, 2601.16275 and 2506.21963 must be read in full before submission. Closest known overlap: species alternation in 2408.15965 and the Z_p^+ rung-shared phases of 2412.11201; 2608.20233 uses the same spacing-dimerized geometry, but for antiblockade dynamics. A keyword search of the local arXiv archive (2024-03 to 2026-09) on 27 Sept found no ground-state study of the laser-driven, single-level vdW dimerized chain. Earlier searches could not read arXiv full text. -->

[FIGURE PENDING: schematic of the dimerized chain and pair pseudo-spin. Panel (a): open chain with strong gaps s = a(1 − ε) and weak gaps w = a(1 + ε), drawn to scale for ε = 0.3, with the blockade radius for R_b/a = 1.5 (g = 11.39). Panel (b): the two states of an occupied pair, τ = +1 and τ = −1, and their exchange by inversion I_s; ideal Z2, Z3 and Z4 patterns of the uniform chain for comparison. Data source: rydberg_chain/geometry.py (chain_positions, strong_weak_couplings). Computation ID: none (schematic); any J_m values printed on the figure come from P4.]

*Figure 1. Dimerized Rydberg chain and pair pseudo-spin. (a) Atoms $2k$ and $2k+1$ form strong pair $k$ with gap $s=a(1-\varepsilon)$, and neighbouring pairs are separated by $w=a(1+\varepsilon)$, so the mean spacing stays $a$ and $g=(R_b/a)^6$ keeps its uniform-chain meaning; the geometry is drawn for $\varepsilon=0.3$ with the blockade radius for $R_b/a=1.5$ ($g=11.39$). (b) An occupied pair carries $\tau=+1$ (left atom excited) or $\tau=-1$ (right atom excited), and inversion $I_s$ about the strong-bond centre exchanges the two. The ideal $\mathbb{Z}_2$, $\mathbb{Z}_3$ and $\mathbb{Z}_4$ patterns of the uniform chain ($\varepsilon=0$) are shown for comparison.*

An excitation on a strong pair sits on either the left or the right atom, defining a pseudo-spin $\tau=n_L-n_R$ that is reversed by inversion about the strong-bond centre. Occupied pairs couple ferromagnetically in $\tau$ through a second difference of $1/r^6$, which vanishes as $\varepsilon\to1$ while the laser-induced transverse field stays finite. This picture predicts a single Ising melting line for the $\mathbb{Z}_2$ crystal. In the $\mathbb{Z}_4$ regime, dimerization lowers the symmetry permuting the four crystal patterns from $D_4$ to $\mathbb{Z}_2\times\mathbb{Z}_2$, which allows the melting to split into two transitions with an intermediate pair-density wave, in which every other strong pair holds one excitation shared between its two atoms. The $\mathbb{Z}_3$ regime keeps its $D_3$ symmetry and should gain no new transition. We also test, as a negative control stated in advance, whether dimerization alone can produce a tricritical point. A classical-limit argument suggests not: tricriticality needs a net attraction between neighbouring occupied pairs with equal $\tau$, but such pairs are always $2a$ apart and repel.

Preliminary pilot calculations are consistent with these predictions. Symmetry-resolved ED of periodic rings of up to 24 atoms at $g=11.39$ ($R_b/a=1.5$) locates the $\mathbb{Z}_2$ line at five dimerizations from $\varepsilon=0$ to 0.7. There the ratio of the lowest inversion-even to the lowest inversion-odd excitation energy lies between 7.93 and 8.00, close to the Ising value 8 and far from the tricritical Ising value 8/3.
<!-- source: research/pilots/2026-09-27/outputs/crit_eps0.0_g11.39_d4.json, crit_eps0.3_g11.39_d3.json, crit_eps0.45_g11.39_d3.json, crit_eps0.6_g11.39_d3.json, crit_eps0.7_g11.39_d3.json. Range quoted for the largest size pair at each eps (fields R_a, R_b of the last "pairs" entry): 7.955/7.958 (eps 0, 16/20 atoms), 8.003/7.998 (0.3), 7.992/7.990 (0.45), 7.939/7.948 (0.6), 7.927/7.938 (0.7; 20/24 atoms) -> 7.927-8.003, the same range as Sec. 3.2.2. Including the smaller size pairs the range is 7.894-8.007. "sizes" in these files count pairs (atoms = 2 x pairs). R = (second k = 0 inversion-even level - E0)/(lowest inversion-odd level - E0), scripts research/pilots/2026-09-27/scripts/pairk0.py and crit.py. Ring sizes: 12-20 atoms at eps = 0 (full four-state pair basis, d = 4); 16-24 atoms at eps >= 0.3 with doubly excited strong pairs excluded (d = 3). Exploratory; to be recomputed with rydberg_chain in P3. -->
A 24-atom ring at $R_b/a=2.7$ and $\Delta/\Omega=5$ keeps the threefold $\mathbb{Z}_3$ ground-state degeneracy to within $10^{-5}\,\hbar\Omega$ at $\varepsilon=0$, 0.3 and 0.6.
<!-- source: research/pilots/2026-09-27/outputs/z3_24_0.0.txt, z3_24_0.3.txt, z3_24_0.6.txt (g = 387.4, i.e. R_b/a = 2.70; three lowest dE printed as 0.00000; basis excludes simultaneous excitation of adjacent atoms only, dim = 103682 = number of 24-site rings with no two adjacent excitations, checked 27 Sept with ced.basis(24, excl=1)), script scripts/ced.py. To be extended to N = 30, 36 in P5. -->
At $R_b/a=3.3$, $\varepsilon=0.5$ and $\Delta/\Omega=2$, the two lowest states of rings of 24 and 36 atoms are both even under inversion and have opposite pair-translation eigenvalues $T_2=\pm1$. Their splitting falls from $0.026$ to $0.0099\,\hbar\Omega$ as the ring grows, while the lowest inversion-odd state stays about $0.25\,\hbar\Omega$ higher: the finite-size signature of a state that breaks pair translation but not inversion, as expected for the pair-density wave.
<!-- source: research/pilots/2026-09-27/outputs/out24_2.0.txt and out36_2.0.txt (g = 1291.5, i.e. R_b/a = 3.30; eps = 0.50; ground state T2 = +1, Is = +1; second state T2 = -1, Is = +1 at dE = 0.02639 (N = 24) and 0.00993 (N = 36); lowest Is = -1 state at 0.25356 and 0.24717), script scripts/ced.py, basis with no two excitations within two sites (dims 9642 and 946781). Same interpretation as Sec. 3.2.3. The earlier draft quoted Delta = 5 (out24_5.0.txt, out36_5.0.txt) as the PDW signature; there the inversion-odd gap is only 0.01743 -> 0.01248 (N x gap 0.42 -> 0.45), i.e. near the tau transition (Sec. 3.2.3), so it is not used. The Delta = 3 gap sequence in docs/RESEARCH_PLAN.md is not saved in outputs/ and is not used here. -->
All of these runs except the uniform ($\varepsilon=0$) ring used restricted bases, and all will be recomputed with the maintained code (computations P3 and P5). The pair-density wave in particular must first be validated in the full Hilbert space and on larger rings (P1 and P2) before it can be claimed.

The paper is organized as follows. The project-direction paragraph that follows states the questions and how they will be decided. The Methods describe exact diagonalization, an exact regrouping of each strong pair into one four-state site, pair-site DMRG, and Qiskit checks on small chains; all calculations are classical simulations. The Results follow the outline: how the physics varies with the coupling constant $R_b/a$, with alternating couplings, and with system size, including the domain wall that boundary fields force into even open chains. The full basis already holds $2^{20}\approx1.05\times10^6$ states at $N=20$, where pair-site DMRG reproduces the ED ground-state energy of an open chain to $4\times10^{-13}\,\hbar\Omega$.
<!-- source: validation/results/scaling_benchmark.json, dmrg_pair[0].E0_minus_ED = -4.16e-13 (open chain, N = 20, eps = 0.3, g = 11.39, Delta = 2.5, chi = 64, no interaction cutoff) and ed[2].dim = 1048576. Rechecked 27 Sept: the file now also holds DMRG runs at N = 50, 100 and 200, all flagged converged = False (neither the energy nor the entropy criterion met; energy variance 3.7e-9, 1.9e-4, 4.3e-4), so no large-N accuracy is quoted here. -->
We model an ideal, fully coherent array without decoherence, state-preparation or detection errors, and discuss at the end of the paper which predictions current hardware can test.

### References for this section

<!-- arXiv identifiers are listed where they were confirmed from repository files or the local arXiv archive; the remaining ones will be added at the reference check. Entries marked [verify] need a bibliographic check before submission. -->

[1] H. Bernien, S. Schwartz, A. Keesling, H. Levine, A. Omran, H. Pichler, S. Choi, A. S. Zibrov, M. Endres, M. Greiner, V. Vuletić, and M. D. Lukin, "Probing many-body dynamics on a 51-atom quantum simulator," Nature **551**, 579 (2017), arXiv:1707.04344.

[2] M. Endres et al., "Atom-by-atom assembly of defect-free one-dimensional cold atom arrays," Science **354**, 1024 (2016).

[3] D. Barredo, S. de Léséleuc, V. Lienhard, T. Lahaye, and A. Browaeys, "An atom-by-atom assembler of defect-free arbitrary two-dimensional atomic arrays," Science **354**, 1021 (2016).

[4] A. Browaeys and T. Lahaye, "Many-body physics with individually controlled Rydberg atoms," Nat. Phys. **16**, 132 (2020).

[5] M. Saffman, T. G. Walker, and K. Mølmer, "Quantum information with Rydberg atoms," Rev. Mod. Phys. **82**, 2313 (2010).

[6] E. Urban, T. A. Johnson, T. Henage, L. Isenhower, D. D. Yavuz, T. G. Walker, and M. Saffman, "Observation of Rydberg blockade between two atoms," Nat. Phys. **5**, 110 (2009).

[7] A. Gaëtan, Y. Miroshnychenko, T. Wilk, A. Chotia, M. Viteau, D. Comparat, P. Pillet, A. Browaeys, and P. Grangier, "Observation of collective excitation of two individual atoms in the Rydberg blockade regime," Nat. Phys. **5**, 115 (2009).

[8] S. R. White, "Density matrix formulation for quantum renormalization groups," Phys. Rev. Lett. **69**, 2863 (1992).

[9] U. Schollwöck, "The density-matrix renormalization group in the age of matrix product states," Ann. Phys. (N.Y.) **326**, 96 (2011).

[10] A. Keesling et al., "Quantum Kibble–Zurek mechanism and critical dynamics on a programmable Rydberg simulator," Nature **568**, 207 (2019).

[11] P. Fendley, K. Sengupta, and S. Sachdev, "Competing density-wave orders in a one-dimensional hard-boson model," Phys. Rev. B **69**, 075106 (2004).

[12] M. Rader and A. M. Läuchli, "Floating phases in one-dimensional Rydberg Ising chains," arXiv:1908.02068 (2019).

[13] N. Chepiga and F. Mila, "Floating phase versus chiral transition in a 1D hard-boson model," Phys. Rev. Lett. **122**, 017205 (2019).

[14] N. Chepiga and F. Mila, "Kibble–Zurek exponent and chiral transition of the period-4 phase of Rydberg chains," Nat. Commun. **12**, 414 (2021). [verify]

[15] K. Slagle, D. Aasen, H. Pichler, R. S. K. Mong, P. Fendley, X. Chen, M. Endres, and J. Alicea, "Microscopic characterization of Ising conformal field theory in Rydberg chains," Phys. Rev. B **104**, 235109 (2021), arXiv:2108.09309.

[16] J. Zhang, S. H. Cantú, F. Liu, A. Bylinskii, B. Braverman, F. Huber, J. Amato-Grill, A. Lukin, N. Gemelke, A. Keesling, S.-T. Wang, Y. Meurice, and S.-W. Tsai, "Probing quantum floating phases in Rydberg atom arrays," arXiv:2401.08087 (2024). [verify journal version] <!-- author list confirmed from a local arXiv listing (Jan 2024) -->

[17] S. Ebadi, T. T. Wang, H. Levine, A. Keesling, G. Semeghini, A. Omran, D. Bluvstein, R. Samajdar, H. Pichler, W. W. Ho, S. Choi, S. Sachdev, M. Greiner, V. Vuletić, and M. D. Lukin, "Quantum phases of matter on a 256-atom programmable quantum simulator," Nature **595**, 227 (2021), arXiv:2012.12281.

[18] P. Scholl, M. Schuler, H. J. Williams, A. A. Eberharter, D. Barredo, K.-N. Schymik, V. Lienhard, L.-P. Henry, T. C. Lang, T. Lahaye, A. M. Läuchli, and A. Browaeys, "Quantum simulation of 2D antiferromagnets with hundreds of Rydberg atoms," Nature **595**, 233 (2021), arXiv:2012.12268.

[19] R. Verresen, M. D. Lukin, and A. Vishwanath, "Prediction of toric code topological order from Rydberg blockade," Phys. Rev. X **11**, 031005 (2021), arXiv:2011.12310.

[20] G. Semeghini, H. Levine, A. Keesling, S. Ebadi, T. T. Wang, D. Bluvstein, R. Verresen, H. Pichler, M. Kalinowski, R. Samajdar, A. Omran, S. Sachdev, A. Vishwanath, M. Greiner, V. Vuletić, and M. D. Lukin, "Probing topological spin liquids on a programmable quantum simulator," Science **374**, 1242 (2021), arXiv:2104.04119.

[21] D. Jaksch, J. I. Cirac, P. Zoller, S. L. Rolston, R. Côté, and M. D. Lukin, "Fast quantum gates for neutral atoms," Phys. Rev. Lett. **85**, 2208 (2000).

[22] M. D. Lukin, M. Fleischhauer, R. Côté, L. M. Duan, D. Jaksch, J. I. Cirac, and P. Zoller, "Dipole blockade and quantum information processing in mesoscopic atomic ensembles," Phys. Rev. Lett. **87**, 037901 (2001).

[23] H. Levine et al., "Parallel implementation of high-fidelity multiqubit gates with neutral atoms," Phys. Rev. Lett. **123**, 170503 (2019).

[24] S. J. Evered et al., "High-fidelity parallel entangling gates on a neutral-atom quantum computer," Nature **622**, 268 (2023).

[25] D. Bluvstein et al., "A quantum processor based on coherent transport of entangled atom arrays," Nature **604**, 451 (2022).

[26] D. Bluvstein et al., "Logical quantum processor based on reconfigurable atom arrays," Nature **626**, 58 (2024).

[27] S. Ebadi et al., "Quantum optimization of maximum independent set using Rydberg atom arrays," Science **376**, 1209 (2022), arXiv:2202.09372.

[28] W. P. Su, J. R. Schrieffer, and A. J. Heeger, "Solitons in polyacetylene," Phys. Rev. Lett. **42**, 1698 (1979).

[29] S. de Léséleuc, V. Lienhard, P. Scholl, D. Barredo, S. Weber, N. Lang, H. P. Büchler, T. Lahaye, and A. Browaeys, "Observation of a symmetry-protected topological phase of interacting bosons with Rydberg atoms," Science **365**, 775 (2019).

[30] Z. Zeybek, R. Mukherjee, and P. Schmelcher, "Quantum phases from competing van der Waals and dipole-dipole interactions of Rydberg atoms," Phys. Rev. Lett. **131**, 203003 (2023), arXiv:2303.17470. [verify volume and article number]

[31] Z. Zeybek, P. Schmelcher, and R. Mukherjee, "Bond-order density wave phases in dimerized extended Bose-Hubbard models," arXiv:2403.06649 (2024).

[32] L.-Y.-N. Liu, S.-Y. Yu, S.-R. Peng, J. Sheng, S. Yi, P. Xu, S.-S. Gong, T. Shi, and J. Cui, "Distinct many-body scars and emergent quantum phases driven by competing interactions in two-species Rydberg arrays," Phys. Rev. B **114**, 115406 (2026), arXiv:2408.15965. <!-- v1 title: "Novel ground states and emergent quantum many-body scars in a two-species Rydberg atom array"; v2 title and journal reference from the local arXiv listing of 18 Aug 2026. -->

[33] N. Chepiga, "Tunable quantum criticality in multicomponent Rydberg arrays," Phys. Rev. Lett. **132**, 076505 (2024), arXiv:2308.12838. [verify journal reference; not yet read in full (P0)]

[34] J. Soto Garcia and N. Chepiga, "Numerical investigation of quantum phases and phase transitions in a two-leg ladder of Rydberg atoms," arXiv:2411.05494 (2024). [verify author name form and journal version: arXiv v2 (Feb 2025) lists the first author as "Jose Soto"]

[35] S.-A. Liao, J. Zhang, and L.-P. Yang, "Phase diagram of Rydberg atoms in a two-leg rectangular ladder," arXiv:2412.11201 (2024).

[36] X. Sun, Y. Le, S. Naus, R. B.-S. Tsai, L. R. B. Picard, S. Murciano, M. Knap, J. Alicea, and M. Endres, "Experimental observation of conformal field theory spectra," arXiv:2601.16275 (2026). [verify journal version; project notes list Nature (2026)]

[37] H. Wang, X. Li, and C. Li, "Tricritical Kibble–Zurek scaling in Rydberg atom ladders," Nat. Commun. **16**, 10584 (2025), arXiv:2505.12979. <!-- journal reference from the local arXiv listing of v2 (2 Dec 2025), DOI 10.1038/s41467-025-65652-9 -->

[38] T. Yan and W. Li, "Nonthermal dynamics of a one-dimensional Rydberg-atom chain with constraint four-body interactions," arXiv:2608.20233 (2026).
