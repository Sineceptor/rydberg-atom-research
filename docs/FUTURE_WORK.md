# Proposed next experiment

**How do alternating neighbour interactions change the excitation profile of a finite Rydberg chain, and can the mean density conceal errors at individual sites when distant interactions are omitted?**

This is a proposed calculation, not a completed result or established novelty claim. It extends the same driven two-level model and keeps every pair interaction consistent with the atom coordinates.

Let a be a reference spacing and g=C₆/(ℏΩa⁶). Choose positions through

$$x_0=0,\qquad x_{i+1}-x_i=a[1+(-1)^i\alpha]^{-1/6},\qquad |\alpha|<1.$$

The two neighbouring couplings become g(1+α) and g(1−α). All other couplings follow their actual distances to the sixth power. This fixes the mean of the two bond types, but not the total chain length or mean over every bond of a finite even chain.

Start with N=8 and 9, g=8, δ=0,2,3 and α=−0.5,−0.25,0,0.25,0.5. The 30 all-pair states describe the physical change. The corresponding nearest-only and second-neighbour models add 60 comparison states. Save coordinates, pair strengths, site occupations, mean density, gap and adjacent connected correlations. Plot physics versus α before approximation errors.

Required checks include recovery of the uniform model at α=0, independent small-matrix agreement, normalisation, low-eigenpair residuals and stability under tighter solver settings. At odd N the ±α chains are related by reflection: checked low-energy eigenvalues should agree and site profiles reverse. At even N the two terminations are different finite chains. Their different counts of strong and weak bonds and total lengths must be reported.

The second-neighbour distance is a[(1+α)⁻¹ᐟ⁶+(1−α)⁻¹ᐟ⁶], which increases with |α|. That directly weakens those interactions. Any improvement in a truncated approximation needs to be checked against this geometric explanation before invoking unexpected many-body behaviour.

Refine only informative regions, then test selected settings at nearby odd and even lengths. State numerical acceptance thresholds and observable stability before interpreting small differences. A null result can establish robustness over a tested domain; it should be retained rather than hidden.

Alternating density interactions do not automatically realise the SSH model or a topological phase. The [Rydberg SSH experiment](https://arxiv.org/abs/1810.13286) used dipolar exchange, a different Hamiltonian term. [Dimerised Rydberg studies with both density and exchange interactions](https://arxiv.org/abs/2303.17470) also require careful model matching.

A later optional sensitivity study could change detuning only at the two ends, keeping the geometry fixed. This would probe boundary response, not establish bulk behaviour. [Boundary-engineering research](https://arxiv.org/abs/2604.18890) provides context but uses a substantially different scale and protocol.
