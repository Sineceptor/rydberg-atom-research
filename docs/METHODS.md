# Methods and code map

These notes describe the implemented finite-chain model, its observables and numerical checks. Qiskit and SciPy have distinct roles.

## From an atom arrangement to a result

1. **Choose the model.** Coordinates determine pair distances. Laser parameters determine the drive Ω and detuning Δ. A basis label such as 0101 specifies which atoms are in the Rydberg state.
2. **Build interactions.** `rydberg_model.interactions` constructs the pair matrix using distance to the power −6. Removing specified couplings is an ideal model control. It is not a simulation of physically losing an atom unless that model is separately defined.
3. **Build the Hamiltonian.** `rydberg_model.hamiltonian` constructs a sparse matrix in the bit basis. Its diagonal entries are −Δ times the number of excited atoms plus the energies of occupied pairs. Flipping a single bit gives a laser matrix element Ω/2 in the implementation's ℏ=1 units.
4. **Solve the intended problem.** Ground-state calculations find the smallest energies and associated vectors. Dynamics apply the matrix exponential to an initial state. A low-energy eigenstate and the state reached by a sudden laser quench are different objects.
5. **Measure the state mathematically.** Squared amplitudes weight each configuration. Sum these probabilities to obtain site occupations and pair probabilities. Reduced density matrices give entropy for a specified division of atoms.
6. **Save a reproducible result.** CSV/JSON files record numerical outputs; figures display them. Managed runs include code hashes and a completion status. Captions and interpretation must still be written and checked by a person.

The bit-basis convention in the reusable model is little-endian: site i is bit i of the integer basis index. A printed binary string can display the opposite order from a left-to-right drawing; the methods and figure labels must state the convention rather than infer site order from the string.

## Equations to explain in the paper

With nᵢ=|1ᵢ⟩⟨1ᵢ|=(I−Zᵢ)/2 and |0⟩ the ground state,

$$H=\frac{\hbar\Omega}{2}\sum_iX_i-\hbar\Delta\sum_i n_i+\sum_{i<j}\frac{C_6}{r_{ij}^6}n_i n_j.$$

This equation uses C₆ in energy×length⁶ units. The code stores the corresponding frequencies with ℏ=1, and usually Ω=1. Report δ=Δ/Ω, g=C₆/(ℏΩa⁶), and τ=Ωt explicitly. The implemented pair sum is i<j and counts each pair once.

For stationary states and the lowest gap,

$$H\phi_k=E_k\phi_k,\qquad G=E_1-E_0.$$

For constant-H dynamics from an explicitly stated initial state,

$$|\psi(t)\rangle=e^{-iHt/\hbar}|\psi(0)\rangle.$$

For site i and average density,

$$\langle n_i\rangle=\sum_s|a_s|^2 n_i(s),\qquad \bar n=\frac1N\sum_i\langle n_i\rangle.$$

For correlations and the corrected staggered-spin diagnostic,

$$C_{ij}=\langle n_i n_j\rangle-\langle n_i\rangle\langle n_j\rangle,\qquad Q=\frac1N\sum_i(-1)^i(2n_i-1),\qquad\langle Q^2\rangle.$$

Q is the quantity called S in parts of the parity package; a separate symbol here avoids confusing it with entropy. Either perfectly alternating binary string has Q²=1, including odd N. A finite independent state can also have a nonzero expectation, so this number alone does not identify a phase.

If entropy is used, define the retained atoms A and use ρ_A=Tr_B|ψ⟩⟨ψ| and S_A=−Tr(ρ_A log₂ρ_A). State that these examples use pure global states, for which this reduced-state entropy measures bipartite entanglement. A Bell-state control is one bit; a product-state control is zero.

## Which file does which job

All names below refer to `study/`.

| File | Role |
|---|---|
| `rydberg_model.py` | Coordinates, distance-derived pair strengths, sparse matrix, dynamics, observables and Qiskit equivalents |
| `extension_model.py` | Finite-chain ground-state helpers and diagnostics used by the interaction-cutoff study |
| `run_august_experiments.py` | Qiskit/SciPy controls; blockade, defects, geometry, entropy, scaling and circuit comparisons |
| `phase_diagram_bridge.py` | Finite-chain ground-state map and size/tail examples |
| `run_extension.py` | Uniform-chain all-pair, first-neighbour and second-neighbour comparisons, with bundled-reference check |
| `parity_diagnostics.py` | Odd/even profiles, corrected order, several low states and configuration probabilities |
| `circuit_controls.py` | Independent tensor/SciPy and split-circuit tests at both signs of detuning |
| `kagome_pilot.py`, `kagome_shape_control.py` | Optional finite patches, clearly separated from the chain question |
| `run_project.py`, `run_support.py` | Reproduction commands and isolated completed/failed run directories |
| `verify_delivery.py` | Checksums of the delivered snapshot, without recomputing physics |

The current `extension_model.py` chain helper is designed for uniform chains. A proposed alternating study must feed a new, tested coordinate/coupling matrix into the model; changing only a plot label or the uniform distance formula is insufficient. Keep that new implementation in a new experiment rather than silently replacing this baseline.

## What the numerical methods mean

For a small control, dense Hermitian diagonalisation returns all eigenpairs. For the larger small-chain scans, the code uses a sparse Hermitian eigensolver to request a few lowest eigenpairs. The 2ᴺ-dimensional state space still grows exponentially; sparse storage does not make arbitrary sizes easy.

Time evolution uses SciPy's action of a matrix exponential. The Qiskit representation independently constructs the same operator from Pauli terms, and small-state evolution provides a comparison. The split-step circuit approximates the evolution through gates and needs a convergence check as the step size decreases. These are classical calculations of quantum models; they did not execute on a physical quantum processor.

The report should include solver settings actually used by the selected run, not generic settings copied from a tutorial. Record atom count, coordinates, boundary convention, drive/detuning, cutoffs, initial state or eigenstate, time/grid choices, tolerance, library versions and the figure's source file. Cite the [SciPy eigsh documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html), [SciPy expm_multiply documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html), and [Qiskit quantum information documentation](https://quantum.cloud.ibm.com/docs/en/api/qiskit/quantum_info) for the methods actually used.

## Validation and limits

Use the hand-built two-atom matrix, g=0 analytic limits, independent matrix/evolution comparisons, probability and energy conservation, symmetry and solver residuals. At close gaps, vary the solver starting vector and inspect more than one low level. The saved JSON/CSV result files contain exact test settings and recorded errors. See [Validation](VALIDATION.md) and [Results and limitations](RESULTS_AND_LIMITS.md).

Passing numerical checks supports implementation accuracy. It does not remove open boundaries, establish a thermodynamic phase, verify an experimental noise model or prove novelty. The methods section should keep these questions separate.
