"""Small Rydberg-array model. Atom 0 is the rightmost bit in a Qiskit label.

We use hbar = 1, |0> = ground and |1> = Rydberg.
H = (omega/2) sum X_i - delta sum n_i + sum_(i<j) V_ij n_i n_j.
Lengths are in nearest-neighbour spacings, and time is in 1/omega units.
"""

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply


def chain(number):
    return np.column_stack((np.arange(number), np.zeros(number)))


def interactions(positions, strength, isolated=None, split=None):
    """V_ij = strength/r_ij**6. All pairs are included once.

    isolated: an atom whose interactions are switched off (laser stays on).
    split: a list of atoms on one side of a completely disconnected cut.
    These disconnections are idealised controls, not missing physical atoms.
    """
    number = len(positions)
    coupling = np.zeros((number, number))
    for i in range(number):
        for j in range(i + 1, number):
            distance = np.linalg.norm(positions[i] - positions[j])
            if distance == 0:
                raise ValueError("Two atoms cannot occupy the same position.")
            value = strength / distance**6
            if i == isolated or j == isolated:
                value = 0
            if split is not None and ((i in split) != (j in split)):
                value = 0
            coupling[i, j] = coupling[j, i] = value
    return coupling


def occupations(number):
    """One row per basis state; one column per atom (0 or 1)."""
    labels = np.arange(2**number)
    return np.column_stack([(labels >> i) & 1 for i in range(number)])


def hamiltonian(coupling, omega=1.0, delta=0.0):
    """Store only nonzero matrix entries so larger examples fit in memory."""
    number = len(coupling)
    occupied = occupations(number)
    labels = np.arange(2**number)
    energy = -delta * occupied.sum(axis=1).astype(float)
    for i in range(number):
        for j in range(i + 1, number):
            energy += coupling[i, j] * occupied[:, i] * occupied[:, j]
    rows, columns = [], []
    for i in range(number):
        # XOR flips atom i, e.g. 000 becomes 001 when i = 0.
        rows.extend(labels)
        columns.extend(labels ^ (1 << i))
    values = np.full(len(rows), omega / 2)
    drive = coo_matrix((values, (rows, columns)), shape=(2**number, 2**number))
    return drive.tocsr() + diags(energy, format="csr")


def evolve(matrix, times):
    """Apply exp(-i H t) to |00...0>. This does not use a Trotter step."""
    times = np.asarray(times)
    if len(times) < 2 or not np.allclose(np.diff(times), times[1] - times[0]):
        raise ValueError("Use at least two equally spaced times.")
    initial = np.zeros(matrix.shape[0], dtype=complex)
    initial[0] = 1
    return expm_multiply(-1j * matrix, initial, start=times[0],
                         stop=times[-1], num=len(times), endpoint=True)


def excitation_probabilities(states, number):
    return np.abs(states)**2 @ occupations(number)


def pair_probability(states, number, i, j):
    occupied = occupations(number)
    return np.abs(states)**2 @ (occupied[:, i] * occupied[:, j])


def reduced_density(state, keep):
    """Group amplitudes by the kept atoms, then sum over the other atoms."""
    number = int(np.log2(len(state)))
    # reshape uses leftmost bit first; atom 0 is the rightmost bit.
    kept_axes = [number - 1 - i for i in keep]
    other_axes = [i for i in range(number) if i not in kept_axes]
    amplitudes = state.reshape([2] * number).transpose(kept_axes + other_axes)
    amplitudes = amplitudes.reshape(2**len(keep), -1)
    return amplitudes @ amplitudes.conj().T


def entropy(state, keep):
    """Entanglement between keep and its complement, for a pure full state."""
    probabilities = np.linalg.eigvalsh(reduced_density(state, keep))
    probabilities = probabilities[probabilities > 1e-14]
    return float(-np.sum(probabilities * np.log2(probabilities)))


def qiskit_hamiltonian(coupling, omega=1.0, delta=0.0):
    """Build the same H independently with Qiskit: n = (I - Z)/2."""
    from qiskit.quantum_info import SparsePauliOp
    number = len(coupling)
    terms = [("I" * number, -delta * number / 2)]

    def label(letters):
        result = ["I"] * number
        for atom, letter in letters.items():
            result[number - 1 - atom] = letter
        return "".join(result)

    for i in range(number):
        terms.append((label({i: "X"}), omega / 2))
        terms.append((label({i: "Z"}), delta / 2))
        for j in range(i + 1, number):
            value = coupling[i, j] / 4
            terms.extend([("I" * number, value),
                          (label({i: "Z"}), -value),
                          (label({j: "Z"}), -value),
                          (label({i: "Z", j: "Z"}), value)])
    return SparsePauliOp.from_list(terms).simplify()


def qiskit_exact_state(coupling, time, omega=1.0, delta=0.0):
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import HamiltonianGate
    from qiskit.quantum_info import Statevector
    number = len(coupling)
    circuit = QuantumCircuit(number)
    matrix = qiskit_hamiltonian(coupling, omega, delta).to_matrix()
    circuit.append(HamiltonianGate(matrix, time), range(number))
    return Statevector.from_instruction(circuit).data


def stepped_circuit(coupling, time, steps, omega=1.0, delta=0.0):
    """Second-order split: half laser step, full energy step, half laser step."""
    from qiskit import QuantumCircuit
    number = len(coupling)
    circuit = QuantumCircuit(number)
    dt = time / steps
    for step in range(steps):
        for i in range(number):
            circuit.rx(omega * dt / 2, i)
        for i in range(number):
            circuit.p(delta * dt, i)
            for j in range(i + 1, number):
                if coupling[i, j] != 0:
                    circuit.cp(-coupling[i, j] * dt, i, j)
        for i in range(number):
            circuit.rx(omega * dt / 2, i)
    return circuit
