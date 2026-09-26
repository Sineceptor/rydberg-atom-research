"""Run the experiments and save every plotted value. No cloud account needed."""
from pathlib import Path
import csv
import json
import platform
import time
import numpy as np
import scipy
import qiskit
import matplotlib
if __name__ == "__main__":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace, entropy as qiskit_entropy
from rydberg_model import (chain, interactions, hamiltonian, evolve, occupations,
    excitation_probabilities, pair_probability, entropy, qiskit_hamiltonian,
    qiskit_exact_state, stepped_circuit)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "results"
FIGURES = ROOT / "figures"
plt.rcParams.update({"font.size": 11, "figure.dpi": 140,
                     "axes.spines.top": False, "axes.spines.right": False})
COLORS = ["#2166ac", "#c47a16", "#76519c", "#57823d"]
plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=COLORS)


def save_csv(name, headers, rows):
    with (DATA / name).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(rows)


def save_figure(name):
    plt.savefig(FIGURES / (name + ".png"), dpi=180, bbox_inches="tight")
    plt.close()


def run(output_root=None):
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "august")
    global DATA, FIGURES
    DATA, FIGURES = Path(output_root) / "results", Path(output_root) / "figures"
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    summary = {"versions": {"python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "qiskit": qiskit.__version__, "matplotlib": matplotlib.__version__},
        "model": "hbar=1; omega=1; delta=0 unless varied; ground initial state; all pairs r^-6"}
    checks = {}

    # 8 August: a current Qiskit baseline, including a Bell-state check.
    times = np.linspace(0, 2 * np.pi, 161)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    baseline_rows = []
    for strength, style in [(0, "-"), (20, "--")]:
        coupling = interactions(chain(2), strength)
        states = np.array([qiskit_exact_state(coupling, t) for t in times])
        sparse_states = evolve(hamiltonian(coupling), times)
        checks[f"qiskit_scipy_state_error_V{strength}"] = float(np.max(np.abs(states-sparse_states)))
        populations = excitation_probabilities(states, 2).sum(axis=1)
        pair = pair_probability(states, 2, 0, 1)
        entanglement = [entropy(state, [0]) for state in states]
        axes[0].plot(times, populations, style, label=f"V / Ω = {strength}")
        axes[1].plot(times, pair, style, label=f"V / Ω = {strength}")
        baseline_rows.extend(zip(np.full(len(times), strength), times, populations, pair, entanglement))
    axes[0].plot(times, np.sin(np.sqrt(2)*times/2)**2, ":", color="black", label="Ideal strong blockade")
    axes[0].set(ylabel="Mean number of excited atoms", ylim=(0, 2.05))
    axes[1].set(ylabel="Probability both atoms are excited", ylim=(0, 1.05))
    for ax in axes:
        ax.set_xlabel("Time Ωt")
        ax.legend(fontsize=9)
    save_figure("01_qiskit_baseline")
    save_csv("qiskit_baseline.csv", ["V_over_omega", "omega_t", "mean_excitations", "P11", "entropy_bits"], baseline_rows)

    bell = QuantumCircuit(2)
    bell.h(0)
    bell.h(1)
    bell.cp(np.pi, 0, 1)
    bell.h(1)
    bell_state = Statevector.from_instruction(bell)
    checks["bell_entropy_bits"] = entropy(bell_state.data, [0])
    summary["bell_probabilities"] = {key: float(value) for key, value in bell_state.probabilities_dict().items() if value > 1e-12}
    (DATA / "bell_circuit.txt").write_text(str(bell.draw(output="text")), encoding="utf-8")
    (DATA / "one_time_step_circuit.txt").write_text(str(stepped_circuit(interactions(chain(2), 5), .1, 1).draw(output="text")), encoding="utf-8")

    # Parameter example: keep the interaction fixed and change detuning.
    rows = []
    plt.figure(figsize=(7, 3.8))
    for delta, style in zip([-2, 0, 2], ["-", "--", ":"]):
        states = evolve(hamiltonian(interactions(chain(2), 5), delta=delta), times)
        total = excitation_probabilities(states, 2).sum(axis=1)
        plt.plot(times, total, style, label=f"Δ / Ω = {delta}")
        rows.extend(zip(np.full(len(times), delta), times, total))
    plt.xlabel("Time Ωt")
    plt.ylabel("Mean number of excited atoms")
    plt.legend()
    save_figure("02_detuning")
    save_csv("detuning.csv", ["delta_over_omega", "omega_t", "mean_excitations"], rows)

    # Reproduce the idea of the tutorial's spin diagrams with seven atoms.
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    for ax, strength in zip(axes, [0, 5]):
        states = evolve(hamiltonian(interactions(chain(7), strength)), times)
        populations = excitation_probabilities(states, 7)
        view = ax.pcolormesh(times, np.arange(7), populations.T, vmin=0, vmax=1, cmap="viridis", shading="auto")
        ax.set(title=f"Seven atoms: V nearest / Ω = {strength}", xlabel="Time Ωt", ylabel="Atom number")
        fig.colorbar(view, ax=ax, label="Excitation probability")
        save_csv(f"spin_diagram_V{strength}.csv", ["omega_t"]+[f"atom_{i}" for i in range(7)], np.column_stack((times, populations)))
    save_figure("03_spin_diagrams")

    # 15 August, part 1: blockade. Compare the entire same time interval.
    scan_times = np.linspace(0, 2*np.pi, 801)
    strengths = [0, .5, 1, 2, 5, 10, 20, 40]
    rows = []
    for strength in strengths:
        coupling = interactions(chain(2), strength)
        states = evolve(hamiltonian(coupling), scan_times)
        pair = pair_probability(states, 2, 0, 1)
        rows.append([strength, np.trapezoid(pair, scan_times)/(2*np.pi), pair.max(), entropy(states[200], [0])])
    summary["blockade"] = rows
    save_csv("blockade_scan.csv", ["V_over_omega", "time_mean_P11", "sampled_max_P11", "entropy_at_pi_over_2_bits"], rows)
    values = np.array(rows)
    plt.figure(figsize=(7, 3.8))
    plt.plot(values[:, 0], values[:, 1], "o-", label="Time mean")
    plt.plot(values[:, 0], values[:, 2], "s--", label="Largest sampled value")
    plt.xlabel("Interaction strength V / Ω")
    plt.ylabel("Probability both atoms are excited")
    plt.ylim(0, 1.05)
    plt.legend()
    save_figure("04_blockade_scan")

    # Part 2: isolate atom 2, then compare with a completely split chain.
    defect_times = np.linspace(0, 4, 161)
    defect_cases = {"Intact": interactions(chain(6), 5),
        "Atom 2 isolated": interactions(chain(6), 5, isolated=2),
        "Two disconnected halves": interactions(chain(6), 5, split=[0, 1, 2])}
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    summary["defects"] = []
    for (name, coupling), style in zip(defect_cases.items(), ["-", "--", ":"]):
        matrix = hamiltonian(coupling)
        states = evolve(matrix, defect_times)
        single = np.array([entropy(state, [2]) for state in states])
        half = np.array([entropy(state, [0, 1, 2]) for state in states])
        axes[0].plot(defect_times, single, style, label=name)
        axes[1].plot(defect_times, half, style, label=name)
        populations = excitation_probabilities(states, 6)
        save_csv("defect_"+name.lower().replace(" ", "_")+".csv",
                 ["omega_t", "S_atom2_bits", "S_012_bits", "P_atom2"],
                 zip(defect_times, single, half, populations[:, 2]))
        summary["defects"].append({"case": name, "S_atom2_at_t2": float(single[80]),
            "S_half_at_t2": float(half[80]), "maximum_S_atom2": float(single.max()), "maximum_S_half": float(half.max())})
        checks[name+"_norm_error"] = float(np.max(np.abs(np.sum(abs(states)**2, axis=1)-1)))
        energies = np.array([np.vdot(state, matrix @ state).real for state in states])
        checks[name+"_energy_drift"] = float(np.max(np.abs(energies-energies[0])))
        if name == "Atom 2 isolated":
            checks["isolated_atom_max_entropy"] = float(single.max())
            checks["isolated_atom_rabi_error"] = float(np.max(np.abs(populations[:, 2]-np.sin(defect_times/2)**2)))
        if name == "Two disconnected halves":
            checks["disconnected_halves_max_entropy"] = float(half.max())
    axes[0].set(ylabel="Atom 2 vs other atoms: entropy (bits)", ylim=(-.02, 1.03))
    axes[1].set(ylabel="Atoms 0, 1, 2 vs 3, 4, 5: entropy (bits)")
    for ax in axes:
        ax.set_xlabel("Time Ωt")
        ax.legend(fontsize=8)
    save_figure("05_defects")

    # Part 3: identical atom count, initial state, laser and minimum spacing.
    geometries = {"Chain": chain(4),
        "Square": np.array([[0, 0], [1, 0], [0, 1], [1, 1]]),
        "Rectangle": np.array([[0, 0], [1, 0], [0, 2], [1, 2]])}
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.8), constrained_layout=True)
    for ax, (name, positions) in zip(axes, geometries.items()):
        ax.scatter(positions[:, 0], positions[:, 1], s=90, color=COLORS[0])
        for i, (x, y) in enumerate(positions):
            ax.annotate(str(i), (x, y), xytext=(5, 6), textcoords="offset points")
        ax.set(title=name, xlabel="x / a", ylabel="y / a", aspect="equal")
        ax.set_xlim(-.4, max(positions[:, 0])+.5)
        ax.set_ylim(-.4, max(positions[:, 1])+.6)
    save_figure("06_geometries")
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    summary["geometry"] = []
    for (name, positions), style in zip(geometries.items(), ["-", "--", ":"]):
        states = evolve(hamiltonian(interactions(positions, 5)), defect_times)
        entropy_01 = np.array([entropy(state, [0, 1]) for state in states])
        entropy_02 = np.array([entropy(state, [0, 2]) for state in states])
        axes[0].plot(defect_times, entropy_01, style, label=name+"; keep 0, 1")
        if name != "Chain":
            axes[1].plot(defect_times, entropy_01, style, label=name+": lower / upper")
            axes[1].plot(defect_times, entropy_02, style, marker="o", markevery=25, markersize=3, label=name+": left / right")
        summary["geometry"].append({"shape": name, "S_01_at_t2": float(entropy_01[80]), "S_02_at_t2": float(entropy_02[80])})
        save_csv("geometry_"+name.lower()+".csv", ["omega_t", "S_01_bits", "S_02_bits"], zip(defect_times, entropy_01, entropy_02))
        if name == "Square":
            checks["square_direction_symmetry_error"] = float(max(abs(entropy_01-entropy_02)))
    for ax in axes:
        ax.set(xlabel="Time Ωt", ylabel="Two atoms vs two atoms: entropy (bits)", ylim=(-.02, 2.02))
        ax.legend(fontsize=8)
    save_figure("07_geometry_entanglement")

    # Part 4: actual timed runs. This is a modest benchmark, not a size limit.
    rows = []
    for number in [4, 6, 8, 10, 12, 14]:
        start = time.perf_counter()
        matrix = hamiltonian(interactions(chain(number), 5))
        build_seconds = time.perf_counter()-start
        start = time.perf_counter()
        states = evolve(matrix, np.linspace(0, 1, 21))
        evolution_seconds = time.perf_counter()-start
        matrix_mib = (matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes)/2**20
        rows.append([number, 2**number, build_seconds, evolution_seconds, matrix_mib,
                     states.nbytes/2**20, float(np.max(abs(np.sum(abs(states)**2, axis=1)-1)))])
    summary["scaling"] = rows
    summary["benchmark_machine"] = {"platform": platform.platform(), "processor": platform.machine(), "repeats": 1,
        "interval": "0 <= omega*t <= 1", "saved_times": 21,
        "scope": "build and evolve only; excludes plotting and entanglement; stored-array memory is not peak process memory"}
    save_csv("scaling.csv", ["atoms", "basis_states", "build_seconds", "evolution_seconds", "sparse_H_MiB", "21_states_MiB", "norm_error"], rows)
    values = np.array(rows)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    axes[0].semilogy(values[:, 0], values[:, 3], "o-")
    axes[0].set(xlabel="Number of atoms", ylabel="Evolution time (s, log scale)")
    axes[1].semilogy(values[:, 0], values[:, 4], "o-", label="Sparse Hamiltonian")
    axes[1].semilogy(values[:, 0], values[:, 5], "s--", label="21 saved states")
    axes[1].set(xlabel="Number of atoms", ylabel="Stored arrays (MiB, log scale)")
    axes[1].legend(fontsize=9)
    save_figure("08_scaling")

    # Independent physics and numerical checks.
    test_coupling = interactions(chain(3), 3)
    matrix = hamiltonian(test_coupling, omega=.7, delta=.4).toarray()
    checks["Hamiltonian_Hermitian_error"] = float(np.max(abs(matrix-matrix.conj().T)))
    checks["Pauli_vs_bit_basis_H_error"] = float(np.max(abs(matrix-qiskit_hamiltonian(test_coupling, .7, .4).to_matrix())))
    expected_two_atom = np.array([[0, .5, .5, 0], [.5, -.4, 0, .5], [.5, 0, -.4, .5], [0, .5, .5, 4.2]])
    checks["explicit_two_atom_H_error"] = float(np.max(abs(hamiltonian(interactions(chain(2), 5), delta=.4).toarray()-expected_two_atom)))
    free = evolve(hamiltonian(interactions(chain(4), 0)), times)
    checks["free_atom_analytic_probability_error"] = float(np.max(abs(excitation_probabilities(free, 4)-np.sin(times[:, None]/2)**2)))
    checks["free_atoms_max_entropy"] = float(max(entropy(state, [0, 1]) for state in free))
    rng = np.random.default_rng(2026)
    state = rng.normal(size=32)+1j*rng.normal(size=32)
    state /= np.linalg.norm(state)
    checks["entropy_vs_Qiskit_error"] = abs(entropy(state, [0, 3])-float(qiskit_entropy(partial_trace(Statevector(state), [1, 2, 4]), base=2)))

    circuit_rows = []
    coupling = interactions(chain(4), 5)
    target = qiskit_exact_state(coupling, 1)
    for steps in [10, 20, 40, 80]:
        circuit = stepped_circuit(coupling, 1, steps)
        approximate = Statevector.from_instruction(circuit).data
        error = max(0., float(1-abs(np.vdot(target, approximate))**2))
        circuit_rows.append([steps, 1/steps, error, circuit.size()])
    summary["circuit_convergence"] = circuit_rows
    save_csv("circuit_convergence.csv", ["steps", "delta_t", "infidelity", "gates"], circuit_rows)
    checks["circuit_error_decreases"] = bool(all(circuit_rows[i+1][2] < circuit_rows[i][2] for i in range(3)))
    fine = np.linspace(0, 2*np.pi, 1601)
    pair = pair_probability(evolve(hamiltonian(interactions(chain(2), 20)), fine), 2, 0, 1)
    checks["blockade_time_mean_grid_change"] = abs(float(np.trapezoid(pair, fine)/(2*np.pi))-summary["blockade"][6][1])
    checks["blockade_sampled_peak_grid_change"] = abs(float(pair.max())-summary["blockade"][6][2])
    assert abs(checks["bell_entropy_bits"]-1) < 1e-10
    assert checks["circuit_error_decreases"]
    for name, value in checks.items():
        if name not in ["bell_entropy_bits", "circuit_error_decreases"]:
            assert abs(value) < 1e-6, (name, value)
    summary["checks"] = checks
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    run()
