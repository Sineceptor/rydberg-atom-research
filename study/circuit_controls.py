"""Nonzero-detuning circuit checks against independent dense SciPy evolution.

The detuning sign is tested directly with an observable phase, and the
second-order split circuit is tested at both delta = -1.3 and +1.3.
No results are written on import; run() uses the managed output wrapper.
"""

from pathlib import Path
import csv
import json
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

from rydberg_model import chain, interactions, stepped_circuit, hamiltonian


def tensor_reference(coupling, omega, delta):
    """Build a dense Hamiltonian independently using Kronecker products."""
    number = len(coupling)
    identity = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    n = np.diag([0., 1.])

    def at(operators):
        result = np.array([[1.]])
        for atom in reversed(range(number)):
            result = np.kron(result, operators.get(atom, identity))
        return result

    matrix = np.zeros((2**number, 2**number))
    for i in range(number):
        matrix += omega/2*at({i: x}) - delta*at({i: n})
        for j in range(i+1, number):
            matrix += coupling[i, j]*at({i: n, j: n})
    return matrix


def phase_aligned_error(reference, trial):
    """Remove an unobservable overall phase before comparing states."""
    overlap = np.vdot(reference, trial)
    aligned = trial * np.exp(-1j*np.angle(overlap))
    return float(np.linalg.norm(aligned-reference))


def csv_write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(output_root=None):
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "circuits")
    root = Path(output_root)
    results, figures = root / "results", root / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    rows, phase_rows, orders = [], [], []
    detunings = [-1.3, 1.3]
    step_counts = [8, 16, 32, 64, 128]
    coupling = interactions(chain(3), strength=5.)
    omega, duration = 1., 1.7
    initial = np.zeros(8, dtype=complex)
    initial[0] = 1.
    maximum_matrix_error = 0.
    for delta in detunings:
        matrix = tensor_reference(coupling, omega, delta)
        matrix_error = float(np.max(abs(matrix-hamiltonian(coupling, omega, delta).toarray())))
        maximum_matrix_error = max(maximum_matrix_error, matrix_error)
        assert matrix_error < 1e-12
        reference = expm(-1j*duration*matrix) @ initial
        assert abs(np.linalg.norm(reference)-1.) < 1e-12
        errors = []
        for steps in step_counts:
            circuit = stepped_circuit(coupling, duration, steps, omega=omega, delta=delta)
            trial = Statevector.from_instruction(circuit).data
            norm_error = float(abs(np.linalg.norm(trial)-1.))
            fidelity = float(abs(np.vdot(reference, trial))**2)
            state_error = phase_aligned_error(reference, trial)
            errors.append(state_error)
            rows.append({"detuning": delta, "atoms": 3, "omega": omega,
                         "nearest_interaction": 5., "duration": duration, "steps": steps,
                         "phase_aligned_state_error": state_error,
                         "infidelity": max(0., 1.-fidelity), "norm_error": norm_error})
            assert norm_error < 1e-11
        observed_order = float(-np.polyfit(np.log(step_counts[-3:]), np.log(errors[-3:]), 1)[0])
        assert np.all(np.diff(errors) < 0), "The state error should decrease on this tested grid."
        assert 1.8 < observed_order < 2.2, "The fine-step state error should scale approximately as steps^-2."
        assert errors[-1] < 1e-4
        orders.append({"detuning": delta, "observed_state_error_order_last_three": observed_order,
                       "finest_step_state_error": errors[-1],
                       "successive_state_error_ratios": [errors[i]/errors[i+1] for i in range(len(errors)-1)]})

        # Laser off, one atom initially |+>: H = -delta*n so the |1> phase is
        # exp(+i delta*t). This exposes the sign that a |0> population cannot.
        phase_time = .7
        preparation = QuantumCircuit(1)
        preparation.h(0)
        preparation.compose(stepped_circuit(np.zeros((1, 1)), phase_time, 1,
                                           omega=0., delta=delta), inplace=True)
        actual = Statevector.from_instruction(preparation).data
        expected = np.array([1., np.exp(1j*delta*phase_time)])/np.sqrt(2)
        reference_phase = expm(-1j*phase_time*tensor_reference(np.zeros((1, 1)), 0., delta)) @ (np.ones(2)/np.sqrt(2))
        phase_error = float(np.max(abs(actual-expected)))
        independent_error = float(np.max(abs(actual-reference_phase)))
        y = np.array([[0., -1j], [1j, 0.]])
        measured_y = float(np.vdot(actual, y @ actual).real)
        expected_y = float(np.sin(delta*phase_time))
        assert phase_error < 1e-12 and independent_error < 1e-12
        assert abs(measured_y-expected_y) < 1e-12
        assert measured_y*delta > 0
        wrong_sign = np.array([1., np.exp(-1j*delta*phase_time)])/np.sqrt(2)
        wrong_sign_error = phase_aligned_error(expected, wrong_sign)
        assert wrong_sign_error > .5, "This phase test must be sensitive to the sign."
        phase_rows.append({"detuning": delta, "duration": phase_time,
                           "expected_y": expected_y, "circuit_y": measured_y,
                           "analytic_amplitude_max_error": phase_error,
                           "scipy_amplitude_max_error": independent_error,
                           "deliberately_reversed_sign_state_error": wrong_sign_error})
    csv_write(results/"circuit_detuning_convergence.csv", rows)
    csv_write(results/"circuit_detuning_phase_sign.csv", phase_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    for delta, color in zip(detunings, ["#a95725", "#006b83"]):
        chosen = [r for r in rows if r["detuning"] == delta]
        axes[0].loglog(step_counts, [r["phase_aligned_state_error"] for r in chosen], "o-", color=color,
                       label=f"Δ / Ω = {delta:+.1f}")
        axes[1].loglog(step_counts, [r["infidelity"] for r in chosen], "o-", color=color)
    guide = rows[0]["phase_aligned_state_error"]*(np.array(step_counts)/step_counts[0])**-2.
    axes[0].loglog(step_counts, guide, "--", color="gray", label="Reference slope: steps⁻²")
    for ax in axes:
        ax.set_xlabel("Number of split steps")
        ax.grid(alpha=.2, which="both")
    axes[0].set_ylabel("State error after removing overall phase")
    axes[1].set_ylabel("1 − fidelity with independent SciPy state")
    axes[0].legend(fontsize=9)
    fig.suptitle("Circuit control at both detuning signs\nThree atoms, all pairs, Vnearest / Ω = 5, Ωt = 1.7", fontsize=12)
    fig.savefig(figures/"circuit_detuning_convergence.png", dpi=180, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)
    summary = {"purpose": "Test both nonzero detuning signs and second-order split convergence.",
               "detunings": detunings, "step_counts": step_counts,
               "maximum_independent_tensor_vs_sparse_matrix_error": maximum_matrix_error,
               "convergence": orders, "phase_sign": phase_rows,
               "maximum_circuit_norm_error": max(r["norm_error"] for r in rows),
               "all_checks_passed": True,
               "scope": "Three-atom closed-system state simulation, not hardware fidelity or ground-state preparation."}
    (results/"circuit_controls_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    run()
