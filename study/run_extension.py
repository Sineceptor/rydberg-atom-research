"""Run the declared finite-chain comparison and its independent controls.

Run: python run_extension.py
The CLI writes raw CSV, checks and figures into a separate managed run folder.
The saved grid is evidence for this scan, not a continuous maximum or a
thermodynamic phase boundary. See VALIDATION_COMPUTE.md after running.
"""
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import platform
import time
import numpy as np
from scipy.linalg import eigh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from extension_model import chain_couplings, ground_point, discarded_tail_bound
from rydberg_model import hamiltonian, qiskit_hamiltonian

ROOT = Path(__file__).resolve().parent
DATA, FIG = ROOT / "results", ROOT / "figures"
SIZES = list(range(4, 13))
STRENGTHS = [2.0, 8.0, 16.0]
DELTAS = np.linspace(-3.0, 7.0, 81)
CUTOFFS = ["nearest1", "nearest2", "all"]
LABELS = {"nearest1": "Nearest neighbours only", "nearest2": "Up to second neighbours", "all": "All pairs"}
COLOURS = {"nearest1": "#d55e00", "nearest2": "#0072b2", "all": "#222222"}


def write_csv(name, rows):
    with (DATA / name).open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def key(number, strength, delta, cutoff):
    return (number, strength, round(float(delta), 12), cutoff)


def comparison(row, reference):
    return {"N": row["N"], "g": row["g"], "delta": row["delta"], "cutoff": row["cutoff"],
            "tail_norm_bound": discarded_tail_bound(row["N"], row["g"], row["cutoff"]),
            "ground_energy_shift_full_minus_cut": reference["ground_energy"]-row["ground_energy"],
            "first_excited_energy_shift_full_minus_cut": reference["ground_energy"]+reference["gap"]-row["ground_energy"]-row["gap"],
            **{f"absolute_{name}_error": abs(row[name]-reference[name]) for name in
               ["density", "alternating_order_squared", "staggered_spin_squared", "gap", "centre_density", "edge_density"]}}


def check_bundled_reference():
    """Fail early if the unchanged-model historical regression cannot pass."""
    # A bundled immutable historical reference makes the regression portable.
    # Missing or changed reference data is an error, never a silent skipped test.
    ref_dir = ROOT / "reference_data"
    metadata = json.loads((ref_dir / "reference_metadata.json").read_text())
    prior = ref_dir / metadata["filename"]
    digest = hashlib.sha256(prior.read_bytes()).hexdigest()
    assert digest == metadata["sha256"], "Bundled historical reference digest mismatch"
    prior_errors = {name: 0. for name in ["density", "alternating_order_squared", "gap"]}
    with prior.open() as file:
        old_rows = list(csv.DictReader(file))
    assert len(old_rows) == metadata["rows"] == 62
    for old in old_rows:
        cutoff = "all" if old["interactions"] == "all pairs" else "nearest1"
        row, _, _ = ground_point(8, 8., float(old["delta_over_omega"]), cutoff)
        for name, old_name in [("density", "density"), ("alternating_order_squared", "alternating_order_squared"), ("gap", "gap_over_omega")]:
            prior_errors[name] = max(prior_errors[name], abs(row[name]-float(old[old_name])))
    result = {"points": len(old_rows), "maximum_errors": prior_errors, "sha256": digest,
        "status": "mandatory bundled reference passed", "reference_file": str(prior.relative_to(ROOT))}
    assert max(prior_errors.values()) < 1e-8, "Historical reference regression failed"
    return result


def independent_checks(rows, saved, reference_check=None):
    checks, control_rows = {}, []
    # The two-atom matrix is written by hand, separately from the bit builder.
    delta, g = 0.7, 3.2
    manual = np.array([[0, .5, .5, 0], [.5, -delta, 0, .5],
                       [.5, 0, -delta, .5], [0, .5, .5, -2*delta+g]])
    actual = hamiltonian(chain_couplings(2, g), delta=delta).toarray()
    checks["handwritten_two_atom_matrix_error"] = float(np.max(abs(manual-actual)))
    # Qiskit constructs the operator as Pauli strings, independently of bit flips.
    errors = []
    for cutoff in CUTOFFS:
        for g, delta in [(2., -1.5), (8., 1.0), (16., 5.0)]:
            coupling = chain_couplings(4, g, cutoff)
            scipy_matrix = hamiltonian(coupling, delta=delta).toarray()
            qiskit_matrix = qiskit_hamiltonian(coupling, delta=delta).to_matrix()
            errors.append(float(np.max(abs(scipy_matrix-qiskit_matrix))))
    checks["qiskit_operator_max_entry_error"] = max(errors)
    # Independent product-state formula for g=0 at several detunings.
    errors, staggered_errors = [], []
    for number in [5, 6, 7]:
        for delta in [-3., -.3, 0., .7, 7.]:
            row, _, _ = ground_point(number, 0, delta, method="sparse")
            exact = (1+delta/np.sqrt(1+delta**2))/2
            exact_s2 = 4*exact*(1-exact)/number + (2*exact-1)**2*(number%2/number)**2
            errors.append(abs(row["density"]-exact))
            staggered_errors.append(abs(row["staggered_spin_squared"]-exact_s2))
    checks["noninteracting_density_error"] = float(max(errors))
    checks["noninteracting_centered_staggered_error"] = float(max(staggered_errors))
    checks["maximum_absolute_eigenpair_residual"] = max(r["residual"] for r in rows)
    checks["maximum_relative_eigenpair_residual"] = max(r["relative_residual"] for r in rows)
    checks["maximum_eigenvector_overlap"] = max(r["orthogonality_error"] for r in rows)
    checks["maximum_normalization_error"] = max(r["normalization_error"] for r in rows)
    # Adding nonnegative diagonal interactions cannot lower E0 (variational principle).
    tail_violation = 0.0
    for number in SIZES:
        for g in STRENGTHS:
            for delta in DELTAS:
                energies = [saved[key(number, g, delta, c)]["ground_energy"] for c in CUTOFFS]
                tail_violation = max(tail_violation, *[-float(x) for x in np.diff(energies)])
    checks["positive_tail_energy_monotonicity_violation"] = tail_violation
    bound_rows, bound_violation, first_energy_violation, norm_error = [], 0.0, 0.0, 0.0
    for number in SIZES:
        for g in STRENGTHS:
            for cutoff in CUTOFFS[:2]:
                bound = discarded_tail_bound(number, g, cutoff)
                # Direct diagonal maximum verifies the combinatorial formula.
                full_h = hamiltonian(chain_couplings(number, g, "all"))
                cut_h = hamiltonian(chain_couplings(number, g, cutoff))
                measured = float((full_h-cut_h).diagonal().max())
                norm_error = max(norm_error, abs(measured-bound))
                largest_energy_shift, largest_gap_change = 0.0, 0.0
                for delta in DELTAS:
                    full = saved[key(number,g,delta,"all")]
                    cut = saved[key(number,g,delta,cutoff)]
                    energy_shifts = [full["ground_energy"]-cut["ground_energy"],
                                     full["ground_energy"]+full["gap"]-cut["ground_energy"]-cut["gap"]]
                    gap_change = abs(full["gap"]-cut["gap"])
                    largest_energy_shift = max(largest_energy_shift, *energy_shifts)
                    largest_gap_change = max(largest_gap_change, gap_change)
                    first_energy_violation = max(first_energy_violation, -energy_shifts[1])
                    bound_violation = max(bound_violation, max(energy_shifts)-bound, gap_change-bound)
                bound_rows.append({"N":number,"g":g,"cutoff":cutoff,"tail_norm_bound":bound,
                                   "largest_sampled_energy_shift":largest_energy_shift,
                                   "largest_sampled_gap_change":largest_gap_change})
    write_csv("discarded_tail_bounds.csv", bound_rows)
    checks["tail_norm_formula_error"] = norm_error
    checks["tail_energy_and_gap_bound_violation"] = bound_violation
    checks["first_excited_energy_monotonicity_violation"] = first_energy_violation
    # E0 is concave in Delta, so -dE0/dDelta=<sum n> is nondecreasing.
    density_violation = 0.0
    for number in SIZES:
        for g in STRENGTHS:
            for cutoff in CUTOFFS:
                values = [saved[key(number, g, d, cutoff)]["density"] for d in DELTAS]
                density_violation = max(density_violation, -float(np.min(np.diff(values))))
    checks["density_detuning_monotonicity_violation"] = density_violation
    # Dense and sparse algorithms on varied points, including each cutoff's
    # lowest coarse-grid gap for N=8 and every g. No continuity assumption.
    points = [(4, 2., -3., "all"), (6, 16., 7., "nearest1"), (8, 8., 0., "nearest2")]
    for g in STRENGTHS:
        for cutoff in CUTOFFS:
            candidates = [r for r in rows if r["N"] == 8 and r["g"] == g and r["cutoff"] == cutoff]
            low = min(candidates, key=lambda r: r["gap"])
            points.append((8, g, low["delta"], cutoff))
    energy_errors, gap_errors, density_errors = [], [], []
    for number, g, delta, cutoff in points:
        sparse, _, _ = ground_point(number, g, delta, cutoff, method="sparse")
        dense, _, _ = ground_point(number, g, delta, cutoff, method="dense")
        energy_errors.append(abs(sparse["ground_energy"]-dense["ground_energy"]))
        gap_errors.append(abs(sparse["gap"]-dense["gap"]))
        density_errors.append(abs(sparse["density"]-dense["density"]))
        control_rows.append({"N": number, "g": g, "delta": delta, "cutoff": cutoff,
                             "dense_gap": dense["gap"], "energy_error": energy_errors[-1],
                             "gap_error": gap_errors[-1], "density_error": density_errors[-1]})
    write_csv("dense_sparse_controls.csv", control_rows)
    checks["dense_sparse_energy_error"] = max(energy_errors)
    checks["dense_sparse_gap_error"] = max(gap_errors)
    checks["dense_sparse_density_error"] = max(density_errors)
    checks["prior_csv_matching"] = reference_check or check_bundled_reference()
    tolerances = {
        "handwritten_two_atom_matrix_error": 1e-12, "qiskit_operator_max_entry_error": 1e-12,
        "noninteracting_density_error": 1e-10, "maximum_absolute_eigenpair_residual": 1e-7,
        "noninteracting_centered_staggered_error": 1e-10,
        "maximum_relative_eigenpair_residual": 1e-9, "maximum_eigenvector_overlap": 1e-10,
        "maximum_normalization_error": 1e-10, "positive_tail_energy_monotonicity_violation": 1e-9,
        "density_detuning_monotonicity_violation": 1e-9, "dense_sparse_energy_error": 1e-9,
        "dense_sparse_gap_error": 1e-9, "dense_sparse_density_error": 1e-9,
        "tail_norm_formula_error": 1e-11, "tail_energy_and_gap_bound_violation": 1e-9,
        "first_excited_energy_monotonicity_violation": 1e-9,
    }
    for name, tolerance in tolerances.items():
        assert checks[name] < tolerance, f"{name}: {checks[name]} >= {tolerance}"
    if "maximum_errors" in checks["prior_csv_matching"]:
        assert max(checks["prior_csv_matching"]["maximum_errors"].values()) < 1e-8
    return checks, tolerances


def make_figures(saved, differences, maxima, profile, selected):
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 190})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for column, number in enumerate([4, 8, 12]):
        for cutoff in CUTOFFS:
            values = [saved[key(number, 8., d, cutoff)]["density"] for d in DELTAS]
            axes[0, column].plot(DELTAS, values, label=LABELS[cutoff], color=COLOURS[cutoff])
            if cutoff != "all":
                errors = [abs(saved[key(number, 8., d, cutoff)]["density"]-saved[key(number, 8., d, "all")]["density"]) for d in DELTAS]
                axes[1, column].plot(DELTAS, errors, color=COLOURS[cutoff])
        axes[0, column].set(title=f"N = {number}", ylabel="Mean excitation density", ylim=(-.02, 1.02))
        axes[1, column].set(xlabel="Detuning Δ / Ω", ylabel="Absolute density difference")
        axes[1, column].grid(alpha=.18)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle("Effect of interaction range at g = V₀ / Ω = 8\nOpen finite chains; ground states on a 0.125 detuning grid")
    fig.savefig(FIG / "01_density_comparison.png", bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.2), constrained_layout=True)
    for ax, g in zip(axes, STRENGTHS):
        for cutoff in CUTOFFS[:2]:
            for parity, marker, style in [(0,"o","-"),(1,"s","--")]:
                selection = sorted([r for r in maxima if r["g"] == g and r["cutoff"] == cutoff and r["N"]%2==parity], key=lambda r:r["N"])
                ax.plot([r["N"] for r in selection], [r["refined_sampled_max_absolute_density_error"] for r in selection], marker+style, color=COLOURS[cutoff], label=LABELS[cutoff]+("; even" if parity==0 else "; odd"))
        ax.set(title=f"g = {g:g}", xlabel="Number of atoms N", ylabel="Largest sampled density difference", xticks=SIZES)
        ax.grid(alpha=.18)
    axes[0].legend(fontsize=8)
    fig.suptitle("Interaction truncation error: separate even and odd chains\nLargest sampled density difference after local half-spacing refinement")
    fig.savefig(FIG / "02_max_error_by_size.png", bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(9.4, 6.8), sharex=True, constrained_layout=True)
    reference = np.array([r["density"] for r in profile if r["cutoff"] == "all"])
    for cutoff in CUTOFFS:
        data = [r for r in profile if r["cutoff"] == cutoff]
        values = np.array([r["density"] for r in data])
        atoms = [r["atom"] for r in data]
        axes[0].plot(atoms, values, "o-", color=COLOURS[cutoff], label=LABELS[cutoff])
        if cutoff != "all":
            axes[1].plot(atoms, values-reference, "o-", color=COLOURS[cutoff])
    axes[0].set(ylabel="Excitation probability at each atom", ylim=(0,1))
    axes[0].legend(fontsize=9)
    axes[1].axhline(0, color="grey", linewidth=.8)
    axes[1].set(ylabel="Signed change from all pairs", xlabel="Atom index (open ends at 0 and 11)", xticks=range(12))
    for ax in axes:
        ax.axvspan(-.3,.3,color="grey",alpha=.10)
        ax.axvspan(10.7,11.3,color="grey",alpha=.10)
        ax.grid(alpha=.15)
    fig.suptitle(f"Spatial profile at a high-error sampled point\nN = 12, g = {selected['g']:g}, Δ / Ω = {selected['delta']:g}")
    fig.savefig(FIG / "03_boundary_profile.png", bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.2), constrained_layout=True)
    for cutoff in CUTOFFS:
        for ax, observable in zip(axes, ["staggered_spin_squared", "gap", "centre_density"]):
            values = [saved[key(12, selected["g"], d, cutoff)][observable] for d in DELTAS]
            ax.plot(DELTAS, values, color=COLOURS[cutoff], label=LABELS[cutoff])
    for ax, title in zip(axes, ["Staggered spin squared ⟨S²⟩", "Lowest energy gap / Ω", "Interior mean density"]):
        ax.set(xlabel="Detuning Δ / Ω", ylabel=title)
        ax.grid(alpha=.18)
    axes[0].legend(fontsize=9)
    fig.suptitle(f"Other diagnostics for N = 12 and g = {selected['g']:g}\nInterior excludes the two endpoint atoms; these curves do not assign bulk phases")
    fig.savefig(FIG / "04_other_diagnostics.png", bbox_inches="tight")
    plt.close(fig)


def run(output_root=None):
    """Compute in an explicit staging folder, or create a fresh managed run."""
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "chains")
    global DATA, FIG
    DATA, FIG = Path(output_root) / "results", Path(output_root) / "figures"
    start = time.perf_counter()
    reference_check = check_bundled_reference()
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rows, saved = [], {}
    for number in SIZES:
        for g in STRENGTHS:
            for delta in DELTAS:
                for cutoff in CUTOFFS:
                    row, _, _ = ground_point(number, g, delta, cutoff)
                    rows.append(row)
                    saved[key(number, g, delta, cutoff)] = row
            print(f"Finished N={number}, g={g:g}; elapsed {time.perf_counter()-start:.1f} s", flush=True)
    write_csv("ground_states_coarse.csv", rows)
    differences = [comparison(r, saved[key(r["N"], r["g"], r["delta"], "all")])
                   for r in rows if r["cutoff"] != "all"]
    write_csv("cutoff_errors_coarse.csv", differences)
    # Refine ±one coarse step around each combination's sampled maximum.
    # This tests local grid sensitivity; it does not bound the continuous max.
    refinements, maxima = [], []
    for number in SIZES:
        for g in STRENGTHS:
            for cutoff in CUTOFFS[:2]:
                coarse = [r for r in differences if r["N"] == number and r["g"] == g and r["cutoff"] == cutoff]
                maximum = max(coarse, key=lambda r:r["absolute_density_error"])
                centre = maximum["delta"]
                local_deltas = np.arange(max(-3., centre-.125), min(7., centre+.125)+.001, .0625)
                local_errors = list(coarse)
                for delta in local_deltas:
                    for mode in CUTOFFS:
                        point_key = key(number, g, delta, mode)
                        if point_key not in saved:
                            row, _, _ = ground_point(number, g, delta, mode)
                            saved[point_key] = row
                            refinements.append(row)
                    local_errors.append(comparison(saved[key(number,g,delta,cutoff)], saved[key(number,g,delta,"all")]))
                refined = max(local_errors, key=lambda r:r["absolute_density_error"])
                maxima.append({"N": number, "g": g, "cutoff": cutoff,
                    "parity": "even" if number%2==0 else "odd",
                    "coarse_max_absolute_density_error": maximum["absolute_density_error"],
                    "coarse_max_delta": centre,
                    "refined_sampled_max_absolute_density_error": refined["absolute_density_error"],
                    "refined_max_delta": refined["delta"],
                    "increase_after_local_refinement": refined["absolute_density_error"]-maximum["absolute_density_error"],
                    "tail_norm_bound": discarded_tail_bound(number, g, cutoff),
                    "coarse_max_absolute_order_error": max(r["absolute_alternating_order_squared_error"] for r in coarse),
                    "coarse_max_absolute_spin_order_error": max(r["absolute_staggered_spin_squared_error"] for r in coarse),
                    "coarse_max_absolute_centre_density_error": max(r["absolute_centre_density_error"] for r in coarse),
                    "coarse_max_centre_density_error_delta": max(coarse,key=lambda r:r["absolute_centre_density_error"])["delta"],
                    "coarse_max_absolute_gap_error": max(r["absolute_gap_error"] for r in coarse)})
    write_csv("ground_states_refinement.csv", refinements)
    write_csv("maximum_errors.csv", maxima)
    high = max([r for r in maxima if r["N"] == 12], key=lambda r:r["refined_sampled_max_absolute_density_error"])
    selected = {"N": 12, "g": high["g"], "delta": high["refined_max_delta"]}
    profile, profile_means = [], []
    for cutoff in CUTOFFS:
        row, sites, _ = ground_point(12, selected["g"], selected["delta"], cutoff)
        profile_means.append(row)
        for atom, density in enumerate(sites):
            profile.append({**selected, "cutoff": cutoff, "atom": atom, "density": float(density)})
    write_csv("boundary_profile.csv", profile)
    write_csv("boundary_profile_means.csv", profile_means)
    checks, tolerances = independent_checks(rows+refinements, saved, reference_check)
    make_figures(saved, differences, maxima, profile, selected)
    versions = {name: importlib.metadata.version(name) for name in ["numpy", "scipy", "matplotlib", "qiskit"]}
    versions["python"] = platform.python_version()
    summary = {"scope": "Finite open chains. No bulk phase boundary, novel phase, or quantum spin liquid is claimed.",
        "model": "H/Omega = (1/2) sum X_i - (Delta/Omega) sum n_i + sum g/r_ij^6 n_i n_j; hbar=Omega=a=1",
        "grid": {"N": SIZES, "g": STRENGTHS, "delta_min": -3., "delta_max": 7., "delta_points": 81, "spacing": .125,
                 "cutoffs": CUTOFFS, "coarse_ground_states": len(rows), "extra_refinement_states": len(refinements)},
        "refinement": "Additional spacing 0.0625 within +/-0.125 of each coarse maximum; only local sampled sensitivity, not a continuous maximum bound.",
        "maxima": maxima, "boundary_profile_setting": selected, "boundary_profile_means": profile_means,
        "smallest_coarse_gap": min(rows, key=lambda r:r["gap"]),
        "checks": checks, "tolerances": tolerances, "versions": versions,
        "elapsed_seconds": time.perf_counter()-start,
        "unchanged_baseline_sha256": hashlib.sha256((ROOT/"rydberg_model.py").read_bytes()).hexdigest()}
    (DATA / "extension_summary.json").write_text(json.dumps(summary, indent=2))
    lines = ["# Computation validation", "", "This record describes the actual saved run, not a claim about untested sizes or bulk phases.", "",
             f"Computed {len(rows):,} coarse-grid ground states and {len(refinements)} additional refinement states.",
             "The all-pair reference is the same finite open chain and laser parameters, with every interaction included.", "",
             "## Numerical checks", "", "| Check | Observed error or violation | Required below |", "|---|---:|---:|"]
    for name, tolerance in tolerances.items():
        lines.append(f"| {name.replace('_', ' ')} | {checks[name]:.6g} | {tolerance:g} |")
    lines += ["", "All stated tolerances passed. The residual includes both lowest eigenvectors; the absolute residual has energy units in Omega=1 units.",
              "Dense/sparse controls cover 12 points, including the minimum coarse-grid gap at N=8 for each strength and cutoff; exact points are in results/dense_sparse_controls.csv.",
              "The handwritten two-atom matrix and Qiskit Pauli construction independently check the bit-basis Hamiltonian.",
              "Free density and the centred-spin second moment are checked at five detunings for each of N=5,6,7 (15 cases). Positive-tail energy and occupation-versus-detuning monotonicity are checked throughout the entire coarse grid.", "",
              "The discarded-interaction norm B = g sum_(r>R) (N-r)/r^6 is verified against the diagonal matrix maximum. Both lowest-energy shifts lie in [0,B], and the absolute gap change is at most B, within tolerance. This bound does not bound density errors.", "",
              "## Previous-result check", "", json.dumps(checks["prior_csv_matching"], indent=2), "",
              "## Sampling limitations", "", summary["refinement"],
              "Maxima use the declared detuning interval; a maximum at its endpoint says nothing about what happens outside it. Centre density excludes only the endpoint atoms and is not a bulk limit.",
              "A nonzero alternating-order value is possible for finite disordered/product states, so no order threshold or phase boundary is assigned.", "",
              "## Provenance and environment", "",
              "rydberg_model.py retains the validated baseline model; its SHA-256 is recorded below.",
              "extension_model.py and run_extension.py implement this separate extension. Results were generated with classical numerical diagonalisation, not quantum hardware.",
              "These calculations and explanatory materials were developed with AI assistance. Numerical validation does not establish independent authorship or novelty.",
              f"Baseline SHA-256: {summary['unchanged_baseline_sha256']}", "", json.dumps(versions, indent=2), "",
              f"Elapsed calculation/validation/plot time: {summary['elapsed_seconds']:.1f} s.", ""]
    (Path(output_root) / "VALIDATION_COMPUTE.md").write_text("\n".join(lines))
    print(json.dumps({"checks":checks,"boundary_profile_setting":selected,"maxima":maxima,"elapsed_seconds":summary["elapsed_seconds"]}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    run()
