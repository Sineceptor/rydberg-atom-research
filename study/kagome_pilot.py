"""A small, open kagome-SITE patch: an alternate direction to discuss.

Run with: python kagome_pilot.py
Uses the existing rydberg_model.py without changing its conventions.
This is finite-array ground-state modelling, not a test for a spin liquid.
"""

from pathlib import Path
import csv
import json
import math
import platform

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rydberg_model import hamiltonian, interactions, occupations

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
SEED = 230926
STRENGTH = 8.0
OMEGA = 1.0
POINT_DELTA = 2.0


def kagome_sites(nx=2, ny=2):
    """Three sites per cell; nearest-neighbour separation is one.

    Atom numbering: cell y, then cell x, then basis site 0, 1, 2.
    The patch is open: no bonds wrap around the edges.
    """
    a1 = np.array([2.0, 0.0])
    a2 = np.array([1.0, np.sqrt(3.0)])
    basis = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2]])
    return np.array([x*a1 + y*a2 + b
                     for y in range(ny) for x in range(nx) for b in basis])


def nearest_bonds(positions):
    return [(i, j) for i in range(len(positions))
            for j in range(i+1, len(positions))
            if np.isclose(np.linalg.norm(positions[i]-positions[j]), 1.0)]


def coupling_models(positions):
    full = interactions(positions, STRENGTH)
    nn = np.zeros_like(full)
    for i, j in nearest_bonds(positions):
        nn[i, j] = nn[j, i] = STRENGTH
    return {"nearest_only": nn, "all_pairs": full}


def independent_tensor_hamiltonian(coupling, omega, delta):
    """Small dense reference using Kronecker products, not bit flips.

    Tensor factor order is atom N-1 down to atom 0 so it matches the
    established rightmost-bit convention.
    """
    n = len(coupling)
    ident = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    occupation = np.diag([0., 1.])

    def operator(factors):
        result = np.array([[1.]])
        for atom in reversed(range(n)):
            result = np.kron(result, factors.get(atom, ident))
        return result

    result = np.zeros((2**n, 2**n))
    for i in range(n):
        result += omega/2 * operator({i: x}) - delta*operator({i: occupation})
        for j in range(i+1, n):
            result += coupling[i, j] * operator({i: occupation, j: occupation})
    return result


def lowest_pair(matrix, seed):
    """Random starts avoid accidentally restricting a symmetry sector."""
    rng = np.random.default_rng(seed)
    energies, vectors = eigsh(matrix, k=2, which="SA", tol=1e-11,
                             ncv=32, maxiter=30000,
                             v0=rng.normal(size=matrix.shape[0]))
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    residuals = np.linalg.norm(matrix @ vectors - vectors*energies, axis=0)
    return energies, vectors, residuals


def observables(state, occupied, bonds):
    weights = np.abs(state)**2
    local = weights @ occupied
    joint = (occupied.T * weights) @ occupied
    connected = joint - np.outer(local, local)
    average_nn_connected = np.mean([connected[i, j] for i, j in bonds])
    return local, connected, float(average_nn_connected)


def write_csv(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def plot_geometry(ax, positions, bonds, values=None, title="", cmap="viridis",
                  vmin=0., vmax=1.):
    for i, j in bonds:
        ax.plot(positions[[i, j], 0], positions[[i, j], 1],
                color="#adb5bd", lw=1.7, zorder=1)
    scatter = ax.scatter(positions[:, 0], positions[:, 1], s=600,
                         c=values if values is not None else "#dcecf6",
                         cmap=cmap if values is not None else None,
                         vmin=vmin if values is not None else None,
                         vmax=vmax if values is not None else None,
                         edgecolors="#264653", zorder=2)
    for i, (x, y) in enumerate(positions):
        ax.text(x, y, str(i), ha="center", va="center", fontsize=10,
                color="black", zorder=3,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.8, pad=1))
    ax.set_title(title, fontsize=12, pad=12)
    ax.set_xlabel("x / nearest-neighbour spacing")
    ax.set_ylabel("y / nearest-neighbour spacing")
    ax.set_aspect("equal")
    ax.set_xlim(-.45, positions[:, 0].max()+.45)
    ax.set_ylim(-.45, positions[:, 1].max()+.45)
    return scatter


def run(output_root=None):
    """Write results only when called, under the caller's chosen output root."""
    global RESULTS, FIGURES
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "kagome")
    destination = Path(output_root)
    RESULTS = destination / "results"
    FIGURES = destination / "figures"
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    positions = kagome_sites()
    bonds = nearest_bonds(positions)
    couplings = coupling_models(positions)
    occupied = occupations(len(positions))
    bond_set = set(bonds)
    triangles = [(i, j, k) for i in range(12) for j in range(i+1, 12)
                 for k in range(j+1, 12)
                 if all(edge in bond_set for edge in [(i, j), (i, k), (j, k)])]
    distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
    assert len(np.unique(np.round(positions, 12), axis=0)) == 12
    assert np.isclose(distances[np.triu_indices(12, 1)].min(), 1)
    assert len(bonds) == 17
    assert len(triangles) == 5
    assert max(np.bincount(np.array(bonds).ravel(), minlength=12)) == 4
    write_csv(RESULTS/"kagome_coordinates.csv",
              [{"atom": i, "x": x, "y": y,
                "nn_degree": sum(i in edge for edge in bonds)}
               for i, (x, y) in enumerate(positions)])
    write_csv(RESULTS/"kagome_bonds.csv",
              [{"atom_i": i, "atom_j": j, "distance": float(distances[i, j]),
                "V_over_omega": float(couplings["all_pairs"][i, j])}
               for i, j in bonds])

    # Two independent reference comparisons, including a six-site open patch.
    checks = {}
    for label, patch in [("triangle", kagome_sites(1, 1)),
                         ("six_site_patch", kagome_sites(2, 1))]:
        small_coupling = interactions(patch, STRENGTH)
        reference = independent_tensor_hamiltonian(small_coupling, OMEGA, 1.3)
        matrix = hamiltonian(small_coupling, OMEGA, 1.3)
        energies, _, residuals = lowest_pair(matrix, SEED)
        expected = eigh(reference, eigvals_only=True)[:2]
        checks[label] = {
            "matrix_max_error": float(np.max(np.abs(reference-matrix.toarray()))),
            "lowest_two_energy_max_error": float(np.max(np.abs(energies-expected))),
            "max_residual": float(max(residuals))}
        assert checks[label]["matrix_max_error"] < 1e-12
        assert checks[label]["lowest_two_energy_max_error"] < 1e-10

    # Reflection about a line at 30 degrees swaps the two primitive vectors.
    reflection = np.array([[.5, np.sqrt(3)/2], [np.sqrt(3)/2, -.5]])
    transformed = positions @ reflection.T
    permutation = np.argmin(np.linalg.norm(transformed[:, None]-positions[None, :], axis=2), axis=1)
    geometric_error = float(np.max(np.abs(transformed-positions[permutation])))
    assert len(set(permutation.tolist())) == 12 and geometric_error < 1e-12
    coupling_symmetry_error = max(float(np.max(np.abs(v-v[np.ix_(permutation, permutation)])))
                                  for v in couplings.values())
    assert coupling_symmetry_error < 1e-12
    checks["reflection"] = {"atom_permutation": permutation.tolist(),
                            "coordinate_max_error": geometric_error,
                            "coupling_max_error": coupling_symmetry_error}

    deltas = np.linspace(-3., 7., 41)
    records, point_states = [], {}
    for model_index, (name, coupling) in enumerate(couplings.items()):
        for delta_index, delta in enumerate(deltas):
            matrix = hamiltonian(coupling, OMEGA, delta)
            energies, vectors, residuals = lowest_pair(matrix, SEED+100*model_index+delta_index)
            local, connected, nn_connected = observables(vectors[:, 0], occupied, bonds)
            records.append({"model": name, "atoms": 12, "delta_over_omega": float(delta),
                            "nearest_V_over_omega": STRENGTH,
                            "ground_energy_over_omega": float(energies[0]),
                            "first_excited_energy_over_omega": float(energies[1]),
                            "gap_over_omega": float(energies[1]-energies[0]),
                            "density": float(np.mean(local)),
                            "nn_connected_mean": nn_connected,
                            "ground_norm_error": float(abs(np.linalg.norm(vectors[:, 0])-1)),
                            "max_eigenpair_residual": float(max(residuals)),
                            "local_reflection_max_error": float(np.max(np.abs(local-local[permutation])))})
            if np.isclose(delta, POINT_DELTA):
                point_states[name] = (vectors[:, 0], local, connected)
        print(f"Completed {name}: 41 detunings", flush=True)
    write_csv(RESULTS/"kagome_scan.csv", records)
    grouped = {name: [r for r in records if r["model"] == name] for name in couplings}
    full, nn = grouped["all_pairs"], grouped["nearest_only"]
    energy_differences = np.array([a["ground_energy_over_omega"]-b["ground_energy_over_omega"]
                                   for a, b in zip(full, nn)])
    density_differences = np.array([a["density"]-b["density"] for a, b in zip(full, nn)])
    checks["scan"] = {
        "points_per_model": len(deltas), "models": 2,
        "max_eigenpair_residual": max(r["max_eigenpair_residual"] for r in records),
        "max_ground_norm_error": max(r["ground_norm_error"] for r in records),
        "max_local_reflection_error": max(r["local_reflection_max_error"] for r in records),
        "minimum_full_minus_nearest_ground_energy": float(min(energy_differences)),
        "minimum_gap": min(r["gap_over_omega"] for r in records)}
    assert checks["scan"]["max_eigenpair_residual"] < 1e-8
    assert checks["scan"]["max_ground_norm_error"] < 1e-10
    assert checks["scan"]["max_local_reflection_error"] < 1e-8
    assert min(energy_differences) >= -1e-9
    assert checks["scan"]["minimum_gap"] > 0
    for row in records:
        assert -1e-12 <= row["density"] <= 1+1e-12

    # Repeat the lowest-gap point from each model with another random start.
    repeats = []
    for name in couplings:
        row = min(grouped[name], key=lambda r: r["gap_over_omega"])
        energies, _, _ = lowest_pair(hamiltonian(couplings[name], OMEGA, row["delta_over_omega"]),
                                     SEED+10000)
        error = float(max(abs(energies - [row["ground_energy_over_omega"],
                                          row["first_excited_energy_over_omega"]])))
        assert error < 1e-8
        repeats.append({"model": name, "delta_over_omega": row["delta_over_omega"],
                        "lowest_two_energy_max_error": error})
    checks["independent_random_start_repeats"] = repeats

    local_rows, corr_rows, basis_rows, triangle_rows = [], [], [], []
    for name, (state, local, connected) in point_states.items():
        for i in range(12):
            local_rows.append({"model": name, "delta_over_omega": POINT_DELTA,
                               "atom": i, "occupation": float(local[i])})
            for j in range(12):
                corr_rows.append({"model": name, "delta_over_omega": POINT_DELTA,
                                  "atom_i": i, "atom_j": j,
                                  "connected_correlation": float(connected[i, j])})
        weights = np.abs(state)**2
        for triangle_index, triangle in enumerate(triangles):
            number_on_triangle = occupied[:, triangle].sum(axis=1)
            triangle_probabilities = [float(weights[number_on_triangle == count].sum())
                                      for count in range(4)]
            assert abs(sum(triangle_probabilities)-1.) < 1e-10
            triangle_rows.append({"model": name, "delta_over_omega": POINT_DELTA,
                                  "triangle": triangle_index,
                                  "atoms": ",".join(map(str, triangle)),
                                  **{f"probability_{count}_excited": probability
                                     for count, probability in enumerate(triangle_probabilities)},
                                  "probability_two_or_more": sum(triangle_probabilities[2:])})
        for rank, basis_index in enumerate(np.argsort(weights)[-10:][::-1], start=1):
            basis_rows.append({"model": name, "delta_over_omega": POINT_DELTA,
                               "rank": rank, "basis_index": int(basis_index),
                               "bits_in_order_atom0_to_atom11": "".join(str(x) for x in occupied[basis_index]),
                               "probability": float(weights[basis_index]),
                               "number_excited": int(sum(occupied[basis_index]))})
    write_csv(RESULTS/"kagome_point_occupations.csv", local_rows)
    write_csv(RESULTS/"kagome_point_correlations.csv", corr_rows)
    write_csv(RESULTS/"kagome_point_top_configurations.csv", basis_rows)
    write_csv(RESULTS/"kagome_triangle_occupations.csv", triangle_rows)
    checks["triangle_probabilities"] = {
        "triangles_per_model": len(triangles),
        "max_normalization_error": max(abs(sum(r[f"probability_{k}_excited"] for k in range(4))-1.)
                                        for r in triangle_rows),
        "probability_two_or_more_range": [min(r["probability_two_or_more"] for r in triangle_rows),
                                           max(r["probability_two_or_more"] for r in triangle_rows)]}

    # A no-interaction reference makes the blockade comparison interpretable.
    # The ground state is a product, so a triangle's count is binomial.
    free_p = (1 + POINT_DELTA / np.sqrt(POINT_DELTA**2 + OMEGA**2))/2
    free_counts = [float(math.comb(3, k)*free_p**k*(1-free_p)**(3-k)) for k in range(4)]
    free_h = independent_tensor_hamiltonian(np.zeros((3, 3)), OMEGA, POINT_DELTA)
    _, free_vectors = eigh(free_h)
    free_weights = np.abs(free_vectors[:, 0])**2
    free_occupied = occupations(3)
    numerical_counts = [float(free_weights[free_occupied.sum(axis=1) == k].sum()) for k in range(4)]
    free_error = float(max(abs(np.array(free_counts)-numerical_counts)))
    assert free_error < 1e-12
    checks["free_triangle_binomial_reference"] = {
        "delta_over_omega": POINT_DELTA, "V_over_omega": 0.,
        "single_atom_excitation_probability": float(free_p),
        "probabilities_0_to_3_excited": free_counts,
        "probability_two_or_more": sum(free_counts[2:]),
        "dense_tensor_probability_max_error": free_error}

    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 140})
    fig, ax = plt.subplots(figsize=(8, 5.7), constrained_layout=True)
    plot_geometry(ax, positions, bonds,
                  title="Open kagome-site patch: 12 atoms, 17 nearest-neighbour bonds")
    fig.savefig(FIGURES/"kagome_geometry.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True, constrained_layout=True)
    for name, color, label in [("nearest_only", "#a65424", "Nearest neighbours only"),
                               ("all_pairs", "#006b83", "All pairs: V = 8 / r⁶")]:
        rows = grouped[name]
        for ax, field in zip(axes, ["density", "gap_over_omega", "nn_connected_mean"]):
            ax.plot(deltas, [r[field] for r in rows], color=color, label=label, lw=2)
            ax.axvline(POINT_DELTA, color="gray", ls=":", lw=1)
            ax.grid(alpha=.18)
    axes[0].set_ylabel("Mean excitation density")
    axes[0].legend(fontsize=10)
    axes[1].set_ylabel("Energy gap / Ω")
    axes[2].set_ylabel("Mean connected correlation\non nearest-neighbour bonds")
    axes[2].set_xlabel("Detuning Δ / Ω")
    fig.suptitle("Twelve-site open patch, Ω = 1, nearest interaction V / Ω = 8\n"
                 "Finite-system ground states; dotted line selects Δ / Ω = 2", fontsize=12)
    fig.savefig(FIGURES/"kagome_ground_state_scan.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    for ax, name, title in zip(axes, couplings, ["Nearest neighbours only", "All pairs"]):
        sc = plot_geometry(ax, positions, bonds, point_states[name][1], title=title, vmax=.8)
    fig.colorbar(sc, ax=axes, shrink=.75, label="Local excitation probability ⟨nᵢ⟩")
    fig.suptitle("Open kagome-site patch at Δ / Ω = 2 and Vnearest / Ω = 8", fontsize=13)
    fig.savefig(FIGURES/"kagome_local_occupations.png", dpi=180)
    plt.close(fig)

    bound = max(np.max(np.abs(c)) for _, _, c in point_states.values())
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), constrained_layout=True)
    for ax, name, title in zip(axes, couplings, ["Nearest neighbours only", "All pairs"]):
        im = ax.imshow(point_states[name][2], vmin=-bound, vmax=bound, cmap="RdBu_r")
        ax.set_title(title)
        ax.set_xlabel("Column: atom j", labelpad=10)
        ax.set_ylabel("Row: atom i", labelpad=10)
        ax.set_xticks(range(12))
        ax.set_yticks(range(12))
    fig.colorbar(im, ax=axes, shrink=.8, label="Cᵢⱼ = ⟨nᵢnⱼ⟩ − ⟨nᵢ⟩⟨nⱼ⟩")
    fig.suptitle("Connected occupation correlations at Δ / Ω = 2\n"
                 "Each square is row i, column j. Diagonal entries are local variances.",
                 fontsize=12)
    # Tight bounding box plus a real white margin preserves the descender in j.
    fig.savefig(FIGURES/"kagome_connected_correlations.png", dpi=180,
                bbox_inches="tight", pad_inches=.20)
    plt.close(fig)

    index = int(np.argmax(abs(density_differences)))
    summary = {
        "purpose": "A bounded alternate finite kagome-site direction, not a quantum-spin-liquid test.",
        "geometry": {"atoms": 12, "cells": [2, 2], "boundary": "open",
                     "primitive_vectors": [[2., 0.], [1., float(np.sqrt(3.))]],
                     "basis": [[0., 0.], [1., 0.], [.5, float(np.sqrt(3.)/2)]],
                     "nearest_bonds": len(bonds), "nearest_triangles": triangles},
        "parameters": {"omega": OMEGA, "nearest_V_over_omega": STRENGTH,
                       "detuning_start": -3., "detuning_stop": 7., "detuning_step": .25,
                       "illustration_delta_over_omega": POINT_DELTA, "seed": SEED},
        "largest_sampled_density_difference": {
            "absolute": float(abs(density_differences[index])),
            "all_pairs_minus_nearest": float(density_differences[index]),
            "delta_over_omega": float(deltas[index]),
            "all_pairs_density": full[index]["density"],
            "nearest_only_density": nn[index]["density"]},
        "illustration": [r for r in records if np.isclose(r["delta_over_omega"], POINT_DELTA)],
        "illustration_redistribution": {
            "atom0_nearest_only": float(point_states["nearest_only"][1][0]),
            "atom0_all_pairs": float(point_states["all_pairs"][1][0]),
            "maximum_absolute_local_occupation_difference": float(np.max(abs(
                point_states["all_pairs"][1] - point_states["nearest_only"][1]))),
            "central_triangle_atoms": [5, 7, 9],
            "central_triangle_empty_probabilities": {
                r["model"]: r["probability_0_excited"]
                for r in triangle_rows if r["atoms"] == "5,7,9"}},
        "checks": checks,
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "all_checks_passed": True}
    (RESULTS/"kagome_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    run()
