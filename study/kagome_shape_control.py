"""Optional appendix: same-N, different-shape open kagome-site patches.

This tests sensitivity to a changed finite patch, not bulk kagome physics.
Importing this module writes no calculation outputs. Call run(output_root).
"""

from pathlib import Path
import json
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kagome_pilot import (kagome_sites, nearest_bonds, coupling_models,
                          lowest_pair, observables, write_csv, plot_geometry)
from rydberg_model import hamiltonian, occupations


def run(output_root=None):
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "kagome-shapes")
    root = Path(output_root)
    results, figures = root / "results", root / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    shapes = {"2x2_patch": kagome_sites(2, 2), "1x4_strip": kagome_sites(1, 4)}
    detunings = [.25, 2., 4.]
    occupied = occupations(12)
    geometry, rows, local_rows, group_rows, coordinates = {}, [], [], [], []
    point_local = {}
    for shape_index, (name, positions) in enumerate(shapes.items()):
        bonds = nearest_bonds(positions)
        degree = np.bincount(np.array(bonds).ravel(), minlength=12)
        bond_set = set(bonds)
        triangles = [list(t) for t in itertools.combinations(range(12), 3)
                     if all(pair in bond_set for pair in itertools.combinations(t, 2))]
        assert len(np.unique(np.round(positions, 12), axis=0)) == 12
        geometry[name] = {"atoms": 12, "open_boundaries": True,
                          "bonds": len(bonds), "triangles": triangles,
                          "degree_counts": {str(d): int(sum(degree == d)) for d in np.unique(degree)},
                          "degree_by_atom": degree.tolist()}
        for i, (x, y) in enumerate(positions):
            coordinates.append({"shape": name, "atom": i, "x": float(x), "y": float(y),
                                "nearest_degree": int(degree[i])})
        for model_index, (model, coupling) in enumerate(coupling_models(positions).items()):
            for delta_index, delta in enumerate(detunings):
                matrix = hamiltonian(coupling, omega=1., delta=delta)
                energy, vectors, residuals = lowest_pair(matrix, 230943+100*shape_index+10*model_index+delta_index)
                local, _, connected = observables(vectors[:, 0], occupied, bonds)
                rows.append({"shape": name, "model": model, "atoms": 12,
                             "omega": 1., "nearest_interaction": 8., "delta": delta,
                             "density": float(np.mean(local)), "ground_energy": float(energy[0]),
                             "gap": float(energy[1]-energy[0]),
                             "mean_nearest_connected": connected,
                             "max_eigenpair_residual": float(max(residuals)),
                             "norm_error": float(abs(np.linalg.norm(vectors[:, 0])-1.))})
                for i, value in enumerate(local):
                    local_rows.append({"shape": name, "model": model, "delta": delta,
                                       "atom": i, "nearest_degree": int(degree[i]),
                                       "excitation_probability": float(value)})
                for d in np.unique(degree):
                    group_rows.append({"shape": name, "model": model, "delta": delta,
                                       "nearest_degree": int(d), "sites_in_group": int(sum(degree == d)),
                                       "mean_excitation_probability": float(np.mean(local[degree == d])),
                                       "minimum": float(min(local[degree == d])),
                                       "maximum": float(max(local[degree == d]))})
                if delta == 2.:
                    point_local[(name, model)] = local
    assert geometry["2x2_patch"]["bonds"] == 17
    assert len(geometry["2x2_patch"]["triangles"]) == 5
    assert geometry["1x4_strip"]["bonds"] == 15
    assert len(geometry["1x4_strip"]["triangles"]) == 4
    assert max(r["max_eigenpair_residual"] for r in rows) < 1e-8
    assert max(r["norm_error"] for r in rows) < 1e-10
    assert all(0 <= r["density"] <= 1 for r in rows)
    shifts = []
    for shape in shapes:
        for delta in detunings:
            chosen = {r["model"]: r for r in rows if r["shape"] == shape and r["delta"] == delta}
            energy_change = chosen["all_pairs"]["ground_energy"]-chosen["nearest_only"]["ground_energy"]
            assert energy_change >= -1e-9
            shifts.append({"shape": shape, "delta": delta,
                           "density_all_minus_nearest": chosen["all_pairs"]["density"]-chosen["nearest_only"]["density"],
                           "ground_energy_all_minus_nearest": energy_change})
    # Summaries compare coordination groups, without treating any group as bulk.
    group_shifts = []
    for shape in shapes:
        degrees = np.array(geometry[shape]["degree_by_atom"])
        local_difference = point_local[(shape, "all_pairs")]-point_local[(shape, "nearest_only")]
        for degree in np.unique(degrees):
            group_shifts.append({"shape": shape, "delta": 2., "nearest_degree": int(degree),
                                 "sites_in_group": int(sum(degrees == degree)),
                                 "mean_signed_local_change": float(np.mean(local_difference[degrees == degree])),
                                 "mean_absolute_local_change": float(np.mean(abs(local_difference[degrees == degree]))),
                                 "maximum_absolute_local_change": float(np.max(abs(local_difference[degrees == degree])))})
    write_csv(results/"kagome_shape_scan.csv", rows)
    write_csv(results/"kagome_shape_coordinates.csv", coordinates)
    write_csv(results/"kagome_shape_local.csv", local_rows)
    write_csv(results/"kagome_shape_degree_groups.csv", group_rows)
    write_csv(results/"kagome_shape_degree_changes.csv", group_shifts)
    write_csv(results/"kagome_shape_density_changes.csv", shifts)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 8), constrained_layout=True)
    for ax, (name, positions) in zip(axes, shapes.items()):
        title = "2 × 2 patch: 17 bonds, 5 triangles" if name == "2x2_patch" else "1 × 4 strip: 15 bonds, 4 triangles"
        sc = plot_geometry(ax, positions, nearest_bonds(positions), point_local[(name, "all_pairs")],
                           title=title, vmax=1.)
    fig.colorbar(sc, ax=axes, shrink=.65, label="Local excitation probability, all-pair model")
    fig.suptitle("Optional shape control: 12 sites in both open arrays\n"
                 "Ω = 1, nearest V / Ω = 8, Δ / Ω = 2; coordination and boundaries change", fontsize=12)
    fig.savefig(figures/"kagome_shape_occupations.png", dpi=180, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.6), constrained_layout=True)
    x = np.array([2, 3, 4])
    for offset, shape, color, label in [(-.16, "2x2_patch", "#006b83", "2 × 2 patch"),
                                       (.16, "1x4_strip", "#ab5729", "1 × 4 strip")]:
        values = {r["nearest_degree"]: r for r in group_shifts if r["shape"] == shape}
        actual_x = np.array(sorted(values))
        ax.bar(actual_x+offset, [100*values[d]["mean_absolute_local_change"] for d in actual_x],
               width=.3, color=color, label=label)
        for d in actual_x:
            ax.text(d+offset, 100*values[d]["mean_absolute_local_change"]+.12,
                    f'n={values[d]["sites_in_group"]}', ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xlabel("Number of nearest neighbours (coordination degree)")
    ax.set_ylabel("Mean absolute change in local probability\n(all pairs versus nearest only), percentage points")
    ax.set_title("At Δ / Ω = 2: where does the interaction-range change appear?\n"
                 "n labels the number of sites; the strip has no degree-4 sites", fontsize=11)
    ax.set_ylim(0, 1.22*max(100*r["mean_absolute_local_change"] for r in group_shifts))
    ax.legend()
    ax.grid(axis="y", alpha=.2)
    fig.savefig(figures/"kagome_shape_degree_changes.png", dpi=180, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)
    summary = {"purpose": "Optional same-N shape sensitivity check; not a bulk or geometry-isolated comparison.",
               "geometry": geometry, "detunings": detunings, "omega": 1., "nearest_V": 8.,
               "density_changes": shifts, "degree_group_changes_at_delta2": group_shifts,
               "max_eigenpair_residual": max(r["max_eigenpair_residual"] for r in rows),
               "max_norm_error": max(r["norm_error"] for r in rows),
               "all_checks_passed": True,
               "limitations": ["Open boundaries in both cases; no wrapping.",
                               "Same atom count does not fix bonds, triangles, aspect ratio or coordination.",
                               "No coordination group is assumed to represent the infinite-system bulk.",
                               "Three detunings do not define a phase diagram."]}
    (results/"kagome_shape_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    run()
