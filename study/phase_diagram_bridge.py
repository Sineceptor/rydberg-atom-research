"""Finite open-chain ground-state maps, not bulk phase boundaries.

Run after the August exercises. This is a new September extension.
Units: omega=1, hbar=1, lattice spacing=1. All pair interactions are included.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rydberg_model import chain, interactions, occupations, hamiltonian, entropy

ROOT = Path(__file__).resolve().parent
DATA, FIG = ROOT / "results", ROOT / "figures"


def ground_point(number, strength, delta, nearest_only=False):
    coupling = interactions(chain(number), strength)
    if nearest_only:
        for i in range(number):
            for j in range(number):
                if abs(i-j) > 1:
                    coupling[i,j] = 0
    H = hamiltonian(coupling, delta=delta)
    # A fixed starting vector makes the iterative calculation reproducible.
    # It contains both reflection symmetries, unlike a symmetric all-ones vector.
    rng = np.random.default_rng(218)
    energies, vectors = eigsh(H, k=2, which="SA", tol=1e-11,
                             v0=rng.normal(size=H.shape[0]))
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    state = vectors[:,0]
    bits = occupations(number)
    weights = abs(state)**2
    density = float(weights @ bits.sum(axis=1) / number)
    # Squaring before averaging detects either alternating pattern.
    # This occupation imbalance equals 1 for even ideal patterns, but not odd.
    # The spin definition below gives 1 for both ideal alternating patterns at any N.
    alternating = 2 * (bits @ ((-1.)**np.arange(number))) / number
    order_squared = float(weights @ alternating**2)
    spin = (2*bits-1) @ ((-1.)**np.arange(number)) / number
    spin_squared = float(weights @ spin**2)
    residual = max(np.linalg.norm(H @ vectors[:,i]-energies[i]*vectors[:,i]) for i in range(2))
    return {"N": number, "V_over_omega": strength, "delta_over_omega": delta,
            "density": density, "alternating_order_squared": order_squared, "staggered_spin_squared": spin_squared,
            "gap_over_omega": float(energies[1]-energies[0]),
            "entropy_bits": entropy(state, list(range(number//2))),
            "energy_over_omega": float(energies[0]), "residual": float(residual)}, state, H


def write_rows(filename, rows):
    with (DATA / filename).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def run(output_root=None):
    if output_root is None:
        from run_support import run_managed
        return run_managed(run, "bridge")
    global DATA, FIG
    DATA, FIG = Path(output_root) / "results", Path(output_root) / "figures"
    DATA.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False})
    detunings = np.linspace(-3, 7, 31)
    strengths = np.linspace(.5, 12, 24)
    rows = []
    for strength in strengths:
        for delta in detunings:
            row, _, _ = ground_point(8, strength, delta)
            rows.append(row)
    write_rows("finite_chain_map_N8.csv", rows)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    for ax, key, title in zip(axes,
            ["density", "staggered_spin_squared", "gap_over_omega"],
            ["Excitation density", "Staggered spin squared", "Lowest energy gap / Ω"]):
        values = np.array([r[key] for r in rows]).reshape(len(strengths), len(detunings))
        p = ax.pcolormesh(detunings, strengths, values, shading="auto", cmap="viridis")
        fig.colorbar(p, ax=ax)
        ax.set(title=title, xlabel="Detuning Δ / Ω", ylabel="Nearest interaction V₀ / Ω")
    fig.suptitle("Eight-atom open chain: finite-system ground-state map")
    fig.savefig(FIG / "09_finite_chain_map.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    slices = []
    dense_error = 0.
    for number in list(range(4,13)):
        for delta in detunings:
            row, state, H = ground_point(number, 8., delta)
            slices.append(row)
            if number == 4:
                exact = eigh(H.toarray(), eigvals_only=True)
                dense_error = max(dense_error, abs(exact[0]-row["energy_over_omega"]),
                                  abs(exact[1]-exact[0]-row["gap_over_omega"]))
    write_rows("size_comparison_V8.csv", slices)
    fig, axes = plt.subplots(1, 3, figsize=(12,3.7), constrained_layout=True)
    for number in list(range(4,13)):
        data = [r for r in slices if r["N"] == number]
        for ax, key in zip(axes, ["density", "staggered_spin_squared", "gap_over_omega"]):
            ax.plot(detunings, [r[key] for r in data], "--" if number % 2 else "-", label=f"N={number}")
    for ax, title in zip(axes,["Excitation density", "Staggered spin squared", "Lowest energy gap / Ω"]):
        ax.set(title=title, xlabel="Detuning Δ / Ω")
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Open chains at g = V₀ / Ω = 8: solid even N, dashed odd N")
    fig.savefig(FIG / "10_size_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    tails=[]
    for delta in detunings:
        for mode in [False,True]:
            row,_,_ = ground_point(8,8.,delta,nearest_only=mode)
            row["interactions"] = "nearest only" if mode else "all pairs"
            tails.append(row)
    write_rows("interaction_tail_comparison.csv", tails)
    # Analytic noninteracting ground-state density and zero cut entropy.
    free, _, _ = ground_point(6,0.,.7)
    analytic = (1+.7/np.sqrt(1+.7**2))/2
    # Hellmann-Feynman: dE0/dDelta = -sum_i <n_i>.
    centre, _, _ = ground_point(6,5.,1.)
    plus, _, _ = ground_point(6,5.,1.+1e-4)
    minus, _, _ = ground_point(6,5.,1.-1e-4)
    derivative = (plus["energy_over_omega"]-minus["energy_over_omega"])/2e-4
    sample = [ground_point(8,8.,x)[0] for x in [-3.,0.,3.,7.]]
    checks = {"maximum_eigenpair_residual": max(r["residual"] for r in rows+slices+tails),
              "dense_sparse_energy_and_gap_error": float(dense_error),
              "free_density_analytic_error": abs(free["density"]-analytic),
              "free_cut_entropy": abs(free["entropy_bits"]),
              "hellmann_feynman_error": abs(derivative+6*centre["density"])}
    assert checks["maximum_eigenpair_residual"] < 1e-7
    assert checks["dense_sparse_energy_and_gap_error"] < 1e-9
    assert checks["free_density_analytic_error"] < 1e-10
    assert checks["free_cut_entropy"] < 1e-9
    assert checks["hellmann_feynman_error"] < 1e-6
    result={"description":"Finite open-chain map; no assigned bulk phase boundaries or QSL claim",
            "grid":{"N":8,"detuning_points":31,"interaction_points":24,"slice_sizes":list(range(4,13))},
            "checks":checks,"sample_points":sample}
    (DATA/"phase_bridge_summary.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result

if __name__ == "__main__":
    run()
