"""Regenerate the publication figures from saved data.

Nothing is recomputed: every plotted value is read from ``study/results/*.csv``
or ``validation/results/scaling_benchmark.json``.  The only computed curves are
analytic reference lines, named in the legends.  The script writes PNG (300 dpi)
and PDF files to ``paper/figures/`` together with the manifest ``FIGURES.md``
(files, sources with SHA-256 prefixes, parameters and draft captions).

    python paper/make_figures.py              # every figure and the manifest
    python paper/make_figures.py fig04 fig06  # selected figures (manifest untouched)
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rydberg_chain import plotstyle as ps  # noqa: E402

STUDY = ROOT / "study" / "results"
BENCHMARK = ROOT / "validation" / "results" / "scaling_benchmark.json"
OUT = ROOT / "paper" / "figures"

FIGURES = {}


def figure(name):
    """Register a figure function under its output stem."""
    def register(function):
        FIGURES[name] = function
        return function
    return register


# ---------------------------------------------------------------- data helpers

def read_csv(path):
    """Rows of a CSV file with numeric fields converted to float."""
    rows = []
    with open(path, newline="") as stream:
        for raw in csv.DictReader(stream):
            row = {}
            for key, value in raw.items():
                try:
                    row[key] = float(value)
                except ValueError:
                    row[key] = value
            rows.append(row)
    return rows


def select(rows, **conditions):
    """Rows whose fields equal the given values."""
    return [r for r in rows if all(r[k] == v for k, v in conditions.items())]


def col(rows, key):
    return np.array([r[key] for r in rows], dtype=float)


def rel(path):
    path = Path(path).resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:12]


SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(x, digits=1):
    """Format x as 'm × 10ⁿ' with Unicode superscripts."""
    mantissa, exponent = f"{x:.{digits}e}".split("e")
    return f"{mantissa} × 10{str(int(exponent)).translate(SUPERSCRIPT)}"


def entry(stem, sources, parameters, caption, replaces=None, skipped=None):
    return {"stem": stem, "sources": [Path(s) for s in sources], "parameters": parameters,
            "caption": caption, "replaces": replaces, "skipped": skipped}


def rb(g):
    """Blockade radius in units of a, R_b/a = g^(1/6)."""
    return g ** (1.0 / 6.0)


PI_TICKS = ([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
            ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])


def parity_handles(colors_labels, marker_color=ps.INK):
    """Legend handles: one line per colour entry, then even/odd marker keys."""
    handles = [Line2D([], [], color=c, label=lab) for c, lab in colors_labels]
    handles += [Line2D([], [], color=marker_color, label=f"{p} N", **ps.PARITY_STYLE[p])
                for p in ("even", "odd")]
    return handles


# --------------------------------------------------------------------- figures

@figure("fig01_two_atom_blockade")
def two_atom_blockade(out):
    base_path, scan_path = STUDY / "qiskit_baseline.csv", STUDY / "blockade_scan.csv"
    base, scan = read_csv(base_path), read_csv(scan_path)
    fig, axes = ps.subplots(1, 3, height=2.3)
    styles = {0.0: (ps.BLUE, "-"), 20.0: (ps.VERMILLION, "--")}
    for V, (color, ls) in styles.items():
        rows = select(base, V_over_omega=V)
        t = col(rows, "omega_t")
        axes[0].plot(t, col(rows, "mean_excitations"), ls, color=color, label=f"$V/\\Omega = {V:g}$")
        axes[1].plot(t, col(rows, "entropy_bits"), ls, color=color, label=f"$V/\\Omega = {V:g}$")
    t = np.linspace(0, 2 * np.pi, 400)
    axes[0].plot(t, np.sin(np.sqrt(2) * t / 2) ** 2, ":", color=ps.INK, lw=1.1,
                 label=r"$\sin^2(\sqrt{2}\,\Omega t/2)$")
    axes[0].set(xlabel=r"$\Omega t$", ylabel=r"$\langle n_0 + n_1\rangle$", ylim=(0, 2.08))
    axes[1].set(xlabel=r"$\Omega t$", ylabel=r"$S_{\mathrm{atom}\,0}$ (bits)", ylim=(-0.02, 1.08))
    for ax in axes[:2]:
        ax.set_xticks(*PI_TICKS)
        ax.set_xlim(0, 2 * np.pi)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3)

    strong = [r for r in scan if r["V_over_omega"] > 0]
    V = col(strong, "V_over_omega")
    axes[2].loglog(V, col(strong, "time_mean_P11"), "o-", color=ps.BLUE, label="time average")
    axes[2].loglog(V, col(strong, "sampled_max_P11"), "s--", color=ps.VERMILLION,
                   markerfacecolor="white", label="largest sampled")
    # Adiabatic elimination of |11> for V >> Omega: P11(t) = (Omega/V)^2/2 sin^2(sqrt2 t/2).
    window = 0.5 * (0.5 - np.sin(2 * np.sqrt(2) * np.pi) / (4 * np.sqrt(2) * np.pi))
    Vline = np.geomspace(2, 60, 50)
    axes[2].loglog(Vline, window / Vline ** 2, "-", color=ps.MUTED, lw=0.8,
                   label=rf"${window:.4f}\,(\Omega/V)^2$")
    axes[2].loglog(Vline, 0.5 / Vline ** 2, "-.", color=ps.MUTED, lw=0.8,
                   label=r"$(\Omega/V)^2/2$")
    axes[2].set(xlabel=r"$V/\Omega$", ylabel=r"$P_{11}$", xlim=(0.4, 60), ylim=(1e-4, 4))
    lines = axes[2].get_lines()
    data_legend = axes[2].legend(handles=lines[:2], loc="upper right")
    axes[2].add_artist(data_legend)
    axes[2].legend(handles=lines[2:], loc="lower left")
    ps.panel_labels(axes)
    ps.save(fig, out / "fig01_two_atom_blockade")

    b20 = select(base, V_over_omega=20.0)
    t20 = col(b20, "omega_t")
    deviation = np.max(np.abs(col(b20, "mean_excitations") - np.sin(np.sqrt(2) * t20 / 2) ** 2))
    s20 = col(b20, "entropy_bits")
    k = int(np.argmax(s20))
    means = {r["V_over_omega"]: r["time_mean_P11"] for r in scan}
    caption = (
        "Two atoms driven resonantly (Δ = 0) from |00⟩ with interaction V = C₆/r⁶. "
        f"(a) Mean number of Rydberg excitations for V/Ω = 0 and 20; for V/Ω = 20 it follows the "
        f"perfect-blockade curve sin²(√2 Ωt/2) (dotted) to within {sci(deviation)}. "
        f"(b) Entanglement entropy of one atom, which stays zero for V = 0 and reaches "
        f"{s20[k]:.4f} bits at Ωt = {t20[k]:.2f} ≈ π/√2 for V/Ω = 20, close to the state (|01⟩ + |10⟩)/√2. "
        f"(c) Time average (circles) and largest sampled value (squares) of the double-excitation "
        f"probability P₁₁ over 0 ≤ Ωt ≤ 2π, with the V ≫ Ω asymptotes {window:.4f}(Ω/V)² and (Ω/V)²/2 (grey); "
        f"the average falls from the free-atom value 3/8 to {sci(means[5.0])} at V/Ω = 5 and "
        f"{sci(means[40.0])} at V/Ω = 40.")
    parameters = [
        "Two atoms, all-pair V = C₆/r⁶ entered as V/Ω; Δ = 0; ħ = Ω = 1; initial state |00⟩.",
        "(a, b): V/Ω ∈ {0, 20}; 161 equally spaced times 0 ≤ Ωt ≤ 2π; states from Qiskit exact "
        "evolution of the Pauli-string Hamiltonian (agrees with SciPy expm_multiply to 1.0 × 10⁻¹⁴, "
        "study/results/summary.json).",
        "(c): V/Ω ∈ {0.5, 1, 2, 5, 10, 20, 40} (V = 0 omitted on the log axis; its values are 0.375 "
        "and 1); 801 equally spaced times 0 ≤ Ωt ≤ 2π; time average by the trapezoid rule.",
        f"Reference lines: P₁₁ ≈ (Ω/V)²/2 · sin²(√2 Ωt/2) for V ≫ Ω, whose window average is "
        f"{window:.5f}(Ω/V)².",
    ]
    return entry("fig01_two_atom_blockade", [base_path, scan_path], parameters, caption,
                 replaces="study/figures/01_qiskit_baseline.png, 04_blockade_scan.png "
                          "(draft Figs. 2–3, paper/sections/03_results_coupling.md)")


@figure("fig02_interaction_switch_off")
def interaction_switch_off(out):
    cases = [("defect_intact.csv", "intact", ps.INK, "-"),
             ("defect_atom_2_isolated.csv", "atom 2 decoupled", ps.VERMILLION, "--"),
             ("defect_two_disconnected_halves.csv", "halves decoupled", ps.BLUE, ":")]
    data = {name: read_csv(STUDY / name) for name, *_ in cases}
    fig, axes = ps.subplots(1, 2, width=ps.DOUBLE_COLUMN * 0.72, height=2.35)
    for name, label, color, ls in cases:
        rows = data[name]
        t = col(rows, "omega_t")
        axes[0].plot(t, col(rows, "S_atom2_bits"), ls, color=color, label=label)
        axes[1].plot(t, col(rows, "S_012_bits"), ls, color=color, label=label)
    axes[0].set(ylabel=r"$S(\{2\}\,|\,\mathrm{rest})$ (bits)", ylim=(-0.03, 1.05))
    axes[1].set(ylabel=r"$S(\{0,1,2\}\,|\,\{3,4,5\})$ (bits)", ylim=(-0.03, 1.35))
    for ax in axes:
        ax.set(xlabel=r"$\Omega t$", xlim=(0, 4))
    axes[0].legend(loc="lower right", bbox_to_anchor=(1.0, 0.06))
    ps.panel_labels(axes)
    ps.save(fig, out / "fig02_interaction_switch_off")

    intact = data["defect_intact.csv"]
    iso = data["defect_atom_2_isolated.csv"]
    cut = data["defect_two_disconnected_halves.csv"]
    iso_rabi = np.max(np.abs(col(iso, "P_atom2") - np.sin(col(iso, "omega_t") / 2) ** 2))
    caption = (
        "Idealised interaction switch-off controls in a six-atom chain (spacing a, V₀/Ω = 5, all pairs "
        "V₀(a/r)⁶, Δ = 0, initial state |000000⟩); couplings are set to zero while every atom and its "
        "drive are kept, so these are not vacancies. "
        "(a) Entanglement entropy of atom 2 with the other five atoms and (b) half-chain entropy "
        "{0,1,2}|{3,4,5}, for the intact chain (solid), for atom 2 with all its interactions removed "
        "(dashed) and for all couplings across the central cut removed (dotted). "
        f"The decoupled atom stays unentangled (S ≤ {sci(np.max(np.abs(col(iso, 'S_atom2_bits'))))} bits) "
        f"and Rabi-oscillates as sin²(Ωt/2), and with the cut removed the half-chain entropy stays below "
        f"{sci(np.max(np.abs(col(cut, 'S_012_bits'))))} bits, against "
        f"{np.max(col(intact, 'S_012_bits')):.2f} bits at Ωt = 4 in the intact chain. "
        f"With only atom 2 decoupled the half-chain entropy still reaches "
        f"{np.max(col(iso, 'S_012_bits')):.3f} bits by Ωt = 4, through the tail V₀/2⁶ = 0.078Ω "
        "between atoms 1 and 3.")
    parameters = [
        "Six atoms at x = 0, 1, ..., 5 (units of a); V_ij = 5/|i − j|⁶ (V₀/Ω = 5, R_b/a = 1.31); Δ = 0; "
        "initial state |000000⟩; exact evolution (SciPy expm_multiply) at 161 times 0 ≤ Ωt ≤ 4.",
        "Atom 2 decoupled: V_2j = 0 for all j, drive on atom 2 kept "
        "(study/rydberg_model.interactions(..., isolated=2)).",
        "Central cut decoupled: V_ij = 0 for i ∈ {0,1,2}, j ∈ {3,4,5} (split=[0, 1, 2]).",
        f"Checks: decoupled-atom Rabi error {sci(iso_rabi)}; norm and energy conserved to < 5 × 10⁻¹⁵ "
        "(study/results/summary.json).",
    ]
    return entry("fig02_interaction_switch_off", [STUDY / c[0] for c in cases], parameters, caption,
                 replaces="study/figures/05_defects.png (draft Fig. 4)")


@figure("fig03_geometry_entanglement")
def geometry_entanglement(out):
    shapes = [("Chain", "chain", np.array([[0, 0], [1, 0], [2, 0], [3, 0]]), ps.BLUE, "-"),
              ("Square", "square", np.array([[0, 0], [1, 0], [0, 1], [1, 1]]), ps.VERMILLION, "--"),
              ("Rectangle", "rectangle", np.array([[0, 0], [1, 0], [0, 2], [1, 2]]), ps.GREEN, ":")]
    paths = {name: STUDY / f"geometry_{name.lower()}.csv" for name, *_ in shapes}
    data = {name: read_csv(p) for name, p in paths.items()}
    fig, axes = ps.subplots(1, 3, height=2.25,
                            gridspec_kw={"width_ratios": [1.15, 1, 1]})
    ax = axes[0]
    offsets = {"Chain": (0.0, 0.0), "Square": (4.3, 0.0), "Rectangle": (6.6, 0.0)}
    for name, label, pos, color, _ in shapes:
        p = pos + np.array(offsets[name])
        for i in range(4):
            for j in range(i + 1, 4):
                if abs(np.linalg.norm(pos[i] - pos[j]) - 1) < 1e-9:
                    ax.plot(*p[[i, j]].T, "-", color=ps.MUTED, lw=0.8, zorder=1)
        ax.plot(*p.T, "o", ms=6.5, color=color, mec="white", mew=1.0, zorder=2)
        for i, (x, y) in enumerate(p):
            ax.annotate(str(i), (x, y), xytext=(3.5, 3.5), textcoords="offset points",
                        fontsize=7, color=ps.INK)
        ax.text(p[:, 0].mean(), -0.75, label, ha="center", va="top", fontsize=7.5)
    ax.plot([0, 1], [2.35, 2.35], "-", color=ps.INK, lw=1.0)
    ax.text(0.5, 2.5, "$a$", ha="center", va="bottom", fontsize=8)
    ax.set(xlim=(-0.5, 8.0), ylim=(-1.5, 2.9), aspect="equal")
    ax.set_anchor("N")
    ax.axis("off")
    for name, label, _, color, ls in shapes:
        rows = data[name]
        t = col(rows, "omega_t")
        axes[1].plot(t, col(rows, "S_01_bits"), ls, color=color, label=label)
        axes[2].plot(t, col(rows, "S_02_bits"), ls, color=color, label=label)
    axes[1].set_title(r"$\{0,1\}\,|\,\{2,3\}$")
    axes[2].set_title(r"$\{0,2\}\,|\,\{1,3\}$")
    for a in axes[1:]:
        a.set(xlabel=r"$\Omega t$", ylabel="entropy (bits)", xlim=(0, 4), ylim=(-0.03, 2.08))
    axes[1].legend(loc="upper left")
    ps.panel_labels(axes)
    ps.save(fig, out / "fig03_geometry_entanglement")

    rect, chain, square = data["Rectangle"], data["Chain"], data["Square"]
    k = int(np.argmax(col(rect, "S_02_bits")))
    t2 = int(np.argmin(np.abs(col(chain, "omega_t") - 2.0)))
    sym = np.max(np.abs(col(square, "S_01_bits") - col(square, "S_02_bits")))
    caption = (
        "Entanglement in four-atom arrays with nearest spacing a, V₀/Ω = 5, all pairs V₀(a/r)⁶, Δ = 0 and "
        "initial state |0000⟩. "
        "(a) Chain, square of side a and a × 2a rectangle; grey lines join atoms at distance a. "
        "(b, c) Entropy of the bipartitions {0,1}|{2,3} and {0,2}|{1,3}; the two square curves coincide "
        "by symmetry. "
        f"Entanglement follows the short, strongly blockaded bonds: in the rectangle the partition "
        f"that cuts both short bonds reaches {rect[k]['S_02_bits']:.2f} bits (Ωt = {rect[k]['omega_t']:.2f}) "
        f"while the one that cuts neither stays below {np.ceil(1e3 * np.max(col(rect, 'S_01_bits'))) / 1e3:.3f} bits, "
        f"and in the chain the interleaved partition, which cuts all three nearest-neighbour bonds, "
        f"carries {chain[t2]['S_02_bits']:.2f} bits at Ωt = 2 against {chain[t2]['S_01_bits']:.2f} bits "
        "for the contiguous one.")
    parameters = [
        "Positions (units of a): chain (0,0),(1,0),(2,0),(3,0); square (0,0),(1,0),(0,1),(1,1); "
        "rectangle (0,0),(1,0),(0,2),(1,2).",
        "V_ij = 5/r_ij⁶ for all pairs (V₀/Ω = 5); Δ = 0; initial state |0000⟩; exact evolution at 161 "
        "times 0 ≤ Ωt ≤ 4.",
        "CSV columns: S_01_bits keeps atoms {0,1}; S_02_bits keeps atoms {0,2}.",
        "Check: the square's S_01 and S_02 columns are identical (symmetry)." if sym == 0 else
        f"Check: square symmetry |S_01 − S_02| ≤ {sym:.1g}.",
    ]
    return entry("fig03_geometry_entanglement", list(paths.values()), parameters, caption,
                 replaces="study/figures/06_geometries.png, 07_geometry_entanglement.png (draft Fig. 5)")


@figure("fig04_cutoff_error_vs_N")
def cutoff_error_vs_N(out):
    path = STUDY / "maximum_errors.csv"
    rows = read_csv(path)
    fig, axes = ps.subplots(1, 3, height=2.55, sharey=True)
    for ax, g in zip(axes, (2.0, 8.0, 16.0)):
        for cutoff in ("nearest1", "nearest2"):
            for parity in ("even", "odd"):
                sel = sorted(select(rows, g=g, cutoff=cutoff, parity=parity), key=lambda r: r["N"])
                ax.semilogy(col(sel, "N"), 100 * col(sel, "refined_sampled_max_absolute_density_error"),
                            color=ps.CUTOFF_COLORS[cutoff], **ps.PARITY_STYLE[parity])
        ax.set_title(f"$g = {g:g}$  ($R_b/a = {rb(g):.2f}$)")
        ax.set(xlabel="$N$", xticks=range(4, 13))
    axes[0].set_ylabel(r"largest $|\Delta\bar n|$ (% points)")
    handles = parity_handles([(ps.CUTOFF_COLORS[c], ps.CUTOFF_LABELS[c]) for c in ("nearest1", "nearest2")])
    fig.legend(handles=handles, loc="outside lower center", ncol=4)
    ps.panel_labels(axes)
    ps.save(fig, out / "fig04_cutoff_error_vs_N")

    nn = select(rows, cutoff="nearest1")
    worst = max(nn, key=lambda r: r["refined_sampled_max_absolute_density_error"])
    n12 = select(nn, N=12.0, g=worst["g"])[0]
    nn2 = max(select(rows, cutoff="nearest2"), key=lambda r: r["refined_sampled_max_absolute_density_error"])
    err = {(r["N"], r["g"]): r["refined_sampled_max_absolute_density_error"] for r in nn}
    odd_peaks = all(err[(n, g)] > err[(n - 1, g)] and err[(n, g)] > err[(n + 1, g)]
                    for n in (5, 7, 9, 11) for g in (2.0, 8.0, 16.0))
    growing = all(np.all(np.diff([err[(n, g)] for n in range(start, 13, 2)]) > 0)
                  for g in (8.0, 16.0) for start in (4, 5))
    assert odd_peaks and growing, "caption statements no longer match maximum_errors.csv"
    caption = (
        "Largest sampled change in the mean excitation density of a uniform open chain when the "
        "interaction g/|i − j|⁶ is truncated to nearest neighbours (vermillion) or to index separations "
        "≤ 2 (blue), measured against the all-pairs ground state; filled circles and solid lines are even N, "
        "open squares and dashed lines odd N. "
        "Each point is the maximum over 81 detunings −3 ≤ Δ/Ω ≤ 7, refined to spacing 0.0625 within ±0.125 "
        "of the coarse maximum. "
        f"Nearest-neighbour truncation shifts the density by up to "
        f"{100 * worst['refined_sampled_max_absolute_density_error']:.2f} percentage points "
        f"(N = {worst['N']:.0f}, g = {worst['g']:g}, against {100 * n12['refined_sampled_max_absolute_density_error']:.2f} "
        f"at N = 12); every odd chain exceeds both even neighbours, and at g = 8 and 16 the error still grows "
        "with N within each parity up to N = 12, so short even-only sequences understate it. "
        f"Keeping second neighbours reduces the error to at most "
        f"{100 * nn2['refined_sampled_max_absolute_density_error']:.2f} points; the dimerized-chain "
        "calculations therefore keep the long-range tail.")
    parameters = [
        "Uniform open chains, N = 4–12; couplings g/|i − j|⁶ with g = V₀/Ω ∈ {2, 8, 16} "
        "(R_b/a = 1.12, 1.41, 1.59).",
        "Models: nearest1 (|i − j| = 1 only), nearest2 (|i − j| ≤ 2), all (every pair, the reference).",
        "Plotted column: refined_sampled_max_absolute_density_error × 100 (percentage points), the maximum "
        "of |n̄_cut − n̄_all| over 81 detunings −3 ≤ Δ/Ω ≤ 7 (spacing 0.125) plus a local grid of spacing "
        "0.0625 within ±0.125 of the coarse maximum. A sampled maximum, not a bound.",
        "Ground states: dense diagonalization for dimension ≤ 64, otherwise SciPy eigsh (tol 10⁻¹¹); "
        "eigenpair residuals ≤ 5.5 × 10⁻¹⁰ (study/results/extension_summary.json).",
    ]
    return entry("fig04_cutoff_error_vs_N", [path], parameters, caption,
                 replaces="study/figures/02_max_error_by_size.png (draft Fig. 13)")


@figure("fig05_finite_chain_map_N8")
def finite_chain_map(out):
    path = STUDY / "finite_chain_map_N8.csv"
    rows = read_csv(path)
    detunings = np.unique(col(rows, "delta_over_omega"))
    strengths = np.unique(col(rows, "V_over_omega"))
    lookup = {(r["V_over_omega"], r["delta_over_omega"]): r for r in rows}

    def grid(key):
        return np.array([[lookup[(g, d)][key] for d in detunings] for g in strengths])

    fig, axes = ps.subplots(1, 3, height=2.3)
    panels = [("density", r"$\bar n$", None, (0, 1)),
              ("staggered_spin_squared", r"$\langle Q^2\rangle$", None, (0, None)),
              ("gap_over_omega", r"$(E_1 - E_0)/\Omega$", "log", None)]
    for ax, (key, label, scale, limits) in zip(axes, panels):
        values = grid(key)
        if scale == "log":
            norm = LogNorm(vmin=values.min(), vmax=values.max())
            mesh = ax.pcolormesh(detunings, strengths, values, shading="nearest", norm=norm,
                                 cmap=ps.SEQUENTIAL, rasterized=True)
        else:
            vmin, vmax = limits
            mesh = ax.pcolormesh(detunings, strengths, values, shading="nearest", vmin=vmin,
                                 vmax=vmax, cmap=ps.SEQUENTIAL, rasterized=True)
        bar = fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.06, aspect=22)
        if scale == "log":
            bar.set_ticks([0.1, 0.2, 0.5, 1, 2, 5], labels=["0.1", "0.2", "0.5", "1", "2", "5"])
            bar.minorticks_off()
        bar.set_label(label)
        bar.outline.set_linewidth(0.5)
        ax.grid(False)
        ax.set(xlabel=r"$\Delta/\Omega$", xticks=[-2, 0, 2, 4, 6])
    axes[0].set_ylabel(r"$g = V_0/\Omega$")
    ps.panel_labels(axes)
    ps.save(fig, out / "fig05_finite_chain_map_N8")

    q2 = max(rows, key=lambda r: r["staggered_spin_squared"])
    gap = min(rows, key=lambda r: r["gap_over_omega"])
    caption = (
        "Ground state of an eight-atom uniform open chain with couplings g/|i − j|⁶ between all pairs, on "
        f"{len(detunings)} detunings −3 ≤ Δ/Ω ≤ 7 and {len(strengths)} couplings 0.5 ≤ g ≤ 12 "
        f"({rb(strengths[0]):.2f} ≤ R_b/a ≤ {rb(strengths[-1]):.2f}). "
        "(a) Mean density, (b) ⟨Q²⟩ with Q = N⁻¹Σᵢ(−1)ⁱ(2nᵢ − 1), and (c) lowest gap E₁ − E₀ on a "
        "logarithmic colour scale. "
        f"Although either perfect alternating pattern has Q² = 1, ⟨Q²⟩ never exceeds "
        f"{q2['staggered_spin_squared']:.3f} (at g = {q2['V_over_omega']:g}, Δ/Ω = {q2['delta_over_omega']:.2f}), "
        "because the open ends of an even chain force a domain wall (a 00 bond) into the ground state; "
        f"the smallest gap is {gap['gap_over_omega']:.3f}Ω (g = {gap['V_over_omega']:g}, "
        f"Δ/Ω = {gap['delta_over_omega']:.2f}). "
        "One finite size with open ends: this map is not a phase diagram.")
    parameters = [
        "N = 8, uniform open chain, all pairs g/|i − j|⁶; ħ = Ω = a = 1.",
        f"Grid: Δ/Ω = linspace(−3, 7, 31) (spacing 1/3), g = linspace(0.5, 12, 24) (spacing 0.5); "
        f"{len(rows)} points; each cell centred on its grid point.",
        "Two lowest eigenpairs by SciPy eigsh (k = 2, tol 10⁻¹¹, fixed random start vector, seed 218); "
        "largest eigenpair residual 3.7 × 10⁻¹⁰ (study/results/phase_bridge_summary.json).",
        "Columns: density, staggered_spin_squared (= ⟨Q²⟩), gap_over_omega.",
    ]
    return entry("fig05_finite_chain_map_N8", [path], parameters, caption,
                 replaces="study/figures/09_finite_chain_map.png (draft Fig. 6)")


@figure("fig06_parity_fixed_setting")
def parity_fixed_setting(out):
    path = STUDY / "parity_fixed_setting_scan.csv"
    rows = read_csv(path)
    fig, axes = ps.subplots(1, 3, height=2.45)
    panels = [("density", r"$\bar n$", "linear"),
              ("staggered_spin_squared", r"$\langle Q^2\rangle$", "linear"),
              ("gap", r"$(E_1 - E_0)/\Omega$", "log")]
    for ax, (key, label, scale) in zip(axes, panels):
        for cutoff in ("nearest1", "all"):
            for parity, remainder in (("even", 0), ("odd", 1)):
                sel = sorted([r for r in select(rows, cutoff=cutoff) if r["N"] % 2 == remainder],
                             key=lambda r: r["N"])
                ax.plot(col(sel, "N"), col(sel, key), color=ps.CUTOFF_COLORS[cutoff],
                        **ps.PARITY_STYLE[parity])
        ax.set_yscale(scale)
        if scale == "log":
            ax.set_yticks([0.05, 0.1, 0.2, 0.5, 1, 2], labels=["0.05", "0.1", "0.2", "0.5", "1", "2"])
            ax.minorticks_off()
        ax.set(xlabel="$N$", ylabel=label, xticks=range(4, 13))
    axes[1].set_ylim(0, 1.02)
    handles = parity_handles([(ps.CUTOFF_COLORS[c], ps.CUTOFF_LABELS[c]) for c in ("nearest1", "all")])
    fig.legend(handles=handles, loc="outside lower center", ncol=4)
    ps.panel_labels(axes)
    ps.save(fig, out / "fig06_parity_fixed_setting")

    full = select(rows, cutoff="all")
    odd = [r for r in full if r["N"] % 2]
    even = {r["N"]: r for r in full if r["N"] % 2 == 0}
    q_odd, gap_odd = col(odd, "staggered_spin_squared"), col(odd, "gap")
    caption = (
        "Even and odd uniform open chains at g = 8 (R_b/a = 1.41) and Δ/Ω = 3 for N = 4–12, with "
        "nearest-neighbour interactions only (vermillion) or all pairs (black); filled circles are even N, "
        "open squares odd N. "
        "(a) Mean density, (b) ⟨Q²⟩ and (c) gap E₁ − E₀ on a logarithmic axis. "
        f"Odd chains, which fit a perfect alternating pattern with both ends excited, keep "
        f"⟨Q²⟩ = {q_odd.min():.3f}–{q_odd.max():.3f} and a gap of {gap_odd.min():.2f}–{gap_odd.max():.2f}Ω "
        f"(all pairs). In even chains the boundary-forced domain wall lowers ⟨Q²⟩ from "
        f"{even[4.0]['staggered_spin_squared']:.2f} to {even[12.0]['staggered_spin_squared']:.2f} and the gap "
        f"from {even[4.0]['gap']:.2f}Ω to {even[12.0]['gap']:.3f}Ω between N = 4 and 12.")
    nn2 = {r["N"]: r for r in select(rows, cutoff="nearest2")}
    nn2_dq = max(abs(nn2[n]["staggered_spin_squared"] - r["staggered_spin_squared"]) for n, r in
                 ((r["N"], r) for r in full))
    parameters = [
        "Uniform open chains N = 4–12; g = 8; Δ/Ω = 3; models nearest1 and all. The nearest2 rows are in "
        f"the CSV but not plotted; they differ from all pairs by at most {nn2_dq:.3f} in ⟨Q²⟩.",
        "Rows are the g = 8, Δ = 3 subset of study/results/ground_states_coarse.csv "
        "(dense diagonalization for dimension ≤ 64, SciPy eigsh tol 10⁻¹¹ otherwise).",
        "Columns: density, staggered_spin_squared (= ⟨Q²⟩), gap.",
        f"Check: E₀(N = 12, all pairs) = {even[12.0]['ground_energy']:.15f}, the saved uniform reference.",
    ]
    return entry("fig06_parity_fixed_setting", [path], parameters, caption,
                 replaces="study/figures/06_fixed_setting_parity.png (draft Fig. 14)")


@figure("fig07_solver_scaling")
def solver_scaling(out):
    stem = "fig07_solver_scaling"
    if not BENCHMARK.exists():
        return entry(stem, [], [], "", skipped=f"{rel(BENCHMARK)} not found")
    bench = json.loads(BENCHMARK.read_text())
    ed = sorted(bench["ed"], key=lambda r: r["n_atoms"])
    lob = sorted(bench["ed_lobpcg"], key=lambda r: r["n_atoms"])
    dm = sorted([r for r in bench["dmrg_pair"] if r["mmax"] == 8], key=lambda r: r["n_atoms"])
    full20 = next(r for r in bench["dmrg_pair"] if r["mmax"] is None)
    p = bench["parameters"]

    fig, axes = ps.subplots(1, 3, height=2.45)
    series = [("ED, Lanczos (eigsh)", ed, ps.BLUE, "-", "total_s"),
              ("ED, LOBPCG", lob, ps.GREEN, ":", "total_s"),
              (f"pair-site DMRG, χ = {p['chi']}", dm, ps.VERMILLION, "--", "seconds")]
    for label, rows, color, ls, tkey in series:
        N = np.array([r["n_atoms"] for r in rows])
        conv = np.array([r.get("converged", True) for r in rows])
        for ax, values in ((axes[0], [r[tkey] for r in rows]),
                           (axes[1], [r["peak_rss_GB"] for r in rows])):
            values = np.array(values)
            ax.plot(N, values, ls, color=color, lw=1.2)
            ax.plot(N[conv], values[conv], "o", color=color)
            ax.plot(N[~conv], values[~conv], "o", color=color, markerfacecolor="white")
    N = np.array([r["n_atoms"] for r in dm])
    conv = np.array([r["converged"] for r in dm])
    var = np.array([abs(r["variance"]) for r in dm])
    axes[2].plot(N, var, "--", color=ps.VERMILLION, lw=1.2)
    axes[2].plot(N[conv], var[conv], "o", color=ps.VERMILLION)
    axes[2].plot(N[~conv], var[~conv], "o", color=ps.VERMILLION, markerfacecolor="white")
    axes[0].set(yscale="log", ylabel="wall time (s)", ylim=(1, 3000))
    axes[1].set(ylabel="peak memory (GB)", ylim=(0, 1.0))
    axes[2].set(yscale="log", ylabel=r"$\langle H^2\rangle - \langle H\rangle^2$")
    for ax in axes:
        ax.set_xscale("log")
        ax.set(xlabel="$N$", xlim=(14, 250))
        ax.set_xticks([20, 50, 100, 200], labels=["20", "50", "100", "200"])
        ax.minorticks_off()
        ax.grid(True, which="major")
    handles = [Line2D([], [], color=c, linestyle=ls, marker="o", label=lab)
               for lab, _, c, ls, _ in series]
    handles += [Line2D([], [], color=ps.INK, marker="o", linestyle="none", label="converged"),
                Line2D([], [], color=ps.INK, marker="o", markerfacecolor="white", linestyle="none",
                       label=f"stopped at {p['dmrg_max_sweeps']} sweeps")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5)
    ps.panel_labels(axes)
    ps.save(fig, out / stem)

    ed20 = next(r for r in ed if r["n_atoms"] == 20)
    shift = next(r for r in dm if r["n_atoms"] == 20)["E0"] - ed20["E0"]
    speed = [e["solve_s"] / lo["solve_s"] for e, lo in zip(ed, lob)]
    agree = max(max(abs(r["E0_minus_eigsh"]), abs(r["E1_minus_eigsh"])) for r in lob)
    unconv = [r for r in dm if not r["converged"]]
    n50, big = unconv[0], unconv[1:]
    assert n50["n_atoms"] == 50 and n50["variance"] < 1e-8 and len(big) == 2
    E_err = 1e-12                      # max_E_err of the DMRG defaults used by the benchmark
    assert n50["last_sweep_dE"] > E_err, "caption says the N = 50 run missed the energy criterion"
    caption = (
        f"Cost of exact diagonalization (ED) and pair-site DMRG for an open dimerized chain at "
        f"ε = {p['eps']}, g = {p['g']} (R_b/a = {rb(p['g']):.2f}) and Δ/Ω = {p['delta']}. "
        "(a) Wall time, (b) peak resident memory and (c) DMRG energy variance; open symbols mark DMRG runs "
        f"stopped at the {p['dmrg_max_sweeps']}-sweep limit before meeting the convergence criteria. "
        f"At N = 20, DMRG reproduces the ED ground-state energy to {sci(abs(full20['E0_minus_ED']))} with all "
        f"couplings and to {sci(abs(shift))} with inter-pair couplings cut beyond 8 pairs, and LOBPCG matches "
        f"the two lowest eigsh energies to {sci(agree)} with a solve {min(speed):.1f} to {max(speed):.1f} times "
        "faster. "
        f"The N = {n50['n_atoms']} run is close to converged but met neither criterion (last-sweep energy change "
        f"{sci(n50['last_sweep_dE'])}, variance {sci(n50['variance'])}); the N = {big[0]['n_atoms']} and "
        f"{big[1]['n_atoms']} runs are not converged "
        f"ground states (variance {sci(big[0]['variance'])} and {sci(big[1]['variance'])}) and give only the "
        "cost per sweep, because the even open chain holds one delocalized domain wall that DMRG moves slowly. "
        "Times come from a shared 4-core machine (load average up to 9) and are upper bounds.")
    parameters = [
        f"Dimerized open chain, gaps 1 ∓ ε alternating (strong first), ε = {p['eps']}, g = {p['g']}, "
        f"Δ/Ω = {p['delta']}; even N (one boundary-forced domain wall).",
        "ED: full 2ᴺ basis, all pairs, rydberg_chain.ed.lowest_states (k = 2) with method 'auto' "
        "(eigsh) or 'lobpcg'; wall time = build + solve.",
        f"DMRG: rydberg_chain.dmrg.pair_chain_ground_state, χ = {p['chi']}, mmax = 8 (inter-pair couplings "
        f"up to 8 pairs), at most {p['dmrg_max_sweeps']} sweeps; the full-range N = 20 run is quoted in the "
        "caption only.",
        "Each case ran in its own subprocess with 2 BLAS threads; memory is that process's peak RSS; "
        f"machine: {bench['machine']['cpu_model']}, {bench['machine']['cpu_count']} cores, shared.",
        "Variance of the full-range N = 20 run is −1.0 × 10⁻¹² (rounding) and is not plotted.",
    ]
    return entry(stem, [BENCHMARK], parameters, caption,
                 replaces="new (Table M1 of paper/sections/02_methods.md gives the same numbers)")


# -------------------------------------------------------------------- manifest

def write_manifest(entries, path):
    lines = [
        "# Publication figures",
        "",
        "Generated by `python paper/make_figures.py` from saved data only; no physics is recomputed. "
        "Each figure is written as a 300 dpi PNG and a vector PDF. Style: `rydberg_chain/plotstyle.py`.",
        "",
        "Conventions: ħ = Ω = a = 1, |0⟩ = ground, |1⟩ = Rydberg, g = V₀/Ω = (R_b/a)⁶, "
        "Q = N⁻¹Σᵢ(−1)ⁱ(2nᵢ − 1). Entropies are in bits. See `docs/NOTATION.md`.",
        "",
        f"Matplotlib {matplotlib.__version__}, NumPy {np.__version__}. Source hashes are SHA-256 prefixes "
        "of the files as read.",
        "",
        "| Figure | Source data | Replaces |",
        "|---|---|---|",
    ]
    for e in entries:
        src = ", ".join(f"`{rel(s)}`" for s in e["sources"]) or "—"
        lines.append(f"| `{e['stem']}` | {src} | {e['replaces'] or (e['skipped'] and 'skipped') or ''} |")
    for e in entries:
        lines += ["", f"## {e['stem']}", ""]
        if e["skipped"]:
            lines += [f"Skipped: {e['skipped']}."]
            continue
        lines += [f"Files: `paper/figures/{e['stem']}.png`, `paper/figures/{e['stem']}.pdf`", ""]
        lines += ["Source data:", ""]
        lines += [f"- `{rel(s)}` (sha256 `{digest(s)}`)" for s in e["sources"]]
        lines += ["", "Parameters:", ""]
        lines += [f"- {p}" for p in e["parameters"]]
        lines += ["", "Draft caption:", "", f"> {e['caption']}"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("names", nargs="*", help="figure stems or prefixes (default: all)")
    parser.add_argument("--out", type=Path, default=OUT, help="output directory")
    args = parser.parse_args(argv)
    chosen = [n for n in FIGURES if not args.names or any(n.startswith(p) for p in args.names)]
    if not chosen:
        parser.error(f"no figure matches {args.names}; choose from {list(FIGURES)}")
    ps.use()
    args.out.mkdir(parents=True, exist_ok=True)
    entries = []
    for name in chosen:
        e = FIGURES[name](args.out)
        entries.append(e)
        print(f"{name}: {'skipped (' + e['skipped'] + ')' if e['skipped'] else 'written'}")
    if len(chosen) == len(FIGURES):
        write_manifest(entries, args.out / "FIGURES.md")
        print(f"manifest: {rel(args.out / 'FIGURES.md')}")


if __name__ == "__main__":
    main()
