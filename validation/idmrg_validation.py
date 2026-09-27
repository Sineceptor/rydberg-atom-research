"""Validation gate for the infinite pair-site iDMRG (rydberg_chain.infinite).

For each parameter point the iDMRG energy per atom is compared with the finite
pair-site DMRG bulk estimate [E(N2) - E(N1)]/(N2 - N1) from two open chains with
a pinned last atom (no domain wall), N1 = 100 and N2 = 200 atoms, using the same
inter-pair cutoff mmax and the same local basis.  A point passes when

  * |e_inf - e_fin| < 1e-6 per atom,
  * the two energy evaluations of the kept iDMRG state agree ('valid'),
  * <n_L - n_R> is equal on both pairs of the unit cell (spread < 1e-6), and
  * the tau = +1 and tau = -1 starts are I_s mirror images: energies within 1e-6,
    |tau_L + tau_R| < 1e-3 (at a critical point the order parameter is set by chi
    and differs between starts at the 1e-4 level).

Near a critical point (iDMRG xi > 20 a) the two-size difference carries an
O(v/(N1 N2)) Casimir bias.  Chains of 400 and 800 atoms are then added and
E(N) = e N + b + c/N + d/N^2 is fitted; such points pass on the fitted e and are
flagged ('energy_fit' in the gate).

Run from the repository root:

    OMP_NUM_THREADS=1 python validation/idmrg_validation.py --workers 2   # all points
    python validation/idmrg_validation.py --table                         # Markdown summary
    python validation/idmrg_validation.py --regate                        # re-evaluate saved rows

Results go to validation/results/idmrg_validation.json.  Points already present are
skipped unless --redo is given.  With --redo --reuse-finite only the infinite part
is recomputed; saved finite energies are reused when chain lengths, chi, mmax and
the local basis match.
"""

import argparse
import hashlib
import json
import logging
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "validation" / "results" / "idmrg_validation.json"

G0, DELTA0 = 11.39, 2.5
EPS_SCAN = [0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
OTHER_POINTS = [(4.83, 1.5, 0.2), (11.39, 4.0, 0.6), (11.39, 1.0, 0.3), (24.0, 3.5, 0.45)]
SETTINGS = {"chi": 64, "mmax": 8, "cell_pairs": 2, "starts": ["L", "R", "sym"], "n_pairs": [50, 100],
            "tol": 1e-6, "order_tol": 1e-6, "mirror_tau_tol": 1e-3,
            "xi_critical": 20.0, "extra_pairs": [200, 400], "extra_max_sweeps": 150}


def points():
    return [(G0, DELTA0, eps) for eps in EPS_SCAN] + OTHER_POINTS


def key(g, delta, eps):
    return f"g={g:g},delta={delta:g},eps={eps:g}"


def truncation_check(g, delta, eps, n_atoms=10):
    """ED of an open n-atom chain: lowest energy per atom, dressed basis minus full basis."""
    import numpy as np
    from rydberg_chain import ed, geometry, infinite as inf

    H = ed.build_hamiltonian(geometry.coupling_matrix(geometry.chain_positions(n_atoms, eps), g), 1.0, delta)
    _, W = inf.dressed_pair_basis(eps, g, delta)
    U = np.ones((1, 1))
    for _ in range(n_atoms // 2):
        U = np.kron(U, W)
    E_full = ed.lowest_states(H, k=1)[0][0]
    E_tr = np.linalg.eigvalsh(U.T @ (H @ U))[0]
    return {"n_atoms": n_atoms, "E_full": float(E_full), "E_truncated": float(E_tr),
            "diff_per_atom": float((E_tr - E_full) / n_atoms),
            "estimate_per_atom": inf.truncation_error_estimate(eps, g, delta)}


def evaluate_gate(row, settings):
    """(Re)compute the finite-size fit, the differences and the gate flags of a row in place."""
    from rydberg_chain import infinite as inf

    obs = row["infinite"]
    energies = {int(N): v["energy"] for N, v in row["finite"].items()}
    row.pop("casimir_coefficient", None)
    row["diff"] = row["e_inf"] - row["e_fin"]
    if len(energies) > 2:
        e, b, c, d = inf.bulk_energy_fit(energies)
        row["finite_fit"] = {"sizes": sorted(energies), "e": e, "b": b, "c": c, "d": d}
        row["e_fin_fit"] = e
        row["diff_fit"] = row["e_inf"] - e
    mirror = obs.get("mirror_check", {})
    gate = {"energy": abs(row["diff"]) < settings["tol"],
            "valid_state": obs["valid"],
            "order_consistent": obs["order_spread"] < settings["order_tol"],
            "mirror_consistent": (mirror.get("dE", 0.0) < settings["tol"]
                                  and mirror.get("dtau", 0.0) < settings["mirror_tau_tol"]),
            "idmrg_converged": obs["converged"],
            "idmrg_energy_converged": obs["energy_converged"],
            "finite_converged": all(v["converged"] for v in row["finite"].values())}
    if "diff_fit" in row:
        gate["energy_fit"] = abs(row["diff_fit"]) < settings["tol"]
    energy_ok = gate["energy"] or gate.get("energy_fit", False)
    gate["pass"] = energy_ok and gate["valid_state"] and gate["order_consistent"] and gate["mirror_consistent"]
    row["gate"] = {k: bool(v) for k, v in gate.items()}
    return row


def run_point(g, delta, eps, settings, cached=None):
    """iDMRG and finite bulk energies at one point; returns a JSON-ready dict.

    cached: a previous row whose finite energies are reused if their settings match.
    """
    logging.disable(logging.WARNING)
    warnings.filterwarnings("ignore")
    from rydberg_chain import effective, infinite as inf

    row = {"g": g, "delta": delta, "eps": eps, "V_s": g / (1 - eps) ** 6, "V_w": g / (1 + eps) ** 6}
    t0 = time.time()
    e_inf, psi, obs = inf.infinite_ground_state(eps, g, delta, chi=settings["chi"], mmax=settings["mmax"],
                                                cell_pairs=settings["cell_pairs"], starts=settings["starts"])
    row["t_infinite"] = time.time() - t0
    row["e_inf"] = e_inf
    row["e_inf_terms"] = inf.energy_from_terms(psi, eps, g, delta, mmax=settings["mmax"])
    row["tail_correction"] = inf.tail_energy_per_atom(eps, g, settings["mmax"], obs["n_left"], obs["n_right"])
    row["infinite"] = obs
    if obs["truncated"]:
        row["truncation_check"] = truncation_check(g, delta, eps)
    tf = effective.tfim_couplings(eps, g, delta)
    row["J1_over_h"] = abs(tf["J1"]) / tf["h"] if tf["h"] == tf["h"] else None

    finite_settings = {"n_pairs": list(settings["n_pairs"]), "chi": settings["chi"],
                       "mmax": settings["mmax"], "truncate": bool(obs["truncated"])}
    reuse = cached is not None and cached.get("finite_settings") == finite_settings
    info = {int(N): v for N, v in cached["finite"].items()} if reuse else {}
    row["t_finite"] = cached["t_finite"] if reuse else 0.0
    row["finite_reused"] = reuse
    sizes = list(settings["n_pairs"])
    if obs["xi"] > settings["xi_critical"]:
        sizes += list(settings["extra_pairs"])
    for n in sizes:
        if 2 * n in info:
            continue
        t0 = time.time()
        max_sweeps = settings["extra_max_sweeps"] if n in settings["extra_pairs"] else 80
        _, _, info[2 * n] = inf.pinned_chain_ground_state(eps, g, delta, n, chi=settings["chi"],
                                                          mmax=settings["mmax"], truncate=obs["truncated"],
                                                          max_sweeps=max_sweeps)
        row["t_finite"] += time.time() - t0
    N1, N2 = 2 * settings["n_pairs"][0], 2 * settings["n_pairs"][1]
    row["e_fin"] = (info[N2]["energy"] - info[N1]["energy"]) / (N2 - N1)
    row["finite"] = info
    row["finite_settings"] = finite_settings
    return evaluate_gate(row, settings)


def _to_json(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _to_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_json(v) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def summarize(results):
    rows = [results[k] for k in sorted(results, key=lambda k: (results[k]["g"], results[k]["delta"],
                                                                 results[k]["eps"]))]
    return {"n_points": len(rows), "n_pass": sum(r["gate"]["pass"] for r in rows),
            "n_pass_two_size": sum(r["gate"]["pass"] and r["gate"]["energy"] for r in rows),
            "max_abs_diff_two_size": max((abs(r["diff"]) for r in rows), default=None),
            "passed_on_fit": [key(r["g"], r["delta"], r["eps"]) for r in rows
                              if r["gate"]["pass"] and not r["gate"]["energy"]],
            "failed": [key(r["g"], r["delta"], r["eps"]) for r in rows if not r["gate"]["pass"]]}


def save(results):
    doc = {"description": __doc__.strip().splitlines()[0],
           "settings": SETTINGS,
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "summary": summarize(results),
           "points": results}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix(".tmp")
    tmp.write_text(json.dumps(_to_json(doc), indent=1))
    tmp.replace(OUT)


def markdown_table(results):
    """Summary table of the saved results (one row per point)."""
    lines = ["| g | Δ | ε | V_s | basis | e_inf (per atom) | e_fin − e_inf | \\|⟨n_L − n_R⟩\\| | density "
             "| ξ/a | S_weak | S_strong | χ | pass |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    rows = sorted(results.values(), key=lambda r: (r["g"] != G0 or r["delta"] != DELTA0, r["g"], r["delta"],
                                                   r["eps"]))
    for r in rows:
        o = r["infinite"]
        diff = f"{-r['diff']:+.1e}"
        if "diff_fit" in r:
            diff += f" (fit: {-r['diff_fit']:+.1e})"
        lines.append(
            f"| {r['g']:g} | {r['delta']:g} | {r['eps']:g} | {r['V_s']:.3g} | "
            f"{'dressed' if o['truncated'] else 'full'} | {r['e_inf']:.10f} | {diff} | "
            f"{o['order']:.5f} | {o['density']:.5f} | {o['xi']:.2f} | {o['entropy_weak'][0]:.4f} | "
            f"{o['entropy_strong'][0]:.4f} | {max(o['chi'])} | {'yes' if r['gate']['pass'] else 'NO'} |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--redo", action="store_true", help="recompute points already saved")
    parser.add_argument("--reuse-finite", action="store_true", help="with --redo: keep saved finite energies")
    parser.add_argument("--regate", action="store_true", help="re-evaluate the gate of the saved rows")
    parser.add_argument("--table", action="store_true", help="print a Markdown table of the saved results")
    parser.add_argument("--only", nargs="*", help="subset of eps values on the g=11.39, delta=2.5 line")
    args = parser.parse_args()
    results = json.loads(OUT.read_text())["points"] if OUT.exists() else {}
    if args.table:
        print(markdown_table(results))
        return
    if args.regate:
        for row in results.values():
            evaluate_gate(row, SETTINGS)
        save(results)
        print(json.dumps(summarize(results), indent=1))
        return
    todo = [p for p in points() if args.redo or key(*p) not in results]
    if args.only:
        wanted = {float(v) for v in args.only}
        todo = [p for p in todo if p[:2] == (G0, DELTA0) and p[2] in wanted]
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, "1")
    print(f"{len(todo)} points to run with {args.workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_point, *p, SETTINGS, results.get(key(*p)) if args.reuse_finite else None): p
                   for p in todo}
        for fut in as_completed(futures):
            row = fut.result()
            results[key(row["g"], row["delta"], row["eps"])] = row
            save(results)
            o = row["infinite"]
            print(f"{key(row['g'], row['delta'], row['eps'])}: e_inf={row['e_inf']:.10f} "
                  f"e_fin={row['e_fin']:.10f} diff={row['diff']:+.2e} order={o['order']:.5f} "
                  f"xi={o['xi']:.2f} pass={row['gate']['pass']} "
                  f"t={row['t_infinite']:.0f}+{row['t_finite']:.0f}s", flush=True)
    print(json.dumps(summarize(results), indent=1))


if __name__ == "__main__":
    main()
