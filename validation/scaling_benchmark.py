"""Timing and memory benchmark for rydberg_chain (ED builder + eigsh, pair-site finite DMRG).

Each case runs in a fresh subprocess so that ru_maxrss is the peak of that case alone.
Usage: python scaling_benchmark.py            (driver; writes the JSON)
       python scaling_benchmark.py --lobpcg   (adds lowest_states(method='lobpcg') timings)
       python scaling_benchmark.py --dmrg     (re-runs only the DMRG cases)
       python scaling_benchmark.py CASE_JSON  (worker; prints one JSON line)
"""
import json
import os
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "results" / "scaling_benchmark.json"
THREADS = {"OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "MKL_NUM_THREADS": "2"}
EPS, G, DELTA, CHI = 0.3, 11.39, 2.5, 64
MAX_SWEEPS = 30
NOTES = ("Each case runs in a separate subprocess; peak_rss_GB is that process's ru_maxrss. "
         "The 4-core machine was shared with other heavy jobs (load average 4-9), so wall times "
         "are upper bounds and vary by ~2x between runs. "
         "ED: build_hamiltonian (all-range couplings) + lowest_states(k=2); 'ed' uses the default "
         "method (eigsh, which='SA', tol=1e-12, above dim 1024), 'ed_lobpcg' uses method='lobpcg' "
         "(Jacobi-preconditioned block LOBPCG, tol=1e-12). "
         "DMRG: pair_chain_ground_state with package defaults (two-site, chi=64, mixer for up to 10 "
         "sweeps, max_E_err=1e-12, max_S_err=1e-8) but max_sweeps=30, started from "
         "default_product_state (one central 00 wall); seconds include MPO construction and the "
         "variance evaluation. 'converged' is False when the energy or the entropy criterion was not met "
         "within max_sweeps; the slow mode is the delocalized domain wall of the even open chain. "
         "Such runs are not converged ground states unless last_sweep_dE and the variance are "
         "small; they always give the cost per sweep.")


def peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # Linux: kB


def worker(case):
    sys.path.insert(0, str(ROOT))
    import numpy as np
    if case["kind"] == "ed":
        from rydberg_chain import ed, geometry, observables as ob
        n = case["n_atoms"]
        t0 = time.perf_counter()
        H = ed.build_hamiltonian(geometry.coupling_matrix(geometry.chain_positions(n, EPS), G), 1.0, DELTA)
        t1 = time.perf_counter()
        method = case.get("method", "auto")
        E, v, res = ed.lowest_states(H, k=2, method=method)
        t2 = time.perf_counter()
        return {"n_atoms": n, "method": method, "dim": 2 ** n, "nnz": int(H.nnz), "build_s": t1 - t0, "solve_s": t2 - t1,
                "total_s": t2 - t0, "E0": float(E[0]), "E1": float(E[1]), "gap": float(E[1] - E[0]),
                "residual_max": float(res.max()), "density": ob.mean_density(v[:, 0], n),
                "peak_rss_GB": peak_gb()}
    if case["kind"] == "study_build":
        sys.path.insert(0, str(ROOT / "study"))
        sys.dont_write_bytecode = True
        import rydberg_model
        from rydberg_chain import geometry
        n = case["n_atoms"]
        V = geometry.coupling_matrix(geometry.chain_positions(n, EPS), G)
        t0 = time.perf_counter()
        H = rydberg_model.hamiltonian(V, 1.0, DELTA)
        return {"n_atoms": n, "build_s": time.perf_counter() - t0, "nnz": int(H.nnz), "peak_rss_GB": peak_gb()}
    if case["kind"] == "dmrg":
        import logging
        logging.disable(logging.WARNING)
        from rydberg_chain import dmrg
        n = case["n_atoms"]
        P = n // 2
        t0 = time.perf_counter()
        M = dmrg.pair_chain_model(EPS, G, DELTA, P, mmax=case["mmax"])
        mpo_chi = int(max(M.H_MPO.chi))
        E, psi, obs = dmrg.pair_chain_ground_state(EPS, G, DELTA, P, chi=CHI, mmax=case["mmax"],
                                                   variance=True, max_sweeps=MAX_SWEEPS)
        t1 = time.perf_counter()
        order = np.abs(obs["pair_order"])
        bulk = order[P // 4: P - P // 4] if P >= 4 else order
        return {"n_atoms": n, "n_pairs": P, "mmax": case["mmax"], "chi": CHI, "mpo_bond_dim": mpo_chi,
                "seconds": t1 - t0, "seconds_per_sweep": (t1 - t0) / obs["n_sweeps"],
                "E0": E, "E_per_atom": E / n, "variance": obs["variance"],
                "n_sweeps": obs["n_sweeps"], "max_sweeps": MAX_SWEEPS, "converged": obs["converged"],
                "last_sweep_dE": obs["last_sweep_dE"], "last_sweep_dS": obs["last_sweep_dS"],
                "max_trunc_err": obs["max_trunc_err"], "final_trunc_err": obs["final_trunc_err"],
                "max_chi": obs["max_chi"],
                "half_chain_entropy_bits": obs["half_chain_entropy"], "density": obs["density"],
                "bulk_abs_pair_order": float(bulk.mean()), "peak_rss_GB": peak_gb()}
    raise ValueError(case)


def run_case(case, timeout):
    env = {**os.environ, **THREADS}
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, __file__, json.dumps(case)], env=env, capture_output=True,
                          text=True, timeout=timeout)
    if proc.returncode != 0:
        return {**case, "error": proc.stderr[-2000:], "wall_s": time.perf_counter() - t0}
    result = json.loads(proc.stdout.strip().splitlines()[-1])
    print(json.dumps(result), flush=True)
    return result


def machine_info():
    import numpy, scipy, tenpy
    info = {"platform": platform.platform(), "python": platform.python_version(),
            "processor": platform.processor() or platform.machine(), "cpu_count": os.cpu_count(),
            "numpy": numpy.__version__, "scipy": scipy.__version__, "tenpy": tenpy.__version__,
            "threads": THREADS, "load_average_at_start": os.getloadavg()}
    try:
        mem = open("/proc/meminfo").read().split("\n")[0]
        info["mem_total"] = mem.split(":")[1].strip()
        model = [l for l in open("/proc/cpuinfo") if l.startswith("model name")]
        info["cpu_model"] = model[0].split(":")[1].strip() if model else None
    except OSError:
        pass
    return info


def main():
    sys.path.insert(0, str(ROOT))
    out = {"created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "machine": machine_info(),
           "parameters": {"eps": EPS, "g": G, "delta": DELTA, "omega": 1.0, "chi": CHI},
           "notes": NOTES,
           "ed": [], "study_builder": [], "dmrg_pair": []}
    out["study_builder"].append(run_case({"kind": "study_build", "n_atoms": 16}, 600))
    for n in (16, 18, 20):
        out["ed"].append(run_case({"kind": "ed", "n_atoms": n}, 1800))
    for n, mmax in DMRG_CASES:
        out["dmrg_pair"].append(run_case({"kind": "dmrg", "n_atoms": n, "mmax": mmax}, 3600))
        OUT.write_text(json.dumps(out, indent=1))
    ed20 = next((r for r in out["ed"] if r.get("n_atoms") == 20 and "E0" in r), None)
    if ed20:
        for r in out["dmrg_pair"]:
            if r.get("n_atoms") == 20 and "E0" in r:
                r["E0_minus_ED"] = r["E0"] - ed20["E0"]
    out["script_source"] = Path(__file__).read_text()
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT)


DMRG_CASES = ((20, None), (20, 8), (50, 8), (100, 8), (200, 8))


def rerun_dmrg():
    """Re-run the DMRG cases with the current package defaults; keep the saved ED results."""
    out = json.loads(OUT.read_text())
    out["dmrg_pair"] = []
    out["notes"] = NOTES
    out["parameters"]["dmrg_max_sweeps"] = MAX_SWEEPS
    out["machine"]["load_average_at_dmrg_pass"] = os.getloadavg()
    for n, mmax in DMRG_CASES:
        out["dmrg_pair"].append(run_case({"kind": "dmrg", "n_atoms": n, "mmax": mmax}, 7200))
        _finish(out)


def _finish(out):
    ed20 = next((r for r in out["ed"] if r.get("n_atoms") == 20 and "E0" in r), None)
    for r in out["dmrg_pair"]:
        if ed20 and r.get("n_atoms") == 20 and "E0" in r:
            r["E0_minus_ED"] = r["E0"] - ed20["E0"]
    out["script_source"] = Path(__file__).read_text()
    OUT.write_text(json.dumps(out, indent=1))


def add_lobpcg():
    """Second pass: time lowest_states(method='lobpcg') and merge into the saved JSON."""
    out = json.loads(OUT.read_text())
    for r in out["ed"]:
        r.setdefault("method", "auto")
    out["notes"] = NOTES
    out["ed_lobpcg"] = [run_case({"kind": "ed", "n_atoms": n, "method": "lobpcg"}, 1800) for n in (16, 18, 20)]
    out["machine"]["load_average_at_lobpcg_pass"] = os.getloadavg()
    for r in out["ed_lobpcg"]:
        ref = next((e for e in out["ed"] if e.get("n_atoms") == r.get("n_atoms")), None)
        if ref and "E0" in r:
            r["E0_minus_eigsh"] = r["E0"] - ref["E0"]
            r["E1_minus_eigsh"] = r["E1"] - ref["E1"]
    out["script_source"] = Path(__file__).read_text()
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--lobpcg":
        add_lobpcg()
    elif len(sys.argv) > 1 and sys.argv[1] == "--dmrg":
        rerun_dmrg()
    elif len(sys.argv) > 1:
        print(json.dumps(worker(json.loads(sys.argv[1]))))
    else:
        main()
