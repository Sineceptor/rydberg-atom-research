"""Effect of the interaction cutoff on the van der Waals crossing of benchmark_uniform.py (part 3).

The vdW cross-check keeps couplings up to fourth neighbours (rmax = 4).  This script
recomputes the (33, 49) crossing of N^(1/8) m_c with rmax = 8 (V_9/V_1 = 2.6e-6) on
the same detuning grid and compares it with the stored rmax = 4 scans.

    python validation/vdw_range_check.py     # about 3 min on 2 cores

Writes validation/results/vdw_range_check.json.
"""

import importlib.util
import json
import multiprocessing
import os
import platform
import time

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("benchmark_uniform", os.path.join(HERE, "benchmark_uniform.py"))
bench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bench)

SIZES = (33, 49)
RMAX = 8
OUT = os.path.join(HERE, "results", "vdw_range_check.json")


def _point(spec):
    n, delta = spec
    return spec, bench.blockade_point(n, delta, bench.VDW_G, rmax=RMAX)


def main():
    t0 = time.time()
    stored = json.load(open(bench.RESULTS))["vdw_check"]
    deltas = stored["deltas"]
    with multiprocessing.get_context("fork").Pool(2) as pool:
        done = dict(pool.map(_point, [(n, d) for n in SIZES for d in deltas]))
    scans = {n: [done[(n, d)] for d in deltas] for n in SIZES}
    _, pairs = bench.analyse_crossings(scans, list(SIZES), deltas)
    ref = next(p for p in stored["pairs"] if tuple(p[:2]) == SIZES)
    result = {
        "script": "validation/vdw_range_check.py", "python": platform.python_version(),
        "g": bench.VDW_G, "sizes": list(SIZES), "deltas": deltas,
        "rmax_stored": stored["rmax"], "rmax": RMAX,
        "crossing_stored": ref[2], "crossing": pairs[0][2], "shift": pairs[0][2] - ref[2],
        "all_converged": all(r["converged"] for scan in scans.values() for r in scan),
        "max_sweeps_used": max(r["n_sweeps"] for scan in scans.values() for r in scan),
        "scans": {str(n): [{k: r[k] for k in ("Delta", "E", "m_c", "n_sweeps", "converged", "max_chi")}
                           for r in scans[n]] for n in SIZES},
        "seconds": time.time() - t0,
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=1)
    print(f"crossing {SIZES}: rmax={stored['rmax']} {ref[2]:.5f}, rmax={RMAX} {pairs[0][2]:.5f}, "
          f"shift {result['shift']:+.1e}; converged={result['all_converged']} -> {OUT}")


if __name__ == "__main__":
    main()
