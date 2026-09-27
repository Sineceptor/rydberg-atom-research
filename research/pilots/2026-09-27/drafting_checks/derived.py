"""Derived numbers for paper/sections/05_results_size.md (27 Sept). Reads saved files only."""
import csv, json, math
import numpy as np
R = "/home/user/rydberg-atom-research/"
def rd(f): return list(csv.DictReader(open(R + f)))

# (b) truncation errors
M = rd("study/results/maximum_errors.csv")
for r in M:
    if r["g"] == "16.0" and r["cutoff"] == "nearest1":
        print("NN g16", r["N"], r["parity"], "refined pp", round(100*float(r["refined_sampled_max_absolute_density_error"]), 4),
              "at delta", r["refined_max_delta"], "coarse pp", round(100*float(r["coarse_max_absolute_density_error"]), 4),
              "gap err", round(float(r["coarse_max_absolute_gap_error"]), 4), "tail bound", round(float(r["tail_norm_bound"]), 4))
nn2 = [(float(r["refined_sampled_max_absolute_density_error"]), r["N"], r["g"]) for r in M if r["cutoff"] == "nearest2"]
print("max nearest2 all:", max(nn2), " max nearest2 N>=8:", max(x for x in nn2 if int(x[1]) >= 8))
# odd larger than both even neighbours?
for g in ("2.0", "8.0", "16.0"):
    e = {int(r["N"]): float(r["refined_sampled_max_absolute_density_error"]) for r in M if r["g"] == g and r["cutoff"] == "nearest1"}
    print("g", g, "odd>even neighbours:", all(e[n] > e[n-1] and e[n] > e[n+1] for n in (5, 7, 9, 11)),
          " even monotone:", all(e[n] < e[n+2] for n in (4, 6, 8, 10)), " odd monotone:", all(e[n] < e[n+2] for n in (5, 7, 9)))
# ratio g=16 / g=8
e8 = {int(r["N"]): float(r["refined_sampled_max_absolute_density_error"]) for r in M if r["g"] == "8.0" and r["cutoff"] == "nearest1"}
e16 = {int(r["N"]): float(r["refined_sampled_max_absolute_density_error"]) for r in M if r["g"] == "16.0" and r["cutoff"] == "nearest1"}
print("ratio g16/g8:", {n: round(e16[n]/e8[n], 2) for n in range(4, 13)})

# slope test: does the NN-only density look like the full density at a shifted detuning?
G = rd("study/results/ground_states_coarse.csv")
for N in (11, 12):
    for g in ("8.0", "16.0"):
        a = sorted((float(r["delta"]), float(r["density"])) for r in G if r["N"] == str(N) and r["g"] == g and r["cutoff"] == "all")
        b = sorted((float(r["delta"]), float(r["density"])) for r in G if r["N"] == str(N) and r["g"] == g and r["cutoff"] == "nearest1")
        d = np.array([x for x, _ in a]); na = np.array([y for _, y in a]); nb = np.array([y for _, y in b])
        slope = np.gradient(na, d)
        i = np.argmax(np.abs(nb - na))
        print(f"N={N} g={g}: max|dn|={abs(nb-na)[i]:.4f} at delta={d[i]}, full slope there={slope[i]:.4f}, implied shift={abs(nb-na)[i]/slope[i]:.4f}, V(2a)={float(g)/64:.4f}, max slope={slope.max():.4f} at {d[np.argmax(slope)]}")

# (c) parity at g=8, delta=3, all pairs
F = rd("study/results/parity_fixed_setting_scan.csv")
gaps = {int(r["N"]): float(r["gap"]) for r in F if r["cutoff"] == "all"}
q2 = {int(r["N"]): float(r["staggered_spin_squared"]) for r in F if r["cutoff"] == "all"}
print("gaps", {n: round(v, 4) for n, v in gaps.items()})
print("Q2", {n: round(v, 4) for n, v in q2.items()})
def box(P):  # level spacing / (2t) for a wall with P-1 positions, hard walls
    return math.cos(math.pi/P) - math.cos(2*math.pi/P)
for N in (6, 8, 10, 12):
    P = N // 2
    print("uniform N", N, "positions", P-1, "t_fit", round(gaps[N]/(2*box(P)), 4), "gap*N^2", round(gaps[N]*N*N, 2))
S = json.load(open(R + "validation/results/scaling_benchmark.json"))
tf = []
for r in S["ed"]:
    N = r["n_atoms"]; P = N//2
    t = r["gap"]/(2*box(P)); tf.append(t)
    print("eps0.3 N", N, "gap", round(r["gap"], 4), "t_fit", round(t, 4))
for N in (50, 100, 200):
    P = N//2
    print("extrapolated wall gap N", N, [f"{2*t*box(P):.2e}" for t in (min(tf), max(tf))])
# leading-order transverse field h0 = (Omega^2/4)[1/Delta + 1/(V_s - Delta)]
def h0(g, eps, D): Vs = g/(1-eps)**6; return 0.25*(1/D + 1/(Vs - D))
print("h0 uniform g8 d3:", round(h0(8, 0, 3), 4), " h0 eps0.3 g11.39 d2.5:", round(h0(11.39, 0.3, 2.5), 4))
print("wall tail energy V(2)-V(3), g=8:", round(8/64 - 8/729, 4))
# standing-wave comparison N=8
B = rd("study/results/parity_bond_probabilities.csv")
for lev in ("0", "1"):
    p = {int(r["left_atom"]): float(r["probability_00"]) for r in B if r["N"] == "8" and r["cutoff"] == "all" and r["level"] == lev and r["g"] == "8.0" and r["delta"] == "3.0"}
    k = int(lev) + 1
    model = np.array([math.sin(k*math.pi*j/4)**2 for j in (1, 2, 3)]); model /= model.sum()
    print("N8 level", lev, "odd bonds 1,3,5:", [round(p[b], 4) for b in (1, 3, 5)], "sum", round(sum(p[b] for b in (1, 3, 5)), 4), "box model:", np.round(model, 3))

# classical (Omega = 0) energy of perfect patterns versus one-wall configurations, N = 8, all pairs
import sys
sys.path.insert(0, R)
from rydberg_chain import geometry
def ecl(bits, V, d):
    occ = [i for i, b in enumerate(bits) if b == "1"]
    return -d*len(occ) + sum(V[i, j] for a, i in enumerate(occ) for j in occ[a+1:])
for eps, g, d in ((0.0, 8.0, 3.0), (0.3, 11.39, 2.5)):
    V = geometry.coupling_matrix(geometry.chain_positions(8, eps), g)
    perfect = ecl("10101010", V, d)
    for s in ("10100101", "10010101"):
        print(f"eps={eps} g={g}: E({s}) - E(perfect) = {ecl(s, V, d) - perfect:.4f}")
# DMRG diagnostics
for r in S["dmrg_pair"]:
    if r["mmax"] == 8:
        print("DMRG N", r["n_atoms"], "sigma_E", round(math.sqrt(max(r["variance"], 0)), 4), "central-half |tau|", round(r["bulk_abs_pair_order"], 3),
              "s/sweep", round(r["seconds_per_sweep"], 1), "GB", round(r["peak_rss_GB"], 3))
print("V(2a) at g=387.4, 1291.5:", round(387.4/64, 2), round(1291.5/64, 2))
