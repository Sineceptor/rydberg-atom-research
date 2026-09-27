# Pilot calculations, 26–27 September 2026

These are exploratory scripts and raw outputs from the direction-finding pilots for the dimerized-chain study. They are kept so that every preliminary number quoted in [RESEARCH_PLAN.md](../../../docs/RESEARCH_PLAN.md) can be traced to the script that produced it.

**Status: exploratory, not package code.** The scripts are not held to the standard of `rydberg_chain/` and `tests/`. Several use restricted bases, and some hard-code parameters. Any number that enters the manuscript must be recomputed with the maintained package, or with a script promoted into `validation/`, and cross-checked against full exact diagonalization. The scripts were written with AI assistance during automated research reviews.

## Contents

| Script | What it computes | Outputs |
|---|---|---|
| `dimer_model.py`, `fast_ed.py` | First pilot. Atom-level ED/DMRG for open dimerized chains; infinite DMRG (atom level, then pair sites) | `iscan_g11.39_D2.5.*` (the **failed** infinite-DMRG scan: unconverged beyond ε ≈ 0.3) |
| `tci_pilot.py` | Open-chain ED, N = 12–16. Peak of d⟨n⟩/dΔ versus N at R_b/a = 1.5 | `tci_pilot.*`. No growth with N up to ε = 0.7. The scan range Δ ≤ 3 does not reach Δ_c for ε ≥ 0.6. |
| `pairk0.py`, `crit.py`, `sfac.py` | Periodic ring in the pair basis, momentum k = 0, resolved by inversion parity. Critical Δ_c(ε) from scaled-gap crossings, gap ratio R, central charge from E0(N) | `crit_eps*_g*_d*.json`, `log_*.txt`, `sf_*.txt` |
| `pairring.py`, `scan.py` | Pair-site ring ED with symmetry-sector penalties; Ising/tricritical test scan | `results.jsonl` |
| `pbc_ed.py`, `eps_scan.py` | Ring ED with at most one excitation per strong pair; Binder-type scans | `eps0*.txt` |
| `ced.py` | Constrained ring ED (no two excitations within two sites). Tests the Z4 → pair-density-wave splitting and the Z3 control | `out24_*.txt`, `out36_*.txt`, `z3_24_*.txt` |
| `dwbox_check.py`, `tfim_check.py`, `quick_checks.py` | Boundary-forced domain wall in even open chains; checks of the effective-model couplings | printed output |
| `diag.py`, `dbg.py`, `runall.sh` | Small diagnostics and a driver | — |

## Findings these pilots support (preliminary)

- **Z2 line (g = 11.39, R_b/a = 1.5), from ring ED:**
  - Δ_c/Ω = 1.0785, 1.326, 1.805, 3.149, 5.611 at ε = 0, 0.3, 0.45, 0.6, 0.7. Values agree to 1e-4 between size pairs.
  - Gap ratios R = 7.87–8.15 (Ising value 8).
  - Leading-order effective Ising estimate misses Δ_c by 17–63% for ε ≤ 0.6. A dressed denominator brings the error to about 0.5% at ε = 0.6–0.7.
- **Z4 regime (R_b/a = 3.3, Δ = 3), constrained ED, N = 16–32:**
  - The inversion-odd gap closes for ε ≤ 0.2 and saturates for ε ≥ 0.4.
  - The two lowest states for ε ≥ 0.4 have T2 = ±1 and are both inversion-even. This is consistent with a pair-density-wave phase.
- **Z3 control (R_b/a = 2.7, Δ = 5, N = 24):** the threefold degeneracy stays exact at ε = 0, 0.3 and 0.6.

Corrections to earlier notes:
- The open-chain interior values of |n_L − n_R| are contaminated by the boundary-forced domain wall.
- The disordered regime with one excitation shared per pair is a crossover, not a separate phase.
