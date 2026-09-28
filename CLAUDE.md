# Working notes for Claude Code sessions in this repository

Student research project: ground-state physics of Rydberg-atom chains, now focused on chains whose spacing alternates between a short gap s = a(1 − ε) and a long gap w = a(1 + ε). The student (Ethan) owns the work and re-derives the key equations himself. AI assistance is disclosed in the manuscript's Acknowledgements and in `docs/PROVENANCE.md`; keep it that way.

## Where things are

- `docs/RESEARCH_PLAN.md`: direction, pilot evidence, literature findings and the computation list P0–P9 with due dates. Read this first.
- `docs/NEXT_STEPS_PROMPTS.md`: ready-to-use prompts for the open tasks (orientation, P1–P9, citations, final audit).
- `docs/DIMERIZED_THEORY.md`: pseudo-spin theory, symmetry arguments (including the Ashkin–Teller labelling of the Z₄ patterns) and corrections log.
- `meeting_notes/`, `private/`, `private_notes/`: gitignored on purpose. Meeting briefs and personal notes live there and never enter the public repository.
- `paper/DRAFT.md`: the master manuscript. `paper/sections/` is superseded; do not edit it. `paper/DRAFT.docx` is an export of `DRAFT.md`.
- `rydberg_chain/`: maintained package (geometry, ED, observables, effective model, finite and infinite DMRG). Tests in `tests/`.
- `validation/`: benchmark and validation scripts with saved results in `validation/results/`.
- `research/pilots/2026-09-27/`: exploratory pilot scripts and outputs. Not package quality; numbers from here must be recomputed with the package before they are final.
- `study/`: the historical finite-chain study. Frozen and checksum-protected.

## Rules

- Never modify anything under `study/`. Check it with `python study/verify_delivery.py`.
- Every number in `paper/DRAFT.md` carries a hidden `<!-- source: ... -->` comment naming the file that produced it. Keep that convention. Mark anything unchecked with `[verify]`.
- Do not write into the shared Google Doc. Produce files (Markdown and .docx) that Ethan reviews and copies over himself.
- This repository is public. Keep personal and college-application material out of it.
- Claims to avoid are listed in `docs/RESEARCH_PLAN.md` ("Not to be claimed"). In particular: no "first", no tricritical point, no new phase for the shared-pair regime at R_b/a = 1.5, and the PDW is presented as the chain counterpart of the rung-resonant crystals of two-leg ladders.

## Commands

```bash
pip install -r requirements.txt
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python -m pytest   # 66 tests, about 1 minute
python paper/make_figures.py                                                  # figures from saved data only
python validation/benchmark_uniform.py --quick                                # quick benchmark (separate output file)
pandoc paper/DRAFT.md -f markdown+tex_math_dollars+raw_html -t docx --resource-path=. -o paper/DRAFT.docx
```

## Conventions

- Units ħ = Ω = a = 1 and g = (R_b/a)⁶. Bit i of a basis state is atom i. Pair k holds atoms 2k (L) and 2k + 1 (R).
- Pseudo-spin τ = n_L − n_R. In the theory notes J_m > 0 means ferromagnetic. `rydberg_chain/effective.py:tfim_couplings` uses the opposite sign (J < 0 is ferromagnetic).
- Even open chains in the Z₂ regime carry a boundary-forced domain wall. For large-N DMRG, use odd chains or pinned ends (last atom detuned by −2Δ).
