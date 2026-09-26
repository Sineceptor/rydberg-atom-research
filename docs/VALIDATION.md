# Validation

The public package was prepared and checked on 26 September 2026. This update cleaned narrative and metadata while preserving numerical calculations. It did not rerun every historical scan.

## Public notebook execution

Both cleaned notebooks ran from fresh kernels with sequential execution counts and no cell errors:

| Notebook | Code cells | Embedded figures |
|---|---:|---:|
| Foundations | 16 | 7 |
| Finite Chain Study | 15 | 6 |

Hashes before and after execution matched for all **88 retained data, figure and reference files**: 61 result files, 25 figures and two regression-reference files. Notebook numerical outputs were regenerated; narrative, notebook metadata and filenames were cleaned. The source notebook checks include hand-built matrices, analytic limits, circuit comparisons, cutoff comparisons and selected parity/kagome points.

## Model and reference checks

The short `tools/check_current_model.py` control checks a hand-built two-atom Hamiltonian, noninteracting ground-state density, dense versus sparse energies/gap and both signs of detuning in split-circuit tests. It writes circuit test outputs only to a temporary folder.

The mandatory regression recomputes all 62 bundled historical N=8 reference rows and requires the maximum density/order/gap difference to be below 10⁻⁸. Its current result is saved in [reference_validation.json](reference_validation.json). The reference CSV and its digest metadata are unchanged.

The two Python edits made for this release replace explanatory string literals in future generated reports. A syntax-tree comparison confirmed that expressions and statements outside string contents were unchanged. `rydberg_model.py`, the Hamiltonian and numerical helper implementations retain their source bytes.

## Existing numerical evidence

`study/results/extension_summary.json` records 6,561 coarse states and 312 local-refinement states with their tolerances and checks. `parity_summary.json` records additional configurations, profiles, independent odd-chain repetitions and dense controls. `circuit_controls_summary.json` records the signed-detuning and convergence checks. Other result files record the foundation and optional kagome studies.

Those larger saved studies were not all repeated as part of public release preparation. Their retained data/figures were compared byte for byte with the corrected baseline. The environment was an already installed Python 3.12 scientific environment; no fresh installation on another platform is claimed.

## Integrity and interpretation

Run `python study/verify_delivery.py` to check the public study snapshot against its new `SHA256SUMS.txt`. Its manifest is specific to this cleaned release and excludes itself and future run folders. Saving an edited notebook will correctly change its digest. This verifies stored bytes, not scientific conclusions.

The checked publication tree contains scientific code, results, figures, cleaned notebooks and public documentation. It starts with a fresh Git history. Personal correspondence, source-document snapshots, private review material, local filesystem provenance and earlier private Git objects are not included. Public commit metadata uses the repository owner's GitHub no-reply identity. The figure files and notebook outputs were reviewed for private content.

Numerical accuracy does not establish a bulk phase, a quantum spin liquid, experimental performance, novelty, independent authorship or peer-reviewed status. See [Results and limitations](RESULTS_AND_LIMITS.md).
