"""Quick independent controls after moving the current package into this repo.

No delivered data or figures are changed. This is not a full scan or a test of
the planned alternating model. Circuit controls write only to a temporary folder.
"""
from pathlib import Path
import json
import sys
import tempfile
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
CURRENT = ROOT / 'study'
sys.path.insert(0, str(CURRENT))
from rydberg_model import hamiltonian, interactions, chain
from extension_model import ground_point
from circuit_controls import run as circuit_controls


def check():
    coupling = interactions(chain(2), 8.0)
    matrix = hamiltonian(coupling, omega=1.0, delta=1.3).toarray()
    reference = np.array([[0, .5, .5, 0], [.5, -1.3, 0, .5],
                          [.5, 0, -1.3, .5], [0, .5, .5, 8-2*1.3]])
    matrix_error = float(np.max(abs(matrix-reference)))
    assert matrix_error < 1e-12
    free, _, _ = ground_point(5, 0.0, 0.7)
    analytic = (1 + .7 / np.sqrt(1 + .7**2)) / 2
    free_error = float(abs(free['density'] - analytic))
    assert free_error < 1e-12
    sparse, _, h = ground_point(7, 8.0, 2.0, method='sparse')
    energy = eigh(h.toarray(), eigvals_only=True)
    dense_error = float(max(abs(sparse['ground_energy']-energy[0]),
                            abs(sparse['gap']-(energy[1]-energy[0]))))
    assert dense_error < 1e-9
    assert sparse['relative_residual'] < 1e-10
    with tempfile.TemporaryDirectory(prefix='rydberg-circuit-check-') as temporary:
        circuit_controls(Path(temporary))
    report = dict(two_atom_matrix_error=matrix_error,
                  free_ground_density_error=free_error,
                  dense_sparse_energy_gap_error=dense_error,
                  relative_eigenpair_residual=sparse['relative_residual'],
                  signed_detuning_circuit_controls='passed')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    check()
