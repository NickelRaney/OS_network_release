"""Export the audited MATLAB histogram fixture (maintainer tool, SciPy needed).

Usage: python tools/extract_onset_reference.py LEGACY_REFERENCE_MAT
Run from the release root. The original experiment remains read only.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat

root = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1])
data = loadmat(source, simplify_cells=True)
roi = loadmat(root / 'data/disordered_lgn.mat', variable_names=['inde'],
              simplify_cells=True)['inde'][0].astype(bool)
cat = data['Cat_select']
cx = np.flatnonzero((cat[:, 0] == 2) & (cat[:, 1] >= 0) & (cat[:, 1] <= 11))
rates = np.stack([values[cx] for values in data['f0_all']])
counts = np.rint(rates / 16.66).astype(np.uint16)
assert np.max(np.abs(counts.astype(float) * 16.66 - rates)) < 1e-12
assert data['T50'].shape == (3, 247)
savemat(root / 'data/historical_onset_reference.mat', {
    'counts': counts,
    'T50': np.asarray(data['T50'], dtype=float),
    't50_best': np.asarray(data['t50_best'], dtype=float),
    'cell_ids': np.flatnonzero(roi)[cx] + 1,
    'Cat': cat[cx], 'trial_count': 60.0,
    'histogram_start_s': .5, 'onset_s': 1.5, 'duration_s': 2.5,
    'legacy_rate_factor': 16.66, 'plot_offset_ms': 4.0,
}, do_compression=True)
print('Exported 3 phases x 247 complex cells x 2000 histogram bins.')
