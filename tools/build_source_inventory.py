"""Record relative source provenance. Usage: script ARCHIVE_ROOT AUDIT_ROOT.

Only the JSON inventory in this release is written; original files are read.
This maintainer utility uses the Python standard library.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
archive, audit = map(Path, sys.argv[1:3])
core = ['V1_4_fast.m', 'LGN3.m', 'LGN4.m', 'NLIF.m', 'ex_gen2.m', 'ex_gen_onset_test.m']
sources = [(archive, name, f'src/{name}', 'runtime source') for name in core]
sources += [(archive, name, None, 'selection evidence') for name in [
    'stat_script.m', 'suddeninput.m', 'suddeninput_multi.m',
    '25_9_28/parameters.txt', '25_9_28/10_13/phase_code.m',
    '2026_1_7_paperfigure/revised/fig5_icell_generate_data.m',
    'native/validation/source_manifest.json',
    'native/validation/optimization_final_v2_vs_matlab_2500ms_seed20261005_20261005_161111_168.json']]
sources += [(archive, '25_9_28/10_13/10_25 version of parameter for new model.mat',
             'data/disordered_lgn.mat', '11 exactly retained variables')]
sources += [(archive, '25_9_28/lowlgn/low lgn full param.mat',
             'data/low_lgn.mat', '11 exactly retained Low-LGN variables')]
base = 'results/onset_reproduction/20261028_old_parameters/'
sources += [(audit, base + name, target, kind) for name, target, kind in [
    ('conclusions.txt', None, 'historical reproduction limits'),
    ('reference/SI_draft_2_numerical_audit_2026-09-19.md', None, 'F1/onset/kernel audit evidence'),
    ('reference/legacy_matlab_reference.mat', 'data/historical_onset_reference.mat', 'compact histogram fixture')]]

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

records = []
for origin, rel, target, kind in sources:
    path = origin / rel
    rec = dict(source=f'{origin.name}/{rel}', source_bytes=path.stat().st_size,
               source_sha256=digest(path), purpose=kind)
    if target:
        rec.update(release_file=target, release_sha256=digest(root / target))
    records.append(rec)
payload = dict(assembled_date='2026-11-23', low_lgn_added_date='2026-11-24',
               scope='Two original MATLAB models: high_cx F1=2.25/1.5 and low_lgn F1=1.3/1.2; later author request includes Low-LGN',
               original_git_head=subprocess.check_output(['git', '-C', str(archive), 'rev-parse', 'HEAD'],text=True).strip(),
               source_state='Working copy had existing tracked modifications and untracked files; hashes define the selected content.',
               files=records)
(root / 'docs/source_inventory.json').write_text(json.dumps(payload, indent=2)+'\n',encoding='utf-8',newline='\n')
print(f'Recorded {len(records)} source files; no machine-specific absolute paths.')
