"""Build a local review ZIP. No upload or publication occurs.

Usage: python tools/package_release.py [DESTINATION.zip]
The destination must not already exist. Generated results and Git internals
are excluded, and archived bytes are checked against the release manifest.
"""
import hashlib
import json
import runpy
import sys
import zipfile
from pathlib import Path

root=Path(__file__).resolve().parents[1]
destination=Path(sys.argv[1]) if len(sys.argv)>1 else root.with_suffix('.zip')
assert not destination.exists(), f'Destination already exists: {destination}'
assert root not in destination.resolve().parents, 'Write the ZIP outside the repository.'
runpy.run_path(str(root/'tools/check_release.py'),run_name='__main__')
manifest=json.loads((root/'validation/release_manifest.json').read_text(encoding='utf-8'))
names=[item['file'] for item in manifest['files']]+['validation/release_manifest.json','outputs/.gitkeep']
with zipfile.ZipFile(destination,'x',compression=zipfile.ZIP_DEFLATED) as archive:
    for name in names:
        archive.write(root/name, f'{root.name}/{name}')
with zipfile.ZipFile(destination) as archive:
    assert archive.testzip() is None
    for item in manifest['files']:
        data=archive.read(f'{root.name}/{item["file"]}')
        assert hashlib.sha256(data).hexdigest()==item['sha256']
print(f'Created and verified {destination.name}: {destination.stat().st_size:,} bytes.')
