"""Check the portable release payload without reading generated experiments.

Run from any directory with Python 3.11+. No third-party packages required.
The manifest is regenerated after publication-document edits are finished.
"""
import hashlib
import json
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1]
excluded = {'.git', 'outputs', '__pycache__'}
files = sorted(p for p in root.rglob('*') if p.is_file()
               and not any(x in excluded for x in p.relative_to(root).parts)
               and p.name != 'release_manifest.json')
problems=[]
for p in files:
    rel=p.relative_to(root).as_posix()
    if p.suffix.lower() in {'.asv','.exe','.mexw64','.xls','.log'}:
        problems.append(f'Unexpected payload: {rel}')
    if p.stat().st_size>50_000_000:
        problems.append(f'Review large file: {rel}')
    if p.suffix in {'.m','.md','.json','.py','.in'}:
        text=p.read_text(encoding='utf-8')
        # Match Windows drive paths and user-home paths, not citations/URLs.
        if re.search(r'(?i)(?<![a-z0-9_])(?:[a-z]:[/\\]|/(?:Users|home)/[a-z0-9_.-]+/)',text):
            problems.append(f'Machine-specific path: {rel}')
assert not problems, '\n'.join(problems)
records=[]
for p in files:
    with p.open('rb') as stream: digest=hashlib.file_digest(stream,'sha256').hexdigest()
    records.append(dict(file=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=digest))
payload=dict(files=records,total_bytes=sum(x['bytes'] for x in records),
             largest_bytes=max(x['bytes'] for x in records),passed=True)
(root/'validation/release_manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
print(f'PASS: {len(files)} files, {payload["total_bytes"]:,} bytes; no absolute local paths or oversized files.')
