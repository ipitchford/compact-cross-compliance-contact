#!/usr/bin/env python3
"""Verify the distributed package's SHA-256 manifest (standard library only)."""
from pathlib import Path
import hashlib
import sys

def main() -> int:
    root=Path(__file__).resolve().parents[1]
    manifest=root/'MANIFEST.sha256'
    if not manifest.is_file():
        print('ERROR: MANIFEST.sha256 is missing',file=sys.stderr);return 2
    count=0;bad=[]
    for line in manifest.read_text().splitlines():
        if not line.strip():continue
        expected,rel=line.split(maxsplit=1)
        rel=rel.lstrip('*');path=(root/rel).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            bad.append(rel+' (missing or unsafe)');continue
        actual=hashlib.sha256(path.read_bytes()).hexdigest();count+=1
        if actual!=expected:bad.append(rel+' (hash mismatch)')
    for entry in bad:print('FAIL '+entry)
    print(f'{count} files checked; {len(bad)} failures')
    return 1 if bad else 0
if __name__=='__main__':raise SystemExit(main())
