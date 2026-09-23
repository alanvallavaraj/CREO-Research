#!/usr/bin/env python3
"""Download the exact 16SGOSeq bacterial-variant FASTA used by the benchmark."""
import hashlib
from pathlib import Path
from urllib.request import urlopen

URL='https://zenodo.org/api/records/15209015/files/bacteria_variants.fasta/content'
SHA='793db5908177dcc3998911cf487d10e5eaf4b17c06a71bf7ca52beafd42c1833'
path=Path(__file__).resolve().parent/'data'/'16SGOSeq_bacteria_variants.fasta'
path.parent.mkdir(exist_ok=True)
if not path.exists():
    with urlopen(URL,timeout=120) as response,path.open('wb') as dest:
        while piece:=response.read(1024*1024):dest.write(piece)
digest=hashlib.sha256(path.read_bytes()).hexdigest()
if digest!=SHA:raise RuntimeError(f'SHA-256 mismatch: {digest}')
print(path,'verified',digest)
