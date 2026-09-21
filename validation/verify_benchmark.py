#!/usr/bin/env python3
"""Verify a user-supplied local CEC2017 archive without modifying original results."""
from pathlib import Path
import hashlib, os, platform, subprocess, sys, tempfile, zipfile

REPO = Path(__file__).resolve().parents[1]
ARCHIVE = REPO / 'cec17_python-master.zip'
EXPECTED_SHA256 = 'a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c'

if not ARCHIVE.exists():
    sys.exit('Missing local benchmark: copy cec17_python-master.zip to the repository root (do not git commit it).')
actual = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
print('Benchmark SHA-256:', actual)
if actual != EXPECTED_SHA256:
    sys.exit('This archive does not match the uploaded benchmark package; do not assume it reproduces the historical experiments.')

with tempfile.TemporaryDirectory(prefix='cec17_validate_') as tmp:
    root = Path(tmp).resolve()
    with zipfile.ZipFile(ARCHIVE) as z:
        for info in z.infolist():
            target = (root / info.filename).resolve()
            if target != root and root not in target.parents:
                sys.exit('Unsafe archive path in CEC benchmark ZIP')
        z.extractall(root)
    package = root / 'cec17_python-master'
    assert (package/'cec17_functions.py').is_file()
    assert (package/'cec17_test_func.c').is_file()
    assert (package/'input_data').is_dir()
    c_path = package/'cec17_test_func.c'
    if sys.platform == 'darwin':
        source = c_path.read_text()
        c_path.write_text(source.replace('#include <malloc.h>', ''))
        command=['clang','-dynamiclib','-O2','-o','cec17_test_func.so','cec17_test_func.c','-lm']
    else:
        command=['gcc','-fPIC','-shared','-O2','-o','cec17_test_func.so','cec17_test_func.c','-lm']
    subprocess.run(command,cwd=package,check=True)
    smoke=r"""from cec17_functions import cec17_test_func
import math
for D in (10,30):
    for f in (1,3,30):
        y=[0.0]
        cec17_test_func([0.0]*D,y,D,1,f)
        if not math.isfinite(y[0]):raise ValueError((D,f,y[0]))
        print('CEC2017: D=%d F%d f(0)=%.12g' % (D,f,y[0]))
"""
    subprocess.run([sys.executable,'-c',smoke],cwd=package,check=True)
print('PASS: supplied benchmark compiled and evaluated six finite test cases. This is not a full optimiser replication.')
