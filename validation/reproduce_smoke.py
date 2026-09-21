#!/usr/bin/env python3
"""Run isolated low-budget copies of the original main and ablation drivers.

Full paper experiment settings and archived data are NEVER changed by this script.
Requires the local ignored cec17_python-master.zip archive in repository root.
"""
from pathlib import Path
import csv, re, shutil, subprocess, sys, tempfile, hashlib
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'cec17_python-master.zip'
if not ARCHIVE.is_file():sys.exit('Copy cec17_python-master.zip into the repository root first.')
print('Benchmark SHA256:',hashlib.sha256(ARCHIVE.read_bytes()).hexdigest())
for study in ('main','ablation'):
    for dimension,fn in ((10,1),(30,3)):
        original=ROOT/(f'experiment_d{dimension}.py' if study=='main' else 'ablation_experiment.py')
        source=original.read_text()
        changes={
            r'^RUNS\s*=\s*30\s*$':'RUNS = 1',
            r'^FUNCTION_IDS\s*=\s*\[fid for fid in range\(1, 31\) if fid (?:not in DELETED_FUNCTION_IDS|!= 2)\]\s*$':f'FUNCTION_IDS = [{fn}]',
            r'^MAX_EVALS_FACTOR\s*=\s*10000.*$':'MAX_EVALS_FACTOR = 20  # small-scale functional test ONLY',
        }
        if study=='ablation':changes[r'^DIMENSIONS\s*=\s*\[30\]\s*$']=f'DIMENSIONS = [{dimension}]'
        for pattern,replacement in changes.items():
            source,n=re.subn(pattern,replacement,source,count=1,flags=re.M)
            if n!=1:raise RuntimeError(f'Smoke adaptation pattern not matched: {pattern} ({study}, D{dimension})')
        with tempfile.TemporaryDirectory(prefix=f'creo_{study}_D{dimension}_') as temporary:
            temp=Path(temporary)
            (temp/'smoke.py').write_text(source)
            shutil.copy2(ARCHIVE,temp/ARCHIVE.name)
            result=subprocess.run([sys.executable,'smoke.py'],cwd=temp,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=240)
            destination=ROOT/'validation'/'smoke_replayed'/f'{study}_D{dimension}'
            destination.mkdir(parents=True,exist_ok=True)
            (destination/'execution.log').write_text(result.stdout)
            (destination/'smoke.py').write_text(source)
            if result.returncode:
                print(result.stdout[-3000:]);raise RuntimeError(f'{study} D{dimension}: process returned {result.returncode}')
            result_folder=temp/(f'cec17_results_D{dimension}' if study=='main' else f'cec17_ablation_small_D{dimension}')
            result_file=result_folder/('all_benchmark_results.csv' if study=='main' else 'all_ablation_results.csv')
            with result_file.open(newline='') as f:rows=list(csv.DictReader(f))
            assert len(rows)==(7 if study=='main' else 3), (study,dimension,len(rows))
            shutil.copy2(result_file,destination/result_file.name)
            print(f'PASS {study}: F{fn}, D{dimension}, 1 run, {20*dimension} evaluations/algorithm, {len(rows)} rows')
print('All reduced-budget smoke checks passed. Original 30-run archives have not been rerun or overwritten.')
