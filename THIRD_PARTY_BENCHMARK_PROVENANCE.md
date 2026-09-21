# CEC2017 dependency: supplied, functionally verified, not redistributable by default

The previously missing local CEC2017 Python/C archive was supplied separately and has now been inspected and tested. To reproduce locally, place **`cec17_python-master.zip` in the root of this repository**. It is intentionally excluded by `.gitignore` and is **NOT** bundled into the public GitHub upload package because the archive contains no standalone license file and redistribution rights have not been established.

- User-supplied archive SHA-256: `a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c`
- Archive ZIP comment includes revision-like identifier: `f60d0000be91ef080deb10b09dce07e4a1c30d9a` (not an independently verified upstream URL).
- Contents: CEC2017 Python wrapper `cec17_functions.py`, C source `cec17_test_func.c`, and function input data.
- Tested: compiled the supplied C source using GCC and evaluated F1, F3 and F30 at zero vectors for D10 and D30; all six outputs were finite. The D10 F1 result agrees with the numerical example in the archive's README at the displayed precision.
- Executed: original main experiment drivers and original ablation experiment driver in **isolated reduced-budget copies** for one function at each dimension. All four smoke runs returned success; results and logs are in `validation/smoke_observed/`. Each used one stochastic run and `20×D` function evaluations per method, not the paper's full 30 runs and `10,000×D` budget.
- **NOT established**: that this user-supplied archive is byte-identical to the version used in the original historical 30-run benchmark result CSVs, as no contemporaneous archive hash was recorded.
- **NOT done**: full 29-function/30-run replication of the two main experiments and two ablation dimensions. Original CSV evidence remains untouched.

## Local verification commands

```bash
python -m pip install -r requirements.txt
python validation/verify_benchmark.py
python validation/reproduce_smoke.py
python analysis/rebuild_archived_results.py
python -m unittest discover -s tests -v
```

To run the **full historical protocol**, use fresh clone and put the benchmark archive at the root before executing `python experiment_d10.py`, `python experiment_d30.py`, and `python ablation_experiment.py` separately with `DIMENSIONS=[10]` and `DIMENSIONS=[30]`. **Do not run the full scripts against the only copy of the published result archives**: they write to the same output directories. This protocol entails roughly 3.48 billion requested objective evaluations for the two main and two ablation dimensions combined, in addition to Python algorithm overhead, so it is not comparable to a smoke test.

For public redistribution, confirm the upstream repository URL, exact revision and licensing terms with the benchmark copyright holder. If permission is granted, add clear third-party attribution and licence files. Otherwise keep the benchmark ZIP external and publish its SHA-256 and acquisition instructions.
