# What was executed for this repository-ready package

**Archived optimisation data:** Preserved from the supplied `CREO-DE-main.zip` byte-for-byte; no numerical optimiser run in this package is represented as a new experiment.

**Reanalysis actually executed:** `python analysis/rebuild_archived_results.py`. It reads and validates four original aggregate CSVs (main D10/D30 and ablation D10/D30). It generates all-seven and main-five ranks, ablation ranks, pairwise wins/ties/losses and a function-wise full-precision mean matrix in `analysis/generated/`. Computed results are included in the paper and in the generated CSV files.

**Unit checks actually executed:** `python -m unittest discover -s tests -v`: 5 tests pass. These use the three original repair functions extracted from the source file without launching the missing third-party benchmark package. They verify clipping, no-op on preferred-region points, shrink-to-preferred behaviour, bounded repair, and original summary structure.

**Archive consistency checked:** Four master aggregate CSVs have the expected 203/203/87/87 rows, all 116 per-function summary table CSVs match their respective master records, and every original experiment source/result/plot file (except the intentionally replaced README and requirements text) matches the supplied archive SHA-256.

**Full CEC2017 experiment NOT rerun:** the original scripts fail at import-time preparation without `cec17_python-master.zip` (not included in the supplied archive); substituting an arbitrary different benchmark implementation would not verify the historical experimental outputs. Total original main benchmark budgets are 29 functions x 30 runs x 7 methods x 10000D evaluations **for each dimension**, plus ablation. Original per-run terminal scores and individual seed assignments are also not archived. It is not scientifically appropriate to invent these from means/SD or assert independent reproduction without them.

**Before submitting an independent reproducibility claim:** Obtain the exact original CEC benchmark archive from the original source, verify SHA-256 and lawful redistribution, rerun the unchanged experiment scripts in a separate clone, preserve seed-specific scores, then compare derived summaries with the archived results. The manuscript explicitly identifies the existing numbers as archived experiments and the new tables as reanalysis.
