# Reproducibility and evidence status (updated after benchmark ZIP upload)

The user-supplied original CREO source scripts and four historical CEC result directories are preserved byte-for-byte in this package. Descriptive rank and win/tie/loss reanalyses were derived from archived CSV summaries, not from new full-budget optimiser runs.

## Confirmed by execution

- A separately supplied local `cec17_python-master.zip` archive has SHA-256 `a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c`.
- The C source from that archive compiled successfully with GCC on the test machine, and its Python wrapper evaluated F1, F3 and F30 at D10/D30 using finite values.
- In isolated temporary directories, the original dimension-specific main and ablation scripts each completed a **one-function, one-run, reduced-budget** execution at D10 and at D30. Main smoke results have seven method rows per dimension; ablation smoke results have three rows per dimension. They are stored only under `validation/smoke_observed/`, NOT under the historical result directories.
- Five pre-existing local source/archive unit tests passed; they do not independently verify the full optimization experiment.
- The paper was recompiled with benchmark provenance and smoke-test status added to Code and Data Availability.

## Not yet established

- No 29-function, 30-run full-budget optimisation experiment was rerun. Each full main dimension uses 29 functions × 30 runs × 7 methods × 10,000D evaluations; the full two-dimensional ablation uses three methods. Together they request approximately **3.48 billion** objective evaluations before Python orchestration costs.
- No contemporaneous SHA-256 hash was saved for the third-party benchmark package used to produce the **historical** CSVs. The subsequently supplied benchmark ZIP may be the original one, but that cannot be proved from the current evidence.
- The original archive contains aggregate summaries, not per-run objective/seed records; exact historical run-level inferential replication is therefore not supported.
- The supplied benchmark ZIP contains no identifiable standalone distribution license. It is a local, git-ignored dependency, not a file to publish on GitHub without confirmed redistribution permission. The upstream URL should also be verified.

To rerun locally, place the supplied benchmark ZIP at repository root and follow `README.md` and `THIRD_PARTY_BENCHMARK_PROVENANCE.md`. Run the full suite only in a fresh clone to avoid overwriting historical archives; log per-run outcomes and the exact benchmark hash for any new publication evidence.
