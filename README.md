> **One-package distribution note:** This archive contains the original numerical experiment code, all four archived result directories, the paper, supplementary analysis/validation, and a LOCAL copy of `cec17_python-master.zip` for convenience. Its third-party redistribution terms have **not** been established. Keep the benchmark archive local and out of a public GitHub repository unless you confirm permission to redistribute it. The `.gitignore` excludes it from ordinary Git commits. For a new public GitHub repository, upload all the other files and folders. Functional smoke checks were run previously, but the complete historical experiment was **not** independently rerun. See [`START_HERE.md`](START_HERE.md) for exact steps.

# CREO-DE: Repair-Guided Evolutionary Optimisation

This repository contains source code and archived outputs for an experimental study of **Repair-Guided Evolutionary Optimisation (CREO)**, including the DE-based hybrid **CREO-DE**. The paper reports a CEC2017-derived numerical experiment and a matched ablation of the repair-guidance term.

> **Scope of the evidence.** The experiments evaluate the official CEC2017 *objective functions* inside an additional, author-defined preferred-region/repair protocol. They are **not** an unmodified CEC2017 unconstrained competition. GA, PSO, DE and GWO use penalty-mode selection/evaluation, while CREO variants use repair-mode evaluation; the strongest test of the guidance term is the within-mode CREO-DE versus CREO-DE-NoRG ablation. The archived outputs do not demonstrate that CREO-DE beats DE overall.


## Benchmark package supplied: functional verification update

The previously missing CEC2017 archive has now been supplied and locally compiled. Its SHA-256 checksum, six function-level evaluations, and four isolated reduced-budget main/ablation execution logs are documented in [`THIRD_PARTY_BENCHMARK_PROVENANCE.md`](THIRD_PARTY_BENCHMARK_PROVENANCE.md) and `validation/smoke_observed/`. Place your own copy of `cec17_python-master.zip` in the repository root locally; `.gitignore` excludes it from GitHub pending confirmed redistribution rights. For independent checks run `python validation/verify_benchmark.py` and `python validation/reproduce_smoke.py`. **These are functional smoke checks, not the full 30-run replication underlying the paper's tables.**

## Contents

| Path | Purpose |
|---|---|
| `experiment_d10.py` | Main benchmark: 29 functions, dimension 10, seven methods. |
| `experiment_d30.py` | Main benchmark: 29 functions, dimension 30, seven methods. |
| `ablation_experiment.py` | DE, CREO-DE-NoRG (beta = 0), CREO-DE; one chosen dimension per run. |
| `CREO-CEC2017.py` | Additional combined numerical experiment driver; not needed for the two reported dimension-specific archives. |
| `CREO.py` | Separate music-oriented CREO demonstration, **not** the driver for the CEC2017 tables. |
| `cec17_results_D10/`, `cec17_results_D30/` | Archived main benchmark CSV summaries and per-function plots. |
| `cec17_ablation_small_D10/`, `cec17_ablation_small_D30/` | Archived ablation CSV summaries and per-function plots. |
| `requirements.txt` | Python runtime dependencies. |
| `results.txt` | Additional results notes. |

## Prerequisites

- Python 3 with `venv` and `pip`; a C compiler (`gcc` on Linux or `clang` on macOS).
- A CEC2017 Python/C benchmark distribution stored as **`cec17_python-master.zip` in the repository root**. The scripts require the extracted distribution to contain `cec17_functions.py`, `cec17_test_func.c`, and the CEC2017 `input_data` used by that implementation. The required package has now been supplied privately and functionally tested. Its SHA-256 is `a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c`. It is **deliberately excluded from the public upload package** pending permission to redistribute third-party code. Place the supplied `cec17_python-master.zip` locally in the repository root before running any benchmark scripts. The original historic experiments did not record an archive hash, so byte-identical provenance remains unverified. See `THIRD_PARTY_BENCHMARK_PROVENANCE.md`.
- These scripts extract and compile a shared C library. Run them on a machine on which you are comfortable compiling the benchmark package. Inspect the third-party C and Python sources before execution.

### Install the Python dependencies

From the repository root, create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate                 # Linux/macOS
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, use `.venv\Scripts\activate` in place of `source ...`; **the current C compiler invocation targets Linux/macOS and needs adaptation before native Windows execution**. `pretty_midi` is needed for `CREO.py` (the separate music demonstration), not for the main numerical benchmark.

`requirements.txt` lists dependencies, not a claim that specific package versions were used in the original archived runs. For a fully pinned experimental environment, use the environment that actually produced the archive to record `python --version`, `python -m pip freeze`, operating system and compiler version.

## Study design reproduced by the scripts

The two main experiment drivers use functions **F1, F3–F30** (29 functions; F2 excluded), `POP_SIZE = 50`, `RUNS = 30`, and an objective-evaluation budget of `10000 × D`: **100,000 evaluations at D10** and **300,000 at D30**. Python and NumPy are seeded with `SEED=42` once at script startup. The pseudorandom state then advances sequentially across methods, runs and functions; run numbers do not identify common independently seeded matched pairs.

The search bounds are `[-100,100]^D`. The auxiliary preferred-region violation score includes departures from `[-80,80]` and a penalty when the mean absolute coordinate exceeds 60. The repair function first clips to benchmark bounds, then repeatedly scales the point by `0.8` up to ten times until the preferred-region violation is at most `1e-12`. A failed ten-step repair is returned even if its remaining violation is nonzero.

In penalty mode the CEC objective is evaluated at the clipped point and `fitness = objective + 1e6 × raw_violation`. In repair mode it is evaluated at the repaired point and `fitness = objective`. Both modes use feasibility-first selection. As a result, the archived `objective` columns describe the *mode-specific evaluated point*, not the penalised fitness and not a common unmodified CEC2017 selection protocol.

CREO-DE uses `F=0.7`, `CR=0.9`, `alpha=0.18`, `beta=0.90` and `gamma=0.02`. The matched CREO-DE-NoRG ablation differs only in setting `beta=0`. Its stored population is the box-clipped trial, with the fully repaired point and displacement retained in the corresponding evaluation record.

## Reproducing Experiment 1 (main comparison)

**Back up the existing `cec17_results_D10/` and `cec17_results_D30/` directories before rerunning**: the scripts write tables and plots under these names and can overwrite archived material. To protect the published archive, clone the repository into a clean separate directory, put the external `cec17_python-master.zip` in that new root, install dependencies, then run:

```bash
python experiment_d10.py
python experiment_d30.py
```

Outputs include `all_benchmark_results.csv` (one per-function, per-algorithm summary row), per-function `tables/` and `plots/`, `overall_friedman_mean_ranks.csv`, `friedman_ranks.csv`, and `wilcoxon_creo_de_vs_others.csv`. The latter contains **unadjusted** Wilcoxon results based on sequential run-indexed samples; it should not be interpreted as a corrected, properly paired multi-problem significance analysis. `overall_algorithm_summary.csv` averages raw objective values across heterogeneous CEC functions and should **not** be used to compare algorithms across functions; prefer function-wise objective values, ranks, and win/tie/loss counts.

The existing per-function plots have a generation/iteration x-axis, not a common objective-evaluation x-axis. The code sets `SAVE_ALL_CURVES=False`; individual-run curves are therefore not exported by default.

## Reproducing Experiment 2 (ablation)

The current `ablation_experiment.py` sets `DIMENSIONS = [30]`. To create each dimension's result separately, run **in a clean working copy**, edit that line, and execute the script once per dimension:

1. Set `DIMENSIONS = [10]`, then run `python ablation_experiment.py`. The output directory will be `cec17_ablation_small_D10/`.
2. Set `DIMENSIONS = [30]`, then run `python ablation_experiment.py` again. The output directory will be `cec17_ablation_small_D30/`.

The reported archived comparison of CREO-DE against CREO-DE-NoRG is **12 wins / 0 ties / 17 losses at D10** and **1 / 0 / 28 at D30**, using the per-function `Mean` in `all_ablation_results.csv` (lower is better). These are descriptive results under the auxiliary repair protocol. The DE comparator runs in penalty mode; it is not a beta-only ablation.

## Reanalyse the archived experimental results (no CEC package required)

This ZIP includes the **original experiment code and all four original result directories unchanged**, plus a new, independent descriptive-analysis script. To check the archived values and regenerate audit tables without rerunning the expensive optimiser:

```bash
python -m pip install -r requirements.txt
python analysis/rebuild_archived_results.py
python -m unittest discover -s tests -v
```

The analysis outputs are stored in `analysis/generated/`. The script reads original full-precision CSV means (never the rounded PDF values), verifies the expected 29 functions and algorithm lists, computes midranks for exact numerical ties, and writes full seven-method and five-method main ranks, three-method ablation ranks, all CREO-DE pairwise win/tie/loss comparisons, and function-wise mean matrices. **These are derivations of archived summaries, not new CEC2017 optimisation runs.** The smoke tests execute the *actual original* repair routines in isolation and verify summary file shapes; passing them is not proof of full-experiment replication.

The manuscript in `paper/` includes the rank tables and the original function-wise results. Do not replace the archived source CSVs with reanalysed results. The user-supplied benchmark archive has now been compiled and smoke-tested. Full replication still requires confirming its historical provenance, a suitably configured C compiler, and enough computing resources to repeat every 30-run, 29-function experiment.

## Results and provenance

The archived main results are in `cec17_results_D10/all_benchmark_results.csv` and `cec17_results_D30/all_benchmark_results.csv`; ablation results are in the corresponding `all_ablation_results.csv` files. These files contain per-function **summary statistics**, not the 30 underlying terminal objective values per method/function, so they do not independently support re-computing a full run-level inferential analysis. Save the actual per-run objective, fitness, violation, seed, algorithm parameters, and benchmark package hash in future exports.

Source-code snapshot represented by the supplied archive: **`5e04fb2`**. Documentation added after that commit belongs to a later commit and should not be described as part of the older snapshot. Before making a permanent data citation or DOI, verify that the public GitHub revision includes exactly the source and output files used by the manuscript.

## Reuse, license and citation

Add the correct **repository license** only after deciding your distribution terms, and **verify that the third-party CEC2017 package's terms permit redistribution** before committing that archive or a compiled derivative. A repository source URL and the exact benchmark-package provenance should be added here by the author. Cite the associated paper once its bibliographic details and DOI have been assigned; do not invent a DOI or claim external replication that has not occurred.

## Publication package

`paper/CREO_REPOSITORY_PACKAGE.tex` is the manuscript source; `paper/CREO_REPOSITORY_PACKAGE.pdf` is the compiled reading copy. These add audited descriptive ranks for all seven archived methods and the ablation variants while retaining the limitations of the augmented evaluation protocol. The `paper/` directory is documentation and does not alter the historical experiment files.

## Uploading this package to GitHub

The repository-ready distribution is a ZIP containing the **files and folders** for the repository root, not a single experimental result. Extract it first, then copy/commit all extracted contents into your existing repository. Preserve the four `cec17_*` folders exactly, and do not upload only the ZIP as a GitHub file. Confirm the new Git commit hash after upload and cite that new snapshot for the revised paper. The benchmark ZIP is intentionally omitted from the public upload while redistribution rights and historical-version equivalence remain unverified. A supplied local copy was successfully compiled and exercised; see the validation records.
