# CREO-Research: Repair-Guided Evolutionary Optimisation

This repository contains the Python implementations, archived numerical results, analysis scripts, verification materials, and manuscript for **Repair-Guided Evolutionary Optimisation (CREO)**, including the Differential Evolution hybrid **CREO-DE**.

**Research scope:** The study evaluates CEC2017 objective functions under an additional, author-defined preferred region and repair protocol. These are **not** unmodified CEC2017 competition results: GA, PSO, DE and GWO use penalty-mode evaluation/selection, while CREO variants use repair-mode evaluation. The most direct test of the repair-guidance term is the CREO-DE versus CREO-DE-NoRG ablation. The archived results do not establish that CREO-DE outperforms DE overall.

## Repository contents

| Path | Contents |
| --- | --- |
| `experiment_d10.py` | Main seven-algorithm experiment, dimension 10. |
| `experiment_d30.py` | Main seven-algorithm experiment, dimension 30. |
| `ablation_experiment.py` | CREO-DE, CREO-DE-NoRG and DE comparison; dimension configured in the script. |
| `CREO-CEC2017.py` | Additional combined numerical experiment script. |
| `CREO.py` | Separate music-oriented CREO demonstration; not the driver for the reported CEC2017 tables. |
| `cec17_results_D10/`, `cec17_results_D30/` | Original main-experiment summary CSV files, tables and convergence plots. |
| `cec17_ablation_small_D10/`, `cec17_ablation_small_D30/` | Original ablation summary CSV files and convergence plots. |
| `analysis/` | Reanalysis of the archived CSV summaries and generated descriptive tables. |
| `validation/` | Benchmark validation scripts, reduced-budget execution logs and smoke-test outputs. |
| `tests/` | Automated checks of selected algorithm functions and archived result structures. |
| `paper/CREO_REPOSITORY_PACKAGE.pdf`, `paper/CREO_REPOSITORY_PACKAGE.tex` | Earlier repository manuscript package, retained for provenance. |
| `paper/CREO_Elsevier/` | Elsevier CAS submission package: source, class assets, highlights, and compiled PDF. |
| `REPRODUCIBILITY_STATUS.md`, `EXPERIMENT_EXECUTION_REPORT.md` | Record what was checked and what remains unreplicated. |
| `THIRD_PARTY_BENCHMARK_PROVENANCE.md` | Provenance and verification details for the external CEC2017 implementation. |
| `requirements.txt` | Python dependencies. |

## CEC2017 benchmark: original source and attribution

The CEC2017 benchmark implementation used by the experiment scripts is the **Python wrapper maintained by Marcelo Lacerda**, available from the following original repository:

**https://github.com/lacerdamarcelo/cec17_python**

That repository contains `cec17_functions.py`, `cec17_test_func.c` and `input_data/`. The wrapper exposes the CEC2017 single-objective, bound-constrained numerical optimisation benchmark functions originally implemented in C. It is **third-party software, not part of the original CREO contribution**.

The upstream repository attributes the original C benchmark and problem definitions to:

> N. H. Awad, M. Z. Ali, J. J. Liang, B. Y. Qu, and P. N. Suganthan, “Problem Definitions and Evaluation Criteria for the CEC 2017 Special Session and Competition on Single Objective Bound Constrained Real-Parameter Numerical Optimization,” Technical Report, Nanyang Technological University, Singapore, November 2016.

### Obtain the benchmark for local reproduction

1. Open https://github.com/lacerdamarcelo/cec17_python.
2. Select **Code → Download ZIP** on the repository's `master` branch.
3. Save the downloaded archive in the **root of this CREO repository**, next to `experiment_d10.py`, under the exact filename **`cec17_python-master.zip`**. Do not extract it yourself: the CREO experiment scripts extract the archive when needed.
4. Ensure that you have a working C compiler: `gcc` on Linux or `clang` on macOS. The original scripts compile the wrapper's C implementation locally.

**Public GitHub upload:** Do **not** commit `cec17_python-master.zip` or its compiled/extracted derivatives unless you have confirmed permission to redistribute them. The supplied `.gitignore` excludes the local archive from normal Git commits; when using GitHub's browser uploader you must exclude it manually. Link to the upstream source instead.

**Archive identity:** A previously supplied and functionally checked local benchmark ZIP had SHA-256:

```text
a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c
```

A fresh download from the upstream `master` branch **may have a different checksum** if the upstream repository changed. `validation/verify_benchmark.py` checks against the recorded hash; a mismatch means that this exact archive was not verified, **not necessarily that the new download is defective**. If you use a different revision, record its commit, checksum and validation outcomes, and do not describe it as the exact historical benchmark version. The benchmark archive used for the original 30-run experiments was not hashed at the time, so byte-identical historical provenance is unconfirmed.

## Installation

The numerical experiments need Python 3, `pip`, a C compiler and the external benchmark ZIP described above. From this repository's root on **macOS or Linux**:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The dependencies listed in `requirements.txt` are NumPy, Pandas, SciPy, Matplotlib and `pretty_midi`. `pretty_midi` is for the separate music demonstration. Dependency versions from the historical experiment environment were not recorded; the current file does not claim to reproduce that environment exactly. The C compilation steps in the current experiment scripts target Linux/macOS and need adaptation for native Windows use.

## Experimental protocol

- **Functions:** CEC2017 F1 and F3–F30 (29 functions; F2 excluded).
- **Dimensions:** 10 and 30.
- **Runs:** 30 per function and algorithm.
- **Population:** 50.
- **Objective-evaluation budget:** `10000 × D` per run: 100,000 at D10 and 300,000 at D30.
- **Seeding:** the scripts seed Python and NumPy once using `SEED = 42`. The random state then advances sequentially across experiments; runs with the same index are not separately seeded matched pairs.
- **CREO-DE settings:** `F = 0.7`, `CR = 0.9`, `alpha = 0.18`, `beta = 0.90`, `gamma = 0.02`. The CREO-DE-NoRG ablation sets `beta = 0`.

The search box is `[-100, 100]^D`. An additional preferred-region violation score measures departures from `[-80, 80]^D` and applies a penalty when the mean absolute coordinate exceeds 60. Repair first clips to the search box and then scales the candidate by `0.8` for up to ten iterations, stopping early when the preferred-region violation is at most `1e-12`.

In **penalty mode**, the CEC objective is evaluated at the box-clipped candidate and the fitness adds `10^6 ×` its violation. In **repair mode**, the CEC objective is evaluated at the repaired point and the fitness is the objective value. Both use feasibility-first selection. Consequently, results from the two modes do not represent identical unmodified benchmark protocols.

## Run the original main experiments

**Use a separate working copy:** the experiment drivers write into the same directory names as the archived results and may overwrite published data.

With the external ZIP in the repository root and dependencies installed, run:

```bash
python experiment_d10.py
python experiment_d30.py
```

The main result summaries are written to:

```text
cec17_results_D10/all_benchmark_results.csv
cec17_results_D30/all_benchmark_results.csv
```

The output directories also contain per-function tables and plots, mean-rank summaries and unadjusted pairwise Wilcoxon result files. The reported Wilcoxon tests are not corrected for multiple comparisons, and the sequential run indices do not establish paired random trials. The function-wise and aggregate-rank tables are more appropriate than averaging raw objective values across different functions. Convergence plots use generation/iteration, not a common objective-evaluation axis.

## Run the original ablation experiments

The existing `ablation_experiment.py` is configured with `DIMENSIONS = [30]`. To reproduce each dimension in a **separate working copy**, set the relevant line and run the script:

```python
DIMENSIONS = [10]
```

```bash
python ablation_experiment.py
```

Then set:

```python
DIMENSIONS = [30]
```

```bash
python ablation_experiment.py
```

The respective summary files are:

```text
cec17_ablation_small_D10/all_ablation_results.csv
cec17_ablation_small_D30/all_ablation_results.csv
```

Based on the archived per-function means, CREO-DE versus CREO-DE-NoRG has **12 wins / 0 ties / 17 losses at D10** and **1 win / 0 ties / 28 losses at D30** (lower objective value is better). These are descriptive outcomes under the additional repair protocol; the DE comparator uses penalty mode and is not a beta-only ablation.

## Reanalyse archived results and run verification checks

The following command recalculates descriptive results **from the stored CSV summaries**, without rerunning the expensive optimisation experiments or requiring the CEC2017 archive:

```bash
python analysis/rebuild_archived_results.py
python -m unittest discover -s tests -v
```

With the locally downloaded benchmark ZIP in place, the following scripts test benchmark compilation/evaluation and run isolated, reduced-budget functional experiments:

```bash
python validation/verify_benchmark.py
python validation/reproduce_smoke.py
```

The observed smoke-test outputs and logs are stored under `validation/smoke_observed/`; they cover one run for selected functions with only `20 × D` evaluations per algorithm. **They do not reproduce the complete 29-function, 30-run historical experiments.** See `REPRODUCIBILITY_STATUS.md` and `EXPERIMENT_EXECUTION_REPORT.md` for verification scope and limitations.

## Data availability and limitations

The four archived result folders preserve original per-function summary CSV files and convergence plots. **They do not contain every individual run's final objective value**, so the complete run-level statistical analysis cannot be independently reconstructed from the published summaries alone. Future experiments should export per-run values, seeds, parameters, feasibility/fitness values and the exact external benchmark revision and checksum.

The original source snapshot was reported as commit `5e04fb2`; the expanded documentation and validation files were prepared later. Record and cite the actual commit or archived DOI corresponding to the final public release rather than assigning the earlier commit to files added subsequently.

## Manuscript and citation

The submission-formatted Elsevier CAS package is available at:

- [Compiled Elsevier manuscript](paper/CREO_Elsevier/main.pdf)
- [Elsevier LaTeX source](paper/CREO_Elsevier/main.tex)
- [Package instructions](paper/CREO_Elsevier/README.md)
- [Research highlights](paper/CREO_Elsevier/highlights.txt)

The package uses Elsevier's CAS single-column class and keeps all required source files at one folder level for Editorial Manager. The section hierarchy, front matter, keywords, highlights, numbered citations, appendix handling, CRediT metadata, cross-references, and compilation were checked. The earlier repository manuscript files remain in `paper/` for provenance.

When describing the numerical benchmark, attribute the Python wrapper to **Marcelo Lacerda** and the CEC2017 problem definitions to **Awad et al. (2016)**, using the upstream source and technical report cited above. Machine-readable citation metadata is provided in [`CITATION.cff`](CITATION.cff). Add the CREO paper's definitive journal citation and DOI after publication; none is asserted here.

## Licence and third-party code

The CREO source and the upstream CEC2017 implementation have separate authorship and distribution considerations. This repository does not grant a licence to redistribute the third-party benchmark. Choose and add a licence for the CREO material only after establishing the appropriate rights and deciding your distribution terms. Obtain the benchmark directly from its upstream repository unless its copyright holder explicitly permits redistribution.
