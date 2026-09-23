# CREO-Research: Repair-Guided Evolutionary Optimisation

This repository contains the research code, archived aggregate results, reproducibility checks, and manuscript package prepared for **Information Sciences** (Elsevier, ISSN 0020-0255) for the CREO family of evolutionary optimisation methods.

> **Submission status:** preprint prepared for journal submission; not peer reviewed or accepted. The repository is currently private.

## Scope and interpretation

The experiments use the CEC 2017 numerical optimisation functions together with an **author-defined preferred-region and repair protocol**. They must not be described as unmodified CEC 2017 competition results.

The main benchmark drivers compare standard optimisers, CREO variants, and CREO-DE. In those drivers, standard baselines use a penalty treatment while repair is enabled for the CREO repair variants. The matched ablation in `ablation_experiment.py` is the appropriate direct test of the region-guidance component: `DE`, `CREO-DE-NoRG`, and `CREO-DE` use the same experimental driver.

## Release status

The original public snapshot is commit `5e04fb2be03071ef601905c2dd3d04390c6f90b3`. The paper-submission release adds documentation, validation, reanalysis, and the Elsevier manuscript without rewriting the archived experiment outputs.

- Archived CSV files, plots, and original experiment scripts are preserved.
- The full 30-run experiments were **not rerun** during release preparation.
- Deterministic integrity tests and lightweight smoke-test tooling are provided.
- The historical benchmark dependency cannot be proven byte-identical because its exact archive was not stored with the original run.

See `REPRODUCIBILITY_STATUS.md`, `EXPERIMENT_EXECUTION_REPORT.md`, and `THIRD_PARTY_BENCHMARK_PROVENANCE.md` before reusing the reported results.

## Repository map

| Path | Purpose |
|---|---|
| `CREO.py` | Core CREO implementation |
| `CREO-CEC2017.py` | Original CEC 2017 benchmark driver |
| `experiment_d10.py`, `experiment_d30.py` | Original D=10 and D=30 experiment drivers |
| `ablation_experiment.py` | Matched DE / CREO-DE-NoRG / CREO-DE ablation |
| `cec17_results_D10/`, `cec17_results_D30/` | Archived main-study aggregates and figures |
| `cec17_ablation_small_D10/`, `cec17_ablation_small_D30/` | Archived ablation aggregates and figures |
| `analysis/` | Rebuilds rank and win/tie/loss summaries from archived CSVs |
| `tests/` | Dataset-integrity and repair-model tests |
| `validation/` | Benchmark verification and lightweight smoke-test tools |
| `paper/Information_Sciences/` | Current Information Sciences `cas-sc` manuscript source, PDF, highlights, and checklist |
| `paper/archive/IEEE_previous/` | Superseded IEEE-formatted manuscript retained only for provenance |

## External CEC 2017 dependency

The original scripts expect `cec17_python-master.zip` in the repository root. The upstream implementation referenced during release preparation is:

<https://github.com/lacerdamarcelo/cec17_python>

Do not commit that third-party ZIP to this repository without confirming its licence and redistribution terms. A revision-time copy used for validation had SHA-256:

```text
a2defd319d4c8dbca1411e3e9a470e524e8165b987ea55a5185e687836a65c6c
```

That checksum records the validation dependency only; it does not establish that the same bytes were used in the historical runs.

## Environment

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Place a legally obtained `cec17_python-master.zip` in the repository root before running benchmark or smoke-test scripts.

## Original experimental protocol

- Functions: CEC 2017 F1 and F3--F30 (F2 is excluded)
- Dimensions: 10 and 30
- Independent runs: 30
- Evaluation budget: `10,000 x dimension`
- Recorded seed: 42
- Ablation population size: 50

Run the original studies only when the required compute time and benchmark dependency are available:

```bash
python experiment_d10.py
python experiment_d30.py
```

`ablation_experiment.py` is configured for D=30 in the original snapshot. For D=10, use a separate copy and change `DIMENSIONS` to `[10]`; do not overwrite archived results unintentionally.

## Reanalysis and checks

The archived aggregate CSVs can be checked and summarised without rerunning the optimisers:

```bash
python analysis/rebuild_archived_results.py
python -m unittest discover -s tests -v
```

To verify a local benchmark archive and execute reduced, temporary smoke configurations:

```bash
python validation/verify_benchmark.py
python validation/reproduce_smoke.py
```

Smoke runs confirm that the pipeline executes; they do not reproduce or replace the published 30-run study.

## Data limitation

The historical archive stores aggregate statistics rather than every independent-run outcome. Exact recomputation of distribution-based tests such as Wilcoxon signed-rank tests is therefore not possible from the public archive alone. Existing statistical-output CSVs are retained as historical artifacts; the reanalysis script limits itself to calculations supported by the archived aggregates.

## Manuscript

The current journal package is in `paper/Information_Sciences/`. It targets [Information Sciences](https://www.elsevier.com/journals/information-sciences/0020-0255) and follows the journal's single-anonymized review format. Compile `main.tex` with the bundled `cas-sc` class files, or use the included `main.pdf` for review.

The earlier IEEE-formatted manuscript is archived in `paper/archive/IEEE_previous/` and must not be submitted to Information Sciences.

## Citation and licence

Machine-readable citation metadata is provided in `CITATION.cff`. Cite the tagged paper-submission release rather than a moving branch.

No open-source licence was present in the original snapshot, so no licence has been inferred or added. Copyright remains with the author(s) unless a licence is added explicitly. Third-party benchmark and Elsevier template files remain subject to their own terms.
