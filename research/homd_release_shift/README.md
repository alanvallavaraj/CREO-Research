# HOMD release-shift CREO experiment

A reference-sequence experiment accompanying [the author manuscript](manuscript.md). **Positive paired result:** full-length 16S assigns oral HMTs more accurately than the V3–V4 region on the same later-release sequences. CREO-DE did not improve the separate novelty-detection outcome over an untuned cosine score.

## Scope and files

| File | Purpose |
|:--|:--|
| homd_data.py | Download and parse the two exact HOMD release pairs |
| fetch_external_data.py | Download and SHA-256-check 16SGOSeq reference FASTA |
| homd_deep_benchmark.py | Temporal, within-release, V3–V4 and genome-derived comparisons |
| homd_equal_budget.py / homd_taxon_balance.py | Additional historical-baseline and taxon-choice sensitivity checks |
| experiment.py | 20-seed, 200-evaluation, four-method optimisation and HMT bootstrap |
| paired_region_experiment.py | Same-ID full-length versus V3–V4 paired comparison and HMT bootstrap |
| supervised_experiment.py | Historical-only supervised confidence-score comparisons |
| make_figures.py | Regenerate three plots from archived results |
| creo_full_experiment_results.json | All fitted weights, objective traces, per-seed metrics, hashes and bootstrap intervals |
| homd_deep_results.json / homd_equal_budget_results.json / homd_taxon_balance_results.json | Secondary-analysis numerical results |
| paired_region_experiment_results.json / supervised_experiment_results.json | Paired positive result and exploratory supervised comparisons |
| references.json | Verified metadata of all 46 manuscript references |
| requirements.txt | Versions used for this run |

The raw public reference files are intentionally excluded. The code requires Python 3.11+ and downloaded public data. Run from this directory:

~~~bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python homd_data.py
python fetch_external_data.py
python homd_deep_benchmark.py
python homd_equal_budget.py
python homd_taxon_balance.py
python experiment.py --seeds 20 --budget 200
python paired_region_experiment.py
python supervised_experiment.py
python make_figures.py
~~~

The full run performs repeated full-length Edlib alignments and can take many minutes. Before interpreting a reproduction, compare SHA-256 of the downloaded HOMD files against the four entries in the archived main result's source_hashes; the external file's hash is checked on download. The external dataset is obtained from https://zenodo.org/records/15209015. The primary HOMD files are from https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/. Inputs and results are archived as of 2026-09-23; later server revisions may change bytes.

## Evaluation safeguards

- All historical v15.23 sequences form the reference pool for the temporal evaluation; every identical full-length sequence is excluded from v16.03 queries before testing.
- Optimisation uses five old-release pseudo-unknown HMT episodes; threshold calibration uses five more old-release episodes. Neither stage uses newer-release HMT labels.
- All optimisers use 200 objective calls per seed and the same simplex repair. CREO without guidance removes only its repair-displacement contribution, leaving attraction and perturbation intact.
- The 20 seeds are search repetitions sharing one temporal biological test, not 20 independent cohorts.
- HMT-cluster bootstrap intervals condition on these releases and fixed fitted weights.
- The in-silico V3–V4 fragments and genome-derived genus data are reference-derived checks. No specimen, assay or clinical diagnosis was evaluated.
- The paired region test fixes query IDs and historical release for both sequence representations. It was designed after the original test was inspected, so the positive result is exploratory and must be replicated externally.

## Main outcomes

On the same 2,236 known-HMT temporal queries, full length was correct in **90.03%** versus **80.23%** with V3–V4. The absolute gain of **9.79 percentage points** has a 95% HMT-cluster-bootstrap interval of **5.84–14.67 points**. The novelty-AUROC gain in this paired subset was +0.0390 with an interval that includes zero.

| Method | Mean later-release AUROC across optimiser seeds |
|:--|--:|
| Untuned 7-mer cosine | **0.8353** |
| Random search | 0.8043 |
| DE | 0.8214 |
| CREO-DE, no repair guidance | 0.8154 |
| CREO-DE, repair guidance | 0.8215 |

The guided-minus-unguided mean seed difference is +0.0061, while guided-minus-raw-cosine is −0.0138. The taxon-cluster bootstrap interval for the representative guided run against raw cosine spans zero (−0.0587 to +0.0407). Additional supervised score models also trailed baseline on temporal novelty AUROC. These results do not justify a novelty-accuracy claim for CREO on this task. See the manuscript for cohort definition and limitations.

**Authorship and submission:** The manuscript is a research draft for the repository owner to check, including affiliation, disclosures and reference formatting; repository storage is not journal publication or peer review.
