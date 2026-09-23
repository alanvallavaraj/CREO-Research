# Information Sciences submission release v1.1.0

This release tailors the current manuscript package for **Information Sciences** (Elsevier, ISSN 0020-0255).

## Journal-specific changes

- Moves the current CAS manuscript package to `paper/Information_Sciences/`.
- Removes the superseded IEEE-formatted manuscript so the repository contains only the current Information Sciences submission package.
- Uses a 178-word abstract to satisfy the submission portal's 200-word limit; retains the maximum of six keywords, 3–5 highlights with an 85-character limit, and single-anonymized review format.
- Adds the required generative-AI-use declaration and journal-specific submission checklist.
- Updates the code/data availability statement with the canonical repository URL and its private-at-submission status.

## Earlier v1.0.0 preparation

This release packages the original CREO research snapshot for manuscript submission and reproducible review.

## Added

- Elsevier `cas-sc` manuscript source, class assets, highlights, and compiled PDF.
- A complete repository guide and corrected Python dependency file.
- Machine-readable citation metadata.
- Reproducibility, experiment-execution, and third-party benchmark provenance notes.
- Aggregate-data reanalysis for mean ranks and matched-ablation win/tie/loss counts.
- Tests for archive structure, reported aggregate comparisons, and repair behaviour.
- Benchmark verification and isolated reduced smoke-test tooling.

## Preserved

- Original snapshot: `5e04fb2be03071ef601905c2dd3d04390c6f90b3`.
- All five original Python research scripts.
- All 219 archived result files checked in the four experiment-result directories.
- No historical numerical output was regenerated or replaced.

## Validation result

- Six deterministic tests passed.
- Four aggregate summary CSVs were rebuilt successfully.
- 224 preservation-scope files were byte-identical to the original source archive.
- The manuscript compiled successfully to a 19-page PDF package (one highlights page plus an 18-page article).

## Important qualification

The full 30-run experiments were not rerun. The exact benchmark ZIP used historically cannot be proven from the old snapshot, and per-run values needed to regenerate distribution-based tests are not present. See `REPRODUCIBILITY_STATUS.md` for the precise boundary of verification.
