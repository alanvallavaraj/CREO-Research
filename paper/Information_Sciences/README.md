# CREO submission package for Information Sciences

This directory contains the submission-formatted manuscript for **Repair-Guided Evolutionary Optimisation (CREO): A Structured Search Paradigm for Continuous Optimisation**, prepared for **Information Sciences** (Elsevier, ISSN 0020-0255).

Journal guide: <https://www.elsevier.com/journals/information-sciences/0020-0255/guide-for-authors>

## Main files

- `main.tex` — Elsevier CAS single-column manuscript source.
- `main.pdf` — compiled manuscript PDF, visually checked after the journal-specific update.
- `highlights.txt` — four submission highlights, each within Elsevier's 85-character limit.
- `SUBMISSION_CHECKLIST.md` — journal-specific format checks and author confirmations still required.
- `cas-sc.cls`, `cas-common.sty`, and `cas-model2-names.bst` — required template files.
- `cas-*.jpeg` — CAS front-matter icons required by the bundled class.

## Compile

Keep all files in this directory and run:

```bash
pdflatex main.tex
pdflatex main.tex
```

No BibTeX run is required because `main.tex` contains a complete `thebibliography` environment.

## Formatting checked

- Information Sciences is a single-anonymized journal, so author details remain in the manuscript.
- The abstract is 243 words, within the journal's 250-word maximum.
- Six English keywords are supplied, within the journal's maximum.
- Four highlights are supplied; each is no more than 85 characters including spaces.
- Elsevier CAS front matter, affiliation, corresponding-author marker, email, and CRediT roles.
- Consistent `section` → `subsection` hierarchy; no unsupported heading levels are used.
- Elsevier keyword and highlights environments.
- Numbered, compressed citations through `natbib`.
- CAS appendix handling and full-width table formatting.
- A generative-AI-use declaration immediately precedes the reference list.
- Successful compilation with resolved citations and cross-references.

## Before journal submission

Complete the author confirmations in `SUBMISSION_CHECKLIST.md`, especially funding, competing interests, ORCID, and repository access. Upload `main.pdf` as the manuscript, upload `highlights.txt` separately as Highlights, and upload the complete source package when Editorial Manager requests LaTeX source files.
