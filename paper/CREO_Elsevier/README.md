# CREO Elsevier CAS submission package

This directory contains the submission-formatted manuscript for **Repair-Guided Evolutionary Optimisation (CREO): A Structured Search Paradigm for Continuous Optimisation**.

## Main files

- `main.tex` — Elsevier CAS single-column manuscript source.
- `main.pdf` — reviewed 19-page output: one highlights page and an 18-page article.
- `highlights.txt` — four submission highlights, each within Elsevier's 85-character limit.
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

- Elsevier CAS front matter, affiliation, corresponding-author marker, email, and CRediT roles.
- Consistent `section` → `subsection` hierarchy; no unsupported heading levels are used.
- Elsevier keyword and highlights environments.
- Numbered, compressed citations through `natbib`.
- CAS appendix handling and full-width table formatting.
- Successful compilation with resolved citations and cross-references.

## Before journal submission

Confirm the selected journal's current Guide for Authors, CRediT roles, funding and acknowledgement statements, competing-interest declaration, and ORCID details. Upload `main.pdf` as the manuscript and the complete contents of this directory when the submission system requests LaTeX source files.
