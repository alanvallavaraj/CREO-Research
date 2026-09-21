# CREO: single-package instructions

This **one ZIP is complete for local setup**. It combines the original CREO source code, existing historical D10 and D30 main/ablation result directories, revised PDF and LaTeX manuscript, analysis and validation scripts, documentation, and the user-provided CEC2017 benchmark archive.

## Create your new GitHub repository

1. Extract the ZIP and open **CREO-Research/**. This is the project root; put its **contents**, not the outer ZIP or another nested `CREO-Research` directory, at the root of a new GitHub repository (for example `CREO-Research`).
2. `cec17_python-master.zip` is third-party source. **Its redistribution licence was not verified.** Keep it locally and **do not commit or upload it to a public repository** unless you obtain permission. It is already excluded by `.gitignore` for normal `git add` usage. If you use GitHub's web uploader, exclude this file manually because `.gitignore` does not control web uploads. You can make a new public repo with everything else and retain this ZIP locally for reproducibility.
3. In GitHub, choose **New repository**, set its name, and create an empty repository **without** initialising it with another README or .gitignore. On your computer, after extracting the ZIP, run:

   ```bash
   cd CREO-Research
   git init
   git add .
   git status
   git commit -m "Publish CREO research code, archived results and manuscript"
   git branch -M main
   git remote add origin https://github.com/YOUR_ACCOUNT/CREO-Research.git
   git push -u origin main
   ```

   Replace YOUR_ACCOUNT with your GitHub user name. Confirm that `git status` does not include `cec17_python-master.zip` before committing.

## Run locally

The benchmark archive is included **for local use only**. The original scripts look for `cec17_python-master.zip` in this directory.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python validation/verify_benchmark.py
python -m unittest discover -s tests -v
python analysis/rebuild_archived_results.py
```

A C compiler is needed to build the third-party benchmark. The two full main experiments and the two full ablation experiments consume substantial CPU time, may overwrite the archived output directories, and have **not** been fully rerun in the verification reported with this package. Use a separate working copy for full re-execution. See `README.md` and `REPRODUCIBILITY_STATUS.md`.

## Where the key files are

- `paper/CREO_REPOSITORY_PACKAGE.pdf` and `.tex` — revised manuscript.
- `experiment_d10.py`, `experiment_d30.py`, `ablation_experiment.py`, `CREO-CEC2017.py`, `CREO.py` — original research code.
- `cec17_results_D10/`, `cec17_results_D30/`, `cec17_ablation_small_D10/`, `cec17_ablation_small_D30/` — archived experiment results and plots.
- `analysis/`, `validation/`, `tests/` — result checking and limited functional verification.
- `cec17_python-master.zip` — third-party benchmark archive (local; do not publish without permission).
