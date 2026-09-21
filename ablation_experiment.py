import os
import sys
import copy
import random
import zipfile
import shutil
import subprocess
from dataclasses import dataclass
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ZIP_PATH = os.path.join(BASE_DIR, "cec17_python-master.zip")
EXTRACT_ROOT = os.path.join(BASE_DIR, "cec17_official_ablation")

DIMENSIONS = [30]
RUNS = 30
POP_SIZE = 50
FUNCTION_IDS = [fid for fid in range(1, 31) if fid != 2]
MAX_EVALS_FACTOR = 10000

ALGORITHMS = [
    "DE",
    "CREO-DE-NoRG",
    "CREO-DE",
]

SEED = 42
LOWER_BOUND = -100.0
UPPER_BOUND = 100.0
SAVE_ALL_CURVES = False

SUPPORTED_DIMENSIONS = {2, 10, 20, 30, 50, 100}
for d in DIMENSIONS:
    if d not in SUPPORTED_DIMENSIONS:
        raise ValueError(f"Unsupported dimension {d}. Supported: {sorted(SUPPORTED_DIMENSIONS)}")

DIMENSION_TAG = "_".join(str(d) for d in DIMENSIONS)
RESULTS_DIR = os.path.join(BASE_DIR, f"cec17_ablation_small_D{DIMENSION_TAG}")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
PLOTS_DIR = os.path.join(RESULTS_DIR, "plots")
CURVES_DIR = os.path.join(RESULTS_DIR, "curves")

random.seed(SEED)
np.random.seed(SEED)

# ============================================================
# PREPARE DIRECTORIES
# ============================================================

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(CURVES_DIR, exist_ok=True)

# ============================================================
# STEP 1: EXTRACT ZIP
# ============================================================

if not os.path.exists(ZIP_PATH):
    raise FileNotFoundError(
        f"Could not find zip file at: {ZIP_PATH}\n"
        "Make sure cec17_python-master.zip is in the same folder as this script."
    )

if os.path.exists(EXTRACT_ROOT):
    shutil.rmtree(EXTRACT_ROOT)
os.makedirs(EXTRACT_ROOT, exist_ok=True)

with zipfile.ZipFile(ZIP_PATH, "r") as zf:
    zf.extractall(EXTRACT_ROOT)

def find_package_dir(root: str) -> str:
    for current_root, _dirs, files in os.walk(root):
        if "cec17_functions.py" in files and "cec17_test_func.c" in files:
            return current_root
    raise FileNotFoundError(
        "Could not locate extracted CEC package folder containing "
        "'cec17_functions.py' and 'cec17_test_func.c'."
    )

PACKAGE_DIR = find_package_dir(EXTRACT_ROOT)
print("Extracted CEC package to:", PACKAGE_DIR)
print("Benchmark functions to run:", FUNCTION_IDS)

# ============================================================
# STEP 2: PATCH C FILE FOR MACOS + COMPILE SHARED LIBRARY
# ============================================================

so_path = os.path.join(PACKAGE_DIR, "cec17_test_func.so")
c_path = os.path.join(PACKAGE_DIR, "cec17_test_func.c")

if not os.path.exists(c_path):
    raise FileNotFoundError(f"Could not find C source file: {c_path}")

if sys.platform == "darwin":
    with open(c_path, "r", encoding="utf-8") as f:
        c_code = f.read()
    if "#include <malloc.h>" in c_code:
        c_code = c_code.replace("#include <malloc.h>", "")
        with open(c_path, "w", encoding="utf-8") as f:
            f.write(c_code)
        print("Patched cec17_test_func.c for macOS: removed <malloc.h>")

if sys.platform == "darwin":
    compile_cmd = [
        "clang",
        "-dynamiclib",
        "-O2",
        "-o", so_path,
        c_path,
        "-lm"
    ]
else:
    compile_cmd = [
        "gcc",
        "-shared",
        "-fPIC",
        "-O2",
        "-o", so_path,
        c_path,
        "-lm"
    ]

print("Compiling shared library...")
print("Compile command:", " ".join(compile_cmd))

result = subprocess.run(compile_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("Compilation failed.")
    print("STDOUT:")
    print(result.stdout)
    print("STDERR:")
    print(result.stderr)
    raise RuntimeError("C compilation failed. See compiler error above.")

print("Compiled:", so_path)

# ============================================================
# STEP 3: IMPORT WRAPPER
# ============================================================

os.chdir(PACKAGE_DIR)
sys.path.insert(0, PACKAGE_DIR)

from cec17_functions import cec17_test_func  # noqa: E402

print("Imported cec17_test_func successfully.")

# ============================================================
# CEC EVALUATION WRAPPER
# ============================================================

def cec_eval(x: np.ndarray, func_id: int) -> float:
    x = np.asarray(x, dtype=float)
    nx = len(x)
    mx = 1
    f = np.zeros(mx, dtype=float)
    cec17_test_func(x, f, nx, mx, func_id)
    return float(f[0])

# ============================================================
# REPAIR / FEASIBILITY MODEL
# ============================================================

@dataclass
class EvalResult:
    x_raw: np.ndarray
    x_used: np.ndarray
    repair_displacement: np.ndarray
    objective: float
    violation: float
    feasible: int
    fitness: float

def clip_to_bounds(x: np.ndarray) -> np.ndarray:
    return np.clip(x, LOWER_BOUND, UPPER_BOUND)

def total_violation(x: np.ndarray) -> float:
    v = 0.0
    soft_low, soft_high = -80.0, 80.0

    below = np.maximum(0.0, soft_low - x)
    above = np.maximum(0.0, x - soft_high)
    v += float(np.sum(below + above))

    mean_abs = np.mean(np.abs(x))
    if mean_abs > 60.0:
        v += float((mean_abs - 60.0) * len(x))

    return float(v)

def repair_solution(x: np.ndarray) -> np.ndarray:
    y = clip_to_bounds(x.copy())
    if total_violation(y) <= 1e-12:
        return y

    for _ in range(10):
        y = 0.80 * y
        y = clip_to_bounds(y)
        if total_violation(y) <= 1e-12:
            break
    return y

def evaluate_candidate(x: np.ndarray, func_id: int, mode: str, penalty_weight: float = 1e6) -> EvalResult:
    x_raw = clip_to_bounds(x.copy())
    x_rep = repair_solution(x_raw)
    disp = x_rep - x_raw

    x_used = x_raw if mode == "penalty" else x_rep
    obj = cec_eval(x_used, func_id)

    raw_violation = total_violation(x_raw)
    used_violation = total_violation(x_used)

    if mode == "penalty":
        fitness = obj + penalty_weight * raw_violation
        feasible = int(raw_violation <= 1e-12)
        violation = raw_violation
    else:
        fitness = obj
        feasible = int(used_violation <= 1e-12)
        violation = used_violation

    return EvalResult(
        x_raw=x_raw,
        x_used=x_used,
        repair_displacement=disp,
        objective=float(obj),
        violation=float(violation),
        feasible=feasible,
        fitness=float(fitness),
    )

def better(a: EvalResult, b: EvalResult) -> bool:
    if a.feasible != b.feasible:
        return a.feasible > b.feasible
    if a.feasible == 1:
        return a.fitness <= b.fitness
    if a.violation != b.violation:
        return a.violation < b.violation
    return a.fitness <= b.fitness

# ============================================================
# POPULATION HELPERS
# ============================================================

def init_population(pop_size: int, dim: int) -> np.ndarray:
    return np.random.uniform(LOWER_BOUND, UPPER_BOUND, size=(pop_size, dim))

def pad_curves(curves: List[np.ndarray]) -> np.ndarray:
    m = max(len(c) for c in curves)
    arr = np.full((len(curves), m), np.nan)
    for i, c in enumerate(curves):
        arr[i, :len(c)] = c
        if len(c) < m:
            arr[i, len(c):] = c[-1]
    return arr

def plot_convergence(curve_dict: Dict[str, np.ndarray], title: str, save_path: str):
    plt.figure(figsize=(10, 6))
    x = np.arange(1, len(next(iter(curve_dict.values()))) + 1)
    for algo, curve in curve_dict.items():
        plt.plot(x, curve, label=algo)
    plt.xlabel("Generation / Iteration")
    plt.ylabel("Best-so-far objective")
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

# ============================================================
# ALGORITHMS
# ============================================================

def run_de(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    curve = [best.objective]

    F, CR = 0.7, 0.9

    while nfe < max_evals:
        new_pop, new_evals = [], []

        for i in range(POP_SIZE):
            idxs = list(range(POP_SIZE))
            idxs.remove(i)
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)

            mutant = pop[r1] + F * (pop[r2] - pop[r3])
            mutant = clip_to_bounds(mutant)

            trial = pop[i].copy()
            jrand = np.random.randint(dim)
            for j in range(dim):
                if np.random.rand() < CR or j == jrand:
                    trial[j] = mutant[j]
            trial = clip_to_bounds(trial)

            trial_eval = evaluate_candidate(trial, func_id, mode="penalty")
            nfe += 1

            if better(trial_eval, evals[i]):
                new_pop.append(trial)
                new_evals.append(trial_eval)
            else:
                new_pop.append(pop[i].copy())
                new_evals.append(evals[i])

            if nfe >= max_evals:
                break

        if new_pop:
            pop[:len(new_pop)] = np.array(new_pop)
            evals[:len(new_evals)] = new_evals

        cur_best = min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness))
        if better(cur_best, best):
            best = copy.deepcopy(cur_best)
        curve.append(best.objective)

    return best, np.array(curve)

def run_creo_de_no_rg(func_id: int, dim: int, max_evals: int):
    """
    Ablation variant:
      - repair-first evaluation
      - DE mutation/crossover
      - best repaired-solution attraction
      - NO repair-guidance term
    """
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="repair") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    gbest = best.x_used.copy()

    F, CR = 0.7, 0.9
    alpha, beta, gamma = 0.18, 0.0, 0.02
    curve = [best.objective]

    while nfe < max_evals:
        new_pop, new_evals = [], []

        for i in range(POP_SIZE):
            idxs = list(range(POP_SIZE))
            idxs.remove(i)
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)

            delta_r = evals[i].repair_displacement

            mutant = (
                pop[r1]
                + F * (pop[r2] - pop[r3])
                + alpha * (gbest - pop[i])
                + beta * delta_r
                + gamma * np.random.randn(dim)
            )
            mutant = clip_to_bounds(mutant)

            trial = pop[i].copy()
            jrand = np.random.randint(dim)
            for j in range(dim):
                if np.random.rand() < CR or j == jrand:
                    trial[j] = mutant[j]
            trial = clip_to_bounds(trial)

            trial_eval = evaluate_candidate(trial, func_id, mode="repair")
            nfe += 1

            if better(trial_eval, evals[i]):
                new_pop.append(trial)
                new_evals.append(trial_eval)
            else:
                new_pop.append(pop[i].copy())
                new_evals.append(evals[i])

            if nfe >= max_evals:
                break

        if new_pop:
            pop[:len(new_pop)] = np.array(new_pop)
            evals[:len(new_evals)] = new_evals

        cur_best = min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness))
        if better(cur_best, best):
            best = copy.deepcopy(cur_best)
            gbest = best.x_used.copy()

        curve.append(best.objective)

    return best, np.array(curve)

def run_creo_de(func_id: int, dim: int, max_evals: int):
    """
    Full CREO-DE:
      - repair-first evaluation
      - DE mutation/crossover
      - best repaired-solution attraction
      - repair-guidance term enabled
    """
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="repair") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    gbest = best.x_used.copy()

    F, CR = 0.7, 0.9
    alpha, beta, gamma = 0.18, 0.90, 0.02
    curve = [best.objective]

    while nfe < max_evals:
        new_pop, new_evals = [], []

        for i in range(POP_SIZE):
            idxs = list(range(POP_SIZE))
            idxs.remove(i)
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)

            delta_r = evals[i].repair_displacement

            mutant = (
                pop[r1]
                + F * (pop[r2] - pop[r3])
                + alpha * (gbest - pop[i])
                + beta * delta_r
                + gamma * np.random.randn(dim)
            )
            mutant = clip_to_bounds(mutant)

            trial = pop[i].copy()
            jrand = np.random.randint(dim)
            for j in range(dim):
                if np.random.rand() < CR or j == jrand:
                    trial[j] = mutant[j]
            trial = clip_to_bounds(trial)

            trial_eval = evaluate_candidate(trial, func_id, mode="repair")
            nfe += 1

            if better(trial_eval, evals[i]):
                new_pop.append(trial)
                new_evals.append(trial_eval)
            else:
                new_pop.append(pop[i].copy())
                new_evals.append(evals[i])

            if nfe >= max_evals:
                break

        if new_pop:
            pop[:len(new_pop)] = np.array(new_pop)
            evals[:len(new_evals)] = new_evals

        cur_best = min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness))
        if better(cur_best, best):
            best = copy.deepcopy(cur_best)
            gbest = best.x_used.copy()

        curve.append(best.objective)

    return best, np.array(curve)

# ============================================================
# DRIVER
# ============================================================

def run_algorithm(name: str, func_id: int, dim: int, max_evals: int):
    if name == "DE":
        return run_de(func_id, dim, max_evals)
    if name == "CREO-DE-NoRG":
        return run_creo_de_no_rg(func_id, dim, max_evals)
    if name == "CREO-DE":
        return run_creo_de(func_id, dim, max_evals)
    raise ValueError(f"Unknown algorithm: {name}")

# ============================================================
# REPORTING HELPERS
# ============================================================

def summarise_runs(run_results: List[EvalResult]) -> Dict[str, float]:
    objs = np.array([r.objective for r in run_results], dtype=float)
    std_val = float(np.std(objs, ddof=1)) if len(objs) > 1 else 0.0
    return {
        "Mean": float(np.mean(objs)),
        "Std": std_val,
        "Best": float(np.min(objs)),
        "Worst": float(np.max(objs)),
    }

def rank_algorithms_from_means(df: pd.DataFrame) -> pd.DataFrame:
    """
    Input df columns expected:
      Function, Dimension, Algorithm, Mean
    Returns aggregate rank table per dimension.
    """
    rows = []
    for dim in sorted(df["Dimension"].unique()):
        sub = df[df["Dimension"] == dim].copy()

        # per function ranking (lower mean is better)
        best_counts = {a: 0 for a in ALGORITHMS}
        worst_counts = {a: 0 for a in ALGORITHMS}
        algo_ranks = {a: [] for a in ALGORITHMS}

        for func in sorted(sub["Function"].unique(), key=lambda x: int(x[1:])):
            sf = sub[sub["Function"] == func].copy()
            sf = sf.sort_values("Mean", ascending=True).reset_index(drop=True)
            sf["Rank"] = np.arange(1, len(sf) + 1)

            best_algo = sf.iloc[0]["Algorithm"]
            worst_algo = sf.iloc[-1]["Algorithm"]
            best_counts[best_algo] += 1
            worst_counts[worst_algo] += 1

            for _, row in sf.iterrows():
                algo_ranks[row["Algorithm"]].append(int(row["Rank"]))

        for algo in ALGORITHMS:
            ranks = np.array(algo_ranks[algo], dtype=float)
            rows.append({
                "Dimension": dim,
                "Algorithm": algo,
                "MeanRank": float(np.mean(ranks)),
                "StdRank": float(np.std(ranks, ddof=1)) if len(ranks) > 1 else 0.0,
                "BestCount": int(best_counts[algo]),
                "WorstCount": int(worst_counts[algo]),
            })

    out = pd.DataFrame(rows)
    out = out.sort_values(["Dimension", "MeanRank", "StdRank"]).reset_index(drop=True)
    return out

# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():
    master_rows = []

    for dim in DIMENSIONS:
        max_evals = MAX_EVALS_FACTOR * dim
        print(f"\n========== DIMENSION {dim} | MAX_EVALS {max_evals} ==========\n")

        for func_id in FUNCTION_IDS:
            print(f"Running F{func_id} | D={dim}")

            per_algo_results: Dict[str, List[EvalResult]] = {a: [] for a in ALGORITHMS}
            per_algo_curves: Dict[str, List[np.ndarray]] = {a: [] for a in ALGORITHMS}

            for algo in ALGORITHMS:
                print(f"  {algo}")
                for run in range(1, RUNS + 1):
                    best, curve = run_algorithm(algo, func_id, dim, max_evals)
                    per_algo_results[algo].append(best)
                    per_algo_curves[algo].append(curve)

                    if SAVE_ALL_CURVES:
                        np.save(
                            os.path.join(CURVES_DIR, f"F{func_id}_D{dim}_{algo}_run{run:02d}.npy"),
                            curve
                        )

            # per-function summary CSV
            table_rows = []
            curve_means = {}

            for algo in ALGORITHMS:
                summary = summarise_runs(per_algo_results[algo])
                row = {
                    "Function": f"F{func_id}",
                    "Dimension": dim,
                    "Algorithm": algo,
                    **summary
                }
                table_rows.append(row)
                master_rows.append(row)

                curve_means[algo] = np.nanmean(pad_curves(per_algo_curves[algo]), axis=0)

            df_table = pd.DataFrame(table_rows).sort_values(["Mean", "Std"])
            df_table.to_csv(
                os.path.join(TABLES_DIR, f"F{func_id}_D{dim}_ablation_summary.csv"),
                index=False
            )

            plot_convergence(
                curve_means,
                title=f"Ablation | CEC2017 F{func_id} | D={dim}",
                save_path=os.path.join(PLOTS_DIR, f"F{func_id}_D{dim}_ablation_convergence.png")
            )

    master_df = pd.DataFrame(master_rows)
    master_path = os.path.join(RESULTS_DIR, "all_ablation_results.csv")
    master_df.to_csv(master_path, index=False)

    rank_df = rank_algorithms_from_means(master_df)
    rank_path = os.path.join(RESULTS_DIR, "ablation_rank_summary.csv")
    rank_df.to_csv(rank_path, index=False)

    print("\nSaved outputs:")
    print(master_path)
    print(rank_path)
    print(TABLES_DIR)
    print(PLOTS_DIR)

if __name__ == "__main__":
    main()