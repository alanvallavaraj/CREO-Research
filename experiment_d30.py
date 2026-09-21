# ============================================================
# OFFICIAL-STYLE CEC2017 + CREO-DE EXPERIMENT (VS CODE VERSION)
# ============================================================
#
# Expected project folder:
#   cec_project/
#       experiment.py
#       cec17_python-master.zip
#
# Outputs:
#   cec_project/cec17_results/
#       all_benchmark_results.csv
#       overall_algorithm_summary.csv
#       wilcoxon_creo_de_vs_others.csv
#       friedman_ranks.csv
#       overall_friedman_mean_ranks.csv
#       tables/
#       plots/
#       curves/
#
# Notes:
# - Uses the local cec17_python-master.zip package
# - Compiles cec17_test_func.c into cec17_test_func.so
# - Runs F1-F30 of CEC2017
# - Adds a repair-oriented search protocol so CREO variants remain meaningful
#   while objective values are always computed by the official CEC evaluator
# ============================================================

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
from scipy.stats import wilcoxon, friedmanchisquare

# ============================================================
# USER SETTINGS
# ============================================================

# Assumes the zip is in the same folder as this script
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

ZIP_PATH = os.path.join(BASE_DIR, "cec17_python-master.zip")
EXTRACT_ROOT = os.path.join(BASE_DIR, "cec17_official")
PACKAGE_DIR = os.path.join(EXTRACT_ROOT, "cec17_python-master")

RESULTS_DIR = os.path.join(BASE_DIR, "cec17_results_D30")
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
PLOTS_DIR = os.path.join(RESULTS_DIR, "plots")
CURVES_DIR = os.path.join(RESULTS_DIR, "curves")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

# ------------------------------------------------------------
# QUICK TEST SETTINGS (recommended first)
# Uncomment these and comment out the full experiment settings
# if you want a quick sanity check.
# ------------------------------------------------------------
# DIMENSIONS = [10]
# RUNS = 2
# POP_SIZE = 20
# FUNCTION_IDS = [1]
# MAX_EVALS_FACTOR = 1000

# ------------------------------------------------------------
# FULL EXPERIMENT SETTINGS
# ------------------------------------------------------------
DIMENSIONS = [30]

SUPPORTED_DIMENSIONS = {2, 10, 20, 30, 50, 100}
for d in DIMENSIONS:
    if d not in SUPPORTED_DIMENSIONS:
        raise ValueError(f"Unsupported dimension {d}")
    
RUNS = 30
POP_SIZE = 50
DELETED_FUNCTION_IDS = {2}
FUNCTION_IDS = [fid for fid in range(1, 31) if fid not in DELETED_FUNCTION_IDS]
MAX_EVALS_FACTOR = 10000            # standard budget = 10000 * D
print("Benchmark functions:", FUNCTION_IDS)


ALGORITHMS = [
    "GA",
    "PSO",
    "DE",
    "GWO",
    "CREO-Repair",
    "CREO-RG",
    "CREO-DE",
]

SAVE_ALL_CURVES = False

# Box bounds commonly used for CEC2017
LOWER_BOUND = -100.0
UPPER_BOUND = 100.0

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

extract_root_abs = os.path.abspath(EXTRACT_ROOT)
package_dir_abs = os.path.abspath(PACKAGE_DIR)

if not os.path.exists(ZIP_PATH):
    raise FileNotFoundError(
        f"Could not find zip file at: {ZIP_PATH}\n"
        f"Make sure cec17_python-master.zip is in the same folder as experiment.py"
    )

if os.path.exists(extract_root_abs):
    shutil.rmtree(extract_root_abs)

os.makedirs(extract_root_abs, exist_ok=True)

with zipfile.ZipFile(ZIP_PATH, "r") as zf:
    zf.extractall(extract_root_abs)

if not os.path.exists(package_dir_abs):
    extracted_items = os.listdir(extract_root_abs)
    raise FileNotFoundError(
        f"Expected package dir not found: {package_dir_abs}\n"
        f"Extracted contents found: {extracted_items}\n"
        f"If the folder name is different, update PACKAGE_DIR in the script."
    )

print("Extracted CEC2017 package to:", package_dir_abs)

# ============================================================
# STEP 2: PATCH C FILE FOR MACOS + COMPILE SHARED LIBRARY
# ============================================================

so_path = os.path.join(package_dir_abs, "cec17_test_func.so")
c_path = os.path.join(package_dir_abs, "cec17_test_func.c")

if not os.path.exists(c_path):
    raise FileNotFoundError(f"Could not find C source file: {c_path}")

# Patch Linux-only malloc.h for macOS
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

result = subprocess.run(
    compile_cmd,
    capture_output=True,
    text=True
)

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
# IMPORTANT: run from package dir so input_data paths work
# ============================================================

os.chdir(package_dir_abs)
sys.path.insert(0, package_dir_abs)

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
# CEC2017 is unconstrained. To keep CREO meaningful as a
# repair-centric framework, we use:
#   - box-feasibility
#   - a soft stability region around the search domain
#
# This does NOT change the official CEC objective evaluator.
# The objective values are always computed by the official CEC code.
# The soft-violation model only influences repair and selection.
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

def tournament_select(pop: np.ndarray, evals: List[EvalResult], k: int = 3) -> np.ndarray:
    idx = np.random.choice(len(pop), size=k, replace=False)
    best_idx = idx[0]
    for j in idx[1:]:
        if better(evals[j], evals[best_idx]):
            best_idx = j
    return pop[best_idx].copy()

def gaussian_mutation(x: np.ndarray, sigma: float = 0.05) -> np.ndarray:
    span = UPPER_BOUND - LOWER_BOUND
    y = x + np.random.randn(len(x)) * sigma * span
    return clip_to_bounds(y)

def one_point_crossover(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    cut = np.random.randint(1, len(a))
    c1 = np.concatenate([a[:cut], b[cut:]])
    c2 = np.concatenate([b[:cut], a[cut:]])
    return c1, c2

# ============================================================
# ALGORITHMS
# ============================================================

def run_ga(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    curve = [best.objective]

    while nfe < max_evals:
        next_pop = []
        elite_idx = min(
            range(len(evals)),
            key=lambda i: (1 - evals[i].feasible, evals[i].violation, evals[i].fitness)
        )
        next_pop.append(pop[elite_idx].copy())

        while len(next_pop) < POP_SIZE:
            p1 = tournament_select(pop, evals)
            p2 = tournament_select(pop, evals)
            c1, c2 = one_point_crossover(p1, p2)
            c1 = gaussian_mutation(c1, sigma=0.05)
            c2 = gaussian_mutation(c2, sigma=0.05)
            next_pop.extend([c1, c2])

        pop = np.array(next_pop[:POP_SIZE])
        evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
        nfe += POP_SIZE

        cur_best = min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness))
        if better(cur_best, best):
            best = copy.deepcopy(cur_best)
        curve.append(best.objective)

    return best, np.array(curve)

def run_pso(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    vel = np.zeros_like(pop)
    vmax = 0.15 * (UPPER_BOUND - LOWER_BOUND)

    evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
    pbest = pop.copy()
    pbest_eval = copy.deepcopy(evals)

    gbest_idx = min(
        range(POP_SIZE),
        key=lambda i: (1 - evals[i].feasible, evals[i].violation, evals[i].fitness)
    )
    gbest = pop[gbest_idx].copy()
    gbest_eval = copy.deepcopy(evals[gbest_idx])

    nfe = POP_SIZE
    curve = [gbest_eval.objective]

    w, c1, c2 = 0.7, 1.5, 1.5

    while nfe < max_evals:
        r1 = np.random.rand(*pop.shape)
        r2 = np.random.rand(*pop.shape)
        vel = w * vel + c1 * r1 * (pbest - pop) + c2 * r2 * (gbest - pop)
        vel = np.clip(vel, -vmax, vmax)
        pop = clip_to_bounds(pop + vel)

        evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
        nfe += POP_SIZE

        for i in range(POP_SIZE):
            if better(evals[i], pbest_eval[i]):
                pbest[i] = pop[i].copy()
                pbest_eval[i] = copy.deepcopy(evals[i])

        best_idx = min(
            range(POP_SIZE),
            key=lambda i: (1 - pbest_eval[i].feasible, pbest_eval[i].violation, pbest_eval[i].fitness)
        )
        if better(pbest_eval[best_idx], gbest_eval):
            gbest = pbest[best_idx].copy()
            gbest_eval = copy.deepcopy(pbest_eval[best_idx])

        curve.append(gbest_eval.objective)

    return gbest_eval, np.array(curve)

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

def run_gwo(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    curve = [best.objective]
    t = 0

    while nfe < max_evals:
        order = sorted(
            range(POP_SIZE),
            key=lambda i: (1 - evals[i].feasible, evals[i].violation, evals[i].fitness)
        )
        alpha = pop[order[0]].copy()
        beta = pop[order[1]].copy()
        delta = pop[order[2]].copy()

        a = 2 - 2 * (t / max(1, max_evals // POP_SIZE))
        t += 1

        new_pop = []
        for i in range(POP_SIZE):
            x = pop[i].copy()

            A1 = 2 * a * np.random.rand(dim) - a
            C1 = 2 * np.random.rand(dim)
            D_alpha = np.abs(C1 * alpha - x)
            X1 = alpha - A1 * D_alpha

            A2 = 2 * a * np.random.rand(dim) - a
            C2 = 2 * np.random.rand(dim)
            D_beta = np.abs(C2 * beta - x)
            X2 = beta - A2 * D_beta

            A3 = 2 * a * np.random.rand(dim) - a
            C3 = 2 * np.random.rand(dim)
            D_delta = np.abs(C3 * delta - x)
            X3 = delta - A3 * D_delta

            x_new = (X1 + X2 + X3) / 3.0
            x_new = clip_to_bounds(x_new)
            new_pop.append(x_new)

        pop = np.array(new_pop)
        evals = [evaluate_candidate(ind, func_id, mode="penalty") for ind in pop]
        nfe += POP_SIZE

        cur_best = min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness))
        if better(cur_best, best):
            best = copy.deepcopy(cur_best)
        curve.append(best.objective)

    return best, np.array(curve)

def run_creo_repair(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="repair") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    gbest = best.x_used.copy()

    alpha, gamma = 0.45, 0.05
    curve = [best.objective]

    while nfe < max_evals:
        new_pop, new_evals = [], []
        for i in range(POP_SIZE):
            eps = np.random.randn(dim)
            x_new = pop[i] + alpha * (gbest - pop[i]) + gamma * eps
            x_new = clip_to_bounds(x_new)

            e_new = evaluate_candidate(x_new, func_id, mode="repair")
            nfe += 1
            new_pop.append(x_new)
            new_evals.append(e_new)

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

def run_creo_rg(func_id: int, dim: int, max_evals: int):
    pop = init_population(POP_SIZE, dim)
    evals = [evaluate_candidate(ind, func_id, mode="repair") for ind in pop]
    nfe = POP_SIZE

    best = copy.deepcopy(min(evals, key=lambda e: (1 - e.feasible, e.violation, e.fitness)))
    gbest = best.x_used.copy()

    alpha, beta, gamma = 0.35, 0.85, 0.04
    curve = [best.objective]

    while nfe < max_evals:
        new_pop, new_evals = [], []
        for i in range(POP_SIZE):
            delta_r = evals[i].repair_displacement
            eps = np.random.randn(dim)

            x_new = pop[i] + alpha * (gbest - pop[i]) + beta * delta_r + gamma * eps
            x_new = clip_to_bounds(x_new)

            e_new = evaluate_candidate(x_new, func_id, mode="repair")
            nfe += 1
            new_pop.append(x_new)
            new_evals.append(e_new)

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
    Hybrid CREO-DE:
      - DE mutation/crossover
      - repair-guided displacement
      - repaired-best attraction
      - feasibility-aware selection
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
    if name == "GA":
        return run_ga(func_id, dim, max_evals)
    if name == "PSO":
        return run_pso(func_id, dim, max_evals)
    if name == "DE":
        return run_de(func_id, dim, max_evals)
    if name == "GWO":
        return run_gwo(func_id, dim, max_evals)
    if name == "CREO-Repair":
        return run_creo_repair(func_id, dim, max_evals)
    if name == "CREO-RG":
        return run_creo_rg(func_id, dim, max_evals)
    if name == "CREO-DE":
        return run_creo_de(func_id, dim, max_evals)
    raise ValueError(f"Unknown algorithm: {name}")

# ============================================================
# REPORTING HELPERS
# ============================================================

def summarise_runs(run_results: List[EvalResult]) -> Dict[str, float]:
    objs = np.array([r.objective for r in run_results], dtype=float)
    feas = np.array([r.feasible for r in run_results], dtype=float)
    viol = np.array([r.violation for r in run_results], dtype=float)

    std_val = float(np.std(objs, ddof=1)) if len(objs) > 1 else 0.0

    return {
        "Mean": float(np.mean(objs)),
        "Std": std_val,
        "Best": float(np.min(objs)),
        "Worst": float(np.max(objs)),
        "FeasibleRate(%)": float(np.mean(feas) * 100.0),
        "MeanViolation": float(np.mean(viol)),
    }

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

def friedman_rank_table(score_matrix: Dict[str, List[float]]) -> pd.DataFrame:
    algos = list(score_matrix.keys())
    arrays = [np.array(score_matrix[a], dtype=float) for a in algos]

    friedmanchisquare(*arrays)

    rows = np.vstack(arrays).T
    ranks_all = []
    for row in rows:
        order = np.argsort(row)   # lower is better
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(row) + 1)
        ranks_all.append(ranks)

    mean_ranks = np.mean(np.array(ranks_all), axis=0)
    df = pd.DataFrame({"Algorithm": algos, "Mean Rank": mean_ranks})
    return df.sort_values("Mean Rank", ascending=True).reset_index(drop=True)

# ============================================================
# MAIN EXPERIMENT
# ============================================================

def main():
    master_rows = []
    pairwise_records = []
    friedman_records = []

    for dim in DIMENSIONS:
        max_evals = MAX_EVALS_FACTOR * dim
        print(f"\n========== DIMENSION {dim} | MAX_EVALS {max_evals} ==========\n")

        for func_id in FUNCTION_IDS:

            # Skip deleted CEC functions
            if func_id == 2:
                print("Skipping deleted function F2")
                continue

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

            table_rows = []
            score_matrix = {}
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

                score_matrix[algo] = [r.objective for r in per_algo_results[algo]]
                curve_means[algo] = np.nanmean(pad_curves(per_algo_curves[algo]), axis=0)

            df_table = pd.DataFrame(table_rows).sort_values(["Mean", "Std"])
            df_table.to_csv(
                os.path.join(TABLES_DIR, f"F{func_id}_D{dim}_summary.csv"),
                index=False
            )

            plot_convergence(
                curve_means,
                title=f"CEC2017 F{func_id} | D={dim}",
                save_path=os.path.join(PLOTS_DIR, f"F{func_id}_D{dim}_convergence.png")
            )

            creo_scores = np.array(score_matrix["CREO-DE"], dtype=float)
            for algo in ALGORITHMS:
                if algo == "CREO-DE":
                    continue
                other_scores = np.array(score_matrix[algo], dtype=float)
                try:
                    stat, p = wilcoxon(
                        creo_scores,
                        other_scores,
                        zero_method="wilcox",
                        alternative="two-sided"
                    )
                except ValueError:
                    stat, p = np.nan, np.nan

                pairwise_records.append({
                    "Function": f"F{func_id}",
                    "Dimension": dim,
                    "Comparison": f"CREO-DE vs {algo}",
                    "Statistic": stat,
                    "p-value": p,
                    "Significant(p<0.05)": int(pd.notna(p) and p < 0.05),
                    "CREO-DE Mean": float(np.mean(creo_scores)),
                    f"{algo} Mean": float(np.mean(other_scores)),
                })

            friedman_df = friedman_rank_table(score_matrix)
            friedman_df.insert(0, "Dimension", dim)
            friedman_df.insert(0, "Function", f"F{func_id}")
            friedman_records.append(friedman_df)

    master_df = pd.DataFrame(master_rows)
    master_df.to_csv(os.path.join(RESULTS_DIR, "all_benchmark_results.csv"), index=False)

    wilcoxon_df = pd.DataFrame(pairwise_records)
    wilcoxon_df.to_csv(os.path.join(RESULTS_DIR, "wilcoxon_creo_de_vs_others.csv"), index=False)

    friedman_all_df = pd.concat(friedman_records, axis=0, ignore_index=True)
    friedman_all_df.to_csv(os.path.join(RESULTS_DIR, "friedman_ranks.csv"), index=False)

    overall_df = (
        master_df.groupby("Algorithm", as_index=False)
        .agg({
            "Mean": "mean",
            "Std": "mean",
            "Best": "mean",
            "Worst": "mean",
            "FeasibleRate(%)": "mean",
            "MeanViolation": "mean"
        })
        .sort_values(["Mean", "MeanViolation"])
        .reset_index(drop=True)
    )
    overall_df.to_csv(os.path.join(RESULTS_DIR, "overall_algorithm_summary.csv"), index=False)

    overall_rank_df = (
        friedman_all_df.groupby("Algorithm", as_index=False)["Mean Rank"]
        .mean()
        .sort_values("Mean Rank", ascending=True)
        .reset_index(drop=True)
    )
    overall_rank_df.to_csv(os.path.join(RESULTS_DIR, "overall_friedman_mean_ranks.csv"), index=False)

    print("\nSaved outputs:")
    print(os.path.join(RESULTS_DIR, "all_benchmark_results.csv"))
    print(os.path.join(RESULTS_DIR, "overall_algorithm_summary.csv"))
    print(os.path.join(RESULTS_DIR, "wilcoxon_creo_de_vs_others.csv"))
    print(os.path.join(RESULTS_DIR, "friedman_ranks.csv"))
    print(os.path.join(RESULTS_DIR, "overall_friedman_mean_ranks.csv"))
    print(os.path.join(TABLES_DIR))
    print(os.path.join(PLOTS_DIR))
    print(os.path.join(CURVES_DIR))

if __name__ == "__main__":
    main()