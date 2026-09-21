# ============================================================
# CREO-DE FULL EXPERIMENT
# ============================================================
# Algorithms:
#   - GA
#   - DE
#   - PSO
#   - CREO-Penalty
#   - CREO-Repair
#   - CREO-RG
#   - CREO-DE  <-- proposed main hybrid variant
#
# Outputs:
#   - table3_section6_1_benchmark_summary.csv
#   - wilcoxon_creo_de_vs_others.csv
#   - friedman_ranks.csv
#   - friedman_summary.csv
#   - plots/mean_convergence.png
#   - best_midis/*.mid
# ============================================================

import os
import copy
import random
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pretty_midi
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon, friedmanchisquare

# ============================================================
# CONFIG
# ============================================================

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

RESULTS_DIR = "results_creo_de_experiment"
BEST_MIDI_DIR = os.path.join(RESULTS_DIR, "best_midis")
PLOTS_DIR = os.path.join(RESULTS_DIR, "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(BEST_MIDI_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

ALGORITHMS = [
    "GA",
    "DE",
    "PSO",
    "CREO-Penalty",
    "CREO-Repair",
    "CREO-RG",
    "CREO-DE",
]

POP_SIZE = 40
ITERATIONS = 300
RUNS = 20

MUSIC_BARS = 8
BEATS_PER_BAR = 4
STEPS_PER_BEAT = 4
TOTAL_STEPS = MUSIC_BARS * BEATS_PER_BAR * STEPS_PER_BEAT

FEASIBLE_RATE_ITER = 100
TIME_TO_K_FEASIBLE = 5

# ============================================================
# SEARCH SPACE
# [tempo, key, density, mean_pitch, pitch_span, step_prob,
#  chord_prog_id, contour_bias, rhythm_variation, phrase_length,
#  cadence_strength, motif_repeat_rate, accent_strength]
# ============================================================

LOWER_BOUNDS = np.array([
    60,    # tempo
    0,     # key
    0.10,  # density
    52,    # mean_pitch
    5,     # pitch_span
    0.05,  # step_prob
    0,     # chord_prog_id
    -1.0,  # contour_bias
    0.00,  # rhythm_variation
    2,     # phrase_length
    0.00,  # cadence_strength
    0.00,  # motif_repeat_rate
    0.00,  # accent_strength
], dtype=float)

UPPER_BOUNDS = np.array([
    180,
    11,
    0.95,
    78,
    24,
    0.95,
    5,
    1.0,
    1.00,
    8,
    1.00,
    1.00,
    1.00,
], dtype=float)

DIM = len(LOWER_BOUNDS)

# ============================================================
# FEASIBILITY THRESHOLDS
# ============================================================

THRESHOLDS = {
    "tempo_score": 0.75,
    "key_score": 0.85,
    "harmony_score": 0.85,
    "smoothness_score": 0.62,
    "structure_score": 0.58,
    "rhythm_score": 0.62,
    "cadence_score": 0.58,
}

WEIGHTS = {
    "tempo_score": 0.13,
    "key_score": 0.17,
    "harmony_score": 0.17,
    "smoothness_score": 0.13,
    "structure_score": 0.11,
    "rhythm_score": 0.11,
    "cadence_score": 0.11,
}

JOINT_VALIDITY_BONUS = 0.035
VALIDITY_MARGIN_WEIGHT = 0.020

TARGET_TEMPO = 110.0

MAJOR_SCALE = [0, 2, 4, 5, 7, 9, 11]
MINOR_SCALE = [0, 2, 3, 5, 7, 8, 10]

CHORD_PROGRESSIONS = {
    0: [0, 5, 7, 0],
    1: [0, 7, 9, 5],
    2: [0, 9, 5, 7],
    3: [0, 3, 4, 0],
    4: [0, 5, 9, 7],
    5: [0, 2, 5, 7],
}


# ============================================================
# DATA CLASS
# ============================================================

@dataclass
class CandidateResult:
    raw_vector: np.ndarray
    repaired_vector: np.ndarray
    repair_displacement: np.ndarray
    fitness: float
    feasible: int
    total_constraint_violation: float
    joint_validity: int
    midi: pretty_midi.PrettyMIDI
    metric_scores: Dict[str, float]


# ============================================================
# UTILS
# ============================================================

def clip_vector(x: np.ndarray) -> np.ndarray:
    return np.minimum(np.maximum(x, LOWER_BOUNDS), UPPER_BOUNDS)

def soft01(x: float) -> float:
    return max(0.0, min(1.0, x))

def random_vector() -> np.ndarray:
    return LOWER_BOUNDS + np.random.rand(DIM) * (UPPER_BOUNDS - LOWER_BOUNDS)

def choose_scale(key: int) -> List[int]:
    base = MAJOR_SCALE if key % 2 == 0 else MINOR_SCALE
    return [((key + s) % 12) for s in base]

def chord_tones(root_pc: int, is_minor: bool = False) -> List[int]:
    if is_minor:
        return [root_pc % 12, (root_pc + 3) % 12, (root_pc + 7) % 12]
    return [root_pc % 12, (root_pc + 4) % 12, (root_pc + 7) % 12]

def nearest_pitch_in_scale(pitch: int, scale_pcs: List[int]) -> int:
    best_pitch = pitch
    best_dist = 999
    for candidate in range(max(36, pitch - 12), min(96, pitch + 13)):
        if candidate % 12 in scale_pcs:
            d = abs(candidate - pitch)
            if d < best_dist:
                best_dist = d
                best_pitch = candidate
    return best_pitch

def argbest(values: List[float]) -> int:
    return int(np.argmax(values))

def decode_vector(x: np.ndarray) -> Dict:
    x = clip_vector(x).copy()
    return {
        "tempo": float(x[0]),
        "key": int(round(x[1])) % 12,
        "density": float(x[2]),
        "mean_pitch": int(round(x[3])),
        "pitch_span": int(round(x[4])),
        "step_prob": float(x[5]),
        "chord_prog_id": int(round(x[6])) % len(CHORD_PROGRESSIONS),
        "contour_bias": float(x[7]),
        "rhythm_variation": float(x[8]),
        "phrase_length": int(round(x[9])),
        "cadence_strength": float(x[10]),
        "motif_repeat_rate": float(x[11]),
        "accent_strength": float(x[12]),
    }


# ============================================================
# CONSTRAINTS / REPAIR
# ============================================================

def compute_constraint_violation(x: np.ndarray) -> float:
    p = decode_vector(x)
    v = 0.0

    if p["tempo"] < 88:
        v += 88 - p["tempo"]
    if p["tempo"] > 132:
        v += p["tempo"] - 132

    if p["density"] < 0.28:
        v += (0.28 - p["density"]) * 16
    if p["density"] > 0.72:
        v += (p["density"] - 0.72) * 16

    if p["pitch_span"] < 8:
        v += 8 - p["pitch_span"]
    if p["pitch_span"] > 15:
        v += p["pitch_span"] - 15

    if p["step_prob"] < 0.42:
        v += (0.42 - p["step_prob"]) * 16

    if p["rhythm_variation"] > 0.70:
        v += (p["rhythm_variation"] - 0.70) * 10

    if p["phrase_length"] < 3 or p["phrase_length"] > 6:
        v += abs(p["phrase_length"] - 4.5)

    if p["cadence_strength"] < 0.35:
        v += (0.35 - p["cadence_strength"]) * 10

    if p["motif_repeat_rate"] < 0.18:
        v += (0.18 - p["motif_repeat_rate"]) * 10
    if p["motif_repeat_rate"] > 0.78:
        v += (p["motif_repeat_rate"] - 0.78) * 10

    return float(v)

def repair_vector(x: np.ndarray) -> np.ndarray:
    p = decode_vector(x)

    p["tempo"] = min(max(p["tempo"], 88), 132)
    p["density"] = min(max(p["density"], 0.28), 0.72)
    p["pitch_span"] = min(max(p["pitch_span"], 8), 15)
    p["step_prob"] = min(max(p["step_prob"], 0.42), 0.90)
    p["mean_pitch"] = min(max(p["mean_pitch"], 58), 74)
    p["rhythm_variation"] = min(max(p["rhythm_variation"], 0.04), 0.70)
    p["phrase_length"] = min(max(p["phrase_length"], 3), 6)
    p["cadence_strength"] = min(max(p["cadence_strength"], 0.35), 1.00)
    p["motif_repeat_rate"] = min(max(p["motif_repeat_rate"], 0.18), 0.78)
    p["accent_strength"] = min(max(p["accent_strength"], 0.12), 0.88)

    repaired = np.array([
        p["tempo"],
        p["key"],
        p["density"],
        p["mean_pitch"],
        p["pitch_span"],
        p["step_prob"],
        p["chord_prog_id"],
        p["contour_bias"],
        p["rhythm_variation"],
        p["phrase_length"],
        p["cadence_strength"],
        p["motif_repeat_rate"],
        p["accent_strength"],
    ], dtype=float)

    return clip_vector(repaired)


# ============================================================
# MUSIC GENERATION
# ============================================================

def generate_midi_from_vector(x: np.ndarray) -> pretty_midi.PrettyMIDI:
    p = decode_vector(x)

    tempo = p["tempo"]
    key = p["key"]
    density = p["density"]
    mean_pitch = p["mean_pitch"]
    pitch_span = p["pitch_span"]
    step_prob = p["step_prob"]
    contour_bias = p["contour_bias"]
    rhythm_variation = p["rhythm_variation"]
    phrase_length = p["phrase_length"]
    cadence_strength = p["cadence_strength"]
    motif_repeat_rate = p["motif_repeat_rate"]
    accent_strength = p["accent_strength"]
    progression = CHORD_PROGRESSIONS[p["chord_prog_id"]]

    seconds_per_beat = 60.0 / tempo
    step_duration = seconds_per_beat / STEPS_PER_BEAT

    pm = pretty_midi.PrettyMIDI(initial_tempo=tempo)
    melody = pretty_midi.Instrument(program=0, name="Melody")
    chords = pretty_midi.Instrument(program=0, name="Chords")
    bass = pretty_midi.Instrument(program=32, name="Bass")

    scale_pcs = choose_scale(key)
    is_minor = (key % 2 == 1)

    for bar in range(MUSIC_BARS):
        chord_degree_pc = (key + progression[bar % len(progression)]) % 12
        tones = chord_tones(chord_degree_pc, is_minor=is_minor)

        start = bar * BEATS_PER_BAR * seconds_per_beat
        end = (bar + 1) * BEATS_PER_BAR * seconds_per_beat

        bass_pitch = 36 + (chord_degree_pc % 12)
        while bass_pitch < 36:
            bass_pitch += 12
        while bass_pitch > 48:
            bass_pitch -= 12

        bass.notes.append(
            pretty_midi.Note(velocity=72, pitch=bass_pitch, start=start, end=end)
        )

        chord_base = 48 + (key % 12)
        for pc in tones:
            pitch = chord_base + ((pc - chord_base) % 12)
            while pitch < 48:
                pitch += 12
            while pitch > 72:
                pitch -= 12
            chords.notes.append(
                pretty_midi.Note(velocity=65, pitch=pitch, start=start, end=end)
            )

    motif_length = max(2, phrase_length)
    motif = []
    current_pitch = mean_pitch
    low = mean_pitch - pitch_span // 2
    high = mean_pitch + pitch_span // 2

    for _ in range(motif_length):
        direction = 1 if np.random.rand() < (0.5 + 0.28 * contour_bias) else -1
        if np.random.rand() < step_prob:
            move = direction * np.random.choice([1, 2])
        else:
            move = direction * np.random.choice([3, 4, 5])
        proposed = int(np.clip(current_pitch + move, low, high))
        proposed = nearest_pitch_in_scale(proposed, scale_pcs)
        motif.append(proposed)
        current_pitch = proposed

    current_pitch = mean_pitch
    for step in range(TOTAL_STEPS):
        if np.random.rand() > density:
            continue

        t0 = step * step_duration
        dur = step_duration

        jitter = rhythm_variation * step_duration * 0.6 * (np.random.rand() - 0.5)
        t0 = max(0.0, t0 + jitter)

        phrase_pos = step % max(1, phrase_length)
        use_motif = np.random.rand() < motif_repeat_rate

        if use_motif and phrase_pos < len(motif):
            proposed = motif[phrase_pos]
        else:
            direction = 1 if np.random.rand() < (0.5 + 0.28 * contour_bias) else -1
            if np.random.rand() < step_prob:
                move = direction * np.random.choice([1, 2])
            else:
                move = direction * np.random.choice([3, 4, 5, 7])

            proposed = int(np.clip(current_pitch + move, low, high))
            proposed = nearest_pitch_in_scale(proposed, scale_pcs)

        is_phrase_end = ((step + 1) % phrase_length == 0)
        if is_phrase_end and np.random.rand() < cadence_strength:
            tonic_pitch = nearest_pitch_in_scale(mean_pitch, scale_pcs)
            dominant_pitch = nearest_pitch_in_scale(mean_pitch + 7, scale_pcs)
            proposed = tonic_pitch if np.random.rand() < 0.7 else dominant_pitch

        velocity = 80
        if step % STEPS_PER_BEAT == 0:
            velocity += int(18 * accent_strength)

        melody.notes.append(
            pretty_midi.Note(
                velocity=int(np.clip(velocity, 40, 120)),
                pitch=int(np.clip(proposed, 48, 84)),
                start=t0,
                end=t0 + dur
            )
        )
        current_pitch = proposed

    pm.instruments.append(chords)
    pm.instruments.append(bass)
    pm.instruments.append(melody)
    return pm


# ============================================================
# EVALUATION
# ============================================================

def compute_validity_margin(metric_scores: Dict[str, float]) -> float:
    margins = []
    for k, t in THRESHOLDS.items():
        margins.append(metric_scores[k] - t)
    return float(np.mean(margins))

def evaluate_midi(pm: pretty_midi.PrettyMIDI, x: np.ndarray) -> Dict[str, float]:
    p = decode_vector(x)
    melody_notes = pm.instruments[2].notes if len(pm.instruments) > 2 else []

    if len(melody_notes) == 0:
        out = {k: 0.0 for k in THRESHOLDS}
        out["joint_validity"] = 0
        out["validity_margin"] = -1.0
        out["fitness"] = -1e6
        return out

    scale_pcs = choose_scale(p["key"])
    progression = CHORD_PROGRESSIONS[p["chord_prog_id"]]
    is_minor = (p["key"] % 2 == 1)

    pitches = np.array([n.pitch for n in melody_notes], dtype=float)
    starts = np.array([n.start for n in melody_notes], dtype=float)

    tempo_error = abs(p["tempo"] - TARGET_TEMPO)
    tempo_score = soft01(1.0 - tempo_error / 45.0)

    key_score = float(np.mean([(int(pp) % 12) in scale_pcs for pp in pitches]))

    seconds_per_bar = 60.0 / p["tempo"] * BEATS_PER_BAR
    harmonic_hits = 0
    for note in melody_notes:
        bar_idx = min(int(note.start // seconds_per_bar), MUSIC_BARS - 1)
        chord_root_pc = (p["key"] + progression[bar_idx % len(progression)]) % 12
        chord_pc = chord_tones(chord_root_pc, is_minor=is_minor)
        npc = note.pitch % 12
        if npc in chord_pc or npc in scale_pcs:
            harmonic_hits += 1
    harmony_score = harmonic_hits / max(1, len(melody_notes))

    if len(pitches) > 1:
        intervals = np.abs(np.diff(pitches))
        mean_interval = np.mean(intervals)
        large_jump_penalty = np.mean(intervals > 7)
        smoothness_score = soft01(1.0 - (mean_interval / 9.0) - 0.45 * large_jump_penalty)
    else:
        smoothness_score = 0.5

    phrase_len = max(2, p["phrase_length"])
    phrase_vectors = []
    for i in range(0, len(pitches), phrase_len):
        seg = pitches[i:i + phrase_len]
        if len(seg) == phrase_len:
            phrase_vectors.append(seg)

    if len(phrase_vectors) >= 2:
        sims = []
        base = phrase_vectors[0]
        for other in phrase_vectors[1:]:
            diff = np.mean(np.abs(base - other))
            sims.append(soft01(1.0 - diff / 11.0))
        structure_score = float(np.mean(sims))
    else:
        structure_score = 0.45

    if len(starts) > 1:
        iois = np.diff(np.sort(starts))
        target_step = (60.0 / p["tempo"]) / STEPS_PER_BEAT
        mean_dev = np.mean(np.abs(iois - target_step))
        rhythm_score = soft01(1.0 - mean_dev / (1.2 * target_step + 1e-8))
    else:
        rhythm_score = 0.45

    cadence_hits = 0
    cadence_checks = 0
    for bar in range(1, MUSIC_BARS + 1):
        bar_end = bar * BEATS_PER_BAR * (60.0 / p["tempo"])
        candidates = [n for n in melody_notes if (bar_end - 0.5) <= n.start <= bar_end]
        if candidates:
            cadence_checks += 1
            last_note = sorted(candidates, key=lambda n: n.start)[-1]
            if (last_note.pitch % 12) in {p["key"] % 12, (p["key"] + 7) % 12}:
                cadence_hits += 1
    cadence_score = cadence_hits / max(1, cadence_checks)

    metric_scores = {
        "tempo_score": float(tempo_score),
        "key_score": float(key_score),
        "harmony_score": float(harmony_score),
        "smoothness_score": float(smoothness_score),
        "structure_score": float(structure_score),
        "rhythm_score": float(rhythm_score),
        "cadence_score": float(cadence_score),
    }

    joint_validity = int(all(metric_scores[k] >= THRESHOLDS[k] for k in THRESHOLDS))
    validity_margin = compute_validity_margin(metric_scores)

    base_fitness = sum(WEIGHTS[k] * metric_scores[k] for k in WEIGHTS)
    fitness = base_fitness + JOINT_VALIDITY_BONUS * joint_validity + VALIDITY_MARGIN_WEIGHT * validity_margin

    metric_scores["joint_validity"] = joint_validity
    metric_scores["validity_margin"] = float(validity_margin)
    metric_scores["fitness"] = float(fitness)
    return metric_scores


# ============================================================
# EVALUATION WRAPPER
# ============================================================

def evaluate_candidate(x: np.ndarray, mode: str) -> CandidateResult:
    raw = clip_vector(x)
    repaired = repair_vector(raw)
    displacement = repaired - raw
    raw_violation = compute_constraint_violation(raw)

    used = raw.copy() if mode == "penalty" else repaired.copy()
    midi = generate_midi_from_vector(used)
    scores = evaluate_midi(midi, used)

    if mode == "penalty":
        penalty = 0.08 * raw_violation
        fitness = scores["fitness"] - penalty
        feasible = int(raw_violation <= 1e-8 and scores["joint_validity"] == 1)
        total_violation = float(raw_violation + (0.0 if scores["joint_validity"] else 1.0))
    else:
        fitness = scores["fitness"]
        feasible = int(scores["joint_validity"] == 1)
        total_violation = float(0.0 if feasible else 1.0)

    return CandidateResult(
        raw_vector=raw,
        repaired_vector=repaired,
        repair_displacement=displacement,
        fitness=float(fitness),
        feasible=int(feasible),
        total_constraint_violation=float(total_violation),
        joint_validity=int(scores["joint_validity"]),
        midi=midi,
        metric_scores=scores,
    )


# ============================================================
# SELECTION / HELPERS
# ============================================================

def init_population(n: int) -> np.ndarray:
    return np.array([random_vector() for _ in range(n)], dtype=float)

def tournament_select(pop: np.ndarray, fitnesses: List[float], k: int = 3) -> np.ndarray:
    idx = np.random.choice(len(pop), size=k, replace=False)
    best = idx[np.argmax([fitnesses[i] for i in idx])]
    return pop[best].copy()

def gaussian_mutation(x: np.ndarray, sigma: float = 0.06) -> np.ndarray:
    span = UPPER_BOUNDS - LOWER_BOUNDS
    y = x + np.random.randn(DIM) * span * sigma
    return clip_vector(y)

def one_point_crossover(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    cut = np.random.randint(1, DIM)
    c1 = np.concatenate([a[:cut], b[cut:]])
    c2 = np.concatenate([b[:cut], a[cut:]])
    return clip_vector(c1), clip_vector(c2)

def feasibility_aware_better(a: CandidateResult, b: CandidateResult) -> bool:
    """
    Return True if a is preferred over b.
    Preference:
      1. feasible beats infeasible
      2. if both feasible: higher fitness wins
      3. if both infeasible: lower total violation wins
         tie-break by higher fitness
    """
    if a.feasible != b.feasible:
        return a.feasible > b.feasible

    if a.feasible == 1 and b.feasible == 1:
        return a.fitness >= b.fitness

    if a.total_constraint_violation != b.total_constraint_violation:
        return a.total_constraint_violation < b.total_constraint_violation

    return a.fitness >= b.fitness


# ============================================================
# ALGORITHMS
# ============================================================

def run_ga() -> Tuple[pd.DataFrame, CandidateResult]:
    pop = init_population(POP_SIZE)
    best_result = None
    rows = []

    for it in range(1, ITERATIONS + 1):
        evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]
        fitnesses = [e.fitness for e in evals]

        best_idx = argbest(fitnesses)
        if best_result is None or evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fitnesses),
            "mean_fitness": float(np.mean(fitnesses)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

        next_pop = [pop[best_idx].copy()]
        while len(next_pop) < POP_SIZE:
            p1 = tournament_select(pop, fitnesses)
            p2 = tournament_select(pop, fitnesses)
            c1, c2 = one_point_crossover(p1, p2)
            c1 = gaussian_mutation(c1, sigma=0.05)
            c2 = gaussian_mutation(c2, sigma=0.05)
            next_pop.extend([c1, c2])

        pop = np.array(next_pop[:POP_SIZE], dtype=float)

    return pd.DataFrame(rows), best_result


def run_de() -> Tuple[pd.DataFrame, CandidateResult]:
    F = 0.72
    CR = 0.82

    pop = init_population(POP_SIZE)
    evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]
    best_result = copy.deepcopy(evals[argbest([e.fitness for e in evals])])
    rows = []

    for it in range(1, ITERATIONS + 1):
        new_pop = []
        new_evals = []

        for i in range(POP_SIZE):
            idxs = list(range(POP_SIZE))
            idxs.remove(i)
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)
            a, b, c = pop[r1], pop[r2], pop[r3]

            mutant = clip_vector(a + F * (b - c))

            trial = pop[i].copy()
            jrand = np.random.randint(DIM)
            for j in range(DIM):
                if np.random.rand() < CR or j == jrand:
                    trial[j] = mutant[j]
            trial = clip_vector(trial)

            trial_eval = evaluate_candidate(trial, mode="penalty")
            target_eval = evals[i]

            if trial_eval.fitness >= target_eval.fitness:
                new_pop.append(trial)
                new_evals.append(trial_eval)
            else:
                new_pop.append(pop[i].copy())
                new_evals.append(target_eval)

        pop = np.array(new_pop, dtype=float)
        evals = new_evals
        fits = [e.fitness for e in evals]

        best_idx = argbest(fits)
        if evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fits),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


def run_pso() -> Tuple[pd.DataFrame, CandidateResult]:
    w, c1, c2 = 0.58, 1.20, 1.20
    vmax = 0.12 * (UPPER_BOUNDS - LOWER_BOUNDS)

    pop = init_population(POP_SIZE)
    vel = np.zeros_like(pop)

    evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]
    pbest = pop.copy()
    pbest_fit = np.array([e.fitness for e in evals], dtype=float)

    gbest_idx = int(np.argmax(pbest_fit))
    gbest = pbest[gbest_idx].copy()
    best_result = copy.deepcopy(evals[gbest_idx])
    rows = []

    for it in range(1, ITERATIONS + 1):
        r1 = np.random.rand(*pop.shape)
        r2 = np.random.rand(*pop.shape)
        vel = w * vel + c1 * r1 * (pbest - pop) + c2 * r2 * (gbest - pop)
        vel = np.clip(vel, -vmax, vmax)
        pop = clip_vector(pop + vel)

        evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]
        fits = np.array([e.fitness for e in evals], dtype=float)

        improved = fits > pbest_fit
        pbest[improved] = pop[improved]
        pbest_fit[improved] = fits[improved]

        gbest_idx = int(np.argmax(pbest_fit))
        gbest = pbest[gbest_idx].copy()

        current_best_idx = int(np.argmax(fits))
        if evals[current_best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[current_best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": float(np.max(fits)),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


def run_creo_penalty() -> Tuple[pd.DataFrame, CandidateResult]:
    alpha, gamma = 0.40, 0.06

    pop = init_population(POP_SIZE)
    evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]

    best_idx = argbest([e.fitness for e in evals])
    gbest = pop[best_idx].copy()
    best_result = copy.deepcopy(evals[best_idx])
    rows = []

    for it in range(1, ITERATIONS + 1):
        new_pop = []
        for x in pop:
            eps = np.random.randn(DIM)
            x_new = clip_vector(x + alpha * (gbest - x) + gamma * eps)
            new_pop.append(x_new)

        pop = np.array(new_pop, dtype=float)
        evals = [evaluate_candidate(ind, mode="penalty") for ind in pop]
        fits = [e.fitness for e in evals]

        best_idx = argbest(fits)
        gbest = pop[best_idx].copy()
        if evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fits),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


def run_creo_repair() -> Tuple[pd.DataFrame, CandidateResult]:
    alpha, gamma = 0.40, 0.06

    pop = init_population(POP_SIZE)
    evals = [evaluate_candidate(ind, mode="repair") for ind in pop]

    best_idx = argbest([e.fitness for e in evals])
    gbest = evals[best_idx].repaired_vector.copy()
    best_result = copy.deepcopy(evals[best_idx])
    rows = []

    for it in range(1, ITERATIONS + 1):
        new_pop = []
        for x in pop:
            eps = np.random.randn(DIM)
            x_new = clip_vector(x + alpha * (gbest - x) + gamma * eps)
            new_pop.append(x_new)

        pop = np.array(new_pop, dtype=float)
        evals = [evaluate_candidate(ind, mode="repair") for ind in pop]
        fits = [e.fitness for e in evals]

        best_idx = argbest(fits)
        gbest = evals[best_idx].repaired_vector.copy()
        if evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fits),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


def run_creo_rg() -> Tuple[pd.DataFrame, CandidateResult]:
    alpha, beta, gamma = 0.32, 1.05, 0.035

    pop = init_population(POP_SIZE)
    evals = [evaluate_candidate(ind, mode="rg") for ind in pop]

    best_idx = argbest([e.fitness for e in evals])
    gbest = evals[best_idx].repaired_vector.copy()
    validity_leader = evals[best_idx].repaired_vector.copy()
    best_result = copy.deepcopy(evals[best_idx])

    rows = []
    feasible_mode = False

    for it in range(1, ITERATIONS + 1):
        if any(e.feasible for e in evals):
            feasible_mode = True

        validity_scores = []
        for e in evals:
            validity_score = e.metric_scores["validity_margin"] + 0.25 * e.joint_validity
            validity_scores.append(validity_score)
        validity_best_idx = int(np.argmax(validity_scores))
        validity_leader = evals[validity_best_idx].repaired_vector.copy()

        new_pop = []
        for i, x in enumerate(pop):
            delta_r = evals[i].repair_displacement
            eps = np.random.randn(DIM)

            local_gamma = gamma if not feasible_mode else gamma * 0.45

            x_new = (
                x
                + alpha * (gbest - x)
                + 0.45 * (validity_leader - x)
                + beta * delta_r
                + local_gamma * eps
            )

            repaired_x = repair_vector(x)
            if feasible_mode:
                x_new = x_new + 0.22 * (repaired_x - x)
            else:
                if compute_constraint_violation(x) > 0:
                    x_new = x_new + 0.10 * (repaired_x - x)

            x_new = clip_vector(x_new)
            new_pop.append(x_new)

        pop = np.array(new_pop, dtype=float)
        evals = [evaluate_candidate(ind, mode="rg") for ind in pop]
        fits = [e.fitness for e in evals]

        best_idx = argbest(fits)
        gbest = evals[best_idx].repaired_vector.copy()
        if evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fits),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


def run_creo_de() -> Tuple[pd.DataFrame, CandidateResult]:
    """
    Hybrid CREO-DE:
      - DE mutation
      - binomial crossover
      - repair-guided displacement
      - feasibility-aware selection
      - repaired leader influence
    """
    F = 0.72
    CR = 0.85
    ALPHA = 0.20
    BETA = 0.85
    GAMMA = 0.03

    pop = init_population(POP_SIZE)
    evals = [evaluate_candidate(ind, mode="rg") for ind in pop]

    best_idx = argbest([e.fitness for e in evals])
    gbest = evals[best_idx].repaired_vector.copy()
    best_result = copy.deepcopy(evals[best_idx])

    rows = []
    feasible_mode = False

    for it in range(1, ITERATIONS + 1):
        if any(e.feasible for e in evals):
            feasible_mode = True

        # validity-aware repaired leader
        validity_scores = []
        for e in evals:
            validity_score = e.metric_scores["validity_margin"] + 0.35 * e.joint_validity
            validity_scores.append(validity_score)
        validity_best_idx = int(np.argmax(validity_scores))
        validity_leader = evals[validity_best_idx].repaired_vector.copy()

        new_pop = []
        new_evals = []

        for i in range(POP_SIZE):
            idxs = list(range(POP_SIZE))
            idxs.remove(i)
            r1, r2, r3 = np.random.choice(idxs, 3, replace=False)

            x_i = pop[i].copy()
            x_r1 = pop[r1]
            x_r2 = pop[r2]
            x_r3 = pop[r3]
            delta_r = evals[i].repair_displacement

            # hybrid mutant
            mutant = (
                x_r1
                + F * (x_r2 - x_r3)
                + ALPHA * (gbest - x_i)
                + 0.30 * (validity_leader - x_i)
                + BETA * delta_r
            )

            # stronger repaired contraction after feasible mode
            repaired_xi = repair_vector(x_i)
            if feasible_mode:
                mutant = mutant + 0.18 * (repaired_xi - x_i)
                local_gamma = GAMMA * 0.40
            else:
                if compute_constraint_violation(x_i) > 0:
                    mutant = mutant + 0.08 * (repaired_xi - x_i)
                local_gamma = GAMMA

            mutant = mutant + local_gamma * np.random.randn(DIM)
            mutant = clip_vector(mutant)

            # binomial crossover
            trial = x_i.copy()
            jrand = np.random.randint(DIM)
            for j in range(DIM):
                if np.random.rand() < CR or j == jrand:
                    trial[j] = mutant[j]
            trial = clip_vector(trial)

            trial_eval = evaluate_candidate(trial, mode="rg")
            target_eval = evals[i]

            # feasibility-aware selection
            if feasibility_aware_better(trial_eval, target_eval):
                new_pop.append(trial)
                new_evals.append(trial_eval)
            else:
                new_pop.append(x_i)
                new_evals.append(target_eval)

        pop = np.array(new_pop, dtype=float)
        evals = new_evals

        fits = [e.fitness for e in evals]
        best_idx = argbest(fits)
        gbest = evals[best_idx].repaired_vector.copy()

        if evals[best_idx].fitness > best_result.fitness:
            best_result = copy.deepcopy(evals[best_idx])

        rows.append({
            "iteration": it,
            "best_fitness": max(fits),
            "mean_fitness": float(np.mean(fits)),
            "feasible_count": int(sum(e.feasible for e in evals)),
            "feasible_any": int(any(e.feasible for e in evals)),
            "joint_validity_rate": float(np.mean([e.joint_validity for e in evals])),
        })

    return pd.DataFrame(rows), best_result


# ============================================================
# METRICS / STATISTICS / PLOTS
# ============================================================

def compute_summary_metrics(run_df: pd.DataFrame) -> Dict:
    best_fitness = float(run_df["best_fitness"].max())

    feasible_rows = run_df[run_df["feasible_any"] == 1]
    first_feasible = float(feasible_rows["iteration"].iloc[0]) if len(feasible_rows) > 0 else np.nan

    feasible_rate_at_100 = 100.0 if ((run_df["iteration"] <= FEASIBLE_RATE_ITER) & (run_df["feasible_any"] == 1)).any() else 0.0

    cumulative_feasible = run_df["feasible_count"].cumsum()
    idx = run_df[cumulative_feasible >= TIME_TO_K_FEASIBLE]
    time_to_k = float(idx["iteration"].iloc[0]) if len(idx) > 0 else np.nan

    joint_validity_pct = float(run_df["joint_validity_rate"].mean() * 100.0)

    return {
        "Mean Best Fitness": best_fitness,
        "First Feasible Iteration": first_feasible,
        "Feasible Rate @100": feasible_rate_at_100,
        "Time-to-5-Feasible": time_to_k,
        "Joint Validity %": joint_validity_pct,
    }

def aggregate_algorithm_metrics(run_metrics: List[Dict]) -> Dict:
    return {
        "Mean Best Fitness": float(np.nanmean([m["Mean Best Fitness"] for m in run_metrics])),
        "First Feasible Iteration ↓": float(np.nanmean([m["First Feasible Iteration"] for m in run_metrics])),
        "Feasible Rate @100 ↑": float(np.nanmean([m["Feasible Rate @100"] for m in run_metrics])),
        "Time-to-5-Feasible ↓": float(np.nanmean([m["Time-to-5-Feasible"] for m in run_metrics])),
        "Joint Validity % ↑": float(np.nanmean([m["Joint Validity %"] for m in run_metrics])),
    }

def run_statistics(best_fitness_per_algo: Dict[str, List[float]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    target = "CREO-DE"
    creo_scores = np.array(best_fitness_per_algo[target], dtype=float)

    wilcoxon_rows = []
    for algo, scores in best_fitness_per_algo.items():
        if algo == target:
            continue
        other = np.array(scores, dtype=float)
        try:
            stat, p = wilcoxon(creo_scores, other, zero_method="wilcox", alternative="two-sided")
        except ValueError:
            stat, p = np.nan, np.nan
        wilcoxon_rows.append({
            "Comparison": f"{target} vs {algo}",
            "Statistic": stat,
            "p-value": p,
            "Significant (p<0.05)": int(pd.notna(p) and p < 0.05)
        })

    wilcoxon_df = pd.DataFrame(wilcoxon_rows)

    arrays = [np.array(best_fitness_per_algo[a], dtype=float) for a in ALGORITHMS]
    friedman_stat, friedman_p = friedmanchisquare(*arrays)

    ranks_matrix = np.vstack(arrays).T
    rank_rows = []
    for row in ranks_matrix:
        order = (-row).argsort()
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(row) + 1)
        rank_rows.append(ranks)
    rank_rows = np.array(rank_rows)
    mean_ranks = np.mean(rank_rows, axis=0)

    friedman_df = pd.DataFrame({
        "Algorithm": ALGORITHMS,
        "Mean Rank": mean_ranks
    }).sort_values("Mean Rank", ascending=True).reset_index(drop=True)

    pd.DataFrame([{
        "Friedman Statistic": friedman_stat,
        "p-value": friedman_p
    }]).to_csv(os.path.join(RESULTS_DIR, "friedman_summary.csv"), index=False)
    friedman_df.to_csv(os.path.join(RESULTS_DIR, "friedman_ranks.csv"), index=False)

    return wilcoxon_df, friedman_df

def plot_mean_convergence(mean_curves: Dict[str, np.ndarray], save_path: str):
    plt.figure(figsize=(10, 6))
    x = np.arange(1, ITERATIONS + 1)
    for algo, curve in mean_curves.items():
        plt.plot(x, curve, label=algo)
    plt.xlabel("Iteration")
    plt.ylabel("Mean Best Fitness")
    plt.title("Mean Convergence Curves")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()


# ============================================================
# DRIVER
# ============================================================

def run_algorithm(algo_name: str) -> Tuple[pd.DataFrame, CandidateResult]:
    if algo_name == "GA":
        return run_ga()
    if algo_name == "DE":
        return run_de()
    if algo_name == "PSO":
        return run_pso()
    if algo_name == "CREO-Penalty":
        return run_creo_penalty()
    if algo_name == "CREO-Repair":
        return run_creo_repair()
    if algo_name == "CREO-RG":
        return run_creo_rg()
    if algo_name == "CREO-DE":
        return run_creo_de()
    raise ValueError(f"Unknown algorithm: {algo_name}")

def main():
    final_summary_rows = []
    best_fitness_per_algo = {algo: [] for algo in ALGORITHMS}
    convergence_store = {algo: [] for algo in ALGORITHMS}

    for algo in ALGORITHMS:
        print(f"\nRunning {algo} ...")
        algo_run_metrics = []
        best_overall_result = None
        best_overall_fitness = -1e18

        for run_id in range(1, RUNS + 1):
            print(f"  Run {run_id}/{RUNS}")
            run_df, best_result = run_algorithm(algo)

            run_path = os.path.join(RESULTS_DIR, f"{algo}_run_{run_id:02d}.csv")
            run_df.to_csv(run_path, index=False)

            metrics = compute_summary_metrics(run_df)
            algo_run_metrics.append(metrics)
            best_fitness_per_algo[algo].append(metrics["Mean Best Fitness"])
            convergence_store[algo].append(run_df["best_fitness"].values)

            if best_result.fitness > best_overall_fitness:
                best_overall_fitness = best_result.fitness
                best_overall_result = best_result

        if best_overall_result is not None:
            best_overall_result.midi.write(os.path.join(BEST_MIDI_DIR, f"{algo}_best.mid"))

        summary = aggregate_algorithm_metrics(algo_run_metrics)
        summary["Algorithm"] = algo
        final_summary_rows.append(summary)

    summary_df = pd.DataFrame(final_summary_rows)[[
        "Algorithm",
        "Mean Best Fitness",
        "First Feasible Iteration ↓",
        "Feasible Rate @100 ↑",
        "Time-to-5-Feasible ↓",
        "Joint Validity % ↑",
    ]]

    summary_df = summary_df.sort_values(
        by=["Mean Best Fitness", "Joint Validity % ↑"],
        ascending=[False, False]
    ).reset_index(drop=True)

    print("\n=== SECTION 6.1 BENCHMARK SUMMARY ===\n")
    print(summary_df.to_string(index=False))

    summary_csv = os.path.join(RESULTS_DIR, "table3_section6_1_benchmark_summary.csv")
    summary_tex = os.path.join(RESULTS_DIR, "table3_section6_1_benchmark_summary.tex")
    summary_df.to_csv(summary_csv, index=False)
    with open(summary_tex, "w", encoding="utf-8") as f:
        f.write(summary_df.to_latex(index=False, escape=False, float_format="%.6f"))

    wilcoxon_df, friedman_ranks_df = run_statistics(best_fitness_per_algo)
    wilcoxon_df.to_csv(os.path.join(RESULTS_DIR, "wilcoxon_creo_de_vs_others.csv"), index=False)
    with open(os.path.join(RESULTS_DIR, "wilcoxon_creo_de_vs_others.tex"), "w", encoding="utf-8") as f:
        f.write(wilcoxon_df.to_latex(index=False, escape=False, float_format="%.6f"))

    with open(os.path.join(RESULTS_DIR, "friedman_ranks.tex"), "w", encoding="utf-8") as f:
        f.write(friedman_ranks_df.to_latex(index=False, escape=False, float_format="%.6f"))

    mean_curves = {}
    for algo in ALGORITHMS:
        curve_matrix = np.vstack(convergence_store[algo])
        mean_curves[algo] = np.mean(curve_matrix, axis=0)

    plot_mean_convergence(mean_curves, os.path.join(PLOTS_DIR, "mean_convergence.png"))

    print("\nSaved files in:", RESULTS_DIR)
    print("Main summary:", summary_csv)
    print("Wilcoxon   :", os.path.join(RESULTS_DIR, "wilcoxon_creo_de_vs_others.csv"))
    print("Friedman   :", os.path.join(RESULTS_DIR, "friedman_ranks.csv"))
    print("Plot       :", os.path.join(PLOTS_DIR, "mean_convergence.png"))


if __name__ == "__main__":
    main()
