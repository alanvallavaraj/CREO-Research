#!/usr/bin/env python3
"""Full HMT release-shift benchmark and CREO-DE ablation.

The outcome is measured after code/data are frozen. All weight tuning and
threshold calibration use v15.23 only; v16.03 is evaluated once. Read README.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from homd_data import DATA, FILES, read_release
from homd_deep_benchmark import (score_batch, current_holdout, standardize_old_new,
                                 cluster_bootstrap_auc)

ROOT = Path(__file__).resolve().parent
METHODS = ('random', 'de', 'creo_no_rg', 'creo_rg')
N_FEATURE = 4
FEATURES = ('cosine_similarity', 'global_edit_similarity_top5',
            'cosine_margin_different_hmt', 'log_best_hmt_reference_count')


def make_features(reference, query):
    s = score_batch(reference, query, align=True)
    counts = Counter(r['hmt'] for r in reference)
    x = np.stack([s['cos_conf'], s['align_conf'],
                  np.clip(s['cos_margin'], 0, 1),
                  np.log1p([counts[t] for t in s['cos_hmt']])], axis=1)
    y = np.array([r['known'] for r in query], dtype=bool)
    return x, y, s['cos_hmt']


def rescale(episodes, test):
    stack = np.concatenate([x for x, y in episodes])
    low = np.quantile(stack, .05, axis=0)
    high = np.quantile(stack, .95, axis=0)
    span = np.maximum(high - low, 1e-8)
    fn = lambda x: np.clip((x-low)/span, 0, 1)
    return [(fn(x),y) for x,y in episodes],fn(test),low.tolist(),high.tolist()


def repair(x):
    """Euclidean projection to simplex, same for every optimiser."""
    x=np.asarray(x, dtype=float)
    u=np.sort(x)[::-1]
    cssv=np.cumsum(u)-1
    idx=np.nonzero(u-cssv/(np.arange(len(u))+1)>0)[0]
    theta=cssv[idx[-1]]/(idx[-1]+1)
    return np.maximum(x-theta,0)


def objective(w, episodes):
    # Mean of five old-release, HMT-held-out task AUROCs.
    return 1.0-np.mean([roc_auc_score(y,x@w) for x,y in episodes])


def optimise(method, episodes, seed, budget=200, pop_size=20):
    rng=np.random.default_rng(seed)
    pop=[np.clip(rng.normal(.25,.4,N_FEATURE),-1,2) for _ in range(pop_size)]
    fitness=[objective(repair(p),episodes) for p in pop]
    best=int(np.argmin(fitness));global_best=repair(pop[best])
    best_val=fitness[best];evals=pop_size
    trace=[(evals,best_val)]
    if method=='random':
        while evals<budget:
            candidate=repair(rng.normal(.25,.4,N_FEATURE))
            val=objective(candidate,episodes);evals+=1
            if val<best_val:global_best,best_val=candidate.copy(),val
            trace.append((evals,best_val))
    else:
        while evals<budget:
            for i in range(pop_size):
                if evals>=budget:break
                r=rng.choice([j for j in range(pop_size) if j!=i],3,replace=False)
                f,cr=.7,.9
                mutant=pop[r[0]]+f*(pop[r[1]]-pop[r[2]])
                if method in ('creo_no_rg','creo_rg'):
                    displacement=repair(pop[i]) - pop[i]
                    mutant += .18*(global_best-pop[i]) + .02*rng.standard_normal(N_FEATURE)
                    if method=='creo_rg':mutant += .90*displacement
                # The same constraint repair is used by all three DE variants.
                crossed=rng.random(N_FEATURE)<cr
                crossed[rng.integers(N_FEATURE)]=True
                raw=np.clip(np.where(crossed,mutant,pop[i]),-1,2)
                candidate=repair(raw)
                val=objective(candidate,episodes);evals+=1
                if val<fitness[i]:pop[i],fitness[i]=raw.copy(),val
                if fitness[i]<best_val:global_best,best_val=repair(pop[i]),fitness[i]
                trace.append((evals,best_val))
    assert evals==budget
    return global_best, best_val, trace


def outcomes(rows, score, predicted, threshold):
    y=np.array([r['known'] for r in rows])
    return {'known_n':int(sum(y)), 'unknown_n':int(sum(~y)),
            'auroc':float(roc_auc_score(y,score)),
            'known_hmt_accuracy':float(np.mean(predicted[y]==np.array([r['hmt'] for r in rows])[y])),
            'known_coverage_old_calibrated':float(np.mean(score[y]>=threshold)),
            'unknown_acceptance_old_calibrated':float(np.mean(score[~y]>=threshold)),
            'unknown_acceptance_at_test_defined_90pct_known_coverage_descriptive':float(np.mean(score[~y]>=np.quantile(score[y],.1)))}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seeds',type=int,default=20)
    parser.add_argument('--budget',type=int,default=200)
    a=parser.parse_args()
    old=read_release(DATA/Path(FILES['old_fasta']).name,DATA/Path(FILES['old_tax']).name,False)
    new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
    episodes=[]
    for seed in range(101,111):
        train,qs=current_holdout(old,seed,max_per_group=1000)
        x,y,_=make_features(train,qs)
        episodes.append((x,y))
        print('old episode',seed,len(train),len(qs),flush=True)
    temporal=standardize_old_new(old,new)
    xt,yt,closest=make_features(old,temporal)
    scaled,xt,low,high=rescale(episodes,xt)
    dev=scaled[:5];cal=scaled[5:]
    results={'design':{'seed_count':a.seeds,'budget_per_run':a.budget,
                       'population':20,'features':FEATURES,
                       'optimization_development_episode_seeds':list(range(101,106)),
                       'calibration_episode_seeds':list(range(106,111)),
                       'heldout_version':'HOMD 16.03',
                       'objective':'mean known-vs-unknown AUROC on five v15.23 pseudo-unknown episodes',
                       'data_policy':'v16.03 labels never enter training/threshold calibration',
                       'warning':'Multiple optimiser seeds share the same temporal test set and are not independent biological replications.'},
             'feature_scale':{'p05':low,'p95':high},'runs':{},
             'source_hashes':{k:hashlib.sha256((DATA/Path(v).name).read_bytes()).hexdigest()
                              for k,v in FILES.items()}}
    # Non-optimised baseline; raw cosine score before clipping to old quantiles.
    scores=score_batch(old,temporal)
    baseline_dev=[x[:,0] for x,y in episodes[:5]]
    baseline_cal=[x[:,0] for x,y in episodes[5:]]
    cal_known=np.concatenate([s[y] for s,(_,y) in zip(baseline_cal,episodes[5:])])
    t0=float(np.quantile(cal_known,.10))
    results['raw_cosine_baseline']={'development_auroc':float(1-objective(np.array([1,0,0,0]),dev)),
                                    'threshold_old':t0,
                                    'temporal':outcomes(temporal,scores['cos_conf'],closest,t0)}
    for method in METHODS:
        runs=[]
        for seed in range(a.seeds):
            w,fitness,trace=optimise(method,dev,seed=seed+7701,budget=a.budget)
            known=np.concatenate([(x@w)[y] for x,y in cal])
            threshold=float(np.quantile(known,.10))
            run={'seed':seed+7701,'weights':w.tolist(),
                 'n_objective_evaluations':a.budget,'development_auroc':1-fitness,
                 'calibration_auroc':float(np.mean([roc_auc_score(y,x@w) for x,y in cal])),
                 'threshold_old':threshold,
                 'temporal':outcomes(temporal,xt@w,closest,threshold),
                 'trace':trace}
            runs.append(run)
        results['runs'][method]=runs
        vals=[r['temporal']['auroc'] for r in runs]
        print(method,'test AUROC mean/range',np.mean(vals),min(vals),max(vals),flush=True)
    # One prespecified representative for a taxon bootstrap: median dev-score run.
    for method in METHODS:
        reps=results['runs'][method]
        mid=sorted(reps,key=lambda r:r['development_auroc'])[len(reps)//2]
        results.setdefault('cluster_bootstrap_vs_raw_cosine',{})[method]=cluster_bootstrap_auc(
            temporal,scores['cos_conf'],xt@np.array(mid['weights']),iters=800,seed=331)
    p=ROOT/'creo_full_experiment_results.json'
    p.write_text(json.dumps(results,indent=2)+'\n')
    print('saved',p,flush=True)

if __name__=='__main__':main()
