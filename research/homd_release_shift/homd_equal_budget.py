#!/usr/bin/env python3
"""Current-release holdouts at 1015-reference budget to control data density."""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from homd_data import DATA,FILES,read_release
from homd_deep_benchmark import current_holdout, score_batch, out_metrics


def main():
    new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
    out={}
    for seed in (11,22,33,44,55):
        full,q=current_holdout(new,seed)
        rng=np.random.default_rng(seed+1000)
        by=defaultdict(list)
        for r in full:by[r['hmt']].append(r)
        # One per known HMT, then randomly allocate the remaining reference budget.
        selected=[by[h][int(rng.integers(len(by[h])))] for h in sorted(by)]
        ids={id(r) for r in selected}
        rest=[r for r in full if id(r) not in ids]
        n=1015-len(selected)
        if n<0:raise ValueError('Budget smaller than known taxa')
        reduced=selected+[rest[i] for i in rng.choice(len(rest),n,replace=False)]
        score=score_batch(reduced,q)
        out[str(seed)]={'train_n':len(reduced),'train_hmts':len(by),
                        'cos':out_metrics(q,score['cos_conf'],score['cos_hmt'])}
        print(seed,out[str(seed)],flush=True)
    Path('homd_equal_budget_results.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
