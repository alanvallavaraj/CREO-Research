#!/usr/bin/env python3
"""Sensitivity: one sampled query per HMT in temporal assessment."""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from homd_data import DATA,FILES,read_release
from homd_deep_benchmark import score_batch,standardize_old_new


def main():
    old=read_release(DATA/Path(FILES['old_fasta']).name,DATA/Path(FILES['old_tax']).name,False)
    new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
    queries=standardize_old_new(old,new)
    scores=score_batch(old,queries)
    group=defaultdict(list)
    for i,r in enumerate(queries):group[r['hmt']].append(i)
    rng=np.random.default_rng(71)
    auc=[];correct=[];false_accept=[]
    for _ in range(300):
        ids=np.array([int(rng.choice(group[h])) for h in sorted(group)])
        known=np.array([queries[i]['known'] for i in ids])
        s=scores['cos_conf'][ids]
        threshold=np.quantile(s[known],.10)
        auc.append(roc_auc_score(known,s))
        correct.append(np.mean([scores['cos_hmt'][i]==queries[i]['hmt'] for i in ids[known]]))
        false_accept.append(np.mean(s[~known]>=threshold))
    out={'definition':'One sequence sampled per HMT, 300 draws, same historical reference and queries; percentiles quantify sequence draw variation, not a confidence interval over independent HMTs.',
         'known_hmts':sum(queries[group[h][0]]['known'] for h in group),
         'unknown_hmts':sum(not queries[group[h][0]]['known'] for h in group),
         'auroc_5_50_95_percentile':[round(float(x),4) for x in np.quantile(auc,[.05,.5,.95])],
         'known_hmt_accuracy_5_50_95_percentile':[round(float(x),4) for x in np.quantile(correct,[.05,.5,.95])],
         'unknown_acceptance_at_test_defined_90pct_known_coverage_5_50_95_percentile':[round(float(x),4) for x in np.quantile(false_accept,[.05,.5,.95])]}
    Path('homd_taxon_balance_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
