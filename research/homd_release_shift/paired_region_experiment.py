#!/usr/bin/env python3
"""Paired, exploratory full-length versus V3-V4 query comparison."""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from homd_data import DATA, FILES, read_release
from homd_deep_benchmark import primer_extract_records,standardize_old_new,score_batch
HERE=Path(__file__).resolve().parent

def bootstrap(rows, correct_full,correct_amp,score_full,score_amp, n=1500,seed=482):
 rng=np.random.default_rng(seed)
 g=defaultdict(list)
 for i,r in enumerate(rows):g[r['hmt']].append(i)
 known=[x for x in g if rows[g[x][0]]['known']]
 new=[x for x in g if not rows[g[x][0]]['known']]
 y=np.array([r['known'] for r in rows],dtype=bool)
 difacc=[];difauc=[]
 for k in range(n):
  ks=rng.choice(known,size=len(known),replace=True)
  ns=rng.choice(new,size=len(new),replace=True)
  ik=np.array([i for name in ks for i in g[name]])
  ix=np.concatenate([ik,np.array([i for name in ns for i in g[name]])])
  difacc.append(np.mean(correct_full[ik])-np.mean(correct_amp[ik]))
  difauc.append(roc_auc_score(y[ix],score_full[ix])-roc_auc_score(y[ix],score_amp[ix]))
 return {'accuracy_full_minus_amp_95pct_hmt_cluster_bootstrap':np.quantile(difacc,[.025,.5,.975]).tolist(),
         'novelty_auc_full_minus_amp_95pct_hmt_cluster_bootstrap':np.quantile(difauc,[.025,.5,.975]).tolist()}
def main():
 old=read_release(DATA/Path(FILES['old_fasta']).name,DATA/Path(FILES['old_tax']).name,False)
 new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
 oldamp=primer_extract_records(old);newamp=primer_extract_records(new)
 qa=standardize_old_new(oldamp,newamp)
 lookup={r['id']:r for r in new}
 qf=[dict(lookup[r['id']],known=r['known']) for r in qa]
 assert len(qf)==len(qa) and all(a['id']==b['id'] for a,b in zip(qf,qa))
 old_sequences={r['seq'] for r in old}
 assert not any(r['seq'] in old_sequences for r in qf)
 sa=score_batch(oldamp,qa);sf=score_batch(old,qf)
 y=np.array([r['known'] for r in qa],bool)
 truth=np.array([r['hmt'] for r in qa])
 cf=sf['cos_hmt']==truth;ca=sa['cos_hmt']==truth
 f=float(np.mean(cf[y]));a=float(np.mean(ca[y]))
 both=np.sum(cf[y]&ca[y]);fonly=np.sum(cf[y]&~ca[y]);aonly=np.sum(~cf[y]&ca[y]);neither=np.sum(~cf[y]&~ca[y])
 result={'exploratory_notice':'Pairing designed after inspecting initial benchmark; needs external validation. Same query IDs and old HMT labels, with full and V3-V4 reference representations.',
   'n_known':int(sum(y)),'n_new':int(sum(~y)),'known_hmts':len(set(truth[y])),'new_hmts':len(set(truth[~y])),
   'full_hmt_accuracy':f,'amp_hmt_accuracy':a,'absolute_accuracy_difference':f-a,
   'paired_accuracy_contingency':{'both_correct':int(both),'full_only_correct':int(fonly),'amp_only_correct':int(aonly),'neither_correct':int(neither)},
   'full_novelty_auroc':float(roc_auc_score(y,sf['cos_conf'])),'amp_novelty_auroc':float(roc_auc_score(y,sa['cos_conf'])),
   'bootstrap':bootstrap(qf,cf,ca,sf['cos_conf'],sa['cos_conf'])}
 (HERE/'paired_region_experiment_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
