#!/usr/bin/env python3
"""Exploratory supervised novelty scoring; v16 is untouched during model fitting.

This study is post hoc relative to the original CREO/linear-score comparison:
selection must be confirmed on a future independent cohort.
"""
import json
from pathlib import Path
import numpy as np
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from experiment import make_features, rescale
from homd_data import DATA,FILES,read_release
from homd_deep_benchmark import current_holdout, standardize_old_new, score_batch, cluster_bootstrap_auc
ROOT=Path(__file__).resolve().parent
def main():
 old=read_release(DATA/Path(FILES['old_fasta']).name,DATA/Path(FILES['old_tax']).name,False)
 new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
 eps=[]
 for seed in range(101,111):
  tr,qs=current_holdout(old,seed,max_per_group=1000)
  x,y,_=make_features(tr,qs)
  eps.append((x,y))
  print('historical',seed,len(y),flush=True)
 temporal=standardize_old_new(old,new)
 xt,yt,_=make_features(old,temporal)
 scaled,xt,low,high=rescale(eps,xt)
 train_x=np.concatenate([x for x,y in scaled[:5]])
 train_y=np.concatenate([y for x,y in scaled[:5]])
 cal_x=np.concatenate([x for x,y in scaled[5:]])
 cal_y=np.concatenate([y for x,y in scaled[5:]])
 methods={
  'logistic_l2':LogisticRegression(C=1,class_weight='balanced',max_iter=1000,random_state=0),
  'random_forest':RandomForestClassifier(n_estimators=300,min_samples_leaf=10,max_features=2,n_jobs=-1,random_state=0),
  'extra_trees':ExtraTreesClassifier(n_estimators=300,min_samples_leaf=10,max_features=2,n_jobs=-1,random_state=0),
  'hist_gradient_boosting':HistGradientBoostingClassifier(max_iter=150,max_leaf_nodes=8,learning_rate=.05,l2_regularization=2,random_state=0),
 }
 base=score_batch(old,temporal)['cos_conf']
 result={'notice':'Exploratory follow-on after inspecting linear CREO experiment; needs independent confirmation. All models train on v15.23 only.', 'train_n':len(train_y),'calibration_n':len(cal_y),'temporal_n':len(yt),'baseline_auroc':roc_auc_score(yt,base),'methods':{}}
 for name,model in methods.items():
  model.fit(train_x,train_y)
  c=model.predict_proba(cal_x)[:,1];t=model.predict_proba(xt)[:,1]
  threshold=float(np.quantile(c[cal_y],.1))
  stats={'calibration_auroc':float(roc_auc_score(cal_y,c)),
   'temporal_auroc':float(roc_auc_score(yt,t)),
   'known_coverage_at_historical_threshold':float(np.mean(t[yt]>=threshold)),
   'new_acceptance_at_historical_threshold':float(np.mean(t[~yt]>=threshold)),
   'threshold':threshold,
   'bootstrap_difference_vs_cosine':cluster_bootstrap_auc(temporal,base,t,iters=800,seed=331)}
  result['methods'][name]=stats
  print(name,json.dumps(stats),flush=True)
 (ROOT/'supervised_experiment_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
