#!/usr/bin/env python3
"""HOMD version-shift benchmark. Run after homd_data.py.

Requires numpy, scipy, scikit-learn, biopython and edlib. No neural network.
The study is exploratory; no clinical or novel-species diagnosis is supported.
"""
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import edlib
import numpy as np
from Bio import SeqIO
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.naive_bayes import MultinomialNB
from sklearn.preprocessing import normalize

from homd_data import DATA, FILES, read_release

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260923)


def score_batch(train, query, length=None, align=False, nb=False):
    """One fit on train, batch inference. Top-hit HMT is a taxonomy baseline.

    Alignment reranks up to five cosine-shortlisted references by global edit
    similarity. For an actual best-hit BLAST baseline use BLAST itself; this is
    an alignment-assisted method, not BLAST/FOMC/QIIME.
    """
    texts = lambda rows: (r['seq'][:length] if length else r['seq'] for r in rows)
    v = CountVectorizer(analyzer='char', ngram_range=(7, 7), binary=True, dtype=np.float32)
    train_binary = v.fit_transform(texts(train))
    tx = normalize(train_binary)
    if nb:
        nbm = MultinomialNB(alpha=1.0).fit(train_binary, [r['hmt'] for r in train])
    hmt = np.array([r['hmt'] for r in train])
    result = defaultdict(list)
    for i in range(0, len(query), 80):
        batch = query[i:i+80]
        q = v.transform(texts(batch))
        sim = (normalize(q) @ tx.T).toarray()
        order = np.argsort(sim, axis=1)[:, ::-1]
        if nb:
            p = nbm.predict_proba(q)
            ix = np.argmax(p, axis=1)
            result['nb_conf'].extend(p[np.arange(len(batch)), ix].tolist())
            result['nb_hmt'].extend(nbm.classes_[ix].tolist())
        for row, ordrow, scores in zip(batch, order, sim):
            best = int(ordrow[0]); best_hmt = hmt[best]
            other = next((int(j) for j in ordrow[1:] if hmt[j] != best_hmt), best)
            bestval, margin = float(scores[best]), float(scores[best]-scores[other])
            result['cos_conf'].append(bestval)
            result['cos_margin'].append(margin)
            result['cos_hmt'].append(best_hmt)
            if align:
                hits = []
                for j in ordrow[:5]:
                    ref = train[int(j)]['seq'][:length] if length else train[int(j)]['seq']
                    seq = row['seq'][:length] if length else row['seq']
                    dist = edlib.align(seq, ref, mode='NW', task='distance')['editDistance']
                    # Normalized edit similarity; gap in length is penalized.
                    hits.append((1.0-dist/max(len(seq),len(ref)), int(j)))
                hits.sort(reverse=True)
                a, j = hits[0]
                other_a = next((s for s, k in hits[1:] if hmt[k] != hmt[j]), 0)
                result['align_conf'].append(a)
                result['align_margin'].append(a-other_a)
                result['align_hmt'].append(hmt[j])
    return {k: np.array(val) for k, val in result.items()}


def out_metrics(rows, scores, predicted, truth_hmt=True):
    known = np.array([r['known'] for r in rows], dtype=bool)
    hmt_true = np.array([r['hmt'] for r in rows])
    au = float(roc_auc_score(known, scores))
    # Interpolated threshold on evaluation set is a ranking summary only.
    threshold = float(np.quantile(scores[known], .10))
    accept = scores >= threshold
    hmt_acc = np.mean(predicted[known] == hmt_true[known]) if truth_hmt else None
    return {'known_n': int(sum(known)), 'unknown_n': int(sum(~known)),
            'known_taxon_accuracy': round(float(hmt_acc), 4) if hmt_acc is not None else None,
            'known_detection_auroc': round(au, 4),
            'unknown_acceptance_at_test_defined_90pct_known_coverage': round(float(np.mean(accept[~known])), 4),
            'known_coverage_at_test_defined_threshold': round(float(np.mean(accept[known])), 4)}


def cluster_bootstrap_auc(rows, score_a, score_b, iters=800, seed=30):
    """Cluster-bootstrap by HMT; taxon weighted rather than sequence weighted."""
    rng = np.random.default_rng(seed)
    groups = defaultdict(list)
    for i, r in enumerate(rows): groups[r['hmt']].append(i)
    known_groups = [g for g in groups if rows[groups[g][0]]['known']]
    unknown_groups = [g for g in groups if not rows[groups[g][0]]['known']]
    y = np.array([r['known'] for r in rows])
    diffs, auc_a, auc_b = [], [], []
    for _ in range(iters):
        gs = rng.choice(known_groups, len(known_groups), replace=True).tolist()
        gs += rng.choice(unknown_groups, len(unknown_groups), replace=True).tolist()
        ids = [i for g in gs for i in groups[g]]
        x = roc_auc_score(y[ids], score_a[ids]); z = roc_auc_score(y[ids], score_b[ids])
        auc_a.append(x);auc_b.append(z);diffs.append(z-x)
    def ci(v):return [round(float(x),4) for x in np.quantile(v,[.025,.5,.975])]
    return {'baseline_auc_95pct_cluster_bootstrap':ci(auc_a),
            'candidate_auc_95pct_cluster_bootstrap':ci(auc_b),
            'candidate_minus_baseline_95pct_cluster_bootstrap':ci(diffs)}


def standardize_old_new(old, new):
    old_hmt={r['hmt'] for r in old};oldseq={r['seq'] for r in old}
    q=[]
    for r in new:
        # Exclude all exact reference reuse, including same-sequence/different-HMT conflicts.
        if r['seq'] in oldseq:
            continue
        q.append(dict(r,known=r['hmt'] in old_hmt))
    return q


PRIMER_EQ=[('N',x) for x in 'ACGT']+[(k,x) for k,chars in
           {'W':'AT','B':'CGT','H':'ACT','V':'ACG'}.items() for x in chars]


def extract_v3v4(s):
    """In-silico 341F/805R V3-V4; allow 3 substitutions per degenerate primer.

    Includes primer sites and is not a wet-lab sequencing simulation.
    """
    endpoints=[]
    for primer in ('CCTACGGGNGGCWGCAG','GGATTAGATACCCBDGTAGTC'):
        match=edlib.align(primer,s,mode='HW',task='locations',k=3,
                          additionalEqualities=PRIMER_EQ)
        if not match['locations']:return None
        endpoints.append(match['locations'][0])
    f,r=endpoints
    if not (200<=f[0]<=700 and 650<=r[0]<=1200 and 300<=r[1]-f[0]<=550):return None
    return s[f[0]:r[1]+1]


def primer_extract_records(rows):
    out=[]
    for r in rows:
        amplicon=extract_v3v4(r['seq'])
        if amplicon:out.append(dict(r,seq=amplicon))
    return out


def current_holdout(new, seed, max_per_group=500):
    rng=np.random.default_rng(seed);groups=defaultdict(list)
    for r in new:groups[r['hmt']].append(r)
    taxa=np.array(sorted(groups));unknown=set(rng.choice(taxa,size=round(.2*len(taxa)),replace=False))
    train, query=[],[]
    for t in taxa:
        group=groups[t]
        if t in unknown:
            query += [dict(x,known=False) for x in group]
        elif len(group)==1:train+=group
        else:
            perm=rng.permutation(len(group));nt=max(1,round(.2*len(group)))
            query += [dict(group[i],known=True) for i in perm[:nt]]
            train += [group[i] for i in perm[nt:]]
    q=[]
    for status in (True,False):
        sub=[x for x in query if x['known']==status]
        chosen=sorted(rng.choice(len(sub),size=min(max_per_group,len(sub)),replace=False))
        q += [sub[i] for i in chosen]
    return train,q


def external_dataset(path, old, new):
    """Genome-derived, names are species not HMT: evaluate named genus only."""
    oldseq={r['seq'] for r in old};newseq={r['seq'] for r in new}
    old_genus={r['genus'] for r in old}
    genome_seen=set();out=[];stats=Counter()
    for s in SeqIO.parse(path,'fasta'):
        seq=str(s.seq).upper();parts=s.description.split('|')
        if len(parts)<3:continue
        tax=parts[2].split(';')
        if len(tax)<7 or not tax[5]:continue
        stats['all']+=1
        if seq in genome_seen:continue
        genome_seen.add(seq);stats['unique']+=1
        if seq in oldseq or seq in newseq: stats['exact_homd_overlap']+=1;continue
        if tax[5] not in old_genus:stats['genus_absent_old']+=1;continue
        if len(seq)<1200 or len(seq)>1800:stats['length_excluded']+=1;continue
        out.append({'seq':seq,'genus':tax[5],'species':tax[6],
                    'hmt':'EXTERNAL','genome':parts[1], 'known':True})
    stats['evaluable_named_genus_n']=len(out)
    return out,dict(stats)


def main():
    old=read_release(DATA/Path(FILES['old_fasta']).name,DATA/Path(FILES['old_tax']).name,False)
    new=read_release(DATA/Path(FILES['new_fasta']).name,DATA/Path(FILES['new_tax']).name,True)
    result={'sources':{'homd':'https://www.homd.org/download/download/all',
                       'genome_derived':'https://zenodo.org/records/15209015'},
            'scope':'curated reference and genome-derived sequences, no patient-level validation',
            'methods':{'7mer':'binary 7-mer cosine nearest reference',
                       'margin':'top cosine minus top cosine from another HMT',
                       'edit':'global edit similarity of top five 7-mer candidates, not exhaustive BLAST',
                       'nb':'sklearn MultinomialNB 7-mer counts (HMT posterior not calibrated)',
                       'cutoff':'test-defined 90% known coverage is descriptive ranking, not deployable calibration'}}
    q=standardize_old_new(old,new)
    print('Temporal',len(q),'known',sum(x['known'] for x in q),flush=True)
    temporal=score_batch(old,q,align=True,nb=True)
    result['temporal']={k:out_metrics(q,temporal[k+'_conf'],temporal[k+'_hmt']) for k in ('cos','align','nb')}
    result['temporal']['cos_margin']=out_metrics(q,temporal['cos_margin'],temporal['cos_hmt'])
    result['temporal']['align_margin']=out_metrics(q,temporal['align_margin'],temporal['align_hmt'])
    result['temporal']['bootstrap_align_vs_cos']=cluster_bootstrap_auc(q,temporal['cos_conf'],temporal['align_conf'])
    print('Temporal results',json.dumps(result['temporal']),flush=True)
    old_amp=primer_extract_records(old)
    new_amp=primer_extract_records(new)
    old_ampseq={r['seq'] for r in old_amp}
    amp_q=standardize_old_new(old_amp,new_amp)
    amp_scores=score_batch(old_amp,amp_q,align=False)
    result['v3v4']={'primer':'341F CCTACGGGNGGCWGCAG, 805R GACTACHVGGGTATCTAATCC',
                    'old_yield':len(old_amp),'new_yield':len(new_amp),
                    'exact_fragment_overlap_excluded':len(new_amp)-len(amp_q),
                    'cos':out_metrics(amp_q,amp_scores['cos_conf'],amp_scores['cos_hmt']),
                    'cos_margin':out_metrics(amp_q,amp_scores['cos_margin'],amp_scores['cos_hmt'])}
    print('V3V4',result['v3v4'],flush=True)
    result['current_repeats']={}
    for seed in (11,22,33,44,55):
        tr,qr=current_holdout(new,seed)
        scores=score_batch(tr,qr,nb=False,align=False)
        result['current_repeats'][str(seed)]={
            'train_n':len(tr),'known_hmts':len({x['hmt'] for x in qr if x['known']}),
            'unknown_hmts':len({x['hmt'] for x in qr if not x['known']}),
            'cos':out_metrics(qr,scores['cos_conf'],scores['cos_hmt']),
            'cos_margin':out_metrics(qr,scores['cos_margin'],scores['cos_hmt'])}
        print('Repeat',seed, result['current_repeats'][str(seed)]['cos'],flush=True)
    ext,inventory=external_dataset(DATA/'16SGOSeq_bacteria_variants.fasta',old,new)
    result['genome_derived_inventory']=inventory
    if ext:
        ext_scores=score_batch(old,ext,align=True)
        result['genome_derived_known_genus']={}
        for meth in ('cos','align'):
            nearest={r['hmt']:r['genus'] for r in old}
            correct=np.array([nearest.get(p)==r['genus'] for r,p in zip(ext,ext_scores[meth+'_hmt'])])
            result['genome_derived_known_genus'][meth]={
                'n':len(ext),'unique_genomes':len({r['genome'] for r in ext}),
                'genus_accuracy':round(float(np.mean(correct)),4),
                'taxon_label_warning':'Genome taxa derive partly from HOMD; genus names may have shifted; not an independent species-level ground truth.'}
    # Same query sequence is never in historical training where evaluated as known.
    result['checks']={'old_sequences':len(old),'new_sequences':len(new),
                      'temporal_known_hmts':len({r['hmt'] for r in q if r['known']}),
                      'temporal_unknown_hmts':len({r['hmt'] for r in q if not r['known']}),
                      'temporal_train_query_exact_overlap':len({r['seq'] for r in old}&{r['seq'] for r in q}),
                      'zenodo_fasta_sha256':hashlib.sha256((DATA/'16SGOSeq_bacteria_variants.fasta').read_bytes()).hexdigest()}
    (OUT/'homd_deep_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('saved',OUT/'homd_deep_results.json',flush=True)


if __name__=='__main__':main()
