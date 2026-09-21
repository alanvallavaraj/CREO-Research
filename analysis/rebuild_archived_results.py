#!/usr/bin/env python3
"""Rebuild descriptive manuscript result tables from the ORIGINAL archived CSVs.

No CEC objective functions are called and no new optimisation runs are claimed.
Run: python analysis/rebuild_archived_results.py
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

METHODS = ("GA", "PSO", "DE", "GWO", "CREO-Repair", "CREO-RG", "CREO-DE")
MAIN_FIVE = ("GA", "PSO", "DE", "GWO", "CREO-DE")
ABLATION = ("DE", "CREO-DE-NoRG", "CREO-DE")
FUNCTIONS = {f"F{i}" for i in range(1,31) if i != 2}


def read_archive(p:Path, d:int, algos:tuple[str,...]):
    if not p.is_file(): raise FileNotFoundError(f"Missing ORIGINAL archived CSV: {p}")
    df=pd.read_csv(p)
    required={'Function','Dimension','Algorithm','Mean','Std','Best','Worst'}
    if not required.issubset(df.columns): raise ValueError(f"Missing columns in {p}: {required-set(df.columns)}")
    if df[['Function','Algorithm']].duplicated().any():raise ValueError(f"Duplicate function/algorithm rows in {p}")
    if not df['Dimension'].eq(d).all():raise ValueError(f"Incorrect dimension in {p}")
    if set(df.Function)!=FUNCTIONS or set(df.Algorithm)!=set(algos):
        raise ValueError(f"Wrong functions/algorithms in {p}: functions={set(df.Function)^FUNCTIONS}, algos={set(df.Algorithm)^set(algos)}")
    if len(df)!=len(FUNCTIONS)*len(algos):raise ValueError(f"Wrong row count in {p}")
    if df[['Mean','Std','Best','Worst']].isna().any().any():raise ValueError(f"Missing numbers in {p}")
    if not ((df['Best']<=df['Mean'])&(df['Mean']<=df['Worst'])).all(): raise ValueError(f"Inconsistent mean/best/worst in {p}")
    return df


def ranks(df:pd.DataFrame, d:int):
    w=df.copy()
    w['Rank']=w.groupby('Function')['Mean'].rank(method='average',ascending=True)
    w['BestCount']=w['Mean'].eq(w.groupby('Function')['Mean'].transform('min')).astype(int)
    w['WorstCount']=w['Mean'].eq(w.groupby('Function')['Mean'].transform('max')).astype(int)
    out=w.groupby('Algorithm').agg(MeanRank=('Rank','mean'),RankSD=('Rank',lambda s:s.std(ddof=1)),BestCount=('BestCount','sum'),WorstCount=('WorstCount','sum')).reset_index()
    out['Dimension']=d
    return out[['Dimension','Algorithm','MeanRank','RankSD','BestCount','WorstCount']].sort_values(['MeanRank','Algorithm']).reset_index(drop=True)


def wtl(piv:pd.DataFrame, one:str, other:str, d:int, label:str):
    a=piv[one].to_numpy();b=piv[other].to_numpy()
    # Strict comparisons of FULL-precision archived means; no rounding-based pseudo-ties.
    return dict(Dimension=d,Dataset=label,Comparison=f'{one} vs {other}',Wins=int((a<b).sum()),Ties=int((a==b).sum()),Losses=int((a>b).sum()))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output-dir',type=Path,default=None)
    args=ap.parse_args();root=args.repo_root.resolve()
    out=(args.output_dir or (root/'analysis'/'generated')).resolve();out.mkdir(parents=True,exist_ok=True)
    all7=[];five=[];abla=[];all_wtl=[];matrices=[]
    for d in (10,30):
        main_df=read_archive(root/f'cec17_results_D{d}'/'all_benchmark_results.csv',d,METHODS)
        abla_df=read_archive(root/f'cec17_ablation_small_D{d}'/'all_ablation_results.csv',d,ABLATION)
        all7.append(ranks(main_df,d));five.append(ranks(main_df[main_df.Algorithm.isin(MAIN_FIVE)],d));abla.append(ranks(abla_df,d))
        mm=main_df.pivot(index='Function',columns='Algorithm',values='Mean').loc[sorted(FUNCTIONS,key=lambda s:int(s[1:]))]
        aa=abla_df.pivot(index='Function',columns='Algorithm',values='Mean').loc[sorted(FUNCTIONS,key=lambda s:int(s[1:]))]
        for other in ('DE','PSO','GA','GWO','CREO-Repair','CREO-RG'):
            all_wtl.append(wtl(mm,'CREO-DE',other,d,'main'))
        for one,other in [('CREO-DE','CREO-DE-NoRG'),('CREO-DE','DE'),('CREO-DE-NoRG','DE')]:
            all_wtl.append(wtl(aa,one,other,d,'ablation'))
        for fn in sorted(FUNCTIONS,key=lambda s:int(s[1:])):
            row={'Function':fn,'Dimension':d}
            for a in METHODS:row[a]=float(mm.loc[fn,a])
            matrices.append(row)
    all7=pd.concat(all7,ignore_index=True);five=pd.concat(five,ignore_index=True);abla=pd.concat(abla,ignore_index=True)
    all7.to_csv(out/'main_seven_algorithm_ranks.csv',index=False)
    five.to_csv(out/'main_five_algorithm_ranks.csv',index=False)
    abla.to_csv(out/'ablation_three_algorithm_ranks.csv',index=False)
    pd.DataFrame(all_wtl).to_csv(out/'pairwise_win_tie_loss.csv',index=False)
    pd.DataFrame(matrices).to_csv(out/'main_function_wise_unrounded_means.csv',index=False)
    # Record traceability without suggesting this rebuild represents new optimiser runs.
    summary=[]
    for d in (10,30):
        summary.extend([f'### Dimension D={d}', 'Main seven-method ranks (per-function means, full precision):', all7[all7.Dimension==d].to_string(index=False),
                        'Ablation ranks:',abla[abla.Dimension==d].to_string(index=False)])
    summary.extend(['', 'All values are REANALYSED from original archived CSVs, not freshly rerun optimisation outcomes.',
                    'Main and ablation experiments are separate runs; do not merge their DE/CREO-DE records.',
                    'The repository does not include the external CEC2017 source archive or individual-run results.'])
    (out/'REANALYSIS_REPORT.txt').write_text('\n'.join(summary)+'\n',encoding='utf-8')
    print('\n'.join(summary))

if __name__=='__main__':main()
