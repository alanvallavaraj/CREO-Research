#!/usr/bin/env python3
"""Generate publication-oriented plots from archived results, no test refitting."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
res=json.loads((HERE/"creo_full_experiment_results.json").read_text())
deep=json.loads((HERE/"homd_deep_results.json").read_text())
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"figure.dpi":150,"savefig.dpi":220})
colors={"random":"#7c7c7c","de":"#4477aa","creo_no_rg":"#cc6677","creo_rg":"#228833"}

fig,axes=plt.subplots(1,2,figsize=(10.2,4.0),layout="constrained")
a,b=axes
names=["Random","DE","CREO no guidance","CREO guided"]
keys=list(colors)
for i,(key,name) in enumerate(zip(keys,names)):
 rows=res["runs"][key]
 values=[r["temporal"]["auroc"] for r in rows]
 jitter=np.random.default_rng(120+i).uniform(-.12,.12,len(values))
 a.scatter(np.full(len(values),i)+jitter,values,s=15,color=colors[key],alpha=.6)
 a.plot([i-.22,i+.22],[np.mean(values)]*2,lw=2.6,color="black")
 b.scatter([r["development_auroc"] for r in rows],values,s=17,color=colors[key],alpha=.7,label=name)
a.axhline(res["raw_cosine_baseline"]["temporal"]["auroc"],lw=1.7,color="black",ls="--",label="Raw cosine")
b.axhline(res["raw_cosine_baseline"]["temporal"]["auroc"],lw=1.3,color="black",ls="--")
b.axvline(res["raw_cosine_baseline"]["development_auroc"],lw=1.3,color="black",ls=":")
a.set_xticks(range(4),names,rotation=18,ha="right")
a.set_ylabel("Temporal novelty AUROC")
a.set_ylim(.74,.855)
a.set_title("A  •  Twenty optimiser seeds per method")
a.legend(fontsize=8,loc="lower left")
b.set_xlabel("Historical development AUROC")
b.set_ylabel("Temporal novelty AUROC")
b.set_ylim(.74,.855)
b.set_title("B  •  Development versus later release")
b.legend(fontsize=7,loc="lower left")
fig.savefig(HERE/"figure_1_optimisation_shift.png")
plt.close(fig)

repeats=deep["current_repeats"]
random_auc=np.array([v["cos"]["known_detection_auroc"] for v in repeats.values()])
temporal=deep["temporal"]["cos"]["known_detection_auroc"]
fig,ax=plt.subplots(figsize=(5.2,3.5),layout="constrained")
ax.scatter(np.arange(1,6),random_auc,color="#4477aa",s=45,label="Within-v16 taxon holdout")
ax.hlines(temporal,.7,5.3,color="#cc6677",lw=2,label="v15 → v16 held-out")
ax.set(xticks=np.arange(1,6),xlabel="Within-release split",ylabel="Known versus new HMT AUROC",ylim=(.78,1.0),xlim=(.7,5.3))
ax.legend(loc="lower left",fontsize=8)
fig.savefig(HERE/"figure_2_split_sensitivity.png")
plt.close(fig)

paired=json.loads((HERE/"paired_region_experiment_results.json").read_text())
fig,ax=plt.subplots(figsize=(5.5,3.6),layout="constrained")
bars=ax.bar(["V3–V4 fragment","Full-length 16S"],
            [paired["amp_hmt_accuracy"],paired["full_hmt_accuracy"]],
            color=["#cc6677","#4477aa"],width=.6)
for bar in bars:
    ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.015,
            f"{bar.get_height():.1%}",ha="center")
ci=paired["bootstrap"]["accuracy_full_minus_amp_95pct_hmt_cluster_bootstrap"]
ax.text(.5,.37,f"Paired gain {paired['absolute_accuracy_difference']:.1%}\n"
        f"95% HMT-cluster interval {ci[0]:.1%}–{ci[2]:.1%}",
        ha="center",va="center",bbox={"facecolor":"white","edgecolor":"#999999"})
ax.set(ylabel="Correct HMT top-hit assignment",ylim=(0,1.02),
       title=f"Same {paired['n_known']:,} known-HMT queries")
fig.savefig(HERE/"figure_3_paired_region_accuracy.png")
plt.close(fig)
print("Created three figures")
