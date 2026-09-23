# Full-length 16S improves oral taxon assignment on a paired HOMD release-shift benchmark

**Author manuscript, 23 September 2026 — for author review, not peer reviewed**  
**Alan Immanuel Benjamin Vallavaraj**  
Correspondence and affiliations: to be completed by the author before submission.

## Abstract

**Background.** We compared full-length 16S with V3–V4 for oral taxon assignment on the same later-release queries, and tested whether CREO-DE improves novelty confidence under reference shift.

**Methods.** HOMD v15.23 provided 1,015 historical reference sequences across 781 human microbial taxa (HMTs). After excluding all exact sequence reuse, v16.03 supplied 5,158 queries from previously represented HMTs and 677 queries from 74 new HMTs. Five taxon-held-out historical episodes defined the objective and five further episodes supplied threshold calibration. Four nonnegative simplex-constrained confidence weights were optimised by random search, differential evolution (DE), CREO-DE without repair guidance, and CREO-DE with guidance, each over 20 seeds and 200 objective evaluations per seed. An unoptimised binary 7-mer cosine score was the comparator.

**Results.** On 2,236 paired known-HMT query sequences, full length attained **90.03%** correct HMT assignments versus **80.23%** for V3–V4 (gain **9.79 percentage points**, 95% HMT-cluster-bootstrap interval **5.84–14.67**). Full length alone resolved 250 cases versus 31 for V3–V4 alone. Paired novelty-AUROC difference was 0.0390 (interval −0.0128 to 0.0930). On the complete full-length temporal set, untuned 7-mer cosine novelty AUROC was **0.8353**. Mean AUROCs were 0.8043 for random search, 0.8214 for DE, 0.8154 for unguided CREO and 0.8215 for guided CREO.

**Conclusions.** Full-length reference sequences improved HMT assignment in this exploratory paired release comparison. CREO and supervised confidence modelling did not improve novelty AUROC. These reference-only results require confirmation with external specimens.

**Keywords:** oral microbiome; 16S rRNA; temporal validation; novelty detection; evolutionary optimisation; CREO.

## 1. Introduction

HOMD and eHOMD supply oral and adjacent-site taxonomies and habitat-specific training collections [1–3]. Other oral catalogues and new full-length datasets reinforce that the reference collection forms part of any classifier's definition [4–6]. Differences among sequencing platforms and regions, evolving taxonomic boundaries and within-mouth variation complicate species-level interpretation [9,11,23,24,43]. A HOMD HMT identifier is a database taxon, not universally a validly named species.

Novelty detection and uncertain assignment have existing solutions, including SSUnique and IDTAXA [8,39]. RDP's naive Bayes, k-mer classifiers, sequence-search software, SpeciateIT, DeepTaxa, and habitat-focused training address related but distinct tasks [3,7,12–17,20,36]. Large benchmarks underscore the dependence of comparative results on reference and test design [10,18]. We make no claim to invent oral open-set 16S classification or to outperform those systems.

This study evaluates whether full-length information improves paired HMT assignment and whether CREO-DE repair guidance improves novelty confidence across releases. DE has several well-studied variants and constraint-handling approaches [27–31,45]; random search is a reasonable equal-budget tuning comparator [46]. Previous CREO archive results in this repository did not consistently favor repair guidance on numerical benchmarks. We included a beta-zero CREO ablation and a temporal evaluation. The paired region experiment was designed after inspecting the original results, so the positive comparison is exploratory and requires external confirmation.

## 2. Materials and methods

### 2.1 Releases and outcomes

We obtained matching full-length 16S FASTA and taxonomy files for HOMD v15.23 and v16.03. Data URLs and hashes are recorded in the code and results. The older set contains 1,015 sequences spanning 781 HMT identifiers, and the newer set contains 6,600 sequences. We defined a newer sequence as *known* if its HMT identifier occurred in v15.23, otherwise *new*. We excluded **every** newer sequence whose exact nucleotide string occurred anywhere in the older FASTA, regardless of label. The resulting temporal test has 5,158 known sequences in 530 HMTs and 677 new sequences in 74 HMTs. No newer-release labels influenced weights or thresholds. This comparison does not disentangle taxonomy revision from sampling and reference-size changes [1–3,25,26,37,38].

### 2.2 Pseudo-unknown development episodes

For each seed from 101 through 110, we randomly set aside about 20% of historical HMT groups as pseudo-unknown and assigned all their sequences to unknown queries. For other HMT groups with at least two sequences, approximately 20% formed known queries and the remainder trained the nearest-reference index. Singleton groups remained in the index. We retained at most 1,000 queries per outcome class in each episode. Seeds 101–105 supplied a mean-AUROC optimisation objective. Seeds 106–110 separately assessed historical generalisation and supplied threshold calibration. These episodes overlap in underlying historical sequences and are not independent biological replicates.

The classifier uses binary character 7-mer counts, cosine similarity, and the HMT of the highest-scoring reference. Confidence features are (i) maximum cosine, (ii) maximum global edit similarity from the five cosine-shortlisted references computed with Edlib [42], (iii) cosine difference to the highest-scoring reference of a different HMT, and (iv) log(1 + best-HMT reference count). This is a shortlisted edit feature, not exhaustive alignment, BLAST or VSEARCH [20]. Each feature was clipped after linear scaling from historical pooled 5th to 95th percentiles. A candidate confidence score is \(s(x)=\sum_{j=1}^{4}w_jx_j\), where \(w_j\geq0\) and \(\sum_jw_j=1\). The untuned comparator uses raw cosine alone.

### 2.3 Equal-budget optimisers

All four searches used population size 20 and exactly 200 objective evaluations for each optimiser seed 7701–7720, including initial evaluations. They evaluated the same historical feature matrices. Every method used the same Euclidean projection \(P\) to the simplex before fitness calculation. Random search sampled Gaussian candidate coordinates. Standard DE used rand/1 mutation \(F=0.7\), binomial crossover \(CR=0.9\), and greedy replacement [45]. Both CREO-DE variants added best-repaired-point attraction \(0.18(g_{\rm best}-x_i)\) and Gaussian perturbation with scale 0.02. The guided variant also added \(0.90[P(x_i)-x_i]\); the ablation set this coefficient to zero. Mutation retains raw off-simplex population vectors, making repair displacement nonzero in applicable candidates. The parameters follow the CREO-DE driver in this repository, but the simplex projection is a task-specific adaptation. General feasibility handling, self-adaptive DE, alternative population strategies and automated tuning are established contexts [27–32,34,46].

For each run, the 10th percentile of *historical calibration* known-query scores, pooled across five held-apart episode seeds, became its frozen acceptance cutoff. The old-calibrated cutoff was applied once to the temporal test. We report AUROC, known coverage, and new-HMT acceptance at this cutoff. A separate test-defined 90%-known-coverage point is strictly descriptive and cannot be used for deployment. Top-hit HMT assignments are always cosine based, so changing weights does not change their HMT accuracy. An exploratory secondary study trained four supervised confidence models (regularised logistic regression, random forest, extremely randomised trees and histogram gradient boosting) on the same old-release development episodes with no new-release labels.

### 2.4 Secondary analyses and uncertainty

The reference-only baseline analysis evaluated cosine, five-candidate global edit, an untuned multinomial naive-Bayes classifier and confidence margins on the same temporal records [36,42]. Five within-v16.03 HMT-withholding splits assessed random-split optimism. An in-silico V3–V4 analysis extracted fragments using 341F/805R degenerate primers with up to three substitutions at each primer, then removed exact historical fragment reuse [43]. These fragments are **not** experimental amplicon reads. For the paired analysis we linked every selected fragment by sequence-record ID to its full-length v16.03 original, then classified the identical 2,742 query IDs against the corresponding historical full-length and extracted-fragment reference sets. We checked that no full-length query exactly matched an old full-length reference. The comparison is restricted to the primer-extractable, nonidentical-fragment subset. A genome-derived 16SGOSeq collection [6] supplied a genus concordance sanity check, but its source material overlaps HOMD reference curation. Read processing tools, sequencing-error filters and broader pipeline comparisons exist outside this reference-only study [19,21,22,40,41,44].

For each search method we selected the run with the median *historical development* AUROC, without using temporal labels to select it. An 800-draw bootstrap separately resampled known and new HMT clusters with replacement, carrying every query sequence for the selected HMT. A separate 1,500-draw HMT-cluster bootstrap resampled the same paired query indices in both representations, yielding a paired accuracy difference among known HMTs and AUROC difference among known and new HMTs. Percentile intervals estimate conditional taxon-cluster variation for these releases and fixed fitted runs. The twenty optimiser seeds share the same biological test and must not be treated as twenty independent test cohorts [33–35].

## 3. Results

### 3.1 Raw baseline and sequence-region checks

Raw cosine yielded temporal novelty AUROC **0.8353**, known-HMT top-hit accuracy **0.8984**, historical-cutoff known coverage **0.9533**, and new-HMT acceptance **0.4801**. Its test-defined 90%-coverage acceptance was **0.4165**, which is a retrospective ranking summary at a different operating point. Shortlisted global-edit AUROC was **0.8116**, untuned naive-Bayes posterior AUROC **0.5472**, and cosine margin AUROC **0.7412**. Those implementations do not stand in for optimised RDP, VSEARCH, BLAST or deep-learning comparators.

Five within-v16.03 pseudo-novelty splits produced cosine AUROCs **0.9289–0.9604**, compared with **0.8353** on the temporal test. The in-silico V3–V4 analysis included 2,236 known and 506 new fragments after exact fragment exclusions, yielding AUROC **0.8477** and known-HMT top-hit accuracy **0.8023**. For the *same* query IDs, full length achieved known-HMT accuracy **0.9003**, a paired gain **0.0979** (HMT-cluster-bootstrap 95% interval **0.0584 to 0.1467**). Both methods were correct on 1,763 known queries; full length alone was correct on 250, V3–V4 alone on 31, and both missed 192. Full-length AUROC on this paired subset was **0.8867**, versus **0.8477** for V3–V4. The AUROC difference interval **−0.0128 to 0.0930** includes zero. The paired subset has a different case mix from the full temporal set and must not be compared directly with its AUROC 0.8353. Newly added HMTs can share identical fragments with historical HMTs, an irreducible label ambiguity for this amplicon.

### 3.2 CREO ablation and optimisation results

| Confidence method | Historical objective AUROC, mean | Historical calibration AUROC, mean | Temporal AUROC, mean (range) | Known coverage at old cutoff, mean | New-HMT acceptance at old cutoff, mean |
|:--|--:|--:|--:|--:|--:|
| Raw cosine, no tuning | 0.9022 | — | **0.8353** | 0.9533 | 0.4801 |
| Random search | 0.9159 | 0.9128 | 0.8043 (0.7583–0.8366) | 0.9584 | 0.4678 |
| Differential evolution | 0.9177 | 0.9161 | 0.8214 (0.7824–0.8355) | 0.9542 | 0.4659 |
| CREO-DE, no guidance | 0.9178 | 0.9162 | 0.8154 (0.7628–0.8346) | 0.9548 | 0.4664 |
| CREO-DE, repair guidance | 0.9173 | 0.9162 | 0.8215 (0.7767–0.8369) | 0.9536 | 0.4663 |

All tuned methods improved the historical objective; none improved mean temporal AUROC over raw cosine. Guided CREO's mean test AUROC was **0.0138 below baseline**. Its mean advantage over the no-guidance ablation was **0.0061 AUROC**, a descriptive difference over optimiser seeds on identical test records. For the median-development-run guided CREO, the taxon-bootstrap AUROC difference versus baseline had 2.5th/median/97.5th percentiles **−0.0587/−0.0027/0.0407**; for unguided CREO the interval was **−0.0454 to 0.0161**. Neither shows reliable superiority. One guided optimiser seed scored 0.8369, but choosing it by temporal outcome would compromise the comparison. Lower false acceptance at the frozen threshold co-occurs with a changed known-coverage rate and does not reverse the AUROC result. The four exploratory supervised models also failed to improve temporal AUROC: logistic regression **0.7844**, random forest **0.7655**, extremely randomised trees **0.7836**, and histogram gradient boosting **0.7885**, although their old-release calibration AUROCs ranged **0.9213–0.9645**.

### 3.3 Dataset dependence

The genome-derived set gave **99.70%** genus concordance across 3,947 eligible sequences from 1,744 genome identifiers; underlying label sourcing prevents treating this as independent clinical or species validation. Repeating a one-sequence-per-HMT draw 300 times on the temporal set gave raw-cosine AUROC 5th/median/95th percentiles **0.8212/0.8360/0.8503**. These quantify within-dataset sequence-choice variation, not prospective population uncertainty. Other cohort-level, technology-level and reference-level differences remain outside this experiment [5,9,10,18,23,25,26,34,35,38].

## 4. Discussion

The paired temporal comparison provides a positive assignment result: full-length references improved HMT accuracy by almost ten percentage points on identical reference queries. The value is a reproducible, paired temporal quantification of an expected region-resolution effect [5,9,11,24], not a claim that longer genes were first discovered to be informative. The positive comparison was added after examining the first data, so confirmation on independent specimens or isolates is essential. Novelty discrimination has a distinct outcome and its paired gain interval includes zero.

The historical novelty objective transported poorly to a later reference release. The result does not prove why the transfer failed: a pseudo-unknown group drawn from the old release may differ from HMTs subsequently added, and the newer release changes taxonomy and representation simultaneously. The repair term was operative under simplex projection, yet the guided algorithm did not establish improvement over raw cosine. Mean performance across seeds is reported to avoid choosing a favorable run using the test labels. Equal objective-evaluation budgets do not imply equal end-to-end CPU costs; precomputed alignments were shared across optimisers.

The sole temporal pair is a proxy for prospective novelty, not independent patient, culture-confirmed isolate, or mock-community testing. Near duplicates can remain after excluding exact sequences. HMT abundance induces sequence weighting and some taxa have few references. Historical episodes share records. Naive Bayes is one untuned configuration; shortlisted alignment is not exhaustive. AUROC obscures operating thresholds: the historical 90% known calibration produced about 95% known coverage after release shift. The V3–V4 analysis is in silico only, and primer extraction changes the eligible case mix. Pairing holds query IDs fixed within that selected subset. No sequencing error, PCR bias, clinical outcome, patient leakage, computation-time or abundance calibration study was performed [5,9,19,41,43].

A more decisive follow-up would freeze algorithms before acquisition of independently verified oral isolates or mock reads; include real primer-defined amplicons, site/subject split, and equal-input comparisons with RDP, IDTAXA, SpeciateIT, DeepTaxa and exhaustive sequence search [7,8,15,36,39]. CREO becomes attractive only if coupled constraints materially shape the task and a beta-only ablation improves independently held-out cohorts at comparable budgets. Its current application does not add demonstrated accuracy to this classifier.

## 5. Data availability and declarations

The scripts, numerical results, 20-run traces, input hashes, pinned dependencies and figures are in this repository under research/homd_release_shift. Public reference FASTA files are downloaded by homd_data.py and fetch_external_data.py and are excluded from version control. Execute homd_deep_benchmark.py, paired_region_experiment.py and supervised_experiment.py for secondary analyses and experiment.py --seeds 20 --budget 200 for the CREO comparison. The external genome-derived download has a hardcoded SHA-256 check; compare the HOMD hashes in the results to verify those downloads before replication. No patient-level data were analysed. Author affiliation, funding, competing interests and intended journal require author verification. This draft has not been peer reviewed or published.

## References

1. Chen, T.; Yu, W.-H.; Izard, J.; Baranova, O. V.; et al. (2010). The Human Oral Microbiome Database: a web accessible resource for investigating oral microbe taxonomic and genomic information. *Database*. https://doi.org/10.1093/database/baq013

2. Escapa, Isabel F.; Chen, Tsute; Huang, Yanmei; Gajare, Prasad; et al. (2018). New Insights into Human Nostril Microbiome from the Expanded Human Oral Microbiome Database (eHOMD): a Resource for the Microbiome of the Human Aerodigestive Tract. *mSystems*. https://doi.org/10.1128/msystems.00187-18

3. F. Escapa, Isabel; Huang, Yanmei; Chen, Tsute; Lin, Maoxuan; et al. (2020). Construction of habitat-specific training sets to achieve species-level assignment in 16S rRNA gene datasets. *Microbiome*. https://doi.org/10.1186/s40168-020-00841-w

4. Griffen, Ann L.; Beall, Clifford J.; Firestone, Noah D.; Gross, Erin L.; et al. (2011). CORE: A Phylogenetically-Curated 16S rDNA Database of the Core Oral Microbiome. *PLoS ONE*. https://doi.org/10.1371/journal.pone.0019051

5. Esberg, Anders; Fries, Niklas; Haworth, Simon; Johansson, Ingegerd (2024). Saliva microbiome profiling by full-gene 16S rRNA Oxford Nanopore Technology versus Illumina MiSeq sequencing. *npj Biofilms and Microbiomes*. https://doi.org/10.1038/s41522-024-00634-1

6. Vázquez-González, Lara; Regueira-Iglesias, Alba; Balsa-Castro, Carlos; Tomás, Inmaculada; et al. (2025). A curated bacterial and archaeal 16S rRNA Gene Oral Sequences dataset. *Scientific Data*. https://doi.org/10.1038/s41597-025-05050-4

7. Salah, Rana; AbdElaal, Khlood R; Ghonaim, Lobna; Awe, Olaitan I; et al. (2026). DeepTaxa: a hybrid CNN-BERT framework for 16S rRNA taxonomic classification. *Bioinformatics Advances*. https://doi.org/10.1093/bioadv/vbag166

8. Murali, Adithya; Bhargava, Aniruddha; Wright, Erik S. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. *Microbiome*. https://doi.org/10.1186/s40168-018-0521-5

9. Wagner, Josef; Coupland, Paul; Browne, Hilary P.; Lawley, Trevor D.; et al. (2016). Evaluation of PacBio sequencing for full-length bacterial 16S rRNA gene classification. *BMC Microbiology*. https://doi.org/10.1186/s12866-016-0891-4

10. Odom, Aubrey R.; Faits, Tyler; Castro-Nallar, Eduardo; Crandall, Keith A.; et al. (2023). Metagenomic profiling pipelines improve taxonomic classification for 16S amplicon sequencing data. *Scientific Reports*. https://doi.org/10.1038/s41598-023-40799-x

11. Martínez-Porchas, Marcel; Villalpando-Canchola, Enrique; Vargas-Albores, Francisco (2016). Significant loss of sensitivity and specificity in the taxonomic classification occurs when short 16S rRNA gene sequences are used. *Heliyon*. https://doi.org/10.1016/j.heliyon.2016.e00170

12. Chaudhary, Nikhil; Sharma, Ashok K.; Agarwal, Piyush; Gupta, Ankit; et al. (2015). 16S Classifier: A Tool for Fast and Accurate Taxonomic Classification of 16S rRNA Hypervariable Regions in Metagenomic Datasets. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0116106

13. Vinje, Hilde; Liland, Kristian Hovde; Almøy, Trygve; Snipen, Lars (2015). Comparing K-mer based methods for improved classification of 16S sequences. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-015-0647-4

14. Ziemski, Michal; Wisanwanichthan, Treepop; Bokulich, Nicholas A.; Kaehler, Benjamin D. (2021). Beating Naive Bayes at Taxonomic Classification of 16S rRNA Gene Sequences. *Frontiers in Microbiology*. https://doi.org/10.3389/fmicb.2021.644487

15. Holm, Johanna B.; Gajer, Pawel; Ravel, Jacques (2024). SpeciateIT and vSpeciateDB: novel, fast, and accurate per sequence 16S rRNA gene taxonomic classification of vaginal microbiota. *BMC Bioinformatics*. https://doi.org/10.1186/s12859-024-05930-3

16. Myer, Phillip R.; McDaneld, Tara G.; Kuehn, Larry A.; Dedonder, Keith D.; et al. (2020). Classification of 16S rRNA reads is improved using a niche-specific database constructed by near-full length sequencing. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0235498

17. Hsieh, Yu-Peng; Hung, Yuan-Mao; Tsai, Mong-Hsun; Lai, Liang-Chuan; et al. (2022). 16S-ITGDB: An Integrated Database for Improving Species Classification of Prokaryotic 16S Ribosomal RNA Sequences. *Frontiers in Bioinformatics*. https://doi.org/10.3389/fbinf.2022.905489

18. Seppey, Mathieu; Benavides, Andres; Berkeley, Matthew R.; Manni, Mosè; et al. (2026). LEMMIv2: benchmarking framework for metagenomic and 16S amplicon profilers with a catalogue of evaluated tools. *Genome Biology*. https://doi.org/10.1186/s13059-026-04089-9

19. Callahan, Benjamin J; McMurdie, Paul J; Rosen, Michael J; Han, Andrew W; et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. https://doi.org/10.1038/nmeth.3869

20. Edgar, Robert C. (2010). Search and clustering orders of magnitude faster than BLAST. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btq461

21. Wood, Derrick E.; Lu, Jennifer; Langmead, Ben (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*. https://doi.org/10.1186/s13059-019-1891-0

22. Caporaso, J Gregory; Kuczynski, Justin; Stombaugh, Jesse; Bittinger, Kyle; et al. (2010). QIIME allows analysis of high-throughput community sequencing data. *Nature Methods*. https://doi.org/10.1038/nmeth.f.303

23. Mukherjee, Chiranjit; Beall, Clifford J.; Griffen, Ann L.; Leys, Eugene J. (2018). High-resolution ISR amplicon sequencing reveals personalized oral microbiome. *Microbiome*. https://doi.org/10.1186/s40168-018-0535-z

24. Hackmann, Timothy J. (2025). Setting new boundaries of 16S rRNA gene identity for prokaryotic taxonomy. *International Journal of Systematic and Evolutionary Microbiology*. https://doi.org/10.1099/ijsem.0.006747

25. Rashidi, Armin; Gem, Hakan; McLean, Jeffrey S.; Kerns, Kristopher; et al. (2024). Multi-cohort shotgun metagenomic analysis of oral and gut microbiota overlap in healthy adults. *Scientific Data*. https://doi.org/10.1038/s41597-024-02916-x

26. Cha, Jun Hyung; Kim, Nayeon; Ma, Junyeong; Lee, Sungho; et al. (2025). A high-quality genomic catalog of the human oral microbiome broadens its phylogeny and clinical insights. *Cell Host & Microbe*. https://doi.org/10.1016/j.chom.2025.10.001

27. Deb, Kalyanmoy (2000). An efficient constraint handling method for genetic algorithms. *Computer Methods in Applied Mechanics and Engineering*. https://doi.org/10.1016/s0045-7825(99)00389-8

28. Brest, Janez; Greiner, Sao; Boskovic, Borko; Mernik, Marjan; et al. (2006). Self-Adapting Control Parameters in Differential Evolution: A Comparative Study on Numerical Benchmark Problems. *IEEE Transactions on Evolutionary Computation*. https://doi.org/10.1109/tevc.2006.872133

29. Qin, A.K.; Huang, V.L.; Suganthan, P.N. (2009). Differential Evolution Algorithm With Strategy Adaptation for Global Numerical Optimization. *IEEE Transactions on Evolutionary Computation*. https://doi.org/10.1109/tevc.2008.927706

30. Liang, J.J.; Qin, A.K.; Suganthan, P.N.; Baskar, S. (2006). Comprehensive learning particle swarm optimizer for global optimization of multimodal functions. *IEEE Transactions on Evolutionary Computation*. https://doi.org/10.1109/tevc.2005.857610

31. Das, Swagatam; Abraham, Ajith; Chakraborty, Uday K.; Konar, Amit (2009). Differential Evolution Using a Neighborhood-Based Mutation Operator. *IEEE Transactions on Evolutionary Computation*. https://doi.org/10.1109/tevc.2008.2009457

32. Vincent, Amala Mary; Jidesh, P. (2023). An improved hyperparameter optimization framework for AutoML systems using evolutionary algorithms. *Scientific Reports*. https://doi.org/10.1038/s41598-023-32027-3

33. DeLong, Elizabeth R.; DeLong, David M.; Clarke-Pearson, Daniel L. (1988). Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. *Biometrics*. https://doi.org/10.2307/2531595

34. Cappellato, Marco; Baruzzo, Giacomo; Di Camillo, Barbara (2022). Investigating differential abundance methods in microbiome data: A benchmark study. *PLOS Computational Biology*. https://doi.org/10.1371/journal.pcbi.1010467

35. Mattiello, Federico; Verbist, Bie; Faust, Karoline; Raes, Jeroen; et al. (2016). A web application for sample size and power calculation in case-control microbiome studies. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btw099

36. Wang, Qiong; Garrity, George M.; Tiedje, James M.; Cole, James R. (2007). Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. *Applied and Environmental Microbiology*. https://doi.org/10.1128/aem.00062-07

37. Quast, Christian; Pruesse, Elmar; Yilmaz, Pelin; Gerken, Jan; et al. (2012). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gks1219

38. McDonald, Daniel; Jiang, Yueyu; Balaban, Metin; Cantrell, Kalen; et al. (2023). Greengenes2 unifies microbial data in a single reference tree. *Nature Biotechnology*. https://doi.org/10.1038/s41587-023-01845-1

39. Lynch, Michael D. J.; Neufeld, Josh D. (2016). SSUnique: Detecting Sequence Novelty in Microbiome Surveys. *mSystems*. https://doi.org/10.1128/msystems.00133-16

40. Schloss, Patrick D.; Westcott, Sarah L.; Ryabin, Thomas; Hall, Justine R.; et al. (2009). Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. *Applied and Environmental Microbiology*. https://doi.org/10.1128/aem.01541-09

41. Bokulich, Nicholas A; Subramanian, Sathish; Faith, Jeremiah J; Gevers, Dirk; et al. (2013). Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencing. *Nature Methods*. https://doi.org/10.1038/nmeth.2276

42. Šošić, Martin; Šikić, Mile (2017). Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btw753

43. Klindworth, Anna; Pruesse, Elmar; Schweer, Timmy; Peplies, Jörg; et al. (2012). Evaluation of general 16S ribosomal RNA gene PCR primers for classical and next-generation sequencing-based diversity studies. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gks808

44. Bolyen, Evan; Rideout, Jai Ram; Dillon, Matthew R.; Bokulich, Nicholas A.; et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*. https://doi.org/10.1038/s41587-019-0209-9

45. Storn, Rainer; Price, Kenneth (1997). Differential Evolution – A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces. *Journal of Global Optimization*. https://doi.org/10.1023/A:1008202821328

46. Bergstra, James; Bengio, Yoshua (2012). Random Search for Hyper-Parameter Optimization. *Journal of Machine Learning Research 13:281–305*. https://jmlr.org/papers/v13/bergstra12a.html
