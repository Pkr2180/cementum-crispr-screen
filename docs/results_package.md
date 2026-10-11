# Results package — every step, with figures and tables

A CRISPR knockout screen for acellular cementum regeneration. This document walks
the workflow step by step, gives the figure and table for each, and closes with
an explicit comparison of **what was designed initially** against **what the data
forced the design to become**.

Final deliverable: **9,631 oligos — 8,127 targeting guides across 861 genes plus
1,504 controls** ([library_complete.csv]({{artifact:f2bc054b-ecf0-44e0-8cb3-eb6082dcbecc}}), [oligo_pool_order.csv]({{artifact:aff3aaea-e46f-45f9-a18d-0739e762f887}})).

---

## The initial design, and what changed

This matters more than any single number, so it comes first. Twelve design
elements were specified before any computation. Seven of them changed, each
because a computed result contradicted the assumption behind it.

| Element | Initially proposed | Final design | What forced the change |
|---|---|---|---|
| Representation | 500–1000× | **200×** | 1/cells is negligible against measured φ = 0.30 |
| Guides per gene | 8 | **10** (9.4 realised) | 8 × 3 gives 57 % power; 10 × 4 gives 80 % |
| Replicates | 3 | **4** | same |
| Library scope | genome-wide + focused | **focused only** | genome-wide decile sort needs 2.8 B cells |
| Readout | ALPL/IBSP ratio, one line | **ratio + parallel osteoblast arm** | tooth vs bone AUROC 0.582 — not separable in one line |
| Primary anchor | transcriptional marker ratio | **PPi/Pi axis** | matrix markers are shared with bone |
| Hit-selection metric | aggregate CFD specificity | **max coding-exon CFD** | 93 % of sites sit at 3 mm; aggregate CFD does not discriminate |
| On-target models | RS3 + DeepSpCas9 + DeepHF | **RS3 (2 scaffolds) + Doench 2014** | deep models need R/basilisk, unavailable on this host |
| Controls | 1,000 NTC + 500 assorted | **1,000 verified NTC + 254 + 250** | NTC must be verified at ≤2 mm; ≤3 mm unattainable |
| L4 model | fine-tuned GEARS/scGPT | **gated on a Perturb-seq arm** | sorted screens yield no perturbation→transcriptome pairs |
| Module balance | as curated | **rebalanced, +39 genes** | curation leaned bone 2.45 : 1 |
| Off-target budget | ≤3 mm, NGG/NAG/NGA | ≤3 mm, NGG/NAG/NGA | unchanged — delivered as specified |

Full table: [table5_initial_vs_final.csv]({{artifact:56bade5a-692f-47cd-8270-f790b4517ee2}})

![Screen architecture]({{artifact:art_c9f60c1b-c8ac-4268-918e-b5690abec195}})

---

## Step 1 — Target space

829 genes assembled from a chromosome-partitioned BioMart query for the kinome
(GO:0004672) plus an eleven-submodule literature curation of the cementum module,
resolved to Ensembl IDs with 1:1 mouse orthologues where they exist (802/829;
*CEMP1* has none).

**Table:** [target_manifest.csv]({{artifact:00111918-6082-4f4f-b01e-2ab86d57a288}}) · [table2_target_space.csv]({{artifact:6736669c-fae4-4183-9d91-50fa927eb7df}})

| Tier | Genes | Guides | Guides/gene |
|---|---|---|---|
| B — cementum module (11 submodules) | 341 | 3,243 | 9.5 |
| C — kinome | 481 | 4,508 | 9.4 |
| B — tooth-programme additions (step 9) | 39 | 376 | 9.6 |

## Step 2 — Genome and annotation

hg38 downloaded from UCSC and indexed to one byte per base: **3,088,269,832 bp**
across 24 primary chromosomes, chr1 verified at 248,956,422. Transcript models
from `ncbiRefSeqCurated` (78,435 transcripts), canonical transcript from
`ncbiRefSeqSelect`. 828 of 829 genes matched.

## Step 3 — Guide enumeration

All NGG protospacers on both strands within coding exons, filtered on cut-site
position (5–65 % of CDS), constitutive-exon fraction ≥0.9, GC 0.30–0.80, no
TTTT, no homopolymer ≥5, no BsmBI/BsaI site, no ambiguous base.

**139,144 candidates across 826 genes** (median 119/gene). 94 genes could not
reach the 90 % constitutive bar — they carry up to 131 RefSeq transcripts — and
went through a recorded per-gene relaxation cascade; 734 passed at the strictest
setting, 10 needed the most permissive.

![Pipeline funnel and target space]({{artifact:art_3fd47487-ed8d-4471-9e7c-1c8120ca1e65}})

**Table:** [table1_pipeline_funnel.csv]({{artifact:f0387ff2-e5f0-44e6-bccd-e5228c0f452c}}) — 94.4 % attrition from candidates to final selection.

## Step 4 — On-target activity

Three scores, rank-aggregated; selection never uses one model alone.

| Model | Family | Spearman vs RS3-Hsu |
|---|---|---|
| Rule Set 3 (Hsu2013 tracr) | gradient-boosted, transcript-aware | — |
| Rule Set 3 (Chen2013 tracr) | same model, improved scaffold | 0.758 |
| Doench 2014 Rule Set 1 | logistic, independent training set | 0.513 |

CRISPRscan was evaluated and **excluded**: trained on T7-transcribed sgRNAs in
zebrafish embryos, which is not the regime of a lentiviral U6 screen.

**Table:** [table3_ontarget_models.csv]({{artifact:94568c70-10b6-4905-812f-f90508661ecf}})

## Step 5 — Off-target enumeration

Exhaustive, not alignment-heuristic: every protospacer expanded into all ~32,551
sequences within 3 substitutions into a sorted 2-bit table, with every
NGG/NAG/NGA-adjacent 20-mer in hg38 streamed against it. A bitmap prefilter on
the top 30 bits rejects ~97 % of sites and gives a 7× speed-up; verified to
return a hit set identical to the unoptimised implementation on chr22.

| Pass | Budget | Guides | Runtime | Sites found |
|---|---|---|---|---|
| A | ≤2 mm | 16,515 | 12.7 min | 111,591 |
| B | ≤3 mm | 9,420 | 82 min | 514,765 |
| Controls | ≤3 mm | 1,408 | 26 min | 108,611 |
| NTC screen | ≤2 mm | 9,000 | 9 min | — |
| Rebalance | ≤3 mm | 858 | 26 min | 59,292 |

**Table:** [table4_offtarget_passes.csv]({{artifact:e3ea6e58-6baf-4e5c-b6ed-8945bb9f432e}}) · per-guide detail in [offtarget_specificity.csv]({{artifact:5da8efbd-6bb1-44f3-9691-cb53b95b40e5}})

![Off-target profile and library composition]({{artifact:art_584d71e2-7e4f-457b-b337-035a0a049e2d}})

**Finding that changed the selection rule.** 93 % of all recovered sites sit at
exactly 3 mismatches, which swamps the aggregate CFD specificity score — median
0.275 in the screened pool versus 0.278 in the selected set. Aggregate CFD does
not discriminate between guides at this budget. Selection used the targeted
question instead: does this guide cut a coding exon elsewhere, and how hard
(max coding CFD > 0.5 flags 3.6 % of guides; > 0.2 flags 17.8 %).

## Step 6 — Power

Dispersion estimated from real replicate counts (GEO GSE226977, 71,090 sgRNAs,
12 samples), not assumed: **φ = 0.299** at endpoint, rising to 0.90 for abundant
guides, endpoint replicate r = 0.43. T0 "replicates" are Poisson (φ ≈ 0, r = 0.90)
and are technical, not biological.

![Power analysis]({{artifact:art_9ad0e766-763e-46a2-9ee4-3a0ce050bd81}})

| Configuration | Power at log₂FC = 1 |
|---|---|
| 8 guides × 3 replicates | 0.57 |
| 8 × 4 | 0.72 |
| **10 × 4** | **0.80** |
| 12 × 4 | 0.87 |

Coverage is not a lever (0.80 / 0.80 / 0.82 at 200× / 500× / 1000×). Dispersion
is: at φ = 0.50 even 12 × 6 reaches only 0.82.

**Tables:** [power_simulation.csv]({{artifact:439bc761-d23f-48e3-8449-a7837d3b933f}}) · [power_sensitivity.csv]({{artifact:9386dd8f-0bc9-4e71-bb6a-f14ba0722237}}) · [power_mde.csv]({{artifact:c3a1e620-3eb7-49d4-aa78-4f24e1144fd0}})

## Step 7 — Scale and feasibility

![Scale and feasibility]({{artifact:art_04645d04-00a6-469f-9a1f-f8cece725f7a}})

| | Focused | Genome-wide equivalent |
|---|---|---|
| Cells to transduce (4 reps) | 354 M | 2.81 B |
| Cells through sorter / rep | 19.5 M | 155 M |
| FACS time / rep | 22 min | 2.9 h |
| Total reads | 83 M | 658 M |

Sorting bin fraction, not library size, drives cell number
(cells = library × coverage ÷ bin fraction). This is the quantitative reason the
design is focused rather than genome-wide.

**Table:** [table6_scale_and_depth.csv]({{artifact:b1dc65ea-7dc8-4a41-b779-435d5b0e6840}}) · grid in [depth_calculator_grid.csv]({{artifact:3daa3aa0-bbbd-43c4-9eb5-b9a3641a2343}})

## Step 8 — Network propagation layer (L3)

Two-layer mean-aggregator GraphSAGE in NumPy, gradients finite-difference
checked to 2.9 × 10⁻⁸. STRING v12 ≥0.7 over the target set (817 nodes, 5,824
edges); GTEx v8 median TPM across 54 tissues as node features.

![Network propagation]({{artifact:art_ac34b01c-55f9-4977-94f8-cdc1364a1650}})

| Model | Held-out AUROC |
|---|---|
| **GraphSAGE (network + expression)** | **0.805 ± 0.088** |
| Expression only (logistic) | 0.704 ± 0.060 |
| Degree only (control) | 0.387 ± 0.043 |

Paired Wilcoxon over 20 splits, p = 2 × 10⁻⁴. The degree control below chance
rules out hub recovery. **Table:** [gnn_validation.csv]({{artifact:9e748b75-bbc2-45c5-b1ce-c5db456f6ed0}})

## Step 9 — Readout specificity, and the rebalancing it triggered

The pre-registered overlap test could not be run: **no published genome-wide
CRISPR screen of osteogenic or mineralization differentiation exists in GEO.**
Re-posed against GO (EBI `goa_human.gaf`, QuickGO closure).

![Readout specificity]({{artifact:art_7953e154-4000-438f-9727-bba790b8422f}})

1. By annotation the programmes are nearly disjoint — Jaccard 0.054, 18.1 % of
   the tooth programme also annotated to bone. **Passes the 70 % threshold.**
2. But they are **not separable computationally**: 0.582 ± 0.084 AUROC, against
   0.805 for the same model on module recovery.
3. And the curated module **leaned bone 2.45 : 1**.

**Table:** [readout_specificity_test.csv]({{artifact:cb3f30d6-5071-423c-b863-7637a03ffeb4}}) · full analysis in [readout_specificity.md]({{artifact:92dea7a9-b764-42e6-b446-126ac67e02e1}})

### Rebalancing (executed)

All 39 tooth-programme genes absent from the manifest were added — *ACP4, AMELY,
AMTN, ANKRD11, AQP1/3/5/6, BCOR, BSG, CAMP, CFTR, CTNNA1, CSF3R, EDA, FAM20A,
FOXC1, FOXI3, GAS1, HDAC2, ID3, ITGB6, KLK5, LAMB1, ODAM, OSR1, RELT, ROGDI,
SCN5A, SCN10A, SERPINE1, SLC4A2, SP6, SSUH2, TMT1A, TSPEAR, TUFT1, WNT6, ZNF22* —
and designed through the identical pipeline: 6,013 candidates, 858 screened at
≤3 mm, 679 passing, **376 guides selected across all 39 genes** (9.6/gene).

![Rebalancing and final pool]({{artifact:art_6c9a0035-78d4-4d4b-8035-d4624d65706a}})

| | Before | After |
|---|---|---|
| Tooth-only genes | 22 | **58** |
| Shared | 10 | 13 |
| Bone-only | 54 | 54 |
| **Bone : tooth ratio** | **2.45 : 1** | **0.93 : 1** |
| GO tooth programme covered | 44 % | **99 %** |

Genes were **added, not deleted**. Removing bone-annotated genes would have
discarded real candidates to improve a ratio; adding the missing tooth genes
corrects the bias without losing discovery power.

## Step 10 — Final library

![Library composition]({{artifact:art_51e4ae88-2393-443f-99ab-229e2764f9b3}})

| Class | Oligos | Genes | Median CFD spec (≤3 mm) |
|---|---|---|---|
| Targeting | 8,127 | 861 | 0.277 |
| Non-targeting control | 1,000 | — | 1.000 |
| Pan-essential control | 254 | 52 | 0.302 |
| Olfactory control | 250 | 97 | 0.227 |
| **Total** | **9,631** | **861** | |

Controls are 15.6 % of the pool. Non-targeting controls were **verified, not
assumed**: 9,000 candidates screened, the 56.2 % with no site at ≤2 mismatches
retained. A zero-site criterion at ≤3 mismatches is unattainable — a random
20-mer averages ~30 such sites.

**Read-back verification of the saved file:** 9,631 rows; all oligos 80 nt; all
protospacers 20 nt; **0** internal BsmBI/BsaI sites; **0** duplicate protospacers;
**0** duplicate guide IDs; **0** non-targeting controls with a genomic match.

**Table:** [table7_library_composition.csv]({{artifact:f05d26f0-7c5f-41c1-ae18-6d5e5d0e3c65}})

---

## Remaining decision

One item, and it is not computational: **confirm the pooled-synthesis flanks
against your vector.** The library uses Broad GPP-style flanks
(`GGAAAGGACGAAACACCG` … `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG`). This is
the only step that cannot be verified from sequence alone and the only one that
silently ruins a synthesis order.

The parallel osteoblast arm (§9) remains a budget decision, not a library one —
the same pool is used in both arms.

## Known limitations

Off-target search covers the primary assembly only, to 3 mismatches, with
NGG/NAG/NGA PAMs. The guide-efficacy distribution in the power model is an
assumption, not a measurement. GTEx contains no dental tissue, so L3 node
features are a proxy. φ derives from one published screen in an unrelated cell
type and should be re-measured in a 1,000-guide pilot before the library is
committed.
