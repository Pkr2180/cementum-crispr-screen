# A CRISPR knockout screen for acellular cementum regeneration

**Design document — built on primary data, not assumptions.**
Target space, guide library, off-target profile, scale, power, and the
predictive layer. Every number below was computed in this project from hg38,
RefSeq, GTEx, STRING, GO, and a published screen's raw counts; sources are named
at each step.

![Screen architecture]({{artifact:art_c9f60c1b-c8ac-4268-918e-b5690abec195}})

---

## 1. Summary, and what changed during the design

Three things changed materially from the architecture I first proposed. Each was
forced by data computed here, and each is worth more than the rest of the
document.

| | Proposed | After computing it | Why |
|---|---|---|---|
| Representation | 500–1000× coverage | **200× is sufficient** | At the dispersion measured in a real screen (φ = 0.30), the 1/cells term is negligible. Coverage above 200× buys nothing. |
| Guides × replicates | 8 × 3 | **10 × 4** | 8 × 3 gives 57 % power at gene-level log₂FC = 1. 10 × 4 gives 80 %. |
| Readout | ALPL/IBSP ratio in one line | **ratio + a parallel osteoblast arm** | Tooth and bone programs are not computationally separable (AUROC 0.582). A single-line marker ratio cannot carry the central claim. |

The library is built and synthesis-ready: **8,127 guides across 861 genes** in a
9,631-oligo pool, every guide screened genome-wide for off-targets to 3
mismatches, with the cementum module rebalanced to 0.93:1 bone:tooth.

## 2. The constraint that shapes the design

Acellular extrinsic fiber cementum is deposited once, during root formation,
by dental-follicle–derived cells after HERS fragmentation. **No cell line makes
it in a dish.** A screen scored on mineralization alone recovers the osteogenic
program and publishes a bone paper under a cementum title.

The design therefore anchors on the pyrophosphate axis, the one axis with
cementum-specific genetic proof: *ALPL* loss (hypophosphatasia) ablates acellular
cementum and causes premature exfoliation of primary teeth; *ANKH* and *ENPP1*
loss produce hypercementosis. The extracellular PPi/Pi ratio — not generic
matrix deposition — is the discriminating variable.

Section 9 tests whether that anchor is sufficient. It is not, on its own.

## 3. Target space

**861 genes**, resolved to Ensembl IDs via BioMart with 1:1 mouse
orthologs where they exist ([target_manifest.csv]({{artifact:00111918-6082-4f4f-b01e-2ab86d57a288}})).

- **Tier B — cementum module, 341 genes** across 11 submodules: PPi/Pi axis,
  root/HERS, PDL fiber identity, Wnt, BMP/TGF-β, growth-factor RTK, remodeling,
  SIBLING/mineralization, mechanotransduction, epigenetic regulators,
  inflammatory.
- **Tier C — kinome, 481 genes** (GO:0004672 protein kinase activity,
  chromosome-partitioned BioMart query). 522 of the 861 genes are kinases once
  overlap with Tier B is counted. Each kinase hit has a chemical probe one step
  behind it.

802 of 829 originally queried genes have 1:1 mouse orthologs. *CEMP1* has none —
relevant if an *in vivo* mouse validation arm is planned.

## 4. Guide design

Enumeration ran against the locally indexed hg38 (3,088,269,832 bp; chr1 length
verified at 248,956,422) using RefSeq Curated transcript structures (78,435
transcripts; 828/829 genes matched).

**Hard filters, applied in order**

1. Cut site inside the canonical (RefSeq Select) CDS, fractional position
   5–65 %, in an exon present in ≥90 % of coding transcripts.
2. No TTTT; GC 30–80 %; no BsmBI/BsaI site; no homopolymer ≥5; no N in the 30-mer.
3. Exactly one perfect genomic match (rejected 564 multi-mapping guides).
4. No 1-mismatch off-target anywhere (rejected 2,274).
5. No coding-exon off-target with CFD > 0.2 at ≤2 mm; CFD > 0.5 at ≤3 mm.
6. **No BsmBI/BsaI site in the assembled oligo, junctions included** — this
   catches 92 guides that a protospacer-only check passes.

Filter 1 is the one most often skipped and the largest driver of real knockout
efficiency. 94 genes could not reach the 90 % constitutive-exon bar (they carry
up to 131 RefSeq transcripts); a per-gene relaxation cascade handled them, and
the tier used is recorded per guide. 734 of 826 genes passed at the strictest
setting.

**Activity scoring — three scores, rank-aggregated**

| Model | Family | Spearman vs RS3-Hsu |
|---|---|---|
| Rule Set 3, Hsu2013 tracr | gradient-boosted, transcript-aware | — |
| Rule Set 3, Chen2013 tracr | same model, improved scaffold | 0.758 |
| Doench 2014 Rule Set 1 | logistic, independent training set | 0.513 |

Selection uses the mean percentile rank, never one model alone. **CRISPRscan was
deliberately excluded** — it is trained on T7-transcribed sgRNAs in zebrafish
embryos, a different regime from a lentiviral U6 screen, and would import
mismatched bias rather than orthogonality.

*Implementation note:* the published Rule Set 3 regressor is pickled under
lightgbm 3.x and raises `TypeError` under lightgbm ≥4. `rs3_helper.py` restores
the regressor default.

## 5. Off-target analysis

**Method.** Exhaustive, not heuristic. Every protospacer is expanded into all
~32,500 sequences within 3 substitutions into a sorted 2-bit table; every
NGG/NAG/NGA-adjacent 20-mer in hg38 is streamed against it with a bitmap
prefilter. Exact for the stated mismatch budget. Verified to return an identical
hit set to the unoptimised implementation on chr22.

Two passes: ≤2 mm on 16,515 guides (12.7 min), ≤3 mm on 9,420 guides (82 min,
514,765 hits). Scored with the official Doench 2016 CFD tables and Hsu 2013 MIT
weights; coding overlap from RefSeq Curated CDS.

![Off-target profile and library composition]({{artifact:art_584d71e2-7e4f-457b-b337-035a0a049e2d}})

**The finding that matters for anyone doing this:** 93 % of all off-target sites
sit at exactly 3 mismatches, so the aggregate CFD specificity score is swamped —
median 0.275 across the screened pool and 0.278 in the selected set. *Aggregate
CFD at ≤3 mm does not discriminate between guides.* What discriminates is the
targeted filter: does this guide cut a coding exon somewhere else, and how hard.
At CFD > 0.5 that flags 3.6 % of guides; at CFD > 0.2, 17.8 %.

## 6. Library composition

[library_complete.csv]({{artifact:f2bc054b-ecf0-44e0-8cb3-eb6082dcbecc}}) — **8,127 targeting guides over 861 genes** (9,631 oligos including controls)

| | |
|---|---|
| Mean guides per gene | 9.44 |
| Genes at the full 10 | 576 |
| Genes spanning ≥2 distinct exons | 762 of 822 original (92.7 %) |
| Genes with <6 guides | 17 |
| Median on-target percentile | 0.900 |
| Oligo length | 80 nt, uniform |
| Internal BsmBI/BsaI sites | 0 (junctions checked) |
| Duplicate protospacers | 0 |

The 17 under-represented genes are paralog pairs where unique targeting is
genuinely hard — *CSNK2A1*/*CSNK2A3*, *TLK1*/*TLK2*, *AMELX* (paralogous with
*AMELY*), *ANXA2* (processed pseudogenes). Treat single-guide results for these
as uncorroborated.

Each guide carries both a pooled-synthesis oligo (Broad GPP-style flanks —
**confirm against your vector before ordering**) and an annealed-oligo pair
(`CACCG…` / `AAAC…C`) for individual lentiCRISPRv2 cloning.

**Controls are built and verified.** [library_complete.csv]({{artifact:c4b6e3fb-55c2-4601-afe5-a19d7dd5d343}}) — the full
**9,631-oligo pool**, 15.6 % controls:

| Class | Oligos | |
|---|---|---|
| Targeting | 8,127 | 861 genes |
| Non-targeting | 1,000 | verified zero genomic sites at ≤2 mm |
| Pan-essential (positive) | 254 | 52 genes |
| Olfactory receptor (neutral) | 250 | 97 genes |

![Library composition]({{artifact:art_51e4ae88-2393-443f-99ab-229e2764f9b3}})

Non-targeting controls were **verified, not assumed**: 9,000 candidates were
screened genome-wide and the 56.2 % with no site at ≤2 mismatches retained. A
zero-site criterion at ≤3 mismatches is unattainable — a random 20-mer is
expected to have ~30 such sites — so ≤2 mismatches is the meaningful bar.

Median CFD specificity is comparable across classes (targeting 0.278,
pan-essential 0.302, olfactory 0.227), so the controls calibrate the null
without being systematically cleaner than the guides they benchmark.

## 7. Scale, depth, and power

Noise model estimated from real replicate data (GEO **GSE226977**, genome-wide
KO screen, 71,090 sgRNAs, 12 samples), not assumed:

- Endpoint biological dispersion **φ = 0.299** (median), rising to 0.90 for
  abundant guides. Endpoint replicate correlation r = 0.43.
- T0 "replicates" are effectively Poisson (φ ≈ 0, r = 0.90) — they are technical,
  not biological, and must not be used to estimate screen noise.

Total count dispersion decomposes as φ_total = φ_bio + 1/cells + 1/reads.

![Power analysis]({{artifact:art_9ad0e766-763e-46a2-9ee4-3a0ce050bd81}})

| Configuration | Power at log₂FC = 1 |
|---|---|
| 8 guides × 3 replicates | 0.57 |
| 8 guides × 4 replicates | 0.72 |
| **10 guides × 4 replicates** | **0.80** |
| 12 guides × 4 replicates | 0.87 |

Minimum detectable effect at the recommended configuration: **log₂FC ≥ 1.0**
(degrades to 1.3 at 8 × 3). Coverage is not a lever: 0.80 / 0.80 / 0.82 at
200× / 500× / 1000×. Dispersion is: at φ = 0.50, even 12 guides × 6 replicates
reaches only 0.82. **Consistent sort gating, matched passage, and single-batch
infection buy more power than any sequencing decision.**

![Scale and feasibility]({{artifact:art_04645d04-00a6-469f-9a1f-f8cece725f7a}})

| | Focused (this library) | Genome-wide equivalent |
|---|---|---|
| Cells to transduce, 4 reps | 354 M | 2.81 B |
| Cells through sorter / rep | 19.5 M | 155 M |
| FACS time / rep (15k events/s) | 22 min | 2.9 h |
| Samples | 17 | 17 |
| Total reads (500/guide) | 83 M | 658 M |

Sorting bin fraction, not library size, drives cell number. A genome-wide
decile-sorted screen needs **2.8 billion cells** — past routine scale. This is
the quantitative reason the design is focused rather than genome-wide.

## 8. Predictive layer

Full specification in [ai_layer_spec.md]({{artifact:ece3552d-6d9e-487f-9f40-15cbc50f5eb5}}). **L1 and L3 are built and validated; L2 and L4 are specified.**

**L3 — network propagation (built).** 2-layer mean-aggregator GraphSAGE in
NumPy, gradients finite-difference checked to 2.9 × 10⁻⁸. Graph: STRING v12,
confidence ≥0.7, 817 nodes, 5,824 edges. Features: GTEx v8 median TPM across 54
tissues plus log degree.

![Network propagation]({{artifact:art_ac34b01c-55f9-4977-94f8-cdc1364a1650}})

Held-out recovery of the PPi/SIBLING mineralization module, 20 random splits:
**AUROC 0.805 ± 0.088** versus 0.704 ± 0.060 for expression alone (paired
Wilcoxon p = 2 × 10⁻⁴). A degree-only control scores 0.387 — below chance — so
this is not hub recovery.

**L4 is gated on an experimental decision.** GEARS and scGPT perturbation
fine-tuning need perturbation → transcriptome pairs. A FACS-sorted pooled screen
produces one scalar per cell and cannot supply them. Either add a direct-capture
Perturb-seq arm on the top ~200 hits (~100k cells), or restrict L4 to zero-shot
Geneformer prioritisation and drop the "predicts unseen perturbations" claim.

## 9. Readout specificity — the falsification test

Full analysis in [readout_specificity.md]({{artifact:92dea7a9-b764-42e6-b446-126ac67e02e1}}).

The pre-registered test was: *if hits would overlap published osteoblast
mineralization screens by >70 %, claim 1 collapses.* **That test could not be
run — there is no published genome-wide CRISPR screen of osteogenic or
mineralization differentiation in GEO.** It was re-posed against GO annotations
(EBI `goa_human.gaf`, QuickGO closure).

![Readout specificity]({{artifact:art_7953e154-4000-438f-9727-bba790b8422f}})

1. By annotation the programs are nearly disjoint — Jaccard 0.054; only 18.1 % of
   the tooth program is also annotated to bone. **Passes the threshold.**
2. But they are **not separable computationally**: the same GraphSAGE that
   reaches 0.805 on module recovery reaches only **0.582 ± 0.084** discriminating
   tooth-only from bone-only genes.
3. And the curated module **leans bone 2.5 : 1** (54 bone-annotated vs 22
   tooth-annotated), covering 77 % of shared genes but only 39 % of tooth-only
   genes — a curation bias pointing the wrong way.

**Verdict: conditional go.** The threshold passes for the wrong reason. Required
changes, in order:

1. **Add a parallel osteoblast-line arm (hFOB1.19 or primary calvarial) and call
   hits differentially.** This substitutes an experimental contrast for the
   computational separation that does not exist. Roughly doubles screen cost.
2. **Move the primary readout from matrix markers to the PPi/Pi axis** — direct
   extracellular-PPi or ALPL enzymatic activity rather than the transcriptional
   reporter ratio.
3. **Rebalance the module before synthesis** — add the 36 tooth-only GO genes not
   currently in the manifest; down-weight generic ossification entries.

Expect ~1 in 5 hits to be generic mineralization regardless. State that prior up
front.

## 10. Bench timeline

| Weeks | |
|---|---|
| 0–6 | Oligo synthesis; reporter knock-ins and clonal iCas9 line; verify single-copy integration |
| 6–9 | Library cloning, ≥50× bacterial coverage; verify representation by NGS of the plasmid pool |
| 9–11 | Virus production and titration; MOI 0.3 in both cell backgrounds |
| 11–13 | Transduction (4 reps × 2 arms), puromycin selection, expansion |
| 13–16 | Differentiation induction, 14–21 d, matched passage, single batch |
| 16–17 | FACS (22 min/rep), gDNA, PCR, sequencing |
| 17–20 | L2 hit calling, L3 propagation, differential calling across arms |
| 20–32 | Validation: arrayed knockouts, chemical probes for kinase hits, *in vivo* defect model |

**Pilot first.** Run a 1,000-guide subset through the full workflow to measure φ
in *your* hands before committing the library. Every power number above is
conditional on φ ≈ 0.30; at φ = 0.5 the design needs rethinking, and that is
cheaper to learn in week 10 than week 20.

## 11. Indicative cost

Order-of-magnitude, US list, excluding labor. **Estimates, not quotes.**

| | Round 1, both arms |
|---|---|
| Oligo pool (~9,250 oligos) | $3–6k |
| Cloning, virus, titration | $3–4k |
| Cell culture (708 M cells) | $8–12k |
| FACS (~3 h instrument) | $1–2k |
| gDNA, PCR, library prep | $2–3k |
| Sequencing (166 M reads) | $3–5k |
| **Total** | **$20–32k** |
| Optional Perturb-seq arm (~100k cells) | +$15–25k |

## 12. Reproducibility

Computed on hg38 with: RefSeq Curated/Select (UCSC), GTEx v8 medians, STRING
v12 (≥0.7), GO via EBI `goa_human.gaf` + QuickGO closure, GSE226977 counts,
`rs3` 0.0.18, Doench 2016 CFD tables, Hsu 2013 MIT weights.

Code: `guide_design.py`, `offtarget.py`, `run_offtarget.py`, `ontarget.py`,
`rs3_helper.py`, `power_sim.py`, `gnn.py`.

**Needed before ordering, in order:**

1. **Confirm the pooled-synthesis flanks against your vector.** The library uses
   Broad GPP-style flanks (`GGAAAGGACGAAACACCG` … `GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG`).
   This is the only step I cannot verify computationally and the only one that
   silently ruins a pool if wrong.
2. **Decide on the parallel osteoblast arm** (§9). It changes the budget, not
   the library.
3. ~~Module rebalancing~~ — **done**. 39 tooth-programme genes added (+376 guides); bone:tooth now 0.93:1, GO tooth coverage 99%.

**Known limitations.** Off-target search covers primary chromosomes only (no
alts/randoms) at ≤3 mismatches with NGG/NAG/NGA PAMs — it will not find a
4-mismatch off-target or one in an alt contig. Guide efficacy in the power model
is a Beta(2, 1.33) assumption, not a measurement. GTEx has no dental tissue, so
L3 node features are a proxy. φ is from one published screen in a different cell
type.
