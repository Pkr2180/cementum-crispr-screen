# Predictive layer specification (L1–L4)

Scope note: **L1 and L3 are built and validated in this project; L2 and L4 are
specified only.** L2 cannot run until real counts exist. L4 needs a GPU, which
this machine does not have, and — see §L4.0 — it also needs an experimental arm
that the sorting screen as designed does not produce.

---

## L1 — Guide activity (BUILT)

| | |
|---|---|
| Models | Rule Set 3 (`rs3` 0.0.18, Hsu2013 and Chen2013 tracrRNA scaffolds) + Doench 2014 Rule Set 1 |
| Input | 30-mer context: 4 nt 5′ + 20 nt protospacer + 3 nt PAM + 3 nt 3′ |
| Output | rank-aggregated percentile across the three scores, plus a disagreement statistic (SD of percentile ranks) |
| Scale run | 139,144 candidate protospacers over 826 genes |
| Observed concordance | Spearman ρ = 0.758 (RS3 Hsu vs RS3 Chen), 0.513 (RS3 Hsu vs Doench 2014), 0.394 (RS3 Chen vs Doench 2014) |

Two implementation notes that will bite anyone reproducing this:

1. The pickled Rule Set 3 regressor was serialised under lightgbm 3.x and raises
   `TypeError` under lightgbm ≥ 4 because `_n_classes` is absent. `rs3_helper.py`
   restores the regressor default before predicting.
2. CRISPRscan (Moreno-Mateos) was deliberately **excluded**. It is trained on
   T7-transcribed sgRNAs injected into zebrafish embryos; a lentiviral U6-driven
   screen in human cells is a different regime and the score would add
   mismatched bias rather than orthogonality.

Selection uses the rank aggregate rather than any single model, so no one
model's training bias determines the library.

## L2 — Hit calling (SPECIFIED)

Ensemble of three estimators over the same count matrix, because each fails
differently:

- **MAGeCK-MLE** — negative-binomial GLM, gene-level beta scores, handles the
  multi-bin sort design directly via a design matrix.
- **BAGEL2** — Bayes factor against reference essential/non-essential sets;
  robust when the effect is a shift in a subpopulation.
- **Chronos** — copy-number–aware; matters because cementoblast and PDL lines
  are frequently aneuploid and CN artefacts otherwise masquerade as depletion.

Consensus rule: a gene is a hit if it passes FDR 0.05 in ≥ 2 of 3 estimators and
the direction agrees. Non-targeting controls (n = 1,000) set the empirical null;
olfactory-receptor controls (n = 250) verify the null is correctly centred.

**Power constraint from the measured noise model** (see `power_analysis.png`):
at the dispersion observed in real screen data (φ ≈ 0.30), the design needs
10 guides/gene × 4 replicates for 80 % power at gene-level log₂FC = 1. Below
that, L2 output is not worth modelling downstream.

## L3 — Network propagation (BUILT)

| | |
|---|---|
| Architecture | 2-layer mean-aggregator GraphSAGE, 32 hidden units, Adam, masked BCE |
| Graph | STRING v12 human, confidence ≥ 0.7, induced over the target set — 817 nodes, 5,824 edges, mean degree 14.3, 74 isolated |
| Node features | GTEx v8 median TPM across 54 tissues (log1p, z-scored per tissue) + log degree |
| Implementation | NumPy with hand-derived gradients, finite-difference checked to 2.9 × 10⁻⁸ |

**Validation.** Trained to recover held-out members of the PPi/SIBLING
mineralization core (46 positives, 60/40 split, 20 random splits):

| Model | held-out AUROC | held-out AP |
|---|---|---|
| GraphSAGE (network + expression) | **0.805 ± 0.088** | 0.390 ± 0.117 |
| Expression only (logistic) | 0.704 ± 0.060 | 0.235 ± 0.054 |
| Degree only (control) | 0.387 ± 0.043 | — |

Paired Wilcoxon over the 20 splits, GraphSAGE vs expression-only: p = 2 × 10⁻⁴.
Base rate 0.058. The degree control sitting *below* chance rules out the
"it just finds hubs" explanation.

This was validated on module recovery, **not** on screen z-scores, because no
screen data exists yet. When counts arrive, the same architecture takes
L2 gene-level effect sizes as the regression target; the module-recovery result
is evidence the graph carries propagatable signal, not a substitute for that fit.

## L4 — Perturbation foundation model (SPECIFIED — and gated)

### L4.0 The constraint to settle first

**GEARS and scGPT perturbation-response fine-tuning require
perturbation → transcriptome pairs. A FACS-sorted pooled screen does not produce
them** — it yields one scalar (bin assignment) per cell, not an expression
profile. Three options, in increasing cost:

1. **Zero-shot only.** Use pre-trained Geneformer for *in-silico* deletion
   ranking. No fine-tuning, no new experiment. Weakest, but free, and it can
   pre-rank the round-2 library.
2. **Add a direct-capture Perturb-seq arm** on the top ~200 hits from round 1
   (CROP-seq / direct-capture 10x, ~100k cells, ~500 cells per perturbation).
   This is what makes a trainable L4 real. **Recommended.**
3. **Genome-scale Perturb-seq.** Not justified at this stage.

Without option 2 or 3, L4 claims in a paper should be confined to zero-shot
prioritisation, and the "predicts unseen perturbations" claim cannot be made.

### Models and compute

| Model | Params | Task | VRAM (fine-tune) | Runtime estimate |
|---|---|---|---|---|
| Geneformer | ~110 M | in-silico deletion, zero-shot or fine-tuned | 24–40 GB full; ~16 GB with LoRA | 1–2 h inference over 100k cells (A100) |
| scGPT | ~50 M | perturbation response | 24 GB | 4–8 h fine-tune on 100k cells |
| GEARS | ~2 M | unseen single + combinatorial perturbations | 8 GB | 1–2 h on a 200-perturbation set |

GEARS is the one that answers the combinatorial question; it is also by far the
cheapest, because it is a GNN over a GO/co-expression graph rather than a
transformer over cells.

### Reference data for fine-tuning

Human periodontal ligament / dental pulp scRNA-seq from GEO, used to adapt the
cell-state prior before perturbation fine-tuning. Accessions to be fixed at run
time by query rather than assumed here.

### Active-learning loop (round 1 → round 2)

Let `μ(g)` and `σ(g)` be the ensemble mean and SD of predicted effect for
perturbation `g` (ensemble over GEARS seeds and L3 folds).

```
acquisition(g) = μ(g) + κ·σ(g)          κ ≈ 1.0   (upper confidence bound)
```

selected under two constraints:

- **Diversity.** Greedy submodular selection penalising STRING-graph proximity
  to already-selected perturbations, so round 2 does not re-screen one complex.
- **Actionability quota.** ≥ 40 % of round-2 slots reserved for perturbations
  with an existing chemical probe, so hits convert to experiments directly.

Round 2 = ~300 perturbations including ~100 pairs, run as a mini-screen at the
same 10 guides × 4 replicates.

**The claim this loop supports, stated precisely:** round-2 hit rate exceeds the
round-1 empirical hit rate at matched FDR. That is a measurable, falsifiable
comparison. It is *not* a claim that the model replaces a screen.

---

## What would make this fail

- φ in the real screen lands near 0.5 rather than 0.30 → L2 underpowered, and
  everything downstream inherits the noise. Mitigation: pilot 1,000-guide run to
  measure φ before committing the full library.
- The sorted readout recovers generic osteogenesis → L3/L4 would faithfully model
  the wrong program. This is tested directly in `readout_specificity.md`.
- No Perturb-seq arm → L4 is zero-shot only, and the headline predictive claim
  must be dropped.
