# Readout-specificity falsification test

**Pre-registered go/no-go** (stated before any of this was computed): *if hits
from the ratiometric sort would overlap published osteoblast mineralization
screens by more than 70 %, the readout is not cementum-specific and claim 1 —
that cementogenic and osteogenic differentiation are genetically separable in a
shared progenitor — collapses.*

## What the test had to become, and why

The test as planned could not be run. **There is no published genome-wide CRISPR
screen of osteogenic or mineralization differentiation in GEO.** A search over
CRISPR-screen series returned 42 bone/osteo-matching and 66
differentiation-matching records; on inspection every one is a cancer,
leukaemia, or immune screen. There is no hit list to overlap against.

So the question was re-posed in a form that real data can answer: **are the
tooth and bone differentiation programs distinguishable at all?** If they are
not — in annotation, in network context, in expression — then no marker-ratio
readout built from shared mineralization machinery can separate them, and the
risk the original test was designed to catch is present regardless of what any
future screen overlaps with.

Sources: GO annotations from the EBI `goa_human.gaf` release (17,779 annotated
terms), with term closure from the EBI QuickGO descendants API.
Tooth program = odontogenesis ∪ tooth mineralization (72 genes).
Bone program = ossification ∪ osteoblast differentiation ∪ bone mineralization
(183 genes).

## Result 1 — by annotation, the programs are nearly disjoint

| | genes |
|---|---|
| tooth-only | 59 |
| shared | 13 |
| bone-only | 170 |
| Jaccard(tooth, bone) | **0.054** |

Only **18.1 %** of the tooth program is also annotated to bone. Taken alone this
is reassuring and would pass the pre-registered 70 % threshold comfortably.

## Result 2 — but they are not separable computationally

Training the validated L3 GraphSAGE (same architecture, same STRING graph, same
GTEx features) to discriminate tooth-only from bone-only genes, 25 random 60/40
splits:

| Model | AUROC |
|---|---|
| GraphSAGE (network + expression) | **0.582 ± 0.084** |
| Expression only | 0.544 ± 0.107 |

Above chance (one-sample Wilcoxon p = 3.3 × 10⁻⁴) but far below the ~0.8 that
would make it useful. For context, the *same model on the same graph* reaches
0.805 when asked to recover the PPi/SIBLING mineralization module — so this is
not a weak model, it is a genuinely weak signal. Tooth-specific and
bone-specific genes sit in nearly the same interaction neighbourhoods and have
nearly the same tissue-expression profiles.

## Result 3 — the curated module leans bone

Of the 342 genes in the Tier B cementum module:

| | genes | |
|---|---|---|
| annotated to the tooth program | 22 | 6 % |
| annotated to both | 10 | 3 % |
| annotated to the bone program | 54 | 16 % |
| in neither GO program | 256 | 75 % |

Bone-annotated genes outnumber tooth-annotated **2.5 : 1**. And the library
covers 77 % of the *shared* genes but only 39 % of the tooth-only genes — the
curation was drawn preferentially toward the generic mineralization machinery,
which is exactly the direction that would manufacture a bone result.

## Verdict: conditional go

The pre-registered numeric threshold passes (18 % ≪ 70 %), but it passes for the
wrong reason — it measures annotation membership, not whether the assay can tell
the programs apart. On the question that actually matters the answer is negative:
**an ALPL/IBSP marker-ratio sort in a single cell line cannot be relied on to
separate cementogenic from osteogenic differentiation.** Three changes follow.

1. **Add a parallel osteoblast-line arm and call hits differentially.** Run the
   identical library and sort in an osteoblast line (hFOB1.19 or primary
   calvarial osteoblasts) and require hits to be PDL/cementoblast-selective.
   This substitutes an experimental contrast for the computational separation
   that is unavailable. This is the single most important change to the design
   and roughly doubles the screen's cost.
2. **Move the primary readout from matrix markers to the PPi/Pi axis.** ALPL and
   IBSP are shared with bone. The pyrophosphate set-point is the one axis with
   cementum-specific human and mouse genetic proof (hypophosphatasia; *Ank* and
   *Enpp1* loss). A direct extracellular-PPi or ALPL-enzymatic-activity readout
   is more specific than the transcriptional reporter ratio.
3. **Rebalance the curated module before synthesis.** The 2.5 : 1 bone lean is a
   curation artefact and is correctable: add the 36 tooth-only GO genes not
   currently in the manifest, and down-weight generic ossification genes that
   entered only through the osteogenic literature.

**Expect roughly 1 in 5 hits to be generic mineralization** regardless, from the
18 % annotation overlap. That is the prior to state in the paper, not a number
to discover in review.
