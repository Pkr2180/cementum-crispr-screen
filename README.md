# A genome-verified CRISPR knockout library for acellular cementum regeneration

Design, genome-wide off-target verification and *in silico* characterisation of a pooled
CRISPR knockout library targeting the cementogenic programme and the human kinome.

**The screen has not been run.** There are no hits, no enrichments and no validated genes in
this repository. Every quantity reported here is a computational property of the design.

---

## The library

| | |
|---|---|
| Oligonucleotides | **9,631** |
| Targeting guides | 8,127 over **861 genes** (mean 9.44 per gene; 576 genes at the full 10) |
| Target space | 380-gene curated cementum module (12 submodules) + 481-gene kinome |
| Controls | 1,000 verified non-targeting · 254 pan-essential (52 genes) · 250 olfactory (97 genes) — 15.6% of the pool |
| Reference genome | GRCh38/hg38, ncbiRefSeqCurated transcript models |
| Off-target verification | exhaustive enumeration to 3 mismatches; **794,259** catalogued sites |
| Oligo length | 80 nt (Hsu 2013 scaffold) / 90 nt (F+E Chen 2013 scaffold) |

Guides were enumerated from 139,144 candidate protospacers and reduced to 7,751 in the main
design pass (94.4% attrition), then extended to 8,127 after a Gene Ontology rebalancing step
that added 39 tooth-programme genes.

## Three design findings

Each contradicts a widely applied but previously untested convention.

**1. Aggregate off-target specificity is uninformative at a three-mismatch budget.**
Within the principal ≤3-mismatch pass, 91.3% of the 514,765 recovered sites lie at exactly
three mismatches, so the summed score is swamped by trivially active sites: median aggregate
CFD specificity was 0.265 across the screened pool and 0.277 across the selected library.
Selection should use a **coding-exon-restricted maximum** instead — a 0.5 threshold on that
statistic flags 3.6% of guides, a 0.2 threshold flags 17.8%.

**2. Library representation is the wrong budget knob.**
With overdispersion measured from published replicate counts (φ = 0.299) rather than assumed,
power at a gene-level log2 fold change of 1.0 is 0.52 (8 guides × 3 replicates), 0.68
(10 × 3), 0.72 (8 × 4) and 0.81 (10 × 4), while raising representation from 200× to 1000×
moves power by about two points. Spend on replicates, not cells.

**3. The cementogenic and osteogenic programmes are not computationally separable.**
The two programmes are nearly disjoint by Gene Ontology annotation (Jaccard 0.054), yet a
GraphSAGE model that recovers cementum module membership at AUROC 0.805 ± 0.088 discriminates
tooth-only from bone-only genes at only **0.582 ± 0.084**. A mineralisation reporter in a
single cell background therefore cannot support a cementum-specific claim; a parallel
osteoblast arm with differential hit calling is required (marginal cost +67%, not a doubling).

## Scale

The focused design needs 350 million cells across four replicates, a 22-minute sort per
replicate, 17 samples and 82 million reads. A genome-wide screen at the same decile gating
would need **2.81 billion cells** and a 2.9-hour sort per replicate — the quantitative basis
for keeping the library focused.

---

## Repository layout

```
manuscript/   main manuscript and supplementary information (Word)
figures/      Figures 1-8 and Supplementary Figures S1-S2 (300 dpi PNG)
tables/       Tables 1-4, Supplementary Tables S1-S5, and the underlying grids
library/      the synthesis-ready pool, target manifest and per-guide off-target records
code/         guide design, off-target enumeration, scoring, power, network model, manuscript build
docs/         methods, results package, design rationale, verified bibliography
```

### Key files

| File | Contents |
|---|---|
| `library/oligo_pool_order.csv` | 9,631 oligos, 80 nt, original Hsu 2013 scaffold |
| `library/oligo_pool_order_FE.csv` | 9,631 oligos, 90 nt, F+E Chen 2013 scaffold |
| `library/oligo_pool_both_scaffolds.csv` | both variants side by side |
| `library/library_complete.csv` | full annotation per guide, including CFD specificity |
| `library/offtarget_specificity.csv` | per-guide off-target summary at ≤3 mismatches |
| `code/verify_flanks.py` | checks which oligo variant matches your vector |

## Before ordering

The 3′ homology flank of each oligo **is** the sgRNA scaffold, so it must match the constant
region already present in your vector. Two variants are provided. Verify which one applies:

```
python code/verify_flanks.py --genbank path/to/your_vector.gb
```

The script searches both strands with circular wrap, reports which oligo file to order, and
exits non-zero on mismatch.

## Reproducing the analysis

```
python code/guide_design.py      # enumerate and filter candidate protospacers
python code/run_offtarget.py     # genome-wide off-target enumeration (CFD / MIT)
python code/ontarget.py          # Doench 2014 Rule Set 1 scoring
python code/power_sim.py         # negative-binomial power grid
python code/gnn.py               # GraphSAGE over STRING + GTEx
python code/assemble.py          # rebuild both Word documents
```

Requires Python 3.11 with NumPy, pandas, SciPy, scikit-learn, LightGBM, Biopython, networkx
and python-docx. Off-target enumeration needs a local copy of GRCh38 and the UCSC
ncbiRefSeqCurated / ncbiRefSeqSelect annotation tables; `run_offtarget.py` uses a bitmap
prefilter on the top 30 bits of the encoded protospacer (~97% rejection, ~7× speed-up) and
was verified against the unoptimised search on chromosome 22.

## Note on simulation

The power analysis is a Monte Carlo calculation, as any prospective power calculation must be.
Its single free parameter, the overdispersion φ, was measured from published replicate screen
counts rather than assumed. No biological result in this repository derives from simulated
data.

## References

All 27 references in the manuscript were verified against PubMed records by DOI; the
machine-readable list is `docs/bib_verified.json`.
