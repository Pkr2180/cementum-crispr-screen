# Methods

## Target gene selection

Two tiers were assembled. Tier C (kinome) was retrieved from Ensembl BioMart
(`hsapiens_gene_ensembl`, Ensembl Genes mart) by filtering on GO:0004672
(protein kinase activity) with `biotype = protein_coding`, partitioned by
chromosome because the unpartitioned query exceeds the service timeout. Tier B
(cementum module) was curated from the literature across eleven submodules:
pyrophosphate/phosphate axis, root and Hertwig's epithelial root sheath,
periodontal-ligament fibre identity, WNT, BMP/TGF-β, growth-factor receptor
tyrosine kinases, bone remodelling, SIBLING and mineralization, mechano-
transduction and hypoxia, epigenetic regulators, and inflammatory mediators.
Symbols were resolved to Ensembl gene identifiers and to mouse orthologues
(`mmusculus_homolog_associated_gene_name`, `mmusculus_homolog_orthology_type`)
in chunks of 30 identifiers. Genes on non-primary scaffolds were discarded.

The union comprised 829 genes (342 Tier B, 487 Tier C; 522 kinases in total),
of which 802 had one-to-one mouse orthologues. *CEMP1* has no annotated mouse
orthologue.

## Genome and transcript annotation

GRCh38/hg38 was obtained from the UCSC Genome Browser
(`goldenPath/hg38/bigZips/hg38.fa.gz`) and encoded to one byte per base
(A/C/G/T/N) for the 24 primary chromosomes; the assembly totalled
3,088,269,832 bp and chromosome 1 measured 248,956,422 bp, matching the GRCh38
primary assembly. Alternate, random, and unplaced contigs were excluded.

Transcript structures were taken from UCSC `ncbiRefSeqCurated` (78,435
transcripts) with the RefSeq Select transcript per gene identified from
`ncbiRefSeqSelect`. Coding exon intervals were derived per transcript by
intersecting exon bounds with the CDS bounds. 828 of 829 target genes matched a
RefSeq record (*ORAI1* did not). Where a gene lacked a RefSeq Select transcript
(3 genes), the transcript with the longest CDS was used as canonical.

## Guide enumeration and hard filtering

All NGG-PAM protospacers on both strands were enumerated within each gene's
coding span. A candidate was retained only if its predicted Cas9 blunt cut site
(3 bp 5′ of the PAM) fell inside the canonical transcript's CDS, at a fractional
CDS position of 0.05–0.65, within an exon present in at least 90 % of that
gene's coding transcripts. Candidates were further required to contain no TTTT
(pol III terminator), to have GC content 0.30–0.80, to contain no BsmBI
(CGTCTC/GAGACG) or BsaI (GGTCTC/GAGACC) recognition site, to contain no
homopolymer run of 5 or more, and to have no ambiguous base in the 30-mer
context.

Ninety-four genes could not satisfy the 90 % constitutive-exon criterion,
predominantly those with large numbers of annotated transcripts (up to 131). For
these a per-gene relaxation cascade was applied, lowering the constitutive-exon
threshold to 0.75, then 0.50 (with the CDS window widened to 0.05–0.80), then 0
(window 0.02–0.90), stopping at the first tier yielding at least 40 candidates;
the tier used is recorded per guide. 734 of 826 genes passed at the strictest
setting and 10 required the most permissive.

This produced 139,144 candidate protospacers across 826 genes (median 119 per
gene). Candidates rejected at the strictest tier by sequence composition were:
homopolymer 7,107; polyT 6,485; GC 5,598; restriction site 2,244.

## On-target activity scoring

Each 30-mer context (4 nt 5′ + 20 nt protospacer + 3 nt PAM + 3 nt 3′) was
scored with Rule Set 3 (`rs3` v0.0.18, sequence mode) under both the Hsu2013 and
Chen2013 tracrRNA scaffolds, and with Doench 2014 Rule Set 1, the latter
re-implemented in Python 3 from the parameter table distributed with CRISPOR.
Scores were converted to percentile ranks and averaged; the standard deviation of
the three ranks was retained as a model-disagreement statistic. Spearman
correlations between models were 0.758 (Rule Set 3 Hsu vs Chen), 0.513 (Rule Set
3 Hsu vs Doench 2014) and 0.394 (Rule Set 3 Chen vs Doench 2014).

CRISPRscan was evaluated and excluded: it is trained on T7-transcribed sgRNAs
microinjected into zebrafish embryos, which is not the regime of a lentivirally
delivered U6-driven screen in human cells.

The Rule Set 3 model is distributed as a LightGBM regressor pickled under
LightGBM 3.x and raises a `TypeError` under LightGBM ≥ 4 because the `_n_classes`
attribute is absent; the regressor default was restored before prediction.

## Genome-wide off-target enumeration

Off-target sites were enumerated exhaustively rather than by alignment
heuristics. For each query protospacer, every sequence within the stated
mismatch budget was generated (1, 61, 1,771 and 32,551 sequences at 0, ≤1, ≤2 and
≤3 substitutions respectively) and encoded as a 40-bit integer (2 bits per base)
into a sorted table. Every 20-mer in hg38 adjacent to an NGG, NAG or NGA PAM, on
both strands and excluding windows containing N, was then encoded identically and
matched against that table by binary search. A bitmap prefilter on the top 30
bits of each code (2²⁷-byte array) rejected approximately 97 % of genomic sites
before the binary search and gave a 7-fold speed-up; the genome was scanned in
20-Mb windows to bound peak memory. The optimised implementation was verified to
return a hit set identical to the unoptimised implementation on chromosome 22.

Two passes were run. Pass A screened the 16,515 highest-ranked candidates (top 20
per gene) at ≤2 mismatches in 12.7 min, returning 111,591 sites (17,497 at 0,
5,791 at 1, and 88,303 at 2 mismatches). Pass B screened 9,420 guides at ≤3
mismatches in 82 min in batches of 1,200, returning 514,765 sites (10,393 at 0,
1,738 at 1, 32,424 at 2 and 470,210 at 3 mismatches).

Each site was scored with the Cutting Frequency Determination score using the
mismatch and PAM tables of Doench et al. (2016) as distributed with CRISPOR, and
with the MIT specificity score using the position weights of Hsu et al. (2013).
Coding-exon overlap was determined against merged RefSeq Curated CDS intervals.

Guides were rejected if they had other than exactly one perfect genomic match
(564 guides), any off-target site at 1 mismatch (2,274), any coding-exon
off-target with CFD > 0.2 at ≤2 mismatches (1,196), or any coding-exon off-target
with CFD > 0.5 at ≤3 mismatches (339 of the deeply screened set).

Aggregate CFD specificity, defined as 1/(1 + Σ CFD over off-target sites), was
computed but not used for selection: 93 % of all recovered sites lie at exactly 3
mismatches, and the aggregate score was indistinguishable between the screened
pool and the selected set (median 0.275 vs 0.278). Selection therefore used the
per-site coding-exon criterion above.

## Library assembly

Up to 10 guides per gene were selected by descending rank-aggregated on-target
score, subject to all filters above and to a constraint that a gene's guides span
at least two distinct canonical-CDS exons where the gene offers more than one.
Pooled oligonucleotides were assembled as a 5′ flank, the 20-nt protospacer and a
3′ scaffold-proximal flank, giving a uniform 80-nt oligo; annealed-oligo pairs
(`CACCG` + protospacer; `AAAC` + reverse complement + `C`) were generated in
parallel for individual cloning. The restriction-site filter was re-applied to
the assembled oligo rather than the protospacer alone, which removed a further 92
guides whose flank–protospacer junction creates a BsmBI site.

The final targeting library comprises 7,751 guides across 822 genes (mean 9.4 per
gene; 540 genes at 10 guides; 762 genes, 92.7 %, spanning ≥2 exons; 17 genes with
fewer than 6 guides, predominantly recent paralogue pairs such as
*CSNK2A1*/*CSNK2A3* and *TLK1*/*TLK2* where unique targeting is constrained). All
oligos were verified to be 80 nt, free of internal BsmBI and BsaI sites, and
non-duplicated.

## Control guides

Non-targeting controls were generated as random 20-mers with GC content
0.40–0.60, no TTTT, no homopolymer run ≥5, and no restriction site in the
assembled oligo, then screened against hg38 at ≤2 mismatches; only sequences with
no genomic site at ≤2 mismatches were retained. A zero-site criterion at ≤3
mismatches is not attainable: with approximately 1.16 × 10⁹ NGG/NAG/NGA-adjacent
20-mers in the genome and 32,551 sequences within 3 substitutions of any given
20-mer, a random 20-mer is expected to have of order 30 such sites.

Of 9,000 generated candidates, 5,062 (56.2 %) had no genomic site at ≤2
mismatches; 1,000 were retained.

Positive controls target 56 pan-essential genes (ribosomal proteins, RNA
polymerase II subunits, proteasome subunits, translation initiation factors,
spliceosome components and core chaperones); neutral controls target olfactory
receptor genes matching `^OR\d+[A-Z]\d*$`. Both were designed through the
identical enumeration, filtering and scoring pipeline used for targeting guides
and screened at ≤3 mismatches. Of 1,408 candidates screened, 662 passed the same
off-target criteria applied to targeting guides; the high attrition among
olfactory receptors reflects their extensive paralogy. The final controls
comprise 254 guides over 52 pan-essential genes and 250 guides over 97 olfactory
receptor genes.

## Module rebalancing

The readout-specificity analysis (below) showed the curated Tier B module was
skewed toward the bone programme (54 bone-annotated versus 22 tooth-annotated
genes, 2.45 : 1) and covered only 44 % of the GO tooth programme. All 39 tooth-
programme genes absent from the manifest were therefore added and designed
through the identical pipeline: 6,013 candidates enumerated, the 858
highest-ranked screened at ≤3 mismatches (59,292 sites), 679 passing the same
off-target criteria, and 376 guides selected across all 39 genes (mean 9.6 per
gene). After this addition the module contains 58 tooth-only, 13 shared and 54
bone-only genes (0.93 : 1) and covers 99 % of the GO tooth programme.

Genes were added rather than removed: deleting bone-annotated genes would
discard genuine candidates to improve a ratio, whereas adding the missing tooth
genes corrects the bias without reducing discovery power.

## Final pool

The complete pool is 9,631 oligonucleotides: 8,127 targeting across 861 genes
(84.4 %), 1,000 non-targeting, 254 pan-essential and 250 olfactory-receptor
(15.6 % controls in total). Median aggregate CFD specificity at ≤3 mismatches was comparable across
classes (targeting 0.278, pan-essential 0.302, olfactory 0.227), confirming the
controls are not systematically cleaner or dirtier than the guides they
calibrate. All 9,255 oligos are 80 nt, free of internal BsmBI and BsaI sites, and
non-duplicated.

## Power analysis

Count dispersion was estimated from published raw sgRNA counts rather than
assumed. Counts for 71,090 sgRNAs across 12 samples of a genome-wide knockout
screen (GEO accession GSE226977) were normalised to counts per million, and
per-sgRNA dispersion was estimated by the method of moments within each
three-replicate arm. Endpoint arms gave a median dispersion of φ = 0.299 (rising
to 0.90 in the most abundant quantile) with a mean inter-replicate Pearson
correlation of 0.43 on log₂ CPM; T0 arms were effectively Poisson (φ ≈ 0,
r = 0.90) and were treated as technical rather than biological replicates.

Total count dispersion was modelled as φ_total = φ_bio + 1/cells + 1/reads, with
φ_bio = 0.288 obtained by subtracting the sequencing and cell-sampling terms of
that experiment from the measured value. Screens were simulated by drawing
negative-binomial counts for control and treated arms, with the per-guide effect
attenuated by a Beta(2, 1.33) efficacy multiplier (mean 0.6) to represent
variable guide activity. Gene-level statistics were the mean log₂ fold change
across guides and replicates; the null distribution was built from 4,000
simulated non-effect genes, p-values taken from its empirical CDF, and power
computed as the fraction of 2,000 effect genes passing Benjamini–Hochberg
control at FDR 0.05. Analytical gradients of the dispersion parameterisation were
verified by sweeping φ_bio across 0.05–0.80.

At the measured dispersion, 10 guides per gene with 4 replicates gives 80 % power
at a gene-level log₂ fold change of 1.0 (8 × 3 gives 57 %; 8 × 4 gives 72 %;
12 × 4 gives 87 %). Power was insensitive to library representation over
200–1000 cells per guide (0.80, 0.80, 0.82 at 10 guides × 4 replicates), because
the 1/cells term is negligible relative to φ_bio.

## Network propagation model

A two-layer mean-aggregator GraphSAGE was implemented in NumPy with
hand-derived gradients, verified against central finite differences to a maximum
relative error of 2.9 × 10⁻⁸. Hidden width 32; Adam, learning rate 0.02; L2
penalty 3 × 10⁻⁴; masked binary cross-entropy; 400–600 epochs.

The graph was the STRING v12 human interaction network restricted to the target
set at combined confidence ≥ 0.7, retrieved through the STRING network API,
giving 817 nodes and 5,824 undirected edges (mean degree 14.3; 74 isolated
nodes). Node features were GTEx v8 median TPM across 54 tissues (log1p, z-scored
per tissue) concatenated with log node degree.

The model was trained to recover held-out members of the pyrophosphate/SIBLING
mineralization module (46 positives), using 60/40 stratified splits repeated 20
times with independent initialisations. Held-out AUROC was 0.805 ± 0.088 for
GraphSAGE versus 0.704 ± 0.060 for L2-regularised logistic regression on the same
features without the graph (paired Wilcoxon over splits, p = 2 × 10⁻⁴); a
degree-only control reached 0.387, below chance, excluding hub recovery as an
explanation. Base positive rate in held-out sets was 0.058.

## Readout-specificity analysis

A pre-registered comparison against published osteoblast mineralization screens
could not be performed: a search of GEO DataSets for CRISPR screen series
matching bone, osteoblast, osteogenic, mineralization or mesenchymal
differentiation terms returned 108 records, none of which is a screen of
osteogenic or mineralization differentiation.

The analysis was therefore re-posed as a separability test on Gene Ontology
annotations. Human annotations were taken from the EBI `goa_human.gaf` release
(17,779 annotated terms), with term closure over `is_a` and `part_of` relations
from the EBI QuickGO descendants endpoint. The tooth programme was defined as
odontogenesis (GO:0042476) ∪ tooth mineralization (GO:0034505), 72 genes; the
bone programme as ossification (GO:0001503) ∪ osteoblast differentiation
(GO:0001649) ∪ bone mineralization (GO:0030282), 183 genes.

The two sets share 13 genes (Jaccard 0.054; 18.1 % of the tooth programme). The
GraphSAGE model above, trained instead to discriminate tooth-only from bone-only
genes over 25 repeated splits, reached AUROC 0.582 ± 0.084 (one-sample Wilcoxon
against 0.5, p = 3.3 × 10⁻⁴), compared with 0.805 for module recovery on the same
graph and features. Within the curated Tier B module, 54 genes are annotated to
the bone programme and 22 to the tooth programme; the library covers 77 % of
genes shared between programmes but 39 % of tooth-only genes.

## Software

Python 3.11; NumPy; pandas; SciPy; scikit-learn 1.9.1; LightGBM 4.7.0; `rs3`
0.0.18; networkx; Matplotlib. Analyses ran on an 8-core x86-64 workstation with
16 GB RAM and no GPU.

## Data availability

Reference data: GRCh38/hg38 and `ncbiRefSeqCurated`/`ncbiRefSeqSelect` (UCSC
Genome Browser); Ensembl BioMart (Ensembl Genes); GTEx v8 median gene
expression; STRING v12; Gene Ontology annotations (EBI `goa_human.gaf`) and
QuickGO; sgRNA counts from GEO GSE226977; CFD scoring tables from the CRISPOR
distribution of Doench et al. (2016).

## Limitations

Off-target enumeration covers the primary assembly only, to 3 mismatches, with
NGG, NAG and NGA PAMs; sites at 4 or more mismatches, on alternate contigs, or
adjacent to other non-canonical PAMs are not recovered. The guide-efficacy
distribution used in the power simulation is an assumption, not a measurement.
GTEx contains no dental tissue, so node features for the network model are a
proxy for the relevant cell states. The dispersion estimate derives from a single
published screen in an unrelated cell type and should be re-measured in a pilot
before the design is committed.
