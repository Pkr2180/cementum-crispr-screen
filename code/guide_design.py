"""SpCas9 guide enumeration over coding exons, hg38.

Genome is stored per chromosome as uint8 arrays: A=0, C=1, G=2, T=3, N/other=4.
All genomic coordinates here are 0-based half-open internally; Ensembl input
coordinates are 1-based inclusive and are converted on load.
"""
import json
import numpy as np

B2C = np.array(list("ACGTN"))
COMP = {0: 3, 1: 2, 2: 1, 3: 0, 4: 4}
_CACHE = {}


def load_chrom(chrom):
    if chrom not in _CACHE:
        _CACHE.clear()  # one chromosome resident at a time
        _CACHE[chrom] = np.load(f"genome/chr/{chrom}.npy")
    return _CACHE[chrom]


def to_str(arr):
    return "".join(B2C[np.clip(arr, 0, 4)])


def revcomp_arr(arr):
    out = arr[::-1].copy()
    m = out < 4
    out[m] = 3 - out[m]
    return out


# ---------------------------------------------------------------- transcripts

def cds_intervals(tx):
    """CDS intervals of one transcript as 0-based half-open [start, end)."""
    cs, ce = tx.get("cds_start"), tx.get("cds_end")
    if cs is None or ce is None:
        return []
    out = []
    for s, e in tx["exons"]:
        a, b = max(s, cs), min(e, ce)
        if a <= b:
            out.append((a - 1, b))
    return sorted(out)


def gene_cds_model(gene):
    """Union CDS intervals, per-base constitutive fraction, and the canonical
    transcript's CDS coordinate map."""
    coding = [t for t in gene["transcripts"]
              if t["biotype"] == "protein_coding" and cds_intervals(t)]
    if not coding:
        return None
    canon = next((t for t in coding if t.get("canonical")), None) or \
        max(coding, key=lambda t: sum(e - s for s, e in cds_intervals(t)))

    ivs = [cds_intervals(t) for t in coding]
    lo = min(s for iv in ivs for s, _ in iv)
    hi = max(e for iv in ivs for _, e in iv)
    cover = np.zeros(hi - lo, dtype=np.int32)
    for iv in ivs:
        for s, e in iv:
            cover[s - lo:e - lo] += 1
    constit = cover / len(ivs)

    canon_iv = cds_intervals(canon)
    strand = gene["strand"]
    ordered = canon_iv if strand == 1 else canon_iv[::-1]
    cds_len = sum(e - s for s, e in ordered)
    return {"lo": lo, "hi": hi, "constit": constit, "n_tx": len(ivs),
            "canon_id": canon["id"], "canon_iv": ordered, "cds_len": cds_len,
            "strand": strand}


def gene_cds_model_rs(gene):
    """CDS model from RefSeq structures (intervals already 0-based half-open)."""
    coding = [t for t in gene["transcripts"] if t["cds"]]
    if not coding:
        return None
    canon = next((t for t in coding if t.get("select")), None) or \
        max(coding, key=lambda t: sum(e - s for s, e in t["cds"]))
    ivs = [t["cds"] for t in coding]
    lo = min(s for iv in ivs for s, _ in iv)
    hi = max(e for iv in ivs for _, e in iv)
    cover = np.zeros(hi - lo, dtype=np.int32)
    for iv in ivs:
        for s, e in iv:
            cover[s - lo:e - lo] += 1
    strand = gene["strand"]
    canon_iv = sorted(canon["cds"]) if strand == 1 else sorted(canon["cds"])[::-1]
    return {"lo": lo, "hi": hi, "constit": cover / len(ivs), "n_tx": len(ivs),
            "canon_id": canon["id"], "canon_iv": canon_iv,
            "cds_len": sum(e - s for s, e in canon_iv), "strand": strand,
            "chrom": gene["chrom"]}


def cds_fraction(model, pos):
    """Fractional position of genomic `pos` along the canonical CDS, or None."""
    off = 0
    for s, e in model["canon_iv"]:
        if s <= pos < e:
            d = (pos - s) if model["strand"] == 1 else (e - 1 - pos)
            return (off + d) / model["cds_len"]
        off += e - s
    return None


# ------------------------------------------------------------------ filtering

BSMBI = ("CGTCTC", "GAGACG")
BSAI = ("GGTCTC", "GAGACC")


def seq_ok(p):
    """Hard sequence filters on the 20nt protospacer."""
    if "TTTT" in p:
        return "polyT"
    gc = (p.count("G") + p.count("C")) / 20.0
    if gc < 0.30 or gc > 0.80:
        return "gc"
    for site in BSMBI + BSAI:
        if site in p:
            return "re_site"
    run = 1
    for i in range(1, 20):
        run = run + 1 if p[i] == p[i - 1] else 1
        if run >= 5:
            return "homopolymer"
    return None


# --------------------------------------------------------------- enumeration

PAM_SET = {"AGG", "TGG", "CGG", "GGG"}


def enumerate_gene(gene, model, min_frac=0.05, max_frac=0.65, min_constit=0.9):
    """Candidate NGG guides whose Cas9 cut site lies in constitutive CDS.

    Returns tuples (protospacer, pam, context30, start, strand, cds_frac,
    constitutive_frac). `start` is the 0-based genomic start of the
    protospacer on the plus strand.
    """
    seq = load_chrom(model["chrom"])
    lo, hi = model["lo"], model["hi"]
    pad = 40
    a = max(0, lo - pad)
    b = min(len(seq), hi + pad)
    sub = seq[a:b]
    constit_ok = model["constit"] >= min_constit
    out = []

    for strand in (1, -1):
        if strand == 1:
            # protospacer [p, p+20), PAM [p+20, p+23) with PAM[1:3] == GG
            hits = np.where((sub[1:-1] == 2) & (sub[2:] == 2))[0]
            cut = a + hits - 4                      # ps - 20 + 16
        else:
            # PAM [q, q+3) = CCN on plus strand, protospacer [q+3, q+23)
            hits = np.where((sub[:-1] == 1) & (sub[1:] == 1))[0]
            cut = a + hits + 6
        if hits.size == 0:
            continue
        m = (cut >= lo) & (cut < hi)
        hits, cut = hits[m], cut[m]
        if hits.size == 0:
            continue
        m = constit_ok[cut - lo]
        hits, cut = hits[m], cut[m]

        for h, c in zip(hits.tolist(), cut.tolist()):
            f = cds_fraction(model, c)
            if f is None or not (min_frac <= f <= max_frac):
                continue
            if strand == 1:
                p = h - 20
                if p < 4 or h + 6 > sub.size:
                    continue
                ctx = sub[p - 4:p + 26]
                if ctx.size != 30 or ctx.max() > 3:
                    continue
                out.append((to_str(sub[p:p + 20]), to_str(sub[h:h + 3]),
                            to_str(ctx), a + p, 1, f,
                            float(model["constit"][c - lo])))
            else:
                if h < 3 or h + 27 > sub.size:
                    continue
                ctx = revcomp_arr(sub[h - 3:h + 27])
                if ctx.size != 30 or ctx.max() > 3:
                    continue
                out.append((to_str(revcomp_arr(sub[h + 3:h + 23])),
                            to_str(revcomp_arr(sub[h:h + 3])), to_str(ctx),
                            a + h + 3, -1, f,
                            float(model["constit"][c - lo])))
    return out


def load_structures(path="handoff/gene_structures.jsonl"):
    genes = {}
    for ln in open(path):
        g = json.loads(ln)
        genes[g["gene_id"]] = g
    return genes
