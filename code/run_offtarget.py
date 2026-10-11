"""Drive the genome-wide off-target search over hg38."""
import numpy as np
import offtarget as ot

CHROMS = ["chr%s" % c for c in list(range(1, 23)) + ["X", "Y"]]


def search(query_codes, max_mm, chroms=CHROMS, verbose=True):
    """Return a dict of arrays: one row per (query, genomic site) hit."""
    vc, vq, vm = ot.expand_variants(query_codes, max_mm)
    order = np.argsort(vc, kind="stable")
    vc, vq, vm = vc[order], vq[order], vm[order]
    if verbose:
        print("  variants: %.1fM" % (vc.size / 1e6), flush=True)

    # Bitmap prefilter on the top 30 bits of each 40-bit code. Rejects ~97% of
    # genomic sites with two cheap gathers, so the binary search (the real cost)
    # only runs on survivors.
    SHIFT = np.uint64(10)
    bm = np.zeros(2 ** 27, dtype=np.uint8)
    bi = (vc >> SHIFT).astype(np.int64)
    np.bitwise_or.at(bm, bi >> 3, (np.uint8(1) << (bi & 7).astype(np.uint8)))
    del bi

    out = {k: [] for k in ("query", "mm", "chrom", "pos", "strand", "pam2", "code")}
    for ci, chrom in enumerate(chroms):
        seq = np.load("genome/chr/%s.npy" % chrom)
        comp = np.where(seq < 4, 3 - seq, 4).astype(np.uint8)
        for codes, lo, strand, pam2 in ot.genome_sites(seq, comp):
            gi = (codes >> SHIFT).astype(np.int64)
            keep = np.nonzero((bm[gi >> 3] >> (gi & 7).astype(np.uint8)) & 1)[0]
            del gi
            if keep.size == 0:
                continue
            sub = codes[keep]
            li = np.searchsorted(vc, sub, "left")
            ri = np.searchsorted(vc, sub, "right")
            cnt = (ri - li).astype(np.int64)
            nz = np.nonzero(cnt)[0]
            if nz.size == 0:
                continue
            idx = keep[nz]
            li, cnt = li[nz], cnt[nz]
            c = cnt
            total = int(c.sum())
            site = np.repeat(idx, c)
            base = np.repeat(np.cumsum(c) - c, c)
            vidx = np.repeat(li, c) + (np.arange(total) - base)
            out["query"].append(vq[vidx])
            out["mm"].append(vm[vidx])
            out["chrom"].append(np.full(total, ci, dtype=np.uint8))
            out["pos"].append(lo[site].astype(np.uint32))
            out["strand"].append(np.full(total, strand, dtype=np.int8))
            out["pam2"].append(pam2[site])
            out["code"].append(codes[site])
        del seq, comp
        if verbose and (ci % 6 == 0 or ci == len(chroms) - 1):
            n = sum(a.size for a in out["query"])
            print("    %s  hits so far %.2fM" % (chrom, n / 1e6), flush=True)
    return {k: (np.concatenate(v) if v else np.array([], dtype=np.int64))
            for k, v in out.items()}
