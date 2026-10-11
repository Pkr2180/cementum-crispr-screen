"""Genome-wide SpCas9 off-target enumeration and CFD/MIT scoring.

Search strategy: for each query protospacer we expand every sequence within
`max_mm` substitutions into a sorted uint64 table (2 bits/base, 20nt = 40 bits),
then stream every PAM-adjacent 20mer in the genome through `searchsorted`.
Exact by construction for the stated mismatch budget; PAMs NGG/NAG/NGA.
"""
import pickle
import numpy as np
from itertools import combinations, product

B2I = {"A": 0, "C": 1, "G": 2, "T": 3}
I2B = "ACGT"
# (second, third) PAM base codes accepted: NGG, NAG, NGA
PAM_PAIRS = [(2, 2), (0, 2), (2, 0)]


def encode(seq):
    c = 0
    for b in seq:
        c = (c << 2) | B2I[b]
    return np.uint64(c)


def decode(code, n=20):
    c = int(code)
    return "".join(I2B[(c >> (2 * (n - 1 - i))) & 3] for i in range(n))


def popcount(x):
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(x).astype(np.uint8)
    v = x.copy()
    out = np.zeros(x.shape, dtype=np.uint8)
    while v.any():
        out += (v & np.uint64(1)).astype(np.uint8)
        v >>= np.uint64(1)
    return out


def n_mismatch(a, b):
    x = np.bitwise_xor(a, b)
    folded = np.bitwise_and(np.bitwise_or(x, np.right_shift(x, np.uint64(1))),
                            np.uint64(0x5555555555555555))
    return popcount(folded)


def expand_variants(codes, max_mm):
    """All codes within <=max_mm substitutions. Returns (codes, query_idx, mm)."""
    n = len(codes)
    out_c, out_q, out_m = [], [], []
    base = np.asarray(codes, dtype=np.uint64)
    out_c.append(base); out_q.append(np.arange(n, dtype=np.uint32))
    out_m.append(np.zeros(n, dtype=np.uint8))
    for mm in range(1, max_mm + 1):
        for pos in combinations(range(20), mm):
            shifts = [np.uint64(2 * (19 - p)) for p in pos]
            masks = [np.uint64(~(np.uint64(3) << s) & np.uint64(0xFFFFFFFFFF))
                     for s in shifts]
            cur_bases = [np.bitwise_and(np.right_shift(base, s), np.uint64(3))
                         for s in shifts]
            for alts in product(range(4), repeat=mm):
                keep = np.ones(n, dtype=bool)
                v = base.copy()
                for s, msk, cb, alt in zip(shifts, masks, cur_bases, alts):
                    keep &= (cb != alt)
                    v = np.bitwise_or(np.bitwise_and(v, msk),
                                      np.left_shift(np.uint64(alt), s))
                if keep.any():
                    out_c.append(v[keep])
                    out_q.append(np.arange(n, dtype=np.uint32)[keep])
                    out_m.append(np.full(keep.sum(), mm, dtype=np.uint8))
    return (np.concatenate(out_c), np.concatenate(out_q), np.concatenate(out_m))


CHUNK = 20_000_000


def genome_sites(seq, comp, chunk=CHUNK):
    """Yield (codes, protospacer_start, strand, pam2code) for all NGG/NAG/NGA.

    Scans the chromosome in windows so that no temporary exceeds ~chunk
    elements; each site is emitted exactly once, assigned to the window
    holding its PAM start.
    """
    L = len(seq)
    for w0 in range(0, L, chunk):
        w1 = min(L, w0 + chunk)
        a = max(0, w0 - 32)
        b = min(L, w1 + 32)
        s = seq[a:b]
        c = comp[a:b]
        n = s.size
        ncum = np.concatenate([[0], np.cumsum(s == 4)])

        for strand in (1, -1):
            ok = np.zeros(n, dtype=bool)
            if strand == 1:
                for b2, b3 in PAM_PAIRS:
                    ok[:n - 2] |= (s[1:n - 1] == b2) & (s[2:n] == b3)
            else:
                for b2, b3 in PAM_PAIRS:
                    ok[:n - 2] |= (c[1:n - 1] == b2) & (c[0:n - 2] == b3)
            p = np.flatnonzero(ok)
            if p.size == 0:
                continue
            # keep PAM starts inside this window, with the full site in range
            gp = p + a
            m = (gp >= w0) & (gp < w1)
            if strand == 1:
                m &= (p >= 20) & (p + 3 <= n)
            else:
                m &= (p + 23 <= n)
            p = p[m]
            if p.size == 0:
                continue
            if strand == 1:
                lo_l = p - 20
                pam2 = (s[p + 1].astype(np.uint8) << 2) | s[p + 2].astype(np.uint8)
            else:
                lo_l = p + 3
                pam2 = (c[p + 1].astype(np.uint8) << 2) | c[p].astype(np.uint8)
            clean = (ncum[lo_l + 20] - ncum[lo_l]) == 0
            lo_l, pam2 = lo_l[clean], pam2[clean]
            if lo_l.size == 0:
                continue
            codes = np.zeros(lo_l.size, dtype=np.uint64)
            buf = np.empty(lo_l.size, dtype=np.int64)
            for k in range(20):
                np.add(lo_l, (k if strand == 1 else 19 - k), out=buf)
                src = s if strand == 1 else c
                codes |= np.left_shift(src[buf].astype(np.uint64),
                                       np.uint64(2 * (19 - k)))
            yield codes, (lo_l + a).astype(np.uint32), strand, pam2


# ------------------------------------------------------------------- scoring

def load_cfd(mm_path="cfd/mismatch_score.pkl", pam_path="cfd/pam_scores.pkl"):
    mm = pickle.load(open(mm_path, "rb"), encoding="latin1")
    pam = pickle.load(open(pam_path, "rb"), encoding="latin1")
    return mm, pam


_RC = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}


def cfd_score(wt, off, pam2, mm_scores, pam_scores):
    """Doench 2016 CFD. wt/off are 20nt DNA; pam2 is the last two PAM bases."""
    score = 1.0
    w = wt.replace("T", "U")
    s = off.replace("T", "U")
    for i, (a, b) in enumerate(zip(w, s)):
        if a == b:
            continue
        score *= mm_scores["r" + a + ":d" + _RC[b] + "," + str(i + 1)]
    return score * pam_scores[pam2]


# Hsu 2013 position weights used by the MIT specificity score (5'->3').
MIT_W = np.array([0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000,
                  0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
                  0.615, 0.804, 0.685, 0.583])


def mit_score(wt, off):
    pos = [i for i, (a, b) in enumerate(zip(wt, off)) if a != b]
    n = len(pos)
    if n == 0:
        return 1.0
    t1 = float(np.prod(1 - MIT_W[pos]))
    d = (max(pos) - min(pos)) / (n - 1) if n > 1 else 19.0
    t2 = 1.0 / (((19 - d) / 19.0) * 4 + 1)
    t3 = 1.0 / (n ** 2)
    return t1 * t2 * t3
