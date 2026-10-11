"""Verify the library's pooled-synthesis flanks against YOUR vector.

This is the one check that cannot be done from the guide sequences alone: the
flanks are homology arms, so they must match the vector's U6 3' end and sgRNA
scaffold 5' end exactly, or Gibson/Golden-Gate assembly will fail silently.

Usage
-----
    python verify_flanks.py my_vector.gbk
    python verify_flanks.py my_vector.fasta
    python verify_flanks.py --seq "ATGC..."

Paste or export the full plasmid sequence from your map (SnapGene, Benchling,
Addgene) in any format; only A/C/G/T characters are read.

Exit code 0 = flanks confirmed, 1 = mismatch (do not order).
"""
import re
import sys

F5 = "GGAAAGGACGAAACACCG"
SCAFFOLDS = {
    "original (Hsu 2013)":
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG",
    "F+E (Chen 2013)":
        "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCG",
}
ORDER_FILE = {"original (Hsu 2013)": "oligo_pool_order.csv",
              "F+E (Chen 2013)": "oligo_pool_order_FE.csv"}


def load(path):
    text = open(path).read()
    if "ORIGIN" in text:                      # GenBank
        text = text.split("ORIGIN", 1)[1]
    text = re.sub(r"^>.*$", "", text, flags=re.M)   # FASTA headers
    return re.sub(r"[^ACGTacgt]", "", text).upper()


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def find(seq, motif):
    """Search both strands of a circular plasmid."""
    doubled = seq + seq                       # handle the origin wrap
    hits = [m.start() % len(seq) for m in re.finditer(motif, doubled)][:1]
    rhits = [m.start() % len(seq) for m in re.finditer(revcomp(motif), doubled)][:1]
    return hits, rhits


def main(argv):
    if len(argv) >= 3 and argv[1] == "--seq":
        seq = re.sub(r"[^ACGTacgt]", "", argv[2]).upper()
    elif len(argv) >= 2:
        seq = load(argv[1])
    else:
        print(__doc__)
        return 1
    print("vector length read: %d bp" % len(seq))

    f, r = find(seq, F5)
    if not (f or r):
        print("FAIL  5' flank %s not found in the vector." % F5)
        print("      Your U6 3' end differs — do not order; send me the map region.")
        return 1
    print("OK    5' flank found at position %d (%s strand)"
          % ((f or r)[0], "+" if f else "-"))

    matched = []
    for name, sc in SCAFFOLDS.items():
        f2, r2 = find(seq, sc)
        if f2 or r2:
            matched.append(name)
            print("OK    scaffold matches %-20s at position %d" % (name, (f2 or r2)[0]))
    if not matched:
        print("FAIL  neither known scaffold found.")
        print("      Your vector uses a different tracrRNA — do not order.")
        for name, sc in SCAFFOLDS.items():
            print("      expected %-20s %s..." % (name, sc[:30]))
        return 1
    if len(matched) > 1:
        print("WARN  more than one scaffold matched; inspect the map manually.")
        return 1

    print("\nCONFIRMED — order %s" % ORDER_FILE[matched[0]])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
