#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio import Align

from common import HGVSRecord


def get_middle(ref_motif, hgvs_motif):
    middle = ["|"] * len(ref_motif)
    for i in range(len(ref_motif)):
        if ref_motif[i] == hgvs_motif[i]:
            middle[i] = " "
    middle = "".join(middle)
    return middle


def is_exact(m: int, g: int) -> str | None:
    if m == 0 and g == 0:
        return "exact"
    return None


def is_mismatched(m: int, g: int) -> str | None:
    if m > 0 and g == 0:
        return "mismatched"
    return None


def is_edited(g: int) -> str | None:
    if g > 0:
        return "edited"
    return None


# %%
# infile = "/tmp/tmp2.tsv"
# ref = "../../../inputs/GRCh38_decoy.fa"
# slop = 30

infile = sys.argv[1]
ref = sys.argv[2]
slop = int(sys.argv[3])

# %%
f = open(infile)
lines = f.readlines()
aligner = Align.PairwiseAligner(mode="global", match_score=0, mismatch_score=-1, gap_score=-1)

with open(ref) as f:
    tmp = list(SeqIO.parse(f, "fasta"))
    fa_records = dict([(x.id, x) for x in tmp])

# %%
for line in lines:
    locus_id, hgvs = line.split("\t")
    hgvs_rec = HGVSRecord(hgvs)

    left_flank = fa_records[hgvs_rec.chrom][
        hgvs_rec.start - 1 - slop : hgvs_rec.start - 1
    ].seq
    ref_motif = fa_records[hgvs_rec.chrom][hgvs_rec.start - 1 : hgvs_rec.end].seq
    hgvs_motif = "".join([seq * n for seq, n in hgvs_rec.units])
    right_flank = fa_records[hgvs_rec.chrom][hgvs_rec.end : hgvs_rec.end + slop].seq

    alignments = aligner.align(ref_motif, hgvs_motif)
    n = len(alignments)
    m = alignments[0].counts().mismatches
    g = alignments[0].counts().gaps

    ref_motif = alignments[0][0, :]
    hgvs_motif = alignments[0][1, :]
    middle = get_middle(ref_motif, hgvs_motif)

    flags: set[str | None] = {None}
    flags.add(is_exact(m, g))
    flags.add(is_mismatched(m, g))
    flags.add(is_edited(g))
    flags.remove(None)

    print(locus_id, n, m, g, list(flags))
    # print(fa_records[hgvs_rec.chrom][hgvs_rec.start - 1 - slop : hgvs_rec.end + slop].seq)
    print(left_flank + ref_motif + right_flank)
    print(left_flank + middle + right_flank)
    print(" " * slop + hgvs_motif)
    print()

    # if "edited" in flags:
    #     print(locus_id, left_flank, str(hgvs_rec), right_flank, sep="\t")
