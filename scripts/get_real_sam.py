from Bio import SeqIO

MAX_REP = 10
SHORT = 3


def apply_vcf(seq: str, chrom: str, start: int, end: int, vcf_h1: list[str], haplotype: int) -> tuple[str, str]:
    i = 0
    result = []
    ecigar = []
    for line in vcf_h1:
        chr_, p, _, ref, alt, _, _, _, _, _ = line.strip().split("\t")
        pos = int(p)

        if chr_ == chrom and start <= pos < end:
            # print(chr_, pos, ref, alt)
            pos = pos - start - 1
            # print(seq[0:pos], seq[pos], seq[pos + 1:])
            result.append(seq[i:pos])
            ecigar.append("M" * len(seq[i:pos]))
            i += len(seq[i:pos])

            result.append(alt)
            i += len(ref)
            change = len(alt) - len(ref)
            if change == 0:
                ecigar.append("X" * len(alt))
            elif change > 0:
                ecigar.append("M" + "I" * change)
            elif change < 0:
                ecigar.append("M" + "D" * abs(change))

            # print("".join(result))
    result.append(seq[i:])
    ecigar.append("M" * len(seq[i:]))
    return ("".join(result), "".join(ecigar))


def visualize_alignment(s1: str, s2: str, ecigar: str) -> tuple[str, str]:
    j = 0
    k = 0
    s1_res = []
    s2_res = []
    for c in ecigar:
        if c == "M" or c == "X":
            s1_res.append(s1[j])
            s2_res.append(s2[k])
            j += 1
            k += 1
        elif c == "D":
            s1_res.append(s1[j])
            j += 1
            s2_res.append("-")
        elif c == "I":
            s1_res.append("-")
            s2_res.append(s2[k])
            k += 1
        else:
            raise NotImplementedError
    return ("".join(s1_res), "".join(s2_res))


def add_slop(s: str, start_slop: int, end_slop: int):
    return s[:start_slop] + " " + s[start_slop:-end_slop] + " " + s[-end_slop:]


def correct_slop(slop: int, ecigar: str) -> tuple[int, int]:
    start_slop = 0
    i = 0
    while start_slop != slop:
        # if ecigar[i] == "M" or ecigar[i] == "X":
        if ecigar[i] != "I":
            start_slop += 1
        i += 1
    start_slop = i

    end_slop = 0
    i = 0
    while end_slop != slop:
        # if ecigar[-i - 1] == "M" or ecigar[-i - 1] == "X":
        if ecigar[-i - 1] != "I":
            end_slop += 1
        i += 1
    end_slop = i
    return start_slop, end_slop


def greedy_sol(seq: str, rep_matrix: list[list[int]]) -> str:
    col = len(seq) - 1
    result = []

    while col > -1:
        column = [row[col] for row in rep_matrix]
        maxi = max(column)
        size = column.index(maxi) + 1
        rep = seq[col + 1 - size:col + 1]
        result.append((rep, maxi))
        col -= maxi * size
    result.reverse()
    # print("".join([f"{seq}[{n}]" for seq, n in result]))

    # remove short repeats (e.g. A[2] -> A[1]A[1])
    result2 = []
    for seq, n in result:
        if len(seq) * n <= SHORT:
            result2.extend([(seq, 1)] * n)
        else:
            result2.append((seq, n))
    # print("".join([f"{seq}[{n}]" for seq, n in result2]))

    # merge nonrepeating neighbors (e.g. A[1]C[1] -> AC[1])
    result3 = []
    buffer = []
    for seq, n in result2:
        if n == 1:
            buffer.append(seq)
        else:
            if len(buffer) != 0:
                tmp = "".join(buffer)
                result3.append((tmp, 1))
                buffer = []
            result3.append((seq, n))
    if len(buffer) != 0:
        tmp = "".join(buffer)
        result3.append((tmp, 1))
        buffer = []

    return "".join([f"{seq}[{n}]" for seq, n in result3])


def convert_to_hgvs(seq: str) -> str:
    seq = seq.replace("-", "")
    # print(seq)
    matrix = []
    for i in range(1, MAX_REP + 1):
        tmp = [0] * (i - 1) + [1] * (len(seq) - i + 1)
        for p in range(2 * i - 1, len(seq)):
            if seq[p + 1 - 2 * i:p + 1 - i] == seq[p + 1 - i:p + 1]:
                tmp[p] = tmp[p - i] + 1
        # print(tmp)
        matrix.append(tmp)

    result = greedy_sol(seq, matrix)
    return result


def adjust_cigar(ref: str, seq: str, ecigar: str, slop: int) -> str:
    # MMMMMMMMMMMMMMMMMMMMMMMMMMIIIIIMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    # MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM IIIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    result = list(ecigar)
    for i in range(slop):
        if result[i] == "I" or result[i] == "D":
            idx1 = result[i:].index("M") if "M" in result[i:] else 1_000_000
            idx2 = result[i:].index("X") if "X" in result[i:] else 1_000_000
            idx = min(idx1, idx2)
            result[i], result[i + idx] = result[i + idx], result[i]

    result.reverse()
    for i in range(slop):
        if result[i] == "I" or result[i] == "D":
            idx1 = result[i:].index("M") if "M" in result[i:] else 1_000_000
            idx2 = result[i:].index("X") if "X" in result[i:] else 1_000_000
            idx = min(idx1, idx2)
            result[i], result[i + idx] = result[i + idx], result[i]

    result.reverse()
    new_ecigar = "".join(result)
    test_ref, test_seq = visualize_alignment(ref, seq, new_ecigar)
    for i in range(len(new_ecigar)):
        if new_ecigar[i] == "M" and test_ref[i] != test_seq[i]:
            print("unable to adjust cigar")
            return ecigar
        if new_ecigar[i] == "X" and test_ref[i] == test_seq[i]:
            print("unable to adjust cigar")
            return ecigar

    return new_ecigar


ref = "./../../../inputs/GRCh38_decoy.fa"
vcf_h1 = "./../../../cache/GIAB_h1.vcf"
vcf_h2 = "./../../../cache/GIAB_h2.vcf"
regions = "./../../../cache/regions.bed"
slop = 30

with open(ref) as f:
    tmp = list(SeqIO.parse(f, "fasta"))
    fa_records = dict([(x.id, x) for x in tmp])

with open(regions) as f:
    lines = f.readlines()

with open(vcf_h1) as f:
    vcf_h1_lines = f.readlines()
    vcf_h1_lines = list(filter(lambda x: not x.startswith("#"), vcf_h1_lines))

with open(vcf_h2) as f:
    vcf_h2_lines = f.readlines()
    vcf_h2_lines = list(filter(lambda x: not x.startswith("#"), vcf_h2_lines))

for line in lines:
    chrom, s, e, desc = line.strip().split("\t")
    start = int(s)
    end = int(e)

    slop_start = start - slop
    slop_end = end + slop
    seq = str(fa_records[chrom][slop_start:slop_end].seq)

    new_seq1, ecigar1 = apply_vcf(seq, chrom, slop_start, slop_end, vcf_h1_lines, 0)
    ref1, seq1 = visualize_alignment(seq, new_seq1, ecigar1)
    start_slop1, end_slop1 = correct_slop(slop, ecigar1)

    adj_ecigar1 = adjust_cigar(seq, new_seq1, ecigar1, slop)
    old_seq1, new_seq1 = visualize_alignment(seq, new_seq1, adj_ecigar1)

    new_seq2, ecigar2 = apply_vcf(seq, chrom, slop_start, slop_end, vcf_h2_lines, 0)
    ref2, seq2 = visualize_alignment(seq, new_seq2, ecigar2)
    start_slop2, end_slop2 = correct_slop(slop, ecigar2)

    adj_ecigar2 = adjust_cigar(seq, new_seq2, ecigar2, slop)
    old_seq2, new_seq2 = visualize_alignment(seq, new_seq2, adj_ecigar2)

    # if slop == start_slop1 == start_slop2 == end_slop1 == end_slop2:
    #     continue

    print("{:>16}:  {}".format("desc", desc))

    print("{:>16}:  {}".format("hg38", add_slop(ref1, start_slop1, end_slop1)))
    print("{:>16}:  {}".format("GIAB_h1", add_slop(seq1, start_slop1, end_slop1)))
    print("{:>16}:  {}".format("", add_slop(ecigar1, start_slop1, end_slop2)))

    print("{:>16}:  {}".format("adj_hg38", add_slop(old_seq1, slop, slop)))
    print("{:>16}:  {}".format("adj_GIAB_h1", add_slop(new_seq1, slop, slop)))
    print("{:>16}:  {}".format("", add_slop(adj_ecigar1, slop, slop)))

    print("{:>16}:  {}".format("hgvs", convert_to_hgvs(new_seq1[start_slop1:-end_slop1])))
    print()

    print("{:>16}:  {}".format("hg38", add_slop(ref2, start_slop2, end_slop2)))
    print("{:>16}:  {}".format("GIAB_h1", add_slop(seq2, start_slop2, end_slop2)))
    print("{:>16}:  {}".format("", add_slop(ecigar2, start_slop2, end_slop2)))

    print("{:>16}:  {}".format("adj_hg38", add_slop(old_seq2, slop, slop)))
    print("{:>16}:  {}".format("adj_GIAB_h2", add_slop(new_seq2, slop, slop)))
    print("{:>16}:  {}".format("", add_slop(adj_ecigar2, slop, slop)))

    print("{:>16}:  {}".format("hgvs", convert_to_hgvs(new_seq2[start_slop2:-end_slop2])))
    print()

    print()
    # exit()
