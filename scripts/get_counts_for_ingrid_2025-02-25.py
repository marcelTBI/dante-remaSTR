from common import HGVSRecord

file = "../../../inputs/polydante/new_set.hgvs.cut.tsv"

with open(file) as f:
    lines = f.readlines()

for line in lines:
    line = line.strip()
    record = HGVSRecord(line)
    print(len(record.units[0][0]))
