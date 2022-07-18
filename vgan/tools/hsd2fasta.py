from Bio import SeqIO, Seq, SeqRecord
import os
import csv

missing_hsds = os.listdir("../data/hsds")

for file in missing_hsds:
    with open("../data/RSRS.fa", "r") as g:
        RSRS = list(SeqIO.parse(g, "fasta"))[0]
        reader = csv.reader(open("missing_hsds/"+file, "r"), delimiter="\t")
        for row in reader:
            RSRS_map = {i:base for i, base in enumerate(str(RSRS.seq))}
            for var in row[2:]:
                if "-" not in var and "." not in var and var[-1] != "D":
                    base = int(var[:-1])
                    derived = var[-1]
                    RSRS_map[base-1] = derived
                elif var[-1] == "D" and "-" not in var:
                    base = int(var[1:-1])
                    RSRS_map[base-1] = ""
                elif "-" in var:
                    base = int(var.split("-")[0])
                    del_end = int(var.split("-")[1][:-1])
                    for b in range(base, del_end + 1):
                        RSRS_map[b-1] = ""
                elif "." in var:
                    base = int(var.split(".")[0])
                    ins_size = int(var.split(".")[1][0])
                    ins_seq = var.split(".")[1][1:]
                    RSRS_map[base-1] = ins_seq

            rec_seq = ""
            for i, base in sorted(RSRS_map.items()):
                rec_seq += RSRS_map[i]

            rec = SeqRecord.SeqRecord(id=file.split(".")[0], description="", seq=Seq.Seq(rec_seq))

            with open("testing"+file.split(".")[0]+".fasta", "w") as h:
                SeqIO.write(rec, h, "fasta")

