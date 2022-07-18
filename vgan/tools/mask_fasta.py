import sys
from Bio import Seq, SeqIO
import random

fastafile = sys.argv[1]
outfile = sys.argv[2]
N = int(sys.argv[3])

with open(fastafile, "r") as f:
    fasta_in = list(SeqIO.parse(f, "fasta"))[0]
    assert(N <= len(fasta_in.seq))
    new_seq = str(fasta_in.seq)

    n = random.randint(0, len(fasta_in.seq))
    if n < len(fasta_in.seq) - N:
        new_seq = new_seq[:n] + "N"*N + new_seq[n+N:]
    else:
        to_junction = len(fasta_in.seq) - n
        from_junction = abs((len(fasta_in.seq)-n)-N)
        new_seq = 'N'*from_junction + new_seq[from_junction:(len(fasta_in.seq) - to_junction)] + 'N'*to_junction

    assert(len(new_seq) == len(fasta_in.seq))
    fasta_in.seq = Seq.Seq(new_seq)

    with open(outfile, "wt") as g:
        SeqIO.write(fasta_in, g, "fasta")
