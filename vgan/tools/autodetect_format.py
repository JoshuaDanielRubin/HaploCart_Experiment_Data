import gzip
import sys
from Bio import SeqIO
import subprocess
import os

def is_fasta(filename):
    with open(filename, "rt") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_zipped_fasta(filename):
    with gzip.open(filename, "rt") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_gam(filename):
    try:
        validate_gam = subprocess.check_output("vg view -a " + filename + " -j > /dev/null 2> /dev/null", shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def autodetect(input_file):

    if not os.path.exists(input_file):
        raise Exception("NONEXISTENT FILE")

    # https://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    is_binary = is_binary_string(open(input_file, 'rb').read(1024))

    if not is_binary:
        # The Vgan team fully supports the rights of the non-binary
        with open(input_file, "rt") as f:
            try:
                records = list(SeqIO.parse(f, "fastq"))
                # Make sure it's bona fide FASTQ
                for record in records:
                    score=record.letter_annotations["phred_quality"]
                if records[0].id.split("/")[0] == records[1].id.split("/")[0]:
                    return "Interleaved"
                else:
                    return "Non-interleaved"
            except:
                if is_fasta(input_file):
                    return "FASTA"
                else:
                    return "ERROR IN FILE"
    else:
        with gzip.open(input_file, "rt") as g:
            if is_gam(input_file):
                return "GAM"
            else:
                try:
                    records = SeqIO.parse(g, "fastq")
                    for record in records:
                        score=record.letter_annotations["phred_quality"]
                    return "FASTQ"
                except:
                    if is_zipped_fasta(input_file):
                       return "FASTA"
                    else:
                        return "ERROR IN FILE"

print(autodetect(sys.argv[1]))
