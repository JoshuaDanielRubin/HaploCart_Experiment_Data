import gzip
import os
import subprocess
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
import pickle
from Bio import SeqIO
import sys
import Levenshtein

def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z

def get_edit_distance(s1, s2, numt=False):
    # Compute edit (Levenshtein) distance between true and predicted haplotypes

    if s2 == "U5a2b1":
        return "NAN"

    if s2 == "C4a1e":
        return "NAN"

    if s2 == "mt-MRCA":
        return "NAN"

    if s2 == "E1b1":
        return "NAN"

    if s2=="B":
        return "NAN"

    with gzip.open("../data/synthetic_fastas/"+s1.split("+")[0]+".fasta.gz", "rt") as f:
        seq1 = next(SeqIO.parse(f, "fasta")).seq
    with gzip.open("../data/synthetic_fastas/"+s2.split("+")[0]+".fasta.gz", "rt") as f:
        seq2 = next(SeqIO.parse(f, "fasta")).seq

    return Levenshtein.distance(seq1, seq2)


thousand_genomes_dict = {"HG00473":"D6a1a","HG01051":"K1a4a1h","HG02666":"L3e4a","HG03112":"L1b1a3",\
                         "NA18510":"L0a1a3","NA19036":"L3b1a1a","NA19661":"D1h1","NA20518":"J2a1a1e"}

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if '"' not in s and "'" not in s:
        return s
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    elif s.startswith(("'", '"')):
        return s[1:]

def score(true, pred, numt=False):
    return get_edit_distance(true, pred, numt=numt)

def get_haplocheck_scores(score_dict, result_filepath, true, bam=False, numt=False):
    with open(result_filepath, "rb") as g:
        next(g)
        for line in g:
            split = [x.decode('utf-8') for x in line.split()]
            if numt:
                file, pred = split[0].split("_")[2], split[1]
            else:
                file, pred = split[0], split[1]
            if bam == False:
                true = file.split("_")[0]
            predscore = score(dequote(true), dequote(pred), numt=numt)
            print(split[0].split("_mem.bam")[0], pred, predscore, "haplocheck")
            if numt:
                score_dict.update({dequote(split[0].split("_mem.bam")[0]):predscore})
            else:
                score_dict.update({file:predscore})
        return score_dict


def get_haplogrep_mask_scores(score_dict, result_dirpath):
    for file in os.listdir(result_dirpath):
        with gzip.open(result_dirpath + file, "r") as g:
            next(g)
            for line in g:
                split = [x.decode('utf-8') for x in line.split()]
                _, pred = split[0], split[1]
                true = file.split("/")[-1].split("_")[0]
                predscore = score(dequote(true), dequote(pred))
                score_dict.update({file:predscore})
    return score_dict

def get_phymer_scores(phymer_dir):
    score_dict = {}
    for file in os.listdir(phymer_dir):
        if "mask" not in file:
            continue
        with open(phymer_dir + file, "r") as g:
            for line in g:
                pm_split = line.split()
                print(pm_split)
                id = pm_split[0]
                pred = pm_split[1].split("-")[0].split(",")[0]
                predscore = score(id, pred)
                score_dict.update({file:predscore})
    return score_dict

def get_haplocart_scores(haplocart_pred_file, bam=False, numt=False):
    haplocart_pred_dict = {}
    with open(haplocart_pred_file, "r") as f:
        for line in f:
            split = line.strip().split()
            if numt:
                file, pred = split[0].split("_")[0], split[1]
            elif bam:
                file, pred = split[0].split("/")[-1], split[1]
            else:
                file, pred = split[0].split("/")[0], split[1]

            if bam == False:
                predscore = score(file.split("_")[0], pred)
                print(file.split("_")[0], pred, predscore)
            else:
                predscore = score(thousand_genomes_dict[split[0].split("_")[0]], pred)
                print(pred, predscore)

            if numt:
                haplocart_pred_dict.update({dequote(split[0].split("_mem.bam")[0]):predscore})

            else:
                haplocart_pred_dict.update({file:predscore})

    return haplocart_pred_dict


def score_bam():

    #haplocheck_score_dict = {}
    #for id in thousand_genome_dict.keys():
    #    haplocheck_score_dict = get_haplocheck_scores(haplocheck_score_dict, "../src/simulations/thousand_genomes/haplocheck_results/"+id+"/"+"haplogroups/haplogroups.txt", thousand_genome_dict[id], bam=True)

    #pickle.dump(haplocheck_score_dict, open("../data/pickles/haplocheck_bam.pk", "wb"))

    haplocart_score_dict = get_haplocart_scores("../data/haplocart_results/bams.txt", bam=True)
    pickle.dump(haplocart_score_dict, open("../data/pickles/haplocart_bams.pk", "wb"))


def score_fastq():

    #haplocheck_score_dict_no_numt = {}
    #haplocheck_score_dict_no_numt = get_haplocheck_scores(haplocheck_score_dict_no_numt, \
    #                         "../src/simulations/thousand_genomes/haplocheck_results/sims/haplogroups/haplogroups.txt", None)
    #with open("../data/pickles/hc_fastq_no_numt.pk", "wb") as g:
    #    pickle.dump(haplocheck_score_dict_no_numt, g)

    #haplocheck_score_dict_with_numt = {}
    #haplocheck_score_dict_with_numt = get_haplocheck_scores(haplocheck_score_dict_with_numt, \
    #                         "../src/simulations/thousand_genomes/haplocheck_results/sims_numts/haplogroups/haplogroups.txt", None, numt=True)
    #with open("../data/pickles/hc_fastq_with_numt.pk", "wb") as h:
    #    pickle.dump(haplocheck_score_dict_with_numt, h)

    haplocart_score_dict_no_numt = get_haplocart_scores("../data/haplocart_results/fastq_no_numt.txt")
    pickle.dump(haplocart_score_dict_no_numt, open("../data/pickles/haplocart_fastq_no_numt.pk", "wb"))

    haplocart_score_dict_with_numt = get_haplocart_scores("../data/haplocart_results/fastq_with_numt.txt", numt=True)
    pickle.dump(haplocart_score_dict_with_numt, open("../data/pickles/haplocart_fastq_with_numt.pk", "wb"))


def score_mask():
    #haplocheck_score_dict = {}
    #haplocheck_score_dict = get_haplogrep_mask_scores(haplocheck_score_dict, \
    #                        "../src/simulations/thousand_genomes/haplocheck_results/mask/")

    #with open("../data/pickles/haplogrep_mask.pk", "wb") as g:
    #    pickle.dump(haplocheck_score_dict, g)

    haplocart_score_dict = get_haplocart_scores("../data/haplocart_results/mask.txt")
    pickle.dump(haplocart_score_dict, open("../data/pickles/haplocart_mask_no_alternate_minimizer.pk", "wb"))

    #phymer_score_dict = get_phymer_scores("../data/phymer_mask/")
    #pickle.dump(phymer_score_dict, open("../data/pickles/phymer_mask.pk", "wb"))


def get_haplogrep_reported_confidence_bam():
    confidence_dict = {}

    for id in thousand_genome_dict.keys():
        file = "../src/simulations/thousand_genomes/haplocheck_results/" + id + "/haplogroups/haplogroups.txt"
        with open(file, "r") as g:
            next(g)
            for line in g:
                split = line.split()
                qual = split[3]
                pred = split[1]
                sample = split[0]
                confidence_dict.update({sample:qual})

    with open("../data/pickles/haplogrep_reported_quals_bam.pk", "wb") as f:
        pickle.dump(confidence_dict, f)


def get_haplogrep_reported_confidence_fastq():
    confidence_dict = {}

    file = "../src/simulations/thousand_genomes/haplocheck_results/sims/haplogroups/haplogroups.txt"
    with open(file, "r") as g:
        next(g)
        for line in g:
            split = line.split()
            qual = split[3]
            pred = split[1]
            sample = split[0]
            print(sample, qual)
            confidence_dict.update({sample:qual})

    with open("../data/pickles/haplogrep_reported_quals_fastq.pk", "wb") as f:
        pickle.dump(confidence_dict, f)


#score_bam()
score_fastq()
#score_mask()
#get_haplogrep_reported_confidence_fastq()


