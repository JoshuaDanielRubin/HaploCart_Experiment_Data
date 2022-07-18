import subprocess
from matplotlib import pyplot as plt
import pickle
import numpy as np

plt.rcParams.update({'font.size': 7})
plt.rc('legend', fontsize=10)
plt.rc('figure', titlesize=12)
plt.rc('axes', titlesize=10)
plt.rc('axes', labelsize=10)
plt.rc('xtick', labelsize=9)
plt.rc('ytick', labelsize=9)

haplogrep_path="../dep/haplogrep/haplogrep"
thousand_genome_dict = {"HG00473":"D6a1a","HG01051":"K1a4a1h","HG02666":"L3e4a","HG03112":"L1b1a3","NA18510":"L0a1a3","NA19036":"L3b1a1a","NA19661":"D1h1","NA20518":"J2a1a1e"}

def is_correct(score):
    if str(score) == "0":
        return 1
    else:
        return 0

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


def make_correctness_plot_fq(haplogrep_dict, haplocart_dict, depthfile, outfile, numt=False):
    depth_dict = {}
    with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file, depth = split[0].split("/")[-1].split("_mem.bam")[0], float(split[1])
            depth_dict.update({file:depth})
    hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_1, hg_1_12, hg_12_14, hg_14_16, hg_16_18, hg_18_2 \
    = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_1, hc_1_12, hc_12_14, hc_14_16, hc_16_18, hc_18_2 \
    = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    for k1,v1 in haplogrep_dict.items():
        v1 = int(v1)
        depth1 = depth_dict[dequote(k1).split("_mem.bam")[0]]
        if depth1 <= 0.1:
            hg_0_01.append(is_correct(v1))
        elif depth1 > 0.1 and depth1 <= 0.2:
            hg_01_02.append(is_correct(v1))
        elif depth1 > 0.2 and depth1 <= 0.3:
            hg_02_03.append(is_correct(v1))
        elif depth1 > 0.3 and depth1 <= 0.4:
            hg_03_04.append(is_correct(v1))
        elif depth1 > 0.4 and depth1 <= 0.5:
            hg_04_05.append(is_correct(v1))
        elif depth1 > 0.5 and depth1 <= 0.6:
            hg_05_06.append(is_correct(v1))
        elif depth1 > 0.6 and depth1 <= 0.7:
            hg_06_07.append(is_correct(v1))
        elif depth1 > 0.7 and depth1 <= 0.8:
            hg_07_08.append(is_correct(v1))
        elif depth1 > 0.8 and depth1 <= 0.9:
            hg_08_09.append(is_correct(v1))
        elif depth1 > 0.9 and depth1 <= 1:
            hg_09_1.append(is_correct(v1))
        elif depth1 > 1 and depth1 <= 1.2:
            hg_1_12.append(is_correct(v1))
        elif depth1 > 1.2 and depth1 <= 1.4:
            hg_12_14.append(is_correct(v1))
        elif depth1 > 1.4 and depth1 <= 1.6:
            hg_14_16.append(is_correct(v1))
        elif depth1 > 1.6 and depth1 <= 1.8:
            hg_16_18.append(is_correct(v1))
        elif depth1 > 1.8 and depth1 <= 2:
            hg_18_2.append(is_correct(v1))


    for k2,v2 in haplocart_dict.items():
        v2 = int(v2)
        if dequote(k2).endswith("tmp"):
            splitted = dequote(k2).split("_")
            depth2 = depth_dict["_".join([splitted[0], "n1000", "l"+splitted[1], "rep"+splitted[2], splitted[3]])]
        else:
            depth2 = depth_dict[dequote(k2)]
        if depth2 <= 0.1:
            hc_0_01.append(is_correct(v2))
        elif depth2 > 0.1 and depth2 <= 0.2:
            hc_01_02.append(is_correct(v2))
        elif depth2 > 0.2 and depth2 <= 0.3:
            hc_02_03.append(is_correct(v2))
        elif depth2 > 0.3 and depth2 <= 0.4:
            hc_03_04.append(is_correct(v2))
        elif depth2 > 0.4 and depth2 <= 0.5:
            hc_04_05.append(is_correct(v2))
        elif depth2 > 0.5 and depth2 <= 0.6:
            hc_05_06.append(is_correct(v2))
        elif depth2 > 0.6 and depth2 <= 0.7:
            hc_06_07.append(is_correct(v2))
        elif depth2 > 0.7 and depth2 <= 0.8:
            hc_07_08.append(is_correct(v2))
        elif depth2 > 0.8 and depth2 <= 0.9:
            hc_08_09.append(is_correct(v2))
        elif depth2 > 0.9 and depth2 <= 1:
            hc_09_1.append(is_correct(v2))
        elif depth2 > 1 and depth2 <= 1.2:
            hc_1_12.append(is_correct(v2))
        elif depth2 > 1.2 and depth2 <= 1.4:
            hc_12_14.append(is_correct(v2))
        elif depth2 > 1.4 and depth2 <= 1.6:
            hc_14_16.append(is_correct(v2))
        elif depth2 > 1.6 and depth2 <= 1.8:
            hc_16_18.append(is_correct(v2))
        elif depth2 > 1.8 and depth2 <= 2:
            hc_18_2.append(is_correct(v2))


    hg_data = [hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_1, \
               hg_1_12, hg_12_14, hg_14_16, hg_16_18, hg_18_2]

    hc_data = [hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_1, \
               hc_1_12, hc_12_14, hc_14_16, hc_16_18, hc_18_2]

    hg_data_prop = [(sum(x) / len(x)) for x in hg_data]
    hc_data_prop = [(sum(x) / len(x)) for x in hc_data]

    plt.title("Proportion of Correct Predictions  \n (Simulated FASTQ)")
    plt.xlabel("Average depth of coverage (x)")
    plt.ylabel("Proportion")

    hg_pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29]
    hc_pos = [0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5]

    hg_bar = plt.bar(hg_pos, hg_data_prop, color="blue", label="Haplogrep2")
    hc_bar = plt.bar(hc_pos, hc_data_prop, color="orange", label="HaploCart")
    plt.bar_label(hg_bar, labels = ["N="+str(len(x)) for x in hg_data])
    plt.bar_label(hc_bar, labels = ["N="+str(len(x)) for x in hc_data])

    plt.xticks([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29], \
               ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1", \
                "1-1.2", "1.2-1.4", "1,4-1.6", "1.6-1.8", "1.8-2.0"], rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

def make_correctness_plot_bam(haplogrep_dict, haplocart_dict, depthfile, outfile):
     bamdepth_dict = {}
     with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file, depth = split[0].split("/")[-1].split("_mem.bam")[0], float(split[1])
            bamdepth_dict.update({file:depth})

     hg_005_01, hg_01_015, hg_015_02, hg_02_025, hg_025_03, hg_03_035, hg_035_04, hg_04_045, hg_045_05, hg_05_055, \
     hg_055_06, hg_06_065, hg_065_07, hg_07_075, hg_075_08, hg_08_085, hg_085_09, hg_09_095, hg_095_1, hg_1_105, hg_105_11, hg_11_115, hg_115_12\
     = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

     hc_005_01, hc_01_015, hc_015_02, hc_02_025, hc_025_03, hc_03_035, hc_035_04, hc_04_045, hc_045_05, hc_05_055, \
     hc_055_06, hc_06_065, hc_065_07, hc_07_075, hc_075_08, hc_08_085, hc_085_09, hc_09_095, hc_095_1, hc_1_105, hc_105_11, hc_11_115, hc_115_12 \
     = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

     for k1,v1 in haplogrep_dict.items():
         v1 = is_correct(int(v1))
         depth1 = bamdepth_dict[dequote(k1).split("_mem.bam")[0]]
         if depth1 > 0.05 and depth1 <= 0.1:
             hg_005_01.append(v1)
         elif depth1 > 0.1 and depth1 <= 0.15:
             hg_01_015.append(v1)
         elif depth1 > 0.15 and depth1 <= 0.2:
             hg_015_02.append(v1)
         elif depth1 > 0.2 and depth1 <= 0.25:
             hg_02_025.append(v1)
         elif depth1 > 0.25 and depth1 <= 0.3:
             hg_025_03.append(v1)
         elif depth1 > 0.3 and depth1 <= 0.35:
             hg_03_035.append(v1)
         elif depth1 > 0.35 and depth1 <= 0.4:
             hg_035_04.append(v1)
         elif depth1 > 0.4 and depth1 <= 0.45:
             hg_04_045.append(v1)
         elif depth1 > 0.45 and depth1 <= 0.5:
             hg_045_05.append(v1)
         elif depth1 > 0.5 and depth1 <= 0.55:
             hg_05_055.append(v1)
         elif depth1 > 0.55 and depth1 <= 0.6:
             hg_055_06.append(v1)
         elif depth1 > 0.6 and depth1 <= 0.65:
             hg_06_065.append(v1)
         elif depth1 > 0.65 and depth1 <= 0.7:
             hg_065_07.append(v1)
         elif depth1 > 0.7 and depth1 <= 0.75:
             hg_07_075.append(v1)
         elif depth1 > 0.75 and depth1 <= 0.8:
             hg_075_08.append(v1)
         elif depth1 > 0.8 and depth1 <= 0.85:
            hg_08_085.append(v1)
         elif depth1 > 0.85 and depth1 <= 0.9:
            hg_085_09.append(v1)
         elif depth1 > 0.9 and depth1 <= 0.95:
            hg_09_095.append(v1)
         elif depth1 > 0.95 and depth1 <= 1:
            hg_095_1.append(v1)
         elif depth1 > 1 and depth1 <= 1.05:
            hg_1_105.append(v1)
         elif depth1 > 1.05 and depth1 <= 1.1:
            hg_105_11.append(v1)
         elif depth1 > 1.1 and depth1 <= 1.15:
            hg_11_115.append(v1)
         elif depth1 > 1.15 and depth1 <= 1.2:
            hg_115_12.append(v1)



     for k2,v2 in haplocart_dict.items():
         v2 = is_correct(int(v2))
         depth2 = bamdepth_dict[dequote(k2).split("_mem.bam")[0]]
         if depth2 > 0.05 and depth2 <= 0.1:
             hc_005_01.append(v2)
         elif depth2 > 0.1 and depth2 <= 0.15:
             hc_01_015.append(v2)
         elif depth2 > 0.15 and depth2 <= 0.2:
             hc_015_02.append(v2)
         elif depth2 > 0.2 and depth2 <= 0.25:
             hc_02_025.append(v2)
         elif depth2 > 0.25 and depth2 <= 0.3:
             hc_025_03.append(v2)
         elif depth2 > 0.3 and depth2 <= 0.35:
             hc_03_035.append(v2)
         elif depth2 > 0.35 and depth2 <= 0.4:
             hc_035_04.append(v2)
         elif depth2 > 0.4 and depth2 <= 0.45:
             hc_04_045.append(v2)
         elif depth2 > 0.45 and depth2 <= 0.5:
             hc_045_05.append(v2)
         elif depth2 > 0.5 and depth2 <= 0.55:
             hc_05_055.append(v2)
         elif depth2 > 0.55 and depth2 <= 0.6:
             hc_055_06.append(v2)
         elif depth2 > 0.6 and depth2 <= 0.65:
             hc_06_065.append(v2)
         elif depth2 > 0.65 and depth2 <= 0.7:
             hc_065_07.append(v2)
         elif depth2 > 0.7 and depth2 <= 0.75:
             hc_07_075.append(v2)
         elif depth2 > 0.75 and depth2 <= 0.8:
             hc_075_08.append(v2)
         elif depth2 > 0.8 and depth2 <= 0.85:
             hc_08_085.append(v2)
         elif depth2 > 0.85 and depth2 <= 0.9:
             hc_085_09.append(v2)
         elif depth2 > 0.9 and depth2 <= 0.95:
             hc_09_095.append(v2)
         elif depth2 > 0.95 and depth2 <= 1:
             hc_095_1.append(v2)
         elif depth2 > 1 and depth2 <= 1.05:
             hc_1_105.append(v2)
         elif depth2 > 1.05 and depth2 <= 1.1:
             hc_105_11.append(v2)
         elif depth2 > 1.1 and depth2 <= 1.15:
             hc_11_115.append(v2)
         elif depth2 > 1.15 and depth2 <= 1.2:
             hc_115_12.append(v2)



     hg_data = [hg_005_01, hg_01_015, hg_015_02, hg_02_025, hg_025_03, hg_03_035, hg_035_04, hg_04_045, hg_045_05, hg_05_055, \
                hg_055_06, hg_06_065, hg_065_07, hg_07_075, hg_075_08, hg_08_085, hg_085_09, hg_09_095, hg_095_1, \
                hg_1_105, hg_105_11, hg_11_115]

     hc_data = [hc_005_01, hc_01_015, hc_015_02, hc_02_025, hc_025_03, hc_03_035, hc_035_04, hc_04_045, hc_045_05, hc_05_055, \
                hc_055_06, hc_06_065, hc_065_07, hc_07_075, hc_075_08, hc_08_085, hc_085_09, hc_09_095, hc_095_1, hc_1_105, \
                hc_105_11, hc_11_115]

     hg_data_prop = [(sum(x) / len(x)) for x in hg_data]
     hc_data_prop = [(sum(x) / len(x)) for x in hc_data]

     hg_pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43]
     hc_pos = [x-0.5 for x in hg_pos]

     plt.figure(figsize=(15, 5))
     plt.title("Proportion of Correct Predictions  \n (Thousand Genomes)")
     plt.xlabel("Average depth of coverage (x)")
     plt.ylabel("Proportion")

     hg_bar = plt.bar(hg_pos, hg_data_prop, color="blue", label="Haplogrep2")
     hc_bar = plt.bar(hc_pos, hc_data_prop, color="orange", label="HaploCart")
     plt.bar_label(hg_bar, labels = ["N="+str(len(x)) for x in hg_data])
     plt.bar_label(hc_bar, labels = ["N="+str(len(x)) for x in hc_data])

     plt.xticks([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43], \
               ["0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.35", \
                "0.35-0.4", "0.4-0.45", "0.45-0.5", "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", \
                "0.7-0.75", "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1", \
                "1-1.05", "1.05-1.1", "1.1-1.15"], rotation=45)

     plt.legend(loc="upper left")
     plt.tight_layout()
     plt.savefig(outfile, dpi=300)
     plt.close()


def plot_bam():
    depth_dict = {}
    with open("../data/1kdepths.txt", "r") as f:
        for line in f:
            split=line.split()
            file, depth = dequote(split[0].split("/")[-1]), float(split[1])
            depth_dict.update({file:depth})

    haplogrep_score_dict = pickle.load(open("../data/pickles/haplogrep_bam.pk", "rb"))
    haplocart_score_dict = pickle.load(open("../data/pickles/haplocart_bams.pk", "rb"))

    make_correctness_plot_bam(haplogrep_score_dict, haplocart_score_dict, "../data/1kdepths.txt", "../data/pngs/bam_correctness.png")


def plot_single_fastq():
    haplogrep_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    haplocart_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_no_numt.pk", "rb"))
    make_correctness_plot_fq(haplogrep_score_dict, haplocart_score_dict, "../data/fastq_no_numt_sim_depths.txt", "../data/pngs/fastq_correctness.png")

plot_bam()
plot_single_fastq()
