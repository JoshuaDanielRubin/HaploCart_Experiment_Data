import subprocess
from matplotlib import pyplot as plt
import pickle
import numpy as np
import matplotlib.patches as mpatches

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
    hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20 \
    = [], [], [], [], [], [], [], [], [], []
    hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20  \
    = [], [], [], [], [], [], [], [], [], []
    for k1,v1 in haplogrep_dict.items():
        v1 = int(v1)
        depth1 = depth_dict[dequote(k1).split("_mem.bam")[0]]
        if depth1 <= 0.2:
            hg_0_02.append(is_correct(v1))
        elif depth1 > 0.2 and depth1 <= 0.4:
            hg_02_04.append(is_correct(v1))
        elif depth1 > 0.4 and depth1 <= 0.6:
            hg_04_06.append(is_correct(v1))
        elif depth1 > 0.6 and depth1 <= 0.8:
            hg_06_08.append(is_correct(v1))
        elif depth1 > 0.8 and depth1 <= 1.0:
            hg_08_10.append(is_correct(v1))
        elif depth1 > 1.0 and depth1 <= 1.2:
            hg_10_12.append(is_correct(v1))
        elif depth1 > 1.2 and depth1 <= 1.4:
            hg_12_14.append(is_correct(v1))
        elif depth1 > 1.4 and depth1 <= 1.6:
            hg_14_16.append(is_correct(v1))
        elif depth1 > 1.6 and depth1 <= 1.8:
            hg_16_18.append(is_correct(v1))
        elif depth1 > 1.8 and depth1 <= 2:
            hg_18_20.append(is_correct(v1))

    for k2,v2 in haplocart_dict.items():
        v2 = int(v2)
        if dequote(k2).endswith("tmp"):
            splitted = dequote(k2).split("_")
            depth2 = depth_dict["_".join([splitted[0], "n1000", "l"+splitted[1], "rep"+splitted[2], splitted[3]])]
        else:
            depth2 = depth_dict[dequote(k2)]
        if depth2 <= 0.2:
            hc_0_02.append(is_correct(v2))
        elif depth2 > 0.2 and depth2 <= 0.4:
            hc_02_04.append(is_correct(v2))
        elif depth2 > 0.4 and depth2 <= 0.6:
            hc_04_06.append(is_correct(v2))
        elif depth2 > 0.6 and depth2 <= 0.8:
            hc_06_08.append(is_correct(v2))
        elif depth2 > 0.8 and depth2 <= 1.0:
            hc_08_10.append(is_correct(v2))
        elif depth2 > 1.0 and depth2 <= 1.2:
            hc_10_12.append(is_correct(v2))
        elif depth2 > 1.2 and depth2 <= 1.4:
            hc_12_14.append(is_correct(v2))
        elif depth2 > 1.4 and depth2 <= 1.6:
            hc_14_16.append(is_correct(v2))
        elif depth2 > 1.6 and depth2 <= 1.8:
            hc_16_18.append(is_correct(v2))
        elif depth2 > 1.8 and depth2 <= 2:
            hc_18_20.append(is_correct(v2))

    total_per_bin = [8678, 8530, 4386, 1243, 1151, 652, 560, 36, 534, 618]

    hc_data = [hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20]
    hg_data = [hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20]

    hg_data_prop = [(sum(x) / len(x)) for x in hg_data]
    hc_data_prop = [(sum(x) / len(x)) for x in hc_data]

    hg_call_rate = [(len(x) / total_per_bin[i]) for i, x in enumerate(hg_data)]
    hc_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hc_data)]

    plt.title("Proportion of Correct Predictions  \n (Simulated FASTQ)")
    plt.xlabel("Average depth of coverage (X)")
    plt.ylabel("Proportion")

    hg_pos = [1,3,5,7,9,11,13,15,17,19]
    hg_callrate_pos = [x - 0.25 for x in hg_pos]
    hc_pos = [x - 0.5 for x in hg_pos]
    hc_callrate_pos = [x - 0.75 for x in hg_pos]

    hc_bar = plt.bar(hc_pos, hc_data_prop, color="orange", label="HaploCart proportion correct", width=0.2)
    hg_bar = plt.bar(hg_pos, hg_data_prop, color="blue", label="Haplogrep2 proportion correct", width=0.2)
    hc_callrate_bar = plt.bar(hc_callrate_pos, hc_call_rate, color="red", label="HaploCart call rate", width=0.2)
    hg_callrate_bar = plt.bar(hg_callrate_pos, hg_call_rate, color="green", label="HaploGrep2 call rate", width=0.2)

    plt.xticks([1,3,5,7,9,11,13,15,17,19], \
               ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0"], rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

def make_correctness_plot_bam(haplogrep_dict, haplocart_dict, depthfile, outfile):
     #print("Size of haplocart dict: ", len(haplocart_dict.keys()))
     #print("Size of haplogrep dict: ", len(haplogrep_dict.keys()))
     bamdepth_dict = {}
     with open(depthfile, "r") as f:
         for line in f:
            split=line.split()
            file, depth = split[0].split("/")[-1].split("_mem.bam")[0], float(split[1])
            bamdepth_dict.update({file:depth})

     hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_10, hg_10_up \
     = [], [], [], [], [], [], [], [], [], [], []

     hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_10, hc_10_up \
     = [], [], [], [], [], [], [], [], [], [], []

     for k1,v1 in haplogrep_dict.items():
         v1 = is_correct(int(v1))
         depth1 = bamdepth_dict[dequote(k1).split("_mem.bam")[0]]
         if depth1 > 0 and depth1 <= 0.1:
             hg_0_01.append(v1)
         elif depth1 > 0.1 and depth1 <= 0.2:
             hg_01_02.append(v1)
         elif depth1 > 0.2 and depth1 <= 0.3:
             hg_02_03.append(v1)
         elif depth1 > 0.3 and depth1 <= 0.4:
             hg_03_04.append(v1)
         elif depth1 > 0.4 and depth1 <= 0.5:
             hg_04_05.append(v1)
         elif depth1 > 0.5 and depth1 <= 0.6:
             hg_05_06.append(v1)
         elif depth1 > 0.6 and depth1 <= 0.7:
             hg_06_07.append(v1)
         elif depth1 > 0.7 and depth1 <= 0.8:
             hg_07_08.append(v1)
         elif depth1 > 0.8 and depth1 <= 0.9:
             hg_08_09.append(v1)
         elif depth1 > 0.9 and depth1 <= 1.0:
             hg_09_10.append(v1)
         elif depth1 > 1.0:
             hg_10_up.append(v1)


     for k2,v2 in haplocart_dict.items():
         v2 = is_correct(int(v2))
         depth2 = bamdepth_dict[dequote(k2).split("_mem.bam")[0]]
         if depth2 > 0 and depth2 <= 0.1:
            hc_0_01.append(v2)
         elif depth2 > 0.1 and depth2 <= 0.2:
            hc_01_02.append(v2)
         elif depth2 > 0.2 and depth2 <= 0.3:
            hc_02_03.append(v2)
         elif depth2 > 0.3 and depth2 <= 0.4:
            hc_03_04.append(v2)
         elif depth2 > 0.4 and depth2 <= 0.5:
            hc_04_05.append(v2)
         elif depth2 > 0.5 and depth2 <= 0.6:
            hc_05_06.append(v2)
         elif depth2 > 0.6 and depth2 <= 0.7:
            hc_06_07.append(v2)
         elif depth2 > 0.7 and depth2 <= 0.8:
            hc_07_08.append(v2)
         elif depth2 > 0.8 and depth2 <= 0.9:
            hc_08_09.append(v2)
         elif depth2 > 0.9 and depth2 <= 1.0:
            hc_09_10.append(v2)
         elif depth2 > 1.0:
            hc_10_up.append(v2)

     total_per_bin = [578, 915, 978, 921, 946, 897, 842, 783, 636, 370, 133]

     hg_data = [hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_10, hg_10_up]
     hc_data = [hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_10,hc_10_up]

     nans = [float('nan'), float('nan')]

     hg_data_prop = [(sum(_) / len(_)) if len(_) > 0 else nans for _ in hg_data]
     hc_data_prop = [(sum(__) / len(__)) if len(__) > 0 else nans for __ in hc_data]

     hg_call_rate = [(len(x) / total_per_bin[i]) for i, x in enumerate(hg_data)]
     hc_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hc_data)]

     #print("hg data lens: ", [len(x) for x in hg_data])
     #print("hc data lens: ", [len(x) for x in hc_data])
     #print("hg call rate: ", hg_call_rate)
     #print("hc call rate: ", hc_call_rate)

     hg_pos = [1,3,5,7,9,11,13,15,17,19,21]
     hg_callrate_pos = [x - 0.5 for x in hg_pos]
     hc_pos = [x-0.25 for x in hg_pos]
     hc_callrate_pos = [x - 0.75 for x in hg_pos]

     plt.figure(figsize=(15, 5))
     plt.title("Proportion of Correct Predictions  \n (Empirical Paired-end FASTQ)")
     plt.xlabel("Average depth of coverage (X)")
     plt.ylabel("Proportion")

     hg_bar = plt.bar(hg_pos, hg_data_prop, color="blue", label="HaploGrep2", width=0.3)
     hc_bar = plt.bar(hc_pos, hc_data_prop, color="orange", label="HaploCart", width=0.3)
     hg_callrate_bar = plt.bar(hg_callrate_pos, hg_call_rate, color="green", label="HaploGrep2 call rate", width=0.3)
     hc_callrate_bar = plt.bar(hc_callrate_pos, hc_call_rate, color="red", label="HaploCart call rate", width=0.3)

     plt.xticks([1,3,5,7,9,11,13,15,17,19,21], \
               ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", \
                "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0", ">=1.0"], rotation=45)

     hc_patch = mpatches.Patch(color='orange', label='HaploCart proportion correct')
     hg_patch = mpatches.Patch(color='blue', label='HaploGrep2 proportion correct')
     hc_callrate_patch = mpatches.Patch(color='red', label='HaploCart call rate')
     hg_callrate_patch = mpatches.Patch(color='green', label='HaploGrep2 call rate')

     plt.legend(handles=[hc_callrate_patch, hg_callrate_patch, hc_patch, hg_patch], loc="upper left")
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
#plot_single_fastq()
