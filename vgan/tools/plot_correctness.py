import subprocess
from matplotlib import pyplot as plt
import pickle
import numpy as np
import matplotlib.patches as mpatches

colorblind_colors = ['#ff7f00', '#377eb8', '#4daf4a', '#f781bf']

plt.rcParams.update({'font.size': 7})
plt.rc('legend', fontsize=16)
plt.rc('figure', titlesize=20)
plt.rc('axes', titlesize=16)
plt.rc('axes', labelsize=16)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

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


def make_correctness_plot_fq(haplogrep_dict, haplocart_dict, haplogrouper_dict, depthfile, outfile):
    depth_dict = {}
    with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file, depth = split[0].split("/")[-1].split("_mem.bam")[0], float(split[1])
            depth_dict.update({file:depth})
    hc_data \
    = [], [], [], [], [], [], [], [], [], []
    hg_data  \
    = [], [], [], [], [], [], [], [], [], []
    for k1,v1 in haplogrep_dict.items():
        v1 = int(v1)
        if "H2a2a1" == dequote(k1).split("_mem.bam")[0].split("_")[0]:
            continue
        depth1 = depth_dict[dequote(k1).split("_mem.bam")[0]]
        if depth1 <= 0.1:
            hg_data[0].append(is_correct(v1))
        elif depth1 > 0.1 and depth1 <= 0.15:
            hg_data[1].append(is_correct(v1))
        elif depth1 > 0.15 and depth1 <= 0.2:
            hg_data[2].append(is_correct(v1))
        elif depth1 > 0.2 and depth1 <= 0.25:
            hg_data[3].append(is_correct(v1))
        elif depth1 > 0.25 and depth1 <= 0.3:
            hg_data[4].append(is_correct(v1))
        elif depth1 > 0.3 and depth1 <= 0.4:
            hg_data[5].append(is_correct(v1))
        elif depth1 > 0.4 and depth1 <= 0.5:
            hg_data[6].append(is_correct(v1))
        elif depth1 > 0.5 and depth1 <= 0.6:
            hg_data[7].append(is_correct(v1))
        elif depth1 > 0.6 and depth1 <= 1:
            hg_data[8].append(is_correct(v1))
        elif depth1 > 1:
            hg_data[9].append(is_correct(v1))

    for k2,v2 in haplocart_dict.items():
        v2 = int(v2)
        if dequote(k2).endswith("tmp"):
            splitted = dequote(k2).split("_")
            if "H2a2a1" == splitted[0]:
                continue
            depth2 = depth_dict["_".join([splitted[0], "n1000", "l"+splitted[1], "rep"+splitted[2], splitted[3]])]
        else:
            depth2 = depth_dict[dequote(k2)]

        if depth2 <= 0.1:
            hc_data[0].append(is_correct(v2))
        elif depth2 > 0.1 and depth2 <= 0.15:
            hc_data[1].append(is_correct(v2))
        elif depth2 > 0.15 and depth2 <= 0.2:
            hc_data[2].append(is_correct(v2))
        elif depth2 > 0.2 and depth2 <= 0.25:
            hc_data[3].append(is_correct(v2))
        elif depth2 > 0.25 and depth2 <= 0.3:
            hc_data[4].append(is_correct(v2))
        elif depth2 > 0.3 and depth2 <= 0.4:
            hc_data[5].append(is_correct(v2))
        elif depth2 > 0.4 and depth2 <= 0.5:
            hc_data[6].append(is_correct(v2))
        elif depth2 > 0.5 and depth2 <= 0.6:
            hc_data[7].append(is_correct(v2))
        elif depth2 > 0.6 and depth2 <= 1:
            hc_data[8].append(is_correct(v2))
        elif depth2 > 1:
            hc_data[9].append(is_correct(v2))

    total_per_bin = [5219, 6290, 6563, 6296, 5631, 5848, 4253, 4874, 5001, 5025]

    hg_data_prop = [(sum(x)) for x in hg_data]
    hc_data_prop = [(sum(x)) for x in hc_data]
    hgrp_data_prop = [(sum(x)) for x in hgrp_data]

    hg_call_rate = [(len(x) / total_per_bin[i]) for i, x in enumerate(hg_data)]
    hc_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hc_data)]
    hgrp_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hgrp_data)]

    plt.title("Total Number of Correct Predictions  \n (Simulated FASTQ)")
    plt.xlabel("Average depth of coverage (X)")
    plt.ylabel("Count")

    hg_pos = [1,2,3,4,5,6,7,8,9,10]
    hg_callrate_pos = [x for x in hg_pos]
    hc_pos = [x - 0.15 for x in hg_pos]
    hc_callrate_pos = [x-0.15 for x in hg_pos]

    fig, ax = plt.subplots(2, figsize=(15, 10), sharex=True)

    hc_callrate_bar = ax[1].bar(hc_callrate_pos, hc_call_rate, color=colorblind_colors[3], label="HaploCart call rate", width=0.2)
    hg_callrate_bar = ax[1].bar(hg_callrate_pos, hg_call_rate, color=colorblind_colors[2], label="HaploGrep2 call rate", width=0.2)
    hc_bar = ax[0].bar(hc_pos, hc_data_prop, color=colorblind_colors[0], label="HaploCart count correct", width=0.2)
    hg_bar = ax[0].bar(hg_pos, hg_data_prop, color=colorblind_colors[1], label="Haplogrep2 count correct", width=0.2)


    hc_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart count correct')
    hg_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2 count correct')
    hc_callrate_patch = mpatches.Patch(color=colorblind_colors[3], label='HaploCart call rate')
    hg_callrate_patch = mpatches.Patch(color=colorblind_colors[2], label='HaploGrep2 call rate')

    ax[0].legend(handles=[hc_patch, hg_patch], loc="upper left")
    ax[1].legend(handles=[hc_callrate_patch, hg_callrate_patch], loc="upper left")

    ax[0].set_ylabel("Count")
    ax[1].set_ylabel("Call rate")
    fig.suptitle("Total Number of Correct Predictions  \n (Simulated Paired-end FASTQ)")
    plt.xlabel("Mean depth of coverage (X)")

    plt.xticks([1,2,3,4,5,6,7,8,9,10], \
               ["0-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", \
                "0.6-1.0", ">1.0"], \
              rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()

def make_correctness_plot_bam(haplogrep_dict, haplocart_dict, haplogrouper_dict, depthfile, outfile):
     #print("Size of haplocart dict: ", len(haplocart_dict.keys()))
     #print("Size of haplogrep dict: ", len(haplogrep_dict.keys()))
     bamdepth_dict = {}
     with open(depthfile, "r") as f:
         for line in f:
            split=line.split()
            file, depth = split[0].split("/")[-1].split("_mem.bam")[0], float(split[1])
            bamdepth_dict.update({file:depth})

     hg_data \
     = [], [], [], [], [], []

     hc_data \
     = [], [], [], [], [], []

     hgrp_data \
     = [], [], [], [], [], []

     for k1,v1 in haplogrep_dict.items():
         v1 = is_correct(int(v1))
         depth1 = bamdepth_dict[dequote(k1).split("_mem.bam")[0]]
         if depth1 > 0 and depth1 <= 0.15:
             hg_data[0].append(v1)
         elif depth1 > 0.15 and depth1 <= 0.3:
             hg_data[1].append(v1)
         elif depth1 > 0.3 and depth1 <= 0.45:
             hg_data[2].append(v1)
         elif depth1 > 0.45 and depth1 <= 0.6:
             hg_data[3].append(v1)
         elif depth1 > 0.6 and depth1 <= 0.75:
             hg_data[4].append(v1)
         elif depth1 > 0.75:
             hg_data[5].append(v1)

     for k2,v2 in haplocart_dict.items():
         v2 = is_correct(int(v2))
         depth2 = bamdepth_dict[dequote(k2).split("_mem.bam")[0]]
         if depth2 > 0 and depth2 <= 0.15:
            hc_data[0].append(v2)
         elif depth2 > 0.1 and depth2 <= 0.3:
            hc_data[1].append(v2)
         elif depth2 > 0.2 and depth2 <= 0.45:
            hc_data[2].append(v2)
         elif depth2 > 0.3 and depth2 <= 0.6:
            hc_data[3].append(v2)
         elif depth2 > 0.4 and depth2 <= 0.75:
            hc_data[4].append(v2)
         elif depth2 > 0.75:
            hc_data[5].append(v2)

     for k3,v3 in haplogrouper_dict.items():
        v3 = is_correct(v3)
        depth3 = bamdepth_dict[k3.replace(".hg.out", ".bam")]
        if depth3 > 0 and depth3 <= 0.15:
           hgrp_data[0].append(v3)
        elif depth3 > 0.1 and depth3 <= 0.3:
           hgrp_data[1].append(v3)
        elif depth3 > 0.2 and depth3 <= 0.45:
           hgrp_data[2].append(v3)
        elif depth3 > 0.3 and depth3 <= 0.6:
           hgrp_data[3].append(v3)
        elif depth3 > 0.4 and depth3 <= 0.75:
           hgrp_data[4].append(v3)
        elif depth3 > 0.75:
           hgrp_data[5].append(v3)

     total_per_bin = [1056, 1415, 1376, 1388, 1235, 1530]

     nans = [float('nan'), float('nan')]

     hg_correct_count = [sum(_) for _ in hg_data]
     hc_correct_count = [sum(_) for _ in hc_data]
     hgrp_correct_count = [sum(_) for _ in hgrp_data]

     hg_call_rate = [(len(x) / total_per_bin[i]) for i, x in enumerate(hg_data)]
     hc_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hc_data)]
     hgrp_call_rate = [(len(y) / total_per_bin[j]) for j, y in enumerate(hgrp_data)]

     print("hc call rate: ", hc_call_rate)

     hg_pos = [1,2,3,4,5,6]
     hg_callrate_pos = [x for x in hg_pos]
     hc_pos = [x-0.25 for x in hg_pos]
     hc_callrate_pos = [x-0.25 for x in hg_pos]
     hgrp_pos = [x+0.25 for x in hg_pos]
     hgrp_callrate_pos = [x+0.25 for x in hg_pos]

     fig, ax = plt.subplots(2, figsize=(15, 10), sharex=True)
     fig.suptitle("Total Number of Correct Predictions  \n (Empirical Paired-end FASTQ)")
     plt.xlabel("Mean depth of coverage (X)")
     ax[0].set_ylabel("Count")
     ax[1].set_ylabel("Call rate")

     hg_bar = ax[0].bar(hg_pos, hg_correct_count, color=colorblind_colors[1], label="HaploGrep2", width=0.25)
     hc_bar = ax[0].bar(hc_pos, hc_correct_count, color=colorblind_colors[0], label="HaploCart", width=0.25)
     hgrp_bar = ax[0].bar(hgrp_pos, hgrp_correct_count, color=colorblind_colors[2], label="HaploGrouper", width=0.25)
     hg_callrate_bar = ax[1].bar(hg_callrate_pos, hg_call_rate, color=colorblind_colors[1], label="HaploGrep2 call rate", width=0.25)
     hc_callrate_bar = ax[1].bar(hc_callrate_pos, hc_call_rate, color=colorblind_colors[0], label="HaploCart call rate", width=0.25)
     hgrp_callrate_bar = ax[1].bar(hgrp_callrate_pos, hgrp_call_rate, color=colorblind_colors[2], label="HaploGrouper call rate", width=0.25)

     plt.xticks([1,2,3,4,5,6], \
               ["0-0.15", "0.15-0.3", "0.3-0.45", "0.45-0.6", "0.6-0.75", ">0.75"], rotation=45)

     hc_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart count correct')
     hg_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2 count correct')
     hgrp_patch = mpatches.Patch(color=colorblind_colors[2], label='HaploGrouper count correct')

     hc_callrate_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart call rate')
     hg_callrate_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2 call rate')
     hgrp_callrate_patch = mpatches.Patch(color=colorblind_colors[2], label='HaploGrouper call rate')

     ax[0].legend(handles=[hc_patch, hg_patch, hgrp_patch], loc="upper left")
     ax[1].legend(handles=[hc_callrate_patch, hg_callrate_patch, hgrp_callrate_patch], loc="upper left")
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
    haplogrouper_score_dict = pickle.load(open("../data/pickles/haplogrouper_empirical.pk", "rb"))

    make_correctness_plot_bam(haplogrep_score_dict, haplocart_score_dict, haplogrouper_score_dict, "../data/1kdepths.txt", "../data/pngs/bam_correctness.png")


def plot_single_fastq():
    haplogrep_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    haplocart_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_no_numt.pk", "rb"))
    make_correctness_plot_fq(haplogrep_score_dict, haplocart_score_dict, "../data/fastq_no_numt_sim_depths.txt", "../data/pngs/fastq_correctness.png")

plot_bam()
#plot_single_fastq()
