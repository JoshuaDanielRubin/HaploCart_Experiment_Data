import math
import matplotlib.patches as mpatches
import subprocess
from matplotlib import pyplot as plt
import pickle
import pdb

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

def safe_log(score):
    if score == 0:
        return 0
    else:
        return math.log(score)

colorblind_colors = ['#ff7f00', '#377eb8']

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

def make_violinplot_bam(haplogrep_dict, haplocart_dict, depthfile):
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

    hg_data = [hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_10, hg_10_up]
    hc_data = [hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_10, hc_10_up]

    pos = [1,3,5,7,9,11,13,15,17,19,21]

    plt.figure(figsize=(15, 5))
    hg_vplot = plt.violinplot(hg_data, positions = pos, showmeans=True)
    hc_vplot = plt.violinplot(hc_data, positions = [x - 0.5 for x in pos], showmeans=True)

    plt.suptitle("Haplogroup Classification Accuracy \n (Empirical Paired-end FASTQ)")
    plt.xlabel("Mean depth of coverage (X)")
    plt.ylabel("Levenshtein distance \n between true and predicted")

    plt.xticks([1,3,5,7,9,11,13,15,17,19,21], \
               ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", \
                "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0", ">=1.0"], rotation=45)

    hg_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2')
    hc_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart')

    plt.legend(handles=[hc_patch, hg_patch], borderpad=0.5, prop={'size': 12})
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.savefig("../data/pngs/bams_downsampled.png", bbox_inches='tight', dpi=300)
    plt.close()


def transform_hg(k, numt):
    splitted = k.split("_")
    if numt == False:
        return "_".join([splitted[0][1:], splitted[2], splitted[3], splitted[4]])
    else:
        return "_".join([splitted[2], splitted[4], splitted[5], splitted[7]])

def transform_hct(k, numt):
    splitted = k.split("_")
    if numt == False:
        return "_".join([splitted[0], 'l'+splitted[1], 'rep'+splitted[2], splitted[3]])
    else:
        return "_".join([splitted[0], splitted[1], splitted[2], splitted[3], splitted[4]])

def make_violinplot_fq(haplogrep_dict, haplocart_dict, depthfile, outfile, numt=False):
    print(len(haplocart_dict.keys()))
    haplogrep_dict = {transform_hg(k, numt):v for k,v in haplogrep_dict.items()}
    haplocart_dict = {transform_hct(k, numt):v for k,v in haplocart_dict.items()}

    fqdepth_dict = {}
    with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file = split[0].split("/")[-1].split("_mem.bam")[0].split("_")
            if numt == False:
                file, depth = "_".join([file[0], file[2], file[3], file[4]]), float(split[1])
            else:
                file, depth = "_".join([file[2], file[4], file[5], file[7]]), float(split[1])
            fqdepth_dict.update({file:depth})

    hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20 \
    = [], [], [], [], [], [], [], [], [], []
    hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20 \
    = [], [], [], [], [], [], [], [], [], []

    for k1,v1 in haplogrep_dict.items():
        splitted = dequote(k1).split("_")
        if "H2a2a1" == splitted[0]:
            continue
        v1 = int(v1)
        depth1 = fqdepth_dict[k1]
        if depth1 <= 0.2:
            hg_0_02.append(v1)
        elif depth1 > 0.2 and depth1 <= 0.4:
            hg_02_04.append(v1)
        elif depth1 > 0.4 and depth1 <= 0.6:
            hg_04_06.append(v1)
        elif depth1 > 0.6 and depth1 <= 0.8:
            hg_06_08.append(v1)
        elif depth1 > 0.8 and depth1 <= 1.0:
           hg_08_10.append(v1)
        elif depth1 > 1.0 and depth1 <= 1.2:
            hg_10_12.append(v1)
        elif depth1 > 1.2 and depth1 <= 1.4:
            hg_12_14.append(v1)
        elif depth1 > 1.4 and depth1 <= 1.6:
            hg_14_16.append(v1)
        elif depth1 > 1.6 and depth1 <= 1.8:
            hg_16_18.append(v1)
        elif depth1 > 1.8 and depth1 <= 2:
            hg_18_20.append(v1)


    for k2,v2 in haplocart_dict.items():
        v2 = int(v2)
        splitted = dequote(k2).split("_")
        if "H2a2a1" == splitted[0]:
            continue
        if numt:
            depth2 = fqdepth_dict["_".join([splitted[0], splitted[1], splitted[2], splitted[4]])]
        else:
            depth2 = fqdepth_dict["_".join([splitted[0], splitted[1], splitted[2], splitted[3]])]
        if depth2 <= 0.2:
            hc_0_02.append(v2)
        elif depth2 > 0.2 and depth2 <= 0.4:
            hc_02_04.append(v2)
        elif depth2 > 0.4 and depth2 <= 0.6:
            hc_04_06.append(v2)
        elif depth2 > 0.6 and depth2 <= 0.8:
            hc_06_08.append(v2)
        elif depth2 > 0.8 and depth2 <= 1.0:
            hc_08_10.append(v2)
        elif depth2 > 1.0 and depth2 <= 1.2:
            hc_10_12.append(v2)
        elif depth2 > 1.2 and depth2 <= 1.4:
            hc_12_14.append(v2)
        elif depth2 > 1.4 and depth2 <= 1.6:
            hc_14_16.append(v2)
        elif depth2 > 1.6 and depth2 <= 1.8:
            hc_16_18.append(v2)
        elif depth2 > 1.8 and depth2 <= 2:
            if safe_log(v2) > 3:
                print(k2)
            hc_18_20.append(v2)

    hc_data = [hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20]
    hg_data = [hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20]

    if numt == False:
        plt.title("Haplogroup Classification Performance \n (Simulated Paired-end FASTQ)")
    else:
        plt.title("Haplogroup Classification Performance \n (Simulated Paired-end FASTQ, With NuMTs)")
    plt.xlabel("Mean depth of coverage (X)")
    plt.ylabel("Log Levenshtein distance \n between true and predicted")
    #plt.ylim(-5, 115)
    nans = [float('nan'), float('nan')]

    hg_pos = [1,3,5,7,9,11,13,15,17,19]
    hc_pos = [x-0.5 for x in hg_pos]

    hg_vplot = plt.violinplot([[safe_log(x) for x in val] or nans for val in hg_data], positions=hg_pos, showmeans=True)
    hc_vplot = plt.violinplot([[safe_log(x) for x in val] or nans for val in hc_data], positions=hc_pos, showmeans=True)

    hg_patch = mpatches.Patch(color="blue", label='HaploCart accuracy distribution')
    hc_patch = mpatches.Patch(color="orange", label='HaploGrep2 accuracy distribution')

    plt.xticks([1,3,5,7,9,11,13,15,17,19], \
               ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0"], rotation=45)

    #for idx, i in enumerate([0,2,4,6,8,10,12,14,16,18]):
    #    plt.text(i, 100, "N="+str(len(hc_data[idx])), rotation=45, fontsize=6, c="orange")
    #    plt.text(i+0.6, 100, "N="+str(len(hg_data[idx])), rotation=45, fontsize=6, c="blue")
    #for idx2, j in enumerate([20,22,24,26,28,30]):
    #    plt.text(j, 60, "N="+str(len(hc_data[idx2+10])), rotation=45, fontsize=6, c="orange")
    #    plt.text(j+0.6, 60, "N="+str(len(hg_data[idx2+10])), rotation=45, fontsize=6, c="blue")

    plt.legend(handles=[hc_patch, hg_patch], borderpad=0.45, prop={'size': 8})
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

    make_violinplot_bam(haplogrep_score_dict, haplocart_score_dict, "../data/1kdepths.txt")


def plot_fastq():

    haplogrep_no_numt_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    haplocart_no_numt_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_no_numt.pk", "rb"))
    haplogrep_with_numt_score_dict = pickle.load(open("../data/pickles/hg_fastq_with_numt.pk", "rb"))
    haplocart_with_numt_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_with_numt.pk", "rb"))

    make_violinplot_fq(haplogrep_no_numt_score_dict, haplocart_no_numt_score_dict, \
                       "../data/fastq_no_numt_sim_depths.txt", "../data/pngs/fastq_no_numt.png", numt=False)

    make_violinplot_fq(haplogrep_with_numt_score_dict, haplocart_with_numt_score_dict, \
                      "../data/fastq_with_numt_sim_depths.txt", "../data/pngs/fastq_with_numt.png", numt=True)


#plot_bam()
plot_fastq()
