import math
import pandas as pd
import matplotlib.patches as mpatches
import subprocess
from matplotlib import pyplot as plt
import pickle
import pdb
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
         'axes.labelsize': 'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

colorblind_colors = ['#ff7f00', '#377eb8', '#009e75']

def safe_log(score):
    if score == 0:
        return 0
    else:
        return math.log(score)

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

def make_violinplot_fq(ax, haplogrep_dict, haplocart_dict, haplogrouper_dict, depthfile, threshold, idx):

    with open("../data/haplocart_results/no_numt_posterior.txt", "rt") as g:
        no_numt_posterior_df = pd.read_csv(g, sep='\t', header=None, names=["Sample", "Pred", "Posterior", "Tree_depth"])

    print("original size: ", no_numt_posterior_df.shape)
    depth_0 = no_numt_posterior_df.loc[no_numt_posterior_df['Tree_depth'] == 0]
    print("depth 0 size: ", depth_0.shape)
    passed_df = depth_0.loc[depth_0['Posterior'] >= threshold]
    passed_samples = list(passed_df['Sample'].unique())
    print("Number passed samples: ", len(passed_samples))
    passed_samples = [transform_hct(x, False) for x in passed_samples]

    haplogrep_dict = {transform_hg(k, False):v for k,v in haplogrep_dict.items()}
    haplocart_dict = {transform_hct(k, False):v for k,v in haplocart_dict.items()}

    fqdepth_dict = {}
    with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file = split[0].split("/")[-1].split("_mem.bam")[0].split("_")
            file, depth = "_".join([file[0], file[2], file[3], file[4]]), float(split[1])
            fqdepth_dict.update({file:depth})

    hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20 \
    = [], [], [], [], [], [], [], [], [], []
    hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20 \
    = [], [], [], [], [], [], [], [], [], []
    hgrp_0_02, hgrp_02_04, hgrp_04_06, hgrp_06_08, hgrp_08_10, hgrp_10_12, hgrp_12_14, hgrp_14_16, hgrp_16_18, hgrp_18_20 \
    = [], [], [], [], [], [], [], [], [], []

    for k1,v1 in haplogrep_dict.items():
        splitted = dequote(k1).split("_")
        if "H2a2a1" == splitted[0] or "L1c2b" == splitted[0]:
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
        if "H2a2a1" == splitted[0] or "L1c2b" == splitted[0]:
            continue

        depth2 = fqdepth_dict["_".join([splitted[0], splitted[1], splitted[2], splitted[3]])]

        if "_".join([splitted[0], splitted[1], splitted[2], splitted[3]]) not in passed_samples:
            continue

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
            hc_18_20.append(v2)

    for k3,v3 in haplogrouper_dict.items():
        v3 = int(v3)
        splitted = k3.split("_")

        if "H2a2a1" == splitted[0] or "L1c2b" == splitted[0]:
            continue
        depth3 = fqdepth_dict["_".join([splitted[0], splitted[2], splitted[3], splitted[4]])]
        if depth3 <= 0.2:
            hgrp_0_02.append(v3)
        elif depth3 > 0.2 and depth3 <= 0.4:
            hgrp_02_04.append(v3)
        elif depth3 > 0.4 and depth3 <= 0.6:
            hgrp_04_06.append(v3)
        elif depth3 > 0.6 and depth3 <= 0.8:
            hgrp_06_08.append(v3)
        elif depth3 > 0.8 and depth3 <= 1.0:
            hgrp_08_10.append(v3)
        elif depth3 > 1.0 and depth3 <= 1.2:
            hgrp_10_12.append(v3)
        elif depth3 > 1.2 and depth3 <= 1.4:
            hgrp_12_14.append(v3)
        elif depth3 > 1.4 and depth3 <= 1.6:
            hgrp_14_16.append(v3)
        elif depth3 > 1.6 and depth3 <= 1.8:
            hgrp_16_18.append(v3)
        elif depth3 > 1.8 and depth3 <= 2:
            if safe_log(v3) > 3:
                print(k3)
            hgrp_18_20.append(v3)

    hc_data = [hc_0_02, hc_02_04, hc_04_06, hc_06_08, hc_08_10, hc_10_12, hc_12_14, hc_14_16, hc_16_18, hc_18_20]
    hg_data = [hg_0_02, hg_02_04, hg_04_06, hg_06_08, hg_08_10, hg_10_12, hg_12_14, hg_14_16, hg_16_18, hg_18_20]
    hgrp_data = [hgrp_0_02, hgrp_02_04, hgrp_04_06, hgrp_06_08, hgrp_08_10, hgrp_10_12, hgrp_12_14, hgrp_14_16, hgrp_16_18, hgrp_18_20]

    plt.suptitle("Haplogroup Classification Performance \n (Simulated Paired-end FASTQ)", fontsize=22)
    ax.set_xlabel("Mean depth of coverage (X)")
    ax.set_ylabel("Log Levenshtein distance \nbetween true and predicted")
    nans = [float('nan'), float('nan')]

    hg_pos = [1,3,5,7,9,11,13,15,17,19]
    hc_pos = [x-0.5 for x in hg_pos]
    hgrp_pos = [x+0.5 for x in hg_pos]

    hg_vplot = ax.violinplot([[safe_log(x) for x in val] or nans for val in hg_data], positions=hg_pos, showmeans=True)
    hc_vplot = ax.violinplot([[safe_log(x) for x in val] or nans for val in hc_data], positions=hc_pos, showmeans=True)
    hgrp_vplot = ax.violinplot([[safe_log(x) for x in val] or nans for val in hgrp_data], positions=hgrp_pos, showmeans=True)

    hg_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart accuracy distribution')
    hc_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2 accuracy distribution')
    hgrp_patch = mpatches.Patch(color=colorblind_colors[2], label='HaploGrouper accuracy distribution')

    ax.set_xticklabels(["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0"], rotation=45)
    ax.set_title("Posterior threshold = " + str(threshold))
    ax.legend(handles=[hc_patch, hg_patch, hgrp_patch], borderpad=0.4, prop={'size': 12})
    return ax

def plot_fastq():
    depth_dict = {}
    with open("../data/fastq_no_numt_sim_depths.txt", "r") as f:
        for line in f:
            split=line.split()
            file, depth = dequote(split[0].split("/")[-1]), float(split[1])
            depth_dict.update({file:depth})

    haplogrep_no_numt_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    haplocart_no_numt_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_no_numt.pk", "rb"))
    haplogrouper_dict = pickle.load(open("../data/pickles/haplogrouper_sim_no_numt.pk", "rb"))

    fig, ax = plt.subplots(3, 1, figsize=(30, 20))
    ax_lst = [ax[0], ax[1], ax[2]]
    for idx, thresh in enumerate([0, 0.3, 0.99]):
        make_violinplot_fq(ax_lst[idx], haplogrep_no_numt_score_dict, haplocart_no_numt_score_dict, haplogrouper_dict, \
                           "../data/fastq_no_numt_sim_depths.txt", thresh, idx)

    plt.setp(ax_lst, xticks=[1,3,5,7,9,11,13,15,17,19], \
                     xticklabels=["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0"])
    plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.92,
                    top=0.92,
                    wspace=0.2,
                    hspace=0.3)
    #plt.tight_layout()
    plt.savefig("../data/pngs/no_numt_threshold.png", dpi=300)

#plot_bam()
plot_fastq()
