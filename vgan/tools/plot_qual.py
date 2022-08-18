
import math
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pickle
from matplotlib import pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small'}
pylab.rcParams.update(params)
hg_patch = mpatches.Patch(color='blue', label='HaploGrep2')
hc_patch = mpatches.Patch(color='orange', label='HaploCart')

thousand_genomes_ids = ["NA19661", "HG00473","HG01051","HG02666","HG03112","NA18510","NA19036","NA20518"]
fastq_ids = ["H2a2a1", "L0a1a1", "Q1", "C1b"]

def transform_hc(k):
    splitted = k.split("_")
    return "_".join([splitted[0], 'l'+splitted[1], 'rep'+splitted[2], splitted[3]])

def transform_hg(k):
    splitted = k.split("_")
    return "_".join([splitted[0][1:], splitted[2], splitted[3], splitted[4]])

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    s = str(s)
    if '"' not in s and "'" not in s:
        return s
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    elif s.startswith(("'", '"')):
        return s[1:]


def plot_bam():
    hg_xs, hg_ys, hg_labels = [], [], []
    hc_xs, hc_ys, hc_labels = [], [], []

    with open("../data/pickles/haplogrep_reported_quals_bam.pk", "rb") as f:
        hg_qual_dict = pickle.load(f)

        with open("../data/pickles/haplogrep_bam.pk", "rb") as g:
            hg_score_dict = pickle.load(g)

        hc_bam_merged = pickle.load(open("../data/pickles/bam_merged.pk", "rb"))
        hc_score_dict = pickle.load(open("../data/pickles/haplocart_bams.pk", "rb"))

        for sample in hg_score_dict.keys():
            score = hg_score_dict[sample]
            reported_qual = hg_qual_dict[sample]
            id = sample.split("_")[0]

            hg_ys.append(float(dequote(reported_qual)))
            hg_xs.append(int(dequote(score)))
            hg_labels.append(id)

        for idx, row in hc_bam_merged.iterrows():
            if row['tree_depth'] > 0:
                continue
            else:
               hc_xs.append(int(hc_score_dict[row['sample']]))
               hc_ys.append(float(row['posterior']))
               hc_labels.append(row['sample'])

        fig,ax = plt.subplots()
        for i, k in enumerate(thousand_genomes_ids):
            hg_idxs = [i for i in range(len(hg_labels)) if dequote(hg_labels[i]).startswith(k)]
            hc_idxs = [i for i in range(len(hc_labels)) if dequote(hc_labels[i]).startswith(k)]
            ax.scatter([hg_xs[j] for j in hg_idxs], [hg_ys[j] for j in hg_idxs], s=8, alpha=0.2, c='blue', label='HaploGrep2' if i==0 else "")
            ax.scatter([hc_xs[j] for j in hc_idxs], [hc_ys[j] for j in hc_idxs], s=8, alpha=0.2, c='orange', label='HaploCart' if i==0 else "")
        plt.suptitle("Reported Confidence on \n Thousand Genomes Project Samples")
        ax.set_xlabel("Levenshtein Distance between true and predicted")
        ax.set_ylabel("Reported Confidence")
        plt.legend(handles=[hg_patch, hc_patch], borderpad=0.45, prop={'size': 8})
        plt.tight_layout()
        plt.xlim(-5, 110)
        plt.savefig("../data/pngs/reported_qual.png", dpi=300)
        plt.close()


def plot_fastq_no_bias():
    depthfile="../data/fastq_no_numt_sim_depths.txt"
    fqdepth_dict = {}
    with open(depthfile, "r") as f:
        for line in f:
            split=line.split()
            file = split[0].split("/")[-1].split("_mem.bam")[0].split("_")
            file, depth = "_".join([file[0], file[2], file[3], file[4]]), float(split[1])
            fqdepth_dict.update({file:depth})


    hg_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    hg_qual_dict = pickle.load(open("../data/pickles/haplogrep_reported_quals_fastq.pk", "rb"))
    hc_fastq_merged = pickle.load(open("../data/pickles/fastq_merged.pk", "rb"))


    hc_data_dict, hg_data_dict = {}, {}
    trend_dict = {}
    for fastq_id in fastq_ids:
        hc_data_dict[fastq_id] = [[] for i in range(10)]
        hg_data_dict[fastq_id] = [[] for i in range(10)]
        trend_dict[fastq_id] = [[], [], [], []]

    for sample in hg_score_dict.keys():
        reported_qual = float(dequote(hg_qual_dict[sample]))
        id = dequote(sample.split("_")[0])
        if id not in fastq_ids:
            continue
        depth = float(fqdepth_dict[transform_hg(sample)])

        trend_dict[id][0].append(depth * 5)
        trend_dict[id][1].append(reported_qual)

        if (depth >= 0 and depth < 0.2):
            hg_data_dict[id][0].append(reported_qual)
        elif (depth >= 0.1 and depth < 0.4):
            hg_data_dict[id][1].append(reported_qual)
        elif (depth >= 0.2 and depth < 0.6):
            hg_data_dict[id][2].append(reported_qual)
        elif (depth >= 0.3 and depth < 0.8):
            hg_data_dict[id][3].append(reported_qual)
        elif (depth >= 0.4 and depth < 1):
            hg_data_dict[id][4].append(reported_qual)
        elif (depth >= 0.5 and depth < 1.2):
            hg_data_dict[id][5].append(reported_qual)
        elif (depth >= 0.6 and depth < 1.4):
            hg_data_dict[id][6].append(reported_qual)
        elif (depth >= 0.7 and depth < 1.6):
            hg_data_dict[id][7].append(reported_qual)
        elif (depth >= 0.8 and depth < 1.8):
            hg_data_dict[id][8].append(reported_qual)
        elif (depth >= 0.9 and depth < 2):
            hg_data_dict[id][9].append(reported_qual)

    for idx, row in hc_fastq_merged.iterrows():
        if int(row['tree_depth']) > 0:
            continue
        depth = float(fqdepth_dict[transform_hc(row['sample'])])
        reported_qual = float(row['posterior'])

        if math.isnan(reported_qual):
            print(row['sample'])
            continue

        id = row['sample'].split("_")[0]
        if id not in fastq_ids:
            continue

        trend_dict[id][2].append(depth * 5)
        trend_dict[id][3].append(reported_qual)

        if (depth >= 0 and depth < 0.2):
            hc_data_dict[id][0].append(reported_qual)
        elif (depth >= 0.1 and depth < 0.4):
            hc_data_dict[id][1].append(reported_qual)
        elif (depth >= 0.2 and depth < 0.6):
            hc_data_dict[id][2].append(reported_qual)
        elif (depth >= 0.3 and depth < 0.8):
            hc_data_dict[id][3].append(reported_qual)
        elif (depth >= 0.4 and depth < 1):
            hc_data_dict[id][4].append(reported_qual)
        elif (depth >= 0.5 and depth < 1.2):
            hc_data_dict[id][5].append(reported_qual)
        elif (depth >= 0.6 and depth < 1.4):
            hc_data_dict[id][6].append(reported_qual)
        elif (depth >= 0.7 and depth < 1.6):
            hc_data_dict[id][7].append(reported_qual)
        elif (depth >= 0.8 and depth < 1.8):
            hc_data_dict[id][8].append(reported_qual)
        elif (depth >= 0.9 and depth < 2):
            hc_data_dict[id][9].append(reported_qual)

    fig,ax = plt.subplots(2, 2, figsize=(12, 12))
    ax_lst = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]

    for i, _ in enumerate(ax_lst):
        _.set_ylabel("Reported Confidence")
        _.set_xlabel("Coverage Depth (X)")
        _.set_title(fastq_ids[i])

    nans = [float('nan'), float('nan')]
    hc_pos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    hg_pos = [x - 0.25 for x in hc_pos]
    xtick_labels = ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1", "1-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2"]

    for i, fastq_id in enumerate(fastq_ids):

        hg_to_plot = [val or nans for val in hg_data_dict[fastq_id]]
        hc_to_plot = [val or nans for val in hc_data_dict[fastq_id]]

        ax_lst[i].violinplot(hg_to_plot, positions=hg_pos)
        ax_lst[i].violinplot(hc_to_plot, positions=hc_pos)
        ax_lst[i].set_xticks([x+0.25 for x in hg_pos])
        ax_lst[i].set_xticklabels(xtick_labels, rotation=45)

        hgp = np.poly1d(np.polyfit(trend_dict[fastq_id][0], trend_dict[fastq_id][1], 3))
        hgt = np.linspace(0, 10, 200)
        hcp = np.poly1d(np.polyfit(trend_dict[fastq_id][2], trend_dict[fastq_id][3], 3))
        hct = np.linspace(0, 10, 200)

        ax_lst[i].plot(hgt, hgp(hgt), ':')
        ax_lst[i].plot(hct, hcp(hct), ':')

    hg_patch = mpatches.Patch(color='blue', label='HaploGrep2 posterior distribution')
    hc_patch = mpatches.Patch(color='orange', label='HaploCart posterior distribution')
    hg_trend_patch = mlines.Line2D([0], [0], color='green', label='HaploGrep2 polynomial regression', linestyle=':', linewidth=1)
    hc_trend_patch = mlines.Line2D([0], [0], color='red', label='HaploCart polynomial regression', linestyle=':', linewidth=1)

    fig.legend(handles=[hg_patch, hc_patch, hg_trend_patch, hc_trend_patch], borderpad=0.45, prop={'size': 6.5})

    plt.suptitle("Coverage Depth vs. Reported Confidence")
    plt.tight_layout()
    plt.savefig("../data/pngs/reported_qual_fastq.png")
    plt.close()

#plot_bam()
plot_fastq_no_bias()

