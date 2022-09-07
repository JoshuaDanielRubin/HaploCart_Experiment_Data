import math
import subprocess
from matplotlib import pyplot as plt
import pickle
import matplotlib.patches as mpatches

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

colorblind_colors = ['#ff7f00', '#377eb8', '#4daf4a']

def safe_log(score):
    if score == 0:
        return 0
    else:
        return math.log(score)

def define_box_properties(plot_name, color_code, label):
    for k, v in plot_name.items():
        plt.setp(plot_name.get(k), color=color_code)

    # use plot function to draw a small line to name the legend.
    plt.plot([], c=color_code, label=label)
    plt.legend()

haplocart_result_path = "../data/haplocart_results/mask.txt"
haplogrep_results_path = "../src/simulations/thousand_genomes/haplocheck_results/mask/"
phymer_mask_path = "../data/phymer_mask/"

with open("../data/pickles/haplocart_mask.pk", "rb") as f:
    haplocart_results = pickle.load(f)

with open("../data/pickles/haplogrep_mask.pk", "rb") as g:
    haplogrep_results = pickle.load(g)

with open("../data/pickles/phymer_mask.pk", "rb") as h:
    phymer_results = pickle.load(h)

hg_0, hg_1000, hg_2000, hg_3000, hg_4000, hg_5000, hg_6000, hg_7000, hg_8000, hg_9000, hg_10000, hg_11000, hg_12000, hg_13000, hg_14000, hg_15000, hg_16000 \
 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

hc_0, hc_1000, hc_2000, hc_3000, hc_4000, hc_5000, hc_6000, hc_7000, hc_8000, hc_9000, hc_10000, hc_11000, hc_12000, hc_13000, hc_14000, hc_15000, hc_16000 \
 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

pm_0, pm_1000, pm_2000, pm_3000, pm_4000, pm_5000, pm_6000, pm_7000, pm_8000, pm_9000, pm_10000, pm_11000, pm_12000, pm_13000, pm_14000, pm_15000, pm_16000 \
 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

hg_scores = [hg_0, hg_1000, hg_2000, hg_3000, hg_4000, hg_5000, hg_6000, hg_7000, hg_8000, hg_9000, hg_10000, hg_11000, hg_12000, \
             hg_13000, hg_14000, hg_15000, hg_16000]

hc_scores = [hc_0, hc_1000, hc_2000, hc_3000, hc_4000, hc_5000, hc_6000, hc_7000, hc_8000, hc_9000, hc_10000, hc_11000, hc_12000, \
             hc_13000, hc_14000, hc_15000, hc_16000]

pm_scores = [pm_0, pm_1000, pm_2000, pm_3000, pm_4000, pm_5000, pm_6000, pm_7000, pm_8000, pm_9000, pm_10000, pm_11000, pm_12000, \
             pm_13000, pm_14000, pm_15000, pm_16000]

for sample, score in haplocart_results.items():
    N = int(sample.split("mask")[-1].split("_")[0])
    hc_scores[int(N/1000)].append(safe_log(int(score)))
    if int(N)<=6000 and int(score) > 30:
        print("sample: ", sample)
        print("N: ", N)
        print("Score: ", score)

for sample, score in haplogrep_results.items():
    N = int(sample.split("mask")[-1].split("_")[0])
    hg_scores[int(N/1000)].append(safe_log(int(score)))

for sample, score in phymer_results.items():
    if score == "NAN":
        continue
    N = int(sample.split("mask")[-1].split("_")[0])
    pm_scores[int(N/1000)].append(safe_log(int(score)))

plt.figure(figsize=(15, 10))
plt.suptitle("Robustness to Missing Bases (Empirical FASTA)")
plt.ylabel("Log Levenshtein distance \n between true and predicted")
plt.xlabel("Number of contiguous missing bases")

pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33]
nans = [float('nan'), float('nan')]

hg_vplot = plt.violinplot([val or nans for val in hg_scores], positions = pos, showmeans=True)
hc_vplot = plt.violinplot([val or nans for val in hc_scores], positions = [x - 0.5 for x in pos], showmeans=True)
pm_vplot = plt.violinplot([val or nans for val in pm_scores], positions = [x + 0.5 for x in pos], showmeans=True)

plt.setp(hc_vplot['bodies'], facecolor=colorblind_colors[0], edgecolor=colorblind_colors[0])
plt.setp(hg_vplot['bodies'], facecolor=colorblind_colors[1], edgecolor=colorblind_colors[1])
plt.setp(pm_vplot['bodies'], facecolor=colorblind_colors[2], edgecolor=colorblind_colors[2])

plt.xticks(pos, [str(x) for x in range(0, 16001, 1000)], rotation=45)

hc_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart')
hg_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2')
pm_patch = mpatches.Patch(color=colorblind_colors[2], label='Phy-Mer')

plt.legend(handles=[hc_patch, hg_patch, pm_patch], loc="upper left")
plt.xlim(2, 34)
plt.tight_layout()
plt.savefig("../data/pngs/mask.png", dpi=300)
plt.close()



