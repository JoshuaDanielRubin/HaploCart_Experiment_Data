import subprocess
from matplotlib import pyplot as plt
import pickle
import matplotlib.patches as mpatches

colorblind_colors = ['#ff7f00', '#377eb8', '#4daf4a']

def is_correct(score):
    if str(score) == "0":
        return 1
    else:
        return 0

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
    hc_scores[int(N/1000)].append(is_correct(int(score)))

for sample, score in haplogrep_results.items():
    id = sample.split("mask")[0].split("_")[0]
    if id == "H2a2a1":
        continue
    N = int(sample.split("mask")[-1].split("_")[0])
    hg_scores[int(N/1000)].append(is_correct(int(score)))

for sample, score in phymer_results.items():
    N = int(sample.split("mask")[-1].split("_")[0])
    pm_scores[int(N/1000)].append(is_correct(int(score)))

plt.figure(figsize=(10, 8))
plt.suptitle("Robustness to Missing Bases (Empirical FASTA)")
plt.ylabel("Total number exactly correct")
plt.xlabel("Number of contiguous missing bases")

hg_correct_count = [(sum(_)) for _ in hg_scores]
hc_correct_count = [(sum(_)) for _ in hc_scores]
pm_correct_count = [(sum(_)) for _ in pm_scores]


hc_pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33]
hg_pos = [x+0.5 for x in hc_pos]
pm_pos = [x+1 for x in hc_pos]

hc_bar = plt.bar(hc_pos, hc_correct_count, color=colorblind_colors[0], label="HaploCart call rate", width=0.5)
hg_bar = plt.bar(hg_pos, hg_correct_count, color=colorblind_colors[1], label="HaploGrep2 call rate", width=0.5)
hg_bar = plt.bar(pm_pos, pm_correct_count, color=colorblind_colors[2], label="HaploGrep2 call rate", width=0.5)

plt.xticks(hc_pos, [str(x) for x in range(0, 16001, 1000)], rotation=45)
plt.xlim(2.5, 34.5)

hc_patch = mpatches.Patch(color=colorblind_colors[0], label='HaploCart')
hg_patch = mpatches.Patch(color=colorblind_colors[1], label='HaploGrep2')
pm_patch = mpatches.Patch(color=colorblind_colors[2], label='Phy-Mer')

plt.legend(handles=[hc_patch, hg_patch, pm_patch], loc="upper left")
plt.tight_layout()
plt.savefig("../data/pngs/mask_correctness.png", dpi=300)
plt.close()



