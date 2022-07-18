import subprocess
from matplotlib import pyplot as plt
import pickle
import matplotlib.patches as mpatches

def define_box_properties(plot_name, color_code, label):
    for k, v in plot_name.items():
        plt.setp(plot_name.get(k), color=color_code)

    # use plot function to draw a small line to name the legend.
    plt.plot([], c=color_code, label=label)
    plt.legend()

haplocart_result_path = "../data/haplocart_results/mask.txt"
haplogrep_results_path = "../src/simulations/thousand_genomes/haplocheck_results/mask/"

with open("../data/pickles/haplocart_mask.pk", "rb") as f:
    haplocart_results = pickle.load(f)

with open("../data/pickles/haplogrep_mask.pk", "rb") as g:
    haplogrep_results = pickle.load(g)

hg_0, hg_1000, hg_2000, hg_3000, hg_4000, hg_5000, hg_6000, hg_7000, hg_8000, hg_9000, hg_10000, hg_11000, hg_12000, hg_13000, hg_14000, hg_15000, hg_16000 \
 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

hj_0, hj_1000, hj_2000, hj_3000, hj_4000, hj_5000, hj_6000, hj_7000, hj_8000, hj_9000, hj_10000, hj_11000, hj_12000, hj_13000, hj_14000, hj_15000, hj_16000 \
 = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

hg_scores = [hg_0, hg_1000, hg_2000, hg_3000, hg_4000, hg_5000, hg_6000, hg_7000, hg_8000, hg_9000, hg_10000, hg_11000, hg_12000, \
             hg_13000, hg_14000, hg_15000, hg_16000]

hj_scores = [hj_0, hj_1000, hj_2000, hj_3000, hj_4000, hj_5000, hj_6000, hj_7000, hj_8000, hj_9000, hj_10000, hj_11000, hj_12000, \
             hj_13000, hj_14000, hj_15000, hj_16000]

for sample, score in haplocart_results.items():
    N = int(sample.split("mask")[-1].split("_")[0])
    hj_scores[int(N/1000)].append(int(score))
    if N==0 and int(score) > 0:
        print("sample: ", sample)
        print("N: ", N)
        print("Score: ", score)

for sample, score in haplogrep_results.items():
    N = int(sample.split("mask")[-1].split("_")[0])
    hg_scores[int(N/1000)].append(int(score))

plt.figure(figsize=(10, 8))
plt.suptitle("Robustness to Missing Bases")
plt.ylabel("Phylogenetic distance between true and predicted")
plt.xlabel("Number of missing bases")

pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33]
nans = [float('nan'), float('nan')]

hg_vplot = plt.violinplot([val or nans for val in hg_scores], positions = pos, showmeans=True)
hj_vplot = plt.violinplot([val or nans for val in hj_scores], positions = [x - 0.5 for x in pos], showmeans=True)

plt.setp(hg_vplot['bodies'], facecolor='blue', edgecolor='blue')
plt.setp(hj_vplot['bodies'], facecolor='orange', edgecolor='orange')

plt.xticks(pos, [str(x) for x in range(0, 16001, 1000)], rotation=45)

hg_patch = mpatches.Patch(color='blue', label='HaploGrep2')
hj_patch = mpatches.Patch(color='orange', label='HaploCart')

plt.legend(handles=[hg_patch, hj_patch], loc="upper left")
plt.tight_layout()
plt.savefig("../data/pngs/mask.png", dpi=300)
plt.close()



