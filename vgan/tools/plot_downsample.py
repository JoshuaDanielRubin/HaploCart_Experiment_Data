import matplotlib.patches as mpatches
import subprocess
from matplotlib import pyplot as plt
import pickle
import pdb

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

    hg_005_01, hg_01_015, hg_015_02, hg_02_025, hg_025_03, hg_03_035, hg_035_04, hg_04_045, hg_045_05, hg_05_055, \
    hg_055_06, hg_06_065, hg_065_07, hg_07_075, hg_075_08, hg_08_085, hg_085_09, hg_09_095, hg_095_1, hg_1_105, hg_105_11, hg_11_115, hg_115_12 \
    = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    hc_005_01, hc_01_015, hc_015_02, hc_02_025, hc_025_03, hc_03_035, hc_035_04, hc_04_045, hc_045_05, hc_05_055, \
    hc_055_06, hc_06_065, hc_065_07, hc_07_075, hc_075_08, hc_08_085, hc_085_09, hc_09_095, hc_095_1, hc_1_105, hc_105_11, hc_11_115, hc_115_12 \
    = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    for k1,v1 in haplogrep_dict.items():
        v1 = int(v1)
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
        v2 = int(v2)
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
               hg_055_06, hg_06_065, hg_065_07, hg_07_075, hg_075_08, hg_08_085, hg_085_09, hg_09_095, hg_095_1, hg_1_105, hg_105_11, hg_11_115, hg_115_12]

    hc_data = [hc_005_01, hc_01_015, hc_015_02, hc_02_025, hc_025_03, hc_03_035, hc_035_04, hc_04_045, hc_045_05, hc_05_055, \
               hc_055_06, hc_06_065, hc_065_07, hc_07_075, hc_075_08, hc_08_085, hc_085_09, hc_09_095, hc_095_1, hc_1_105, hc_105_11, hc_11_115, hc_115_12]

    pos = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45]

    plt.figure(figsize=(15, 5))
    hg_vplot = plt.violinplot(hg_data, positions = pos, showmeans=True)
    hc_vplot = plt.violinplot(hc_data, positions = [x - 0.5 for x in pos], showmeans=True)

    plt.suptitle("Haplogroup Classification Performance \n (Thousand Genomes BAM files)")
    plt.xlabel("Mean depth of coverage (x)")
    plt.ylabel("Levenshtein distance between true and predicted")

    plt.xticks([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45], \
               ["0.05-0.1", "0.1-0.15", "0.15-0.2", "0.2-0.25", "0.25-0.3", "0.3-0.35", \
                "0.35-0.4", "0.4-0.45", "0.45-0.5", "0.5-0.55", "0.55-0.6", "0.6-0.65", "0.65-0.7", \
                "0.7-0.75", "0.75-0.8", "0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-1", "1-1.05", "1.05-1.1", "1.1-1.15", "1.15-1.2"], rotation=45)

    hg_patch = mpatches.Patch(color='blue', label='HaploGrep2')
    hc_patch = mpatches.Patch(color='orange', label='HaploCart')

    plt.legend(handles=[hg_patch, hc_patch], borderpad=0.45, prop={'size': 9})
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


    hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, hg_09_1, hg_1_12, hg_12_14, hg_14_16, hg_16_18, hg_18_2, hg_2_22 = \
    [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, hc_09_1, hc_1_12, hc_12_14, hc_14_16, hc_16_18, hc_18_2, hc_2_22  = \
    [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    for k1,v1 in haplogrep_dict.items():
        v1 = int(v1)
        depth1 = fqdepth_dict[k1]
        if depth1 <= 0.1:
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
        elif depth1 > 0.9 and depth1 <= 1:
            hg_09_1.append(v1)
        elif depth1 > 1 and depth1 <= 1.2:
            hg_1_12.append(v1)
        elif depth1 > 1.2 and depth1 <= 1.4:
            hg_12_14.append(v1)
        elif depth1 > 1.4 and depth1 <= 1.6:
            hg_14_16.append(v1)
        elif depth1 > 1.6 and depth1 <= 1.8:
            hg_16_18.append(v1)
        elif depth1 > 1.8 and depth1 <= 2:
            hg_18_2.append(v1)
        elif depth1 > 2 and depth1 <= 2.2:
            hg_2_22.append(v1)

    for k2,v2 in haplocart_dict.items():
        v2 = int(v2)
        splitted = dequote(k2).split("_")
        if numt:
            depth2 = fqdepth_dict["_".join([splitted[0], splitted[1], splitted[2], splitted[4]])]
        else:
            depth2 = fqdepth_dict["_".join([splitted[0], splitted[1], splitted[2], splitted[3]])]
        if depth2 <= 0.1:
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
        elif depth2 > 0.9 and depth2 <= 1:
            hc_09_1.append(v2)
        elif depth2 > 1 and depth2 <= 1.2:
            hc_1_12.append(v2)
        elif depth2 > 1.2 and depth2 <= 1.4:
            hc_12_14.append(v2)
        elif depth2 > 1.4 and depth2 <= 1.6:
            hc_14_16.append(v2)
        elif depth2 > 1.6 and depth2 <= 1.8:
            hc_16_18.append(v2)
        elif depth2 > 1.8 and depth2 <= 2:
            hc_18_2.append(v2)
        elif depth2 > 2 and depth2 <= 2.2:
            hc_2_22.append(v2)

    hg_data = [hg_0_01, hg_01_02, hg_02_03, hg_03_04, hg_04_05, hg_05_06, hg_06_07, hg_07_08, hg_08_09, \
               hg_09_1, hg_1_12, hg_12_14, hg_14_16, hg_16_18, hg_18_2, hg_2_22]

    hc_data = [hc_0_01, hc_01_02, hc_02_03, hc_03_04, hc_04_05, hc_05_06, hc_06_07, hc_07_08, hc_08_09, \
               hc_09_1, hc_1_12, hc_12_14, hc_14_16, hc_16_18, hc_18_2, hc_2_22]


    if numt == False:
        plt.title("Haplogroup Classification Performance \n (Simulated Paired-end FASTQ)")
    else:
        plt.title("Haplogroup Classification Performance \n (Simulated Paired-end FASTQ, With NuMTs)")
    plt.xlabel("Mean depth of coverage (x)")
    plt.ylabel("Levenshtein distance between true and predicted")
    plt.ylim(-5, 115)
    nans = [float('nan'), float('nan')]
    hg_vplot = plt.violinplot([val or nans for val in hg_data], positions=[1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31], showmeans=True)
    hc_vplot = plt.violinplot([val or nans for val in hc_data], positions=[0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,\
                                                                           24.5,26.5,28.5,30.5], showmeans=True)

    hg_patch = mpatches.Patch(color='blue', label='HaploGrep2')
    hc_patch = mpatches.Patch(color='orange', label='HaploCart')

    plt.xticks([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31], \
               ["0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1", "1-1.2", \
                "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0", "2.0-2.2"], rotation=45)

    for idx, i in enumerate([0,2,4,6,8,10,12,14,16,18]):
        plt.text(i, 100, "N="+str(len(hc_data[idx])), rotation=45, fontsize=6, c="orange")
        plt.text(i+0.5, 100, "N="+str(len(hg_data[idx])), rotation=45, fontsize=6, c="blue")
    for idx2, j in enumerate([20,22,24,26,28,30]):
        plt.text(j, 60, "N="+str(len(hc_data[idx2+10])), rotation=45, fontsize=6, c="orange")
        plt.text(j+0.5, 60, "N="+str(len(hg_data[idx2+10])), rotation=45, fontsize=6, c="blue")
    plt.legend(handles=[hg_patch, hc_patch], borderpad=0.45, prop={'size': 8})
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
    depth_dict = {}
    with open("../data/fastq_no_numt_sim_depths.txt", "r") as f:
        for line in f:
            split=line.split()
            file, depth = dequote(split[0].split("/")[-1]), float(split[1])
            depth_dict.update({file:depth})

    haplogrep_no_numt_score_dict = pickle.load(open("../data/pickles/hg_fastq_no_numt.pk", "rb"))
    haplocart_no_numt_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_no_numt.pk", "rb"))
    haplogrep_with_numt_score_dict = pickle.load(open("../data/pickles/hg_fastq_with_numt.pk", "rb"))
    haplocart_with_numt_score_dict = pickle.load(open("../data/pickles/haplocart_fastq_with_numt.pk", "rb"))

    make_violinplot_fq(haplogrep_no_numt_score_dict, haplocart_no_numt_score_dict, \
                       "../data/fastq_no_numt_sim_depths.txt", "../data/pngs/fastq_no_numt.png", numt=False)

    make_violinplot_fq(haplogrep_with_numt_score_dict, haplocart_with_numt_score_dict, \
                      "../data/fastq_with_numt_sim_depths.txt", "../data/pngs/fastq_with_numt.png", numt=True)


plot_bam()
plot_fastq()
