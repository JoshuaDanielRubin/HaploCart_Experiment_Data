import seaborn as sns
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

fasta_dict = {"HV4b":"KP340180", "H2a2a1":"JX154035", "L0a1a1":"MK295855", "L1c4b":"MN894773", "B2b3a":"MW057682",
              "J2a1a1a1":"MZ190830", "L2a1a3c":"KR135866", "T2e1a1b1":"JN828512", "L2a1j":"KR135846",
              "Q1":"MN849793", "C1b":"MN894713", "Z1a":"MG660559", "D1":"KP172430",
               "A2":"MZ387838", "Y1b":"GU123044", "P":"MN849673", "L3b1":"KT819256",
               "U2e1b1":"KT698031", "F1a1":"MH553920", "I2b":"MN516596", "S1a":"DQ404440", "V3":"MN516629",
               "X3a":"JQ245804", "E1a1b":"EF061151"}

thousand_genomes_dict = {"HG00473":"D6a1a","HG01051":"K1a4a1h","HG02666":"L3e4a","HG03112":"L1b1a3",\
                         "NA18510":"L0a1a3","NA19036":"L3b1a1a","NA19661":"D1h1","NA20518":"J2a1a1e"}

color_map = cm.get_cmap('BuPu')

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14
HUGE_SIZE=17

def plot_posterior(fastq = False, numt=False, bam=False):
    plt.rc('font', size=HUGE_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=HUGE_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=HUGE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=HUGE_SIZE)  # fontsize of the figure title

    if bam:
        result_file = "../data/haplocart_results/bams.txt"
        posterior_file = '../data/haplocart_results/bam_posterior.txt'
        title = 'Posterior Probabilities for Thousand Genome Samples'
        outfile = "../data/pngs/bam_posterior.png"
    elif fastq == True and numt == False:
        result_file = "../data/haplocart_results/fastq_no_numt.txt"
        posterior_file = '../data/haplocart_results/no_numt_posterior.txt'
        title = 'Posterior Probabilities for Simulated FASTQ'
        outfile = '../data/pngs/fastq_posterior.png'
    elif fastq and numt:
        result_file = "../data/haplocart_results/fastq_with_numt.txt"
        posterior_file = '../data/haplocart_results/with_numt_posterior.txt'
        outfile = '../data/pngs/fastq_posterior_with_numt.png'

    with open(result_file, "rt") as f:
        with open(posterior_file, "rt") as g:
            result_df = pd.read_csv(f, sep='\t', header=None)
            if fastq:
                result_df.columns=['sample', 'pred', 'n_reads']
            else:
                result_df.columns=['sample', 'pred', 'threads', 'time', 'n_reads']
            posterior_df = pd.read_csv(g, sep='\t', header=None)
            posterior_df.columns=['sample', 'clade', 'posterior', 'tree_depth']
            posterior_df['ground_truth'] = posterior_df.apply (lambda row: row['sample'].split("_")[0], axis=1)
            print(posterior_df['ground_truth'].unique())
            if fastq:
                posterior_df['target_depth'] = posterior_df.apply (lambda row: float(row['sample'].split("_s")[-1].split("_interleave_tmp")[0]), axis=1)
            elif bam:
                posterior_df['target_depth'] = posterior_df.apply (lambda row: float(row['sample'].split("_")[-1].split(".", 1)[-1].split(".bam")[0]), axis=1)
            merged_df = pd.merge(posterior_df, result_df, how='outer', on='sample')
            if bam:
                pickle.dump(merged_df, open("../data/pickles/bam_merged.pk", "wb"))
            elif fastq and (numt == False):
                pickle.dump(merged_df, open("../data/pickles/fastq_merged.pk", "wb"))
            merged_dfs = [v for k, v in merged_df.groupby('ground_truth')]

            if bam:
                for i in range(8):
                    plt.figure(figsize=(10, 15))
                    plt.title("Posterior probability of confidence in assignment\n for sample " + merged_dfs[i]['ground_truth'].iloc[0] \
                                       + " (haplogroup " + thousand_genomes_dict[merged_dfs[i]['ground_truth'].iloc[0]] + ")")

                    plt.xlabel('Coverage depth (X)')
                    plt.ylabel('Posterior probability of sample \n being subtended by the given clade')
                    merged_dfs[i]['tree_depth'] = merged_dfs[i]['tree_depth'].apply(pd.to_numeric)
                    print(i)
                    merged_dfs[i] = merged_dfs[i].sample(frac=1).reset_index(drop=True)
                    for label, grp in merged_dfs[i].groupby('tree_depth'):
                        sns.lineplot(color=color_map(1/(float(label)+1.0001)), x=grp.target_depth, y=grp.posterior, alpha=0.8, linewidth=7, \
                                     sort=True, orient='x', label=label, estimator=np.mean, errorbar=('ci', 0))
                    plt.tight_layout()
                    plt.ylim(0.25, 1.02)
                    leg = plt.legend(title="Tree depth")
                    for legobj in leg.legendHandles:
                        legobj.set_linewidth(2.0)
                    plt.savefig("../data/pngs/bam_posteriors/bam_posterior_"+str(i)+".png", dpi=300)
                    plt.clf()

            elif fastq:
                for i in range(24):
                    plt.figure(figsize=(15, 10))
                    plt.title("Posterior probability of confidence in assignment\n for sample " + fasta_dict[merged_dfs[0]['ground_truth'].iloc[0]] \
                                   + " (haplogroup " + merged_dfs[0]['ground_truth'].iloc[0] + ")")

                    plt.xlabel('Coverage Depth (X)')
                    plt.ylabel('Posterior probability of sample \n being subtended by the given clade')
                    merged_dfs[i]['tree_depth'] = merged_dfs[i]['tree_depth'].apply(pd.to_numeric)
                    print(i)
                    merged_dfs[i] = merged_dfs[i].sample(frac=1).reset_index(drop=True)
                    for label, grp in merged_dfs[i].groupby('tree_depth'):
                        sns.lineplot(color=color_map(1/(float(label)+1.0001)), x=grp.target_depth, y=grp.posterior, alpha=0.8, \
                        linewidth=7, sort=True, orient='x', errorbar=('ci', 0), label=label, estimator=np.mean)
                    plt.ylim(0.25, 1.02)
                    leg = plt.legend(title="Tree depth")
                    for legobj in leg.legendHandles:
                        legobj.set_linewidth(2.0)
                    if numt == False:
                        plt.savefig("../data/pngs/fastq_posteriors/fastq_posterior_no_numt_"+str(i)+".png")
                    elif numt:
                        plt.savefig("../data/pngs/fastq_posteriors/fastq_posterior_with_numt_"+str(i)+".png")
                    plt.clf()


plot_posterior(bam=True)
#plot_posterior(fastq=True, numt=False)
#plot_posterior(fastq=True, numt=True)

