import seaborn as sns
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

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
                result_df.columns=['sample', 'pred', 'threads', 'ms', 'n_reads']
            else:
                result_df.columns=['sample', 'pred', 'threads', 'ms', 'n_reads']
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
                fig, ax = plt.subplots(2,4,figsize=(32, 27))
                ax_dict = {0:ax[0][0], 1:ax[0][1], 2:ax[0][2], 3:ax[0][3],  4:ax[1][0], 5:ax[1][1],
                           6:ax[1][2], 7:ax[1][3]}
                ax[0][0].set_title(merged_dfs[0]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][1].set_title(merged_dfs[1]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][2].set_title(merged_dfs[2]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][3].set_title(merged_dfs[3]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][0].set_title(merged_dfs[4]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][1].set_title(merged_dfs[5]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][2].set_title(merged_dfs[6]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][3].set_title(merged_dfs[7]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")

                ax[0][0].set_xlabel('Coverage depth (X)')
                ax[0][1].set_xlabel('Coverage depth (X)')
                ax[0][2].set_xlabel('Coverage depth (X)')
                ax[0][3].set_xlabel('Coverage depth (X)')
                ax[1][0].set_xlabel('Coverage depth (X)')
                ax[1][1].set_xlabel('Coverage depth (X)')
                ax[1][2].set_xlabel('Coverage depth (X)')
                ax[1][3].set_xlabel('Coverage depth (X)')

                ax[0][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')

            elif fastq:
                fig, ax = plt.subplots(5, 5, figsize=(64, 54))
                ax_dict = {0:ax[0][0], 1:ax[0][1], 2:ax[0][2], 3:ax[0][3], 4:ax[0][4], 5:ax[1][0], 6:ax[1][1],
                           7:ax[1][2], 8:ax[1][3], 9:ax[1][4], 10:ax[2][0], 11:ax[2][1], 12:ax[2][2],
                           13:ax[2][3], 14:ax[2][4], 15:ax[3][0], 16:ax[3][1], 17:ax[3][2], 18:ax[3][3], 19:ax[3][4],
                           20:ax[4][0], 21:ax[4][1], 22:ax[4][2], 23:ax[4][3], 24:ax[4][4]}
                ax[0][0].set_title(merged_dfs[0]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][1].set_title(merged_dfs[1]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][2].set_title(merged_dfs[2]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][3].set_title(merged_dfs[3]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[0][4].set_title(merged_dfs[4]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][0].set_title(merged_dfs[5]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][1].set_title(merged_dfs[6]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][2].set_title(merged_dfs[7]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][3].set_title(merged_dfs[8]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[1][4].set_title(merged_dfs[9]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[2][0].set_title(merged_dfs[10]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[2][1].set_title(merged_dfs[11]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[2][2].set_title(merged_dfs[12]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[2][3].set_title(merged_dfs[13]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[2][4].set_title(merged_dfs[14]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[3][0].set_title(merged_dfs[15]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[3][1].set_title(merged_dfs[16]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[3][2].set_title(merged_dfs[17]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[3][3].set_title(merged_dfs[18]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[3][4].set_title(merged_dfs[19]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[4][0].set_title(merged_dfs[20]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[4][1].set_title(merged_dfs[21]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[4][2].set_title(merged_dfs[22]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[4][3].set_title(merged_dfs[23]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")
                ax[4][4].set_title(merged_dfs[24]['ground_truth'].iloc[0] + " posterior confidence estimates (mean over replicates) ")

                ax[0][0].set_xlabel('Coverage Depth (X)')
                ax[0][1].set_xlabel('Coverage Depth (X)')
                ax[0][2].set_xlabel('Coverage Depth (X)')
                ax[0][3].set_xlabel('Coverage Depth (X)')
                ax[0][4].set_xlabel('Coverage Depth (X)')
                ax[1][0].set_xlabel('Coverage Depth (X)')
                ax[1][1].set_xlabel('Coverage Depth (X)')
                ax[1][2].set_xlabel('Coverage Depth (X)')
                ax[1][3].set_xlabel('Coverage Depth (X)')
                ax[1][4].set_xlabel('Coverage Depth (X)')
                ax[2][0].set_xlabel('Coverage Depth (X)')
                ax[2][1].set_xlabel('Coverage Depth (X)')
                ax[2][2].set_xlabel('Coverage Depth (X)')
                ax[2][3].set_xlabel('Coverage Depth (X)')
                ax[2][4].set_xlabel('Coverage Depth (X)')
                ax[3][0].set_xlabel('Coverage Depth (X)')
                ax[3][1].set_xlabel('Coverage Depth (X)')
                ax[3][2].set_xlabel('Coverage Depth (X)')
                ax[3][3].set_xlabel('Coverage Depth (X)')
                ax[3][4].set_xlabel('Coverage Depth (X)')
                ax[4][0].set_xlabel('Coverage Depth (X)')
                ax[4][1].set_xlabel('Coverage Depth (X)')
                ax[4][2].set_xlabel('Coverage Depth (X)')
                ax[4][3].set_xlabel('Coverage Depth (X)')
                ax[4][4].set_xlabel('Coverage Depth (X)')

                ax[0][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[0][4].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[1][4].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[2][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[2][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[2][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[2][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[2][4].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[3][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[3][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[3][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[3][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[3][4].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[4][0].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[4][1].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[4][2].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[4][3].set_ylabel('Posterior probability of sample \n being subtended by the given clade')
                ax[4][4].set_ylabel('Posterior probability of sample \n being subtended by the given clade')


            for i, df in enumerate(merged_dfs):
                df['tree_depth'] = df['tree_depth'].apply(pd.to_numeric)
                print(i)
                df = df.sample(frac=1).reset_index(drop=True)
                for label, grp in df.groupby('tree_depth'):
                    if bam:
                        sns.lineplot(color=color_map(1/(float(label)+1.0001)), x=grp.target_depth, y=grp.posterior, alpha=0.8, linewidth=7, ax=ax_dict[i], \
                                     sort=True, ci=0, label=label, estimator=np.mean)
                    elif fastq:
                        sns.lineplot(color=color_map(1/(float(label)+1.0001)), x=grp.target_depth, y=grp.posterior, alpha=0.8, linewidth=7, ax=ax_dict[i], sort=True, ci=0, label=label, estimator=np.mean)


            for ax_ in ax_dict.values():
                leg = ax_.legend(title="Tree depth")
                for legobj in leg.legendHandles:
                    legobj.set_linewidth(2.0)
                ax_.set_ylim(0.25, 1.02)
            plt.tight_layout()
            plt.savefig(outfile, dpi=300)


#plot_posterior(bam=True)
plot_posterior(fastq=True, numt=False)
plot_posterior(fastq=True, numt=True)

