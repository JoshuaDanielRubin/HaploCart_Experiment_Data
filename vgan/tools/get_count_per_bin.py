import pandas as pd


def bam():

    bam_depth_file = "../data/1kdepths.txt"
    bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

    with open(bam_depth_file, "r") as f:
        df = pd.read_csv(f, sep=' ', names=["sample", "depth"])
        df["depth"] = df["depth"].astype(float)
        print(df["depth"].value_counts(bins=bins))


def fastq_no_numt():
    bins=[0, 0.2, 0.4, 0.6, 0.8,1.0,1.2, 1.4, 1.6, 1.8, 2.0]
    fastq_no_numt_depth_file = "../data/fastq_no_numt_sim_depths.txt"

    with open(fastq_no_numt_depth_file, "r") as f:
        df = pd.read_csv(f, sep=' ', names=["sample", "depth"])
        df["depth"] = df["depth"].astype(float)
        print(df["depth"].value_counts(bins=bins))


fastq_no_numt()
