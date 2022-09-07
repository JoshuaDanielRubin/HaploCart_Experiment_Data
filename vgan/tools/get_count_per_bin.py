import pandas as pd

def bam():
    bam_depth_file = "../data/1kdepths.txt"
    bins=[0, 0.15, 0.3, 0.45, 0.6, 0.75, 10.0]

    with open(bam_depth_file, "r") as f:
        df = pd.read_csv(f, sep=' ', names=["sample", "depth"])
        df["depth"] = df["depth"].astype(float)
        print(df["depth"].value_counts(bins=bins))


def fastq_no_numt():
    bins=[0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,1,10]
    fastq_no_numt_depth_file = "../data/fastq_no_numt_sim_depths.txt"

    with open(fastq_no_numt_depth_file, "r") as f:
        df = pd.read_csv(f, sep=' ', names=["sample", "depth"])
        df["depth"] = df["depth"].astype(float)
        print(df["depth"].value_counts(bins=bins))

def fastq_with_numt():
    bins=[0,0.1,0.2,0.3,0.4,0.5,10.0]
    fastq_with_numt_depth_file = "../data/fastq_with_numt_sim_depths.txt"

    with open(fastq_with_numt_depth_file, "r") as f:
        df = pd.read_csv(f, sep=' ', names=["sample", "depth"])
        df["depth"] = df["depth"].astype(float)
        print(df["depth"].value_counts(bins=bins))


bam()
#fastq_no_numt()
#fastq_with_numt()
