import time
import subprocess
from subprocess import STDOUT, call, PIPE
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

#IDS =["H2a2a1"]
IDS = ["L0a1a1", "HV4b", "C1b", "L1c4b", "B2b3a", "J2a1a1a1", "L2a1a3c", "T2e1a1b1", "L2a1j", "L1c2b", "Q1", \
      "X3a", "E1a1b", "V3", "S1a", "I2b", "F1a1", "U2e1b1", "L3b1", "P", "D1", "A2", "Z1a", "Y1b"]
THOUSAND_GENOME_IDS = ["NA19661", "HG00473","HG01051","HG02666","HG03112","NA18510","NA19036","NA20518"]
LENGTHS=["50", "100"]
FQRATES=["0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.1", "0.2", "0.3"]
BAMTARGETS = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
MASKN=range(0, 16001, 1000)
REPLICATES = [str(x) for x in range(100)]

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s


def fastq_no_numt(tup):
    (id, length, rep, rate) = tup
    return "nice -n 19 /home/ctools/interleave_fastq/interleavefastq.sh " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+id+"_n1000_l" + length+"_r1_rep"+rep+"_s"+rate+".fq.gz " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+id+"_n1000_l" + length+"_r2_rep"+rep+"_s"+rate+".fq.gz " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+id+"_n1000_l" + length+"_rep"+rep+"_s"+rate+ \
                                    ".fq.gz | sed -r '/^\s*$/d' > /home/projects/mito_haplotype/vgan/src/vgan_tmpdir/" + \
                                    id + "_" + length + "_" + rep + "_s" + rate + "_interleave_tmp && " + \
                                    " /home/projects/mito_haplotype/vgan/bin/vgan haplocart -t 1 -p -o /home/projects/mito_haplotype/vgan/data/haplocart_results/no_numt_dir/" + id + "_n1000_l" + length + "_rep" + rep + "_s"+rate + \
                                    " -p -pf /home/projects/mito_haplotype/vgan/data/haplocart_results/no_numt_dir/posterior_" + id + "_n1000_l" + length + "_rep" + rep + "_s"+rate +" -s " + \
                                    id + "_n1000_l" + length + "_rep" + rep + "_s"+rate \
                                    + " -fq1 /home/projects/mito_haplotype/vgan/src/vgan_tmpdir/" + \
                                    id + "_" + length + "_" + rep + "_s" + rate + "_interleave_tmp -i -q\n"


def fastq_with_numt(tup):
    (id, length, rep, rate) = tup
    return "/home/ctools/interleave_fastq/interleavefastq.sh " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+".fq.gz " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+"_r1.fq.gz " + \
                                    "/home/projects/mito_haplotype/vgan/src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+ \
                                    "_r2.fq.gz | sed -r '/^\s*$/d' > /home/projects/mito_haplotype/vgan/src/vgan_tmpdir/" + \
                                    id + "_l" + length + "_rep" + rep + "_nr200_s" + rate + "_interleave_tmp && " + \
                                    "/home/projects/mito_haplotype/vgan/bin/vgan haplocart -p -o /home/projects/mito_haplotype/vgan/data/haplocart_results/with_numt_dir/numtS_and_" + id + "_n1000_l" + length + "_rep" + rep + "_nr200_s"+rate + \
                                    " -s numtS_and_" + id + "_n1000_l" + length + "_rep" + rep + "_nr200_s"+rate + \
                                    " -pf /home/projects/mito_haplotype/vgan/data/haplocart_results/with_numt_dir/posterior_numtS_and_" + id + \
                                    "_n1000_l" + length + "_rep" + rep + "_nr200_s"+rate +" -fq1 /home/projects/mito_haplotype/vgan/src/vgan_tmpdir/" + id + "_l" + length + \
                                    "_rep" + rep + "_nr200_s" + rate + "_interleave_tmp -i -q\n"


def mask(tup):
    (id, rep, N) = tup
    return "nice -n 19 /home/projects/mito_haplotype/vgan/tools/../bin/vgan haplocart -t 1 -f " + \
           "/home/projects/mito_haplotype/vgan/tools/../src/simulations/mask/" + id + "_mask" + str(N) + "_rep" + str(rep) + \
           ".fa -s " + id + "_mask" + str(N) + "_rep" + str(rep) + " -o /home/projects/mito_haplotype/vgan/tools/../data/haplocart_results/mask_dir/"+ \
           id + "_mask" + str(N) + "_rep" + str(rep)+ " -q -p -pf /home/projects/mito_haplotype/vgan/data/haplocart_results/mask_dir/posterior_"+ id + "_mask" + str(N) + "_rep" + str(rep) + "\n"


def bam(tup):
    (id, target, replicate, file) = tup
    return subprocess.Popen(["samtools sort -@ 30 ../src/simulations/thousand_genomes/"+id+"/" + file + \
            " | samtools bam2fq -@ 20 > ../src/vgan_tmpdir/" + id + "_" + target + "_" + replicate + \
            "; ./../bin/vgan haplocart -p -pf ../data/haplocart_results/bam_posterior.txt -t 30 -i -fq1 ../src/vgan_tmpdir/"+ id + "_" + target + "_" + replicate + " -s " + \
            file.split("/")[-1] + " -o ../data/haplocart_results/bams.txt; rm ../src/vgan_tmpdir/" + id + "_" + target + "_" + replicate + \
            "; rm ../src/vgan_tmpdir/" + file.split("/")[-1] + ".gam; rm ../src/vgan_tmpdir/" + file.split("/")[-1] + "_sorted.gam"] ,shell=True)

def run():

    # Consensus FASTA
    for id in IDS:
        subprocess.Popen(["./../bin/vgan haplocart -q -t -1 -o ../data/haplocart_results/consensus.txt -s " + \
                        id + "_consensus " + "-f ../src/simulations/"+id+".fa"], shell=True)

    # No NuMT FASTQ

    #fastq_no_numt_cmds = []
    #for id in IDS:
    #    for length in LENGTHS:
    #        for rep in REPLICATES:
    #            for rate in FQRATES:
    #                fastq_no_numt_cmds.append(fastq_no_numt((id, length, rep, rate)))
    #with open("no_numt_cmds", "wt") as f:
    #    for cmd in fastq_no_numt_cmds:
    #        f.write(cmd)

   # With NuMTs

    #fastq_with_numt_cmds = []
    #for id in IDS:
    #    for length in LENGTHS:
    #        for rep in REPLICATES:
    #            for rate in FQRATES:
    #                fastq_with_numt_cmds.append(fastq_with_numt((id, length, rep, rate)))

    #with open("with_numt_cmds", "wt") as f:
    #    for cmd in fastq_with_numt_cmds:
    #        f.write(cmd)

   # Thousand Genome BAMS

    #counter = 0
    #for replicate in REPLICATES:
    #    for id in THOUSAND_GENOME_IDS:
    #        for target in BAMTARGETS:
    #            file = id+"_"+replicate+"."+target+".bam"
    #            p = bam((id, target, replicate, file))
    #            counter += 1
    #            if counter % 15 == 0:
    #                p.wait()

    # Masking

    #counter = 0
    #mask_cmds = []
    #for id in IDS:
    #    for rep in REPLICATES:
    #        for N in MASKN:
    #            p = mask((id, rep, N))
    #            mask_cmds.append(p)

    #with open("mask_cmds", "wt") as f:
    #    for cmd in mask_cmds:
    #        f.write(cmd)

if __name__ == "__main__":
    run()


