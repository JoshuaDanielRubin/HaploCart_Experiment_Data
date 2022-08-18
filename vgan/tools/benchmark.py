import time
import subprocess
from subprocess import STDOUT, call, PIPE
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

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
    subprocess.run("/home/ctools/interleave_fastq/interleavefastq.sh " + \
                                    "../src/simulations/"+id+"_n1000_l" + length+"_r1_rep"+rep+"_s"+rate+".fq.gz " + \
                                    "../src/simulations/"+id+"_n1000_l" + length+"_r2_rep"+rep+"_s"+rate+".fq.gz " + \
                                    "../src/simulations/"+id+"_n1000_l" + length+"_rep"+rep+"_s"+rate+ \
                                    ".fq.gz | sed -r '/^\s*$/d' > ../src/vgan_tmpdir/" + \
                                    id + "_" + length + "_" + rep + "_s" + rate + "_interleave_tmp; " + \
                                    "./../bin/vgan haplocart -t 30 -p -o ../data/haplocart_results/fastq_no_numt.txt" + \
                                    " -pf ../data/haplocart_results/no_numt_posterior.txt -s " + id + "_n1000_l" + length + "_rep" + rep + "_s"+rate \
                                    + " -fq1 ../src/vgan_tmpdir/" + \
                                    id + "_" + length + "_" + rep + "_s" + rate + "_interleave_tmp -i", shell=True)


def fastq_with_numt(tup):
    (id, length, rep, rate) = tup
    subprocess.run(["/home/ctools/interleave_fastq/interleavefastq.sh " + \
                                    "../src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+".fq.gz " + \
                                    "../src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+"_r1.fq.gz " + \
                                    "../src/simulations/"+"numtS_and_"+id+"_n1000_l" + length+"_rep"+rep+"_nr200_s"+rate+ \
                                    "_r2.fq.gz | sed -r '/^\s*$/d' > ../src/vgan_tmpdir/" + \
                                    id + "_l" + length + "_rep" + rep + "_nr200_s" + rate + "_interleave_tmp; " + \
                                    "./../bin/vgan haplocart -p -o ../data/haplocart_results/fastq_with_numt.txt" + \
                                    " -s numtS_and_" + id + "_n1000_l" + length + "_rep" + rep + "_nr200_s"+rate + \
                                    " -pf ../data/haplocart_results/with_numt_posterior.txt -fq1 ../src/vgan_tmpdir/" + id + "_l" + length + \
                                    "_rep" + rep + "_nr200_s" + rate + "_interleave_tmp -i"], shell=True)


def mask(tup):
    (id, rep, N) = tup
    return subprocess.Popen(["./../bin/vgan haplocart -t 35 -f " + "../src/simulations/mask/" + id + "_mask" + str(N) + "_rep" + str(rep) + \
    ".fa -s " + id + "_mask" + str(N) + "_rep" + str(rep) + " -o ../data/haplocart_results/mask.txt"], shell=True)


def bam(tup):
    (id, target, replicate, file) = tup
    return subprocess.Popen(["samtools sort -@ 30 ../src/simulations/thousand_genomes/"+id+"/" + file + \
            " | samtools bam2fq -@ 20 > ../src/vgan_tmpdir/" + id + "_" + target + "_" + replicate + \
            "; ./../bin/vgan haplocart -p -pf ../data/haplocart_results/bam_posterior.txt -t 30 -i -fq1 ../src/vgan_tmpdir/"+ id + "_" + target + "_" + replicate + " -s " + \
            file.split("/")[-1] + " -o ../data/haplocart_results/bams.txt; rm ../src/vgan_tmpdir/" + id + "_" + target + "_" + replicate + \
            "; rm ../src/vgan_tmpdir/" + file.split("/")[-1] + ".gam; rm ../src/vgan_tmpdir/" + file.split("/")[-1] + "_sorted.gam"] ,shell=True)

def run():

    # Consensus FASTA
    #for id in IDS:
    #    subprocess.Popen(["./../src/vgan haplojuice -v -o ../data/haplojuice_results/consensus.txt -s " + \
    #                    id + "_consensus " + "-f ../src/simulations/"+id+".fa"], shell=True)

    # No NuMT FASTQ

    #fastq_no_numt_cmds = []
    #for id in IDS:
    #    for length in LENGTHS:
    #        for rep in REPLICATES:
    #            for rate in FQRATES:
    #                fastq_no_numt_cmds.append((id, length, rep, rate))
    #p = Pool(25)
    #p.map(fastq_no_numt, fastq_no_numt_cmds)

   # With NuMTs

    #fastq_with_numt_cmds = []
    #for id in IDS:
    #    for length in LENGTHS:
    #        for rep in REPLICATES:
    #            for rate in FQRATES:
    #                fastq_with_numt_cmds.append((id, length, rep, rate))

    #p = Pool(20)
    #p.map(fastq_with_numt, fastq_with_numt_cmds)

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
    counter = 0
    for id in IDS:
        for rep in REPLICATES:
            for N in MASKN:
                p = mask((id, rep, N))
                counter += 1
                if counter % 20 == 0:
                    p.wait()


if __name__ == "__main__":
    run()



