import subprocess

THOUSAND_GENOME_IDS = ["HG00473","HG01051","HG02666","HG03112","NA18510","NA19036","NA19661","NA20518"]
RATES = ["0.5", "0.3", "0.1", "0.05", "0.03", "0.01", "0.005", "0.003", "0.001", "0.0005", "0.0003", "0.0001", "0.00005"]

for id in THOUSAND_GENOME_IDS:
    for rate in RATES:
        f = open(id+"/"+id+"_"+rate+".bam", "wb")
        subprocess.run(["samtools", "view", "-b", id+"/"+id+".bam", "-s", rate], stdout=f)
