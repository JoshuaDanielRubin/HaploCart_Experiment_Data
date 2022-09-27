import subprocess

IDS = ["H2a2a1", "L0a1a1", "HV4b", "C1b", "L1c4b", "B2b3a", "J2a1a1a1", "L2a1a3c", "T2e1a1b1", "L2a1j", "L1c2b", "Q1", \
      "X3a", "E1a1b", "V3", "S1a", "I2b", "F1a1", "U2e1b1", "L3b1", "P", "D1", "A2", "Z1a", "Y1b"]

for id in IDS:
    subprocess.call("nice python /home/ctools/phy-mer-master/Phy-Mer.py /home/ctools/phy-mer-master/PhyloTree_b16_k12.txt ../src/simulations/backup_fastas/" + id + ".fa", shell=True)
