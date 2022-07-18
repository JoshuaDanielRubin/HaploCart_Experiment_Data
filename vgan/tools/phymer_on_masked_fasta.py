import subprocess
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures

phymer="nice python /home/ctools/phy-mer-master/Phy-Mer.py /home/ctools/phy-mer-master/PhyloTree_b16_k12.txt"

IDS = ["L0a1a1", "HV4b", "L1c4b", "B2b3a", "H2a2a1", "J2a1a1a1", "L2a1a3c", "T2e1a1b1", "L2a1j", "L1c2b", "Q1"]
MASKN=range(0, 16001, 1000)
REPLICATES = [str(x) for x in range(100)]

def mask(tup):
    (id,N,rep) = tup
    subprocess.run(phymer+" ../src/simulations/mask/" + id + "_mask" + str(N) + "_rep" + str(rep) + \
                     ".fa > ../data/phymer_mask/"+id+"_mask"+str(N)+"_rep"+rep, shell=True)

cmds = []
for id in IDS:
    for N in MASKN:
        for rep in REPLICATES:
            cmds.append((id,N,rep))

futures = []
with ThreadPoolExecutor(max_workers=50) as executor:
    for cmd in cmds:
        futures.append(executor.submit(mask, cmd))
    for future in concurrent.futures.as_completed(futures):
        print(future.result())
