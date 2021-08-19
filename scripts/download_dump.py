#needs to be run from data/seqs

import subprocess as sp
import os

Runs = []
with open("../../SRR_Acc_List.txt", "r") as filin:
    Runs = filin.read().splitlines()


for r in Runs:
    if not os.path.exists(r + "/" + r + ".sra"):
        if os.path.exists(r + "/" + r + ".sra.lock"):
           sp.run(["rm", "-r", r])
        sp.run(["prefetch", "-v", r])
    if not os.path.exists(r + ".fastq.gz"):
        sp.run(["fastq-dump", "--gzip", r])
        sp.run(["vdb-validate", r])

# Maybe delete folders with the .sra only inside afterwards

for r in Runs:
    sp.run(["rm", "-r", r])
    
