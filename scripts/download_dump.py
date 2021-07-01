#needs to be run from data/seqs

import subprocess as sp

Runs = []
with open("../../SRR_Acc_List.txt", "r") as filin:
    Runs = filin.read().splitlines()


for r in Runs:
    sp.run(["prefetch", "-v", r])

for r in Runs:
    sp.run(["fastq-dump", "--split-files", "--gzip", r])

for r in Runs:
    sp.run(["vdb-validate", r])
    
