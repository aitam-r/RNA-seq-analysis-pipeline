import subprocess as sp
import os

if not os.path.exists("./data/seqs"):
    os.mkdir("./data/seqs")

# Should then be run from project root, absolutely
os.chdir("./data/seqs")

Runs = []
with open("../../SRR_Acc_List.txt", "r") as filin:
    Runs = filin.read().splitlines()


for r in Runs:
    if not os.path.exists(r) and not (os.path.exists(r + "_2.fastq.gz") or os.path.exists(r + "_1.fastq.gz")):
        print("Prefetching...")
        sp.run(["prefetch", "-v", r])
    if not (os.path.exists(r + "_2.fastq.gz") and os.path.exists(r + "_1.fastq.gz")):
        print("fastq_dumping and validating...")
        sp.run(["fastq-dump", "--gzip", "--split-files", r])
        sp.run(["vdb-validate", r])

# Maybe delete folders with the .sra only inside afterwards
for r in Runs:
    if os.path.exists(r):
        sp.run(["rm", "-r", r])
    
