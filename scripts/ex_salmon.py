import subprocess as sp

Runs = []
with open("SraAccList.txt", "r") as filin:
    for ligne in filin:
        Runs.append(ligne[:-1])

Runs.pop()

#a looper si ça marche
r2 = "data/" + Runs[1] + "_1.fastq"
dest = "quants/" + Runs[1]
sp.run(["salmon", "quant", "-i", "salmon__partial_sa_index__default/default", "-l", "A", "-r", r2, "--validateMappings", "--gcBias", "-p", "4", "-o", dest])
