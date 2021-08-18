with open("SRR_Acc_List.txt", "r") as filin:
    DATASETS = filin.read().splitlines()
    
DIRS = ["seqs", "trim"]

rule all:
  input:
      expand("data/quants/{dataset}/quant.sf", dataset=DATASETS),
      expand("analysis/{dirs}/multiqc_report.html", dirs=DIRS)

rule multiqc:
    input:
        expand("analysis/{{dirs}}/{dataset}_{idx}_fastqc.zip", dataset=DATASETS, idx=["1","2"])
    output:
        "analysis/{dirs}/multiqc_report.html"
    params:
        dirs="{dirs}"
    shell:
        "multiqc ./analysis/{params.dirs} -o analysis/{params.dirs}" 

rule fastqc:
    input:
        "data/{dirs}/{sample}.fastq.gz"
    output:
        "analysis/{dirs}/{sample}_fastqc.zip"
    params:
        dirs="{dirs}"
    threads: 1
    shell:
        "fastqc -o analysis/{params.dirs} -t {threads} {input}"

rule cutadapt:
    input:
        r1="data/seqs/{sample}_1.fastq.gz",
        r2="data/seqs/{sample}_2.fastq.gz"
    output:
        o1="data/trim/{sample}_1.fastq.gz",
        o2="data/trim/{sample}_2.fastq.gz"
    threads: 1
    params:
        a1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        a2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        ml="20" #minimum length
    shell:
        "cutadapt -a {params.a1} -A {params.a2} -m {params.ml} -o "
        "{output.o1} -p {output.o2} --length-tag 'length=' "
        "--cores={threads} {input.r1} {input.r2}" 


rule salmon_quant:
    input:
        r1 = "data/trim/{sample}_1.fastq.gz",
        r2 = "data/trim/{sample}_2.fastq.gz",
        index = "data/index/salmon_partial_sa_index__default/default"
    output:
        "data/quants/{sample}/quant.sf"
    params:
        dir = "data/quants/{sample}"
    threads: 1
    shell:
        "salmon quant -i {input.index} -l A -p {threads} --validateMappings "
        "--gcBias --numGibbsSamples 20 -o {params.dir} "
        "-1 {input.r1} -2 {input.r2}"
