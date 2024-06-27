import os

genes = []
with open("genes_HLA.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("../data/HLA-{gene}/reads_3.fa", gene=genes),
        expand("../data/HLA-{gene}/reads_5.fa", gene=genes),
        expand("../data/HLA-{gene}/reads_10.fa", gene=genes)

rule generate_3:
    input:
        reads_file = "../data/HLA-{gene}/reads.fa"
    output:
        results_file = "../data/HLA-{gene}/reads_3.fa"
    shell:
        "python fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file}"
        
rule generate_5:
    input:
        reads_file = "../data/HLA-{gene}/reads.fa"
    output:
        results_file = "../data/HLA-{gene}/reads_5.fa"
    shell:
        "python fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file}"

rule generate_10:
    input:
        reads_file = "../data/HLA-{gene}/reads.fa"
    output:
        results_file = "../data/HLA-{gene}/reads_10.fa"
    shell:
        "python fasta_pert.py {input.reads_file} --p 0.10 > {output.results_file}"