import os

genes = []
with open("genes_HLA.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("output/HLA/genes/{gene}/reads_{err}.fa", 
        gene=genes,
        err=[3, 5, 10]
        ),
        expand("output/HLA/genes/{gene}/reads_rec.fa",
        gene=genes,
        )

rule generate_3:
    input:
        reads_file = "output/HLA/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA/genes/{gene}/reads_3.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file}"
        
rule generate_5:
    input:
        reads_file = "output/HLA/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA/genes/{gene}/reads_5.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file}"

rule generate_10:
    input:
        reads_file = "output/HLA/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA/genes/{gene}/reads_10.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.10 > {output.results_file}"

rule generate_rec:
    input:
        graph = "output/HLA/genes/{gene}/graph.gfa"
    output:
        results_file = "output/HLA/genes/{gene}/reads_rec.fa"
    shell:
        "python scripts/generate_rec_reads.py {input.graph} 10 > {output.results_file}"