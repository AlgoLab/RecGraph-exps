import os

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("output/HLA_full/genes/{gene}/{rec}/reads_{err}.fa", 
        gene=genes,
        rec=[0,1,2],
        err=[0,3, 5]
        ),

rule generate_0:
    input:
        reads_file = "output/HLA_full/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/0/reads_0.fa"
    shell:
        "cp {input.reads_file}  {output.results_file}"

rule generate_3:
    input:
        reads_file = "output/HLA_full/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/0/reads_3.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file}"
        
rule generate_5:
    input:
        reads_file = "output/HLA_full/genes/{gene}/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/0/reads_5.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file}"

rule generate_2_rec_exact:
    input:
        graph = "output/HLA_full/genes/{gene}/graph.gfa"
    output:
        results_file = "output/HLA_full/genes/{gene}/2/reads_0.fa"
    shell:
        "python scripts/generate_rec_reads.py {input.graph} 10 > {output.results_file}"

rule generate_2_rec_3:
    input:
        reads_file = "output/HLA_full/genes/{gene}/2/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/2/reads_3.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file}"

rule generate_2_rec_5:
    input:
        reads_file = "output/HLA_full/genes/{gene}/2/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/2/reads_5.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file}"

rule generate_1_rec_exact:
    input:
        graph = "output/HLA_full/genes/{gene}/graph.gfa"
    output:
        results_file = "output/HLA_full/genes/{gene}/1/reads_0.fa"
    shell:
        "python scripts/recomb_reads.py {input.graph} 10 > {output.results_file}"

rule generate_1_rec_3:
    input:
        reads_file = "output/HLA_full/genes/{gene}/1/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/1/reads_3.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.03 > {output.results_file}"

rule generate_1_rec_5:
    input:
        reads_file = "output/HLA_full/genes/{gene}/1/reads_0.fa"
    output:
        results_file = "output/HLA_full/genes/{gene}/1/reads_5.fa"
    shell:
        "python scripts/fasta_pert.py {input.reads_file} --p 0.05 > {output.results_file}"  



