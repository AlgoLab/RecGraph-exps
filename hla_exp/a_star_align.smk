import os
from os.path import join as pjoin
import glob

OUTPUT_DIR = "alignments/"
LOG_DIR = "logs/"
RECGRAPH_BIN = (
    config["recgraph"]
    if "recgraph" in config
    else "/data/RecGraph/target/release/recgraph"
)
genes = []
with open("genes_HLA.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("alignments/0/{gene}.gaf", gene=genes),
        expand("alignments/3/{gene}.gaf", gene=genes),
        expand("alignments/5/{gene}.gaf", gene=genes),
        expand("alignments/10/{gene}.gaf", gene=genes),

rule align_exact:
    input:
        read = "../data/HLA-{gene}/reads.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "0/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "0/{gene}.gaf"
    shell:
        "{RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 30 -r 16 > {output.gaf} 2> {log}"

rule align_3:
    input:
        read = "../data/HLA-{gene}/reads_3.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "3/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "3/{gene}.gaf"
    shell:
        "{RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 8 -r 16 > {output.gaf} 2> {log}"

rule align_5:
    input:
        read = "../data/HLA-{gene}/reads_5.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "5/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "5/{gene}.gaf"
    shell:
        "{RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 8 -r 16 > {output.gaf} 2> {log}"

rule align_10:
    input:
        read = "../data/HLA-{gene}/reads_10.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "10/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "10/{gene}.gaf"
    shell:
        "{RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 6 -r 16 > {output.gaf} 2> {log}"


