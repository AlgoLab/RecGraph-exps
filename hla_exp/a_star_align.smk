import os
from os.path import join as pjoin
import glob

OUTPUT_DIR = "alignments/"
LOG_DIR = "logs/"

RECGRAPH_BIN = config["recgraph"]
OLD_RECGRAPH_BIN = config["old_recgraph"]
genes = []
with open("genes_HLA.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

rule all:
    input:
        expand("alignments/0/{gene}.gaf", gene=genes),
        expand("alignments/3/{gene}.gaf", gene=genes),
        expand("alignments/5/{gene}.gaf", gene=genes),
        expand("alignments/10/{gene}.gaf", gene=genes),
        expand("alignments/old_0/{gene}.gaf", gene=genes),
        expand("alignments/old_3/{gene}.gaf", gene=genes),
        expand("alignments/old_5/{gene}.gaf", gene=genes),
        expand("alignments/old_10/{gene}.gaf", gene=genes),
        expand("alignments/graph_aligner_0/{gene}.gaf", gene=genes),
        expand("alignments/graph_aligner_3/{gene}.gaf", gene=genes),
        expand("alignments/graph_aligner_5/{gene}.gaf", gene=genes),
        expand("alignments/graph_aligner_10/{gene}.gaf", gene=genes),

rule align_exact:
    input:
        read = "../data/HLA-{gene}/reads.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "0/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "0/{gene}.gaf"
    shell:
        " /usr/bin/time -v {RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 30 -r 16 > {output.gaf} 2> {log}"

rule align_3:
    input:
        read = "../data/HLA-{gene}/reads_3.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "3/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "3/{gene}.gaf"
    shell:
        " /usr/bin/time -v {RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 8 -r 16 > {output.gaf} 2> {log}"

rule align_5:
    input:
        read = "../data/HLA-{gene}/reads_5.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "5/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "5/{gene}.gaf"
    shell:
        " /usr/bin/time -v {RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 8 -r 16 > {output.gaf} 2> {log}"

rule align_10:
    input:
        read = "../data/HLA-{gene}/reads_10.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "10/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "10/{gene}.gaf"
    shell:
        " /usr/bin/time -v {RECGRAPH_BIN} -q {input.read} -g {input.graph} -s 6 -r 16 > {output.gaf} 2> {log}"

rule align_exact_old_recgraph:
	input:
		read = "../data/HLA-{gene}/reads.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
	log:
        LOG_DIR + "old_0/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "old_0/{gene}.gaf"
	shell:
		"/usr/bin/time -v {OLD_RECGRAPH_BIN} -q {input.read} -g {input.graph} -k 1 -M 1 -X 0 -O 0 -E 0 -r 16 -d 0 > {output.cigar} 2> {log}"

rule align_3_old_recgraph:
    input:
        read = "../data/HLA-{gene}/reads_3.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "old_3/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "old_3/{gene}.gaf"
    shell:
        "/usr/bin/time -v {OLD_RECGRAPH_BIN} -q {input.read} -g {input.graph} -k 1 -M 1 -X 0 -O 0 -E 0 -r 16 -d 0 > {output.cigar} 2> {log}"

rule align_5_old_recgraph:
    input:
        read = "../data/HLA-{gene}/reads_5.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "old_5/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "old_5/{gene}.gaf"
    shell:
        "/usr/bin/time -v {OLD_RECGRAPH_BIN} -q {input.read} -g {input.graph} -k 1 -M 1 -X 0 -O 0 -E 0 -r 16 -d 0 > {output.cigar} 2> {log}"

rule align_10_old_recgraph:
    input:
        read = "../data/HLA-{gene}/reads_10.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "old_10/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "old_10/{gene}.gaf"
    shell:
        "/usr/bin/time -v {OLD_RECGRAPH_BIN} -q {input.read} -g {input.graph} -k 1 -M 1 -X 0 -O 0 -E 0 -r 16 -d 0 > {output.cigar} 2> {log}"

rule graph_aligner_exact:
    input:
        read = "../data/HLA-{gene}/reads.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "graph_aligner_0/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "graph_aligner_0/{gene}.gaf"
    conda:
        "envs/graphaligner.yaml"

    threads: 1

    shell:
        "/usr/bin/time -v  GraphAligner -g {input.graph} -f {input.read} -t 1 -a {output.gaf} -x vg 2> {log}"

rule graph_aligner_3:
    input:
        read = "../data/HLA-{gene}/reads_3.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "graph_aligner_3/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "graph_aligner_3/{gene}.gaf"
    conda:
        "envs/graphaligner.yaml"

    threads: 1

    shell:
        "/usr/bin/time -v  GraphAligner -g {input.graph} -f {input.read} -t 1 -a {output.gaf} -x vg 2> {log}"

rule graph_aligner_5:
    input:
        read = "../data/HLA-{gene}/reads_5.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "graph_aligner_5/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "graph_aligner_5/{gene}.gaf"
    conda:
        "envs/graphaligner.yaml"

    threads: 1

    shell:
        "/usr/bin/time -v  GraphAligner -g {input.graph} -f {input.read} -t 1 -a {output.gaf} -x vg 2> {log}"

rule graph_aligner_10:
    input:
        read = "../data/HLA-{gene}/reads_10.fa",
        graph = "../data/HLA-{gene}/graph.gfa"
    log:
        LOG_DIR + "graph_aligner_10/{gene}.log"
    output:
        gaf = OUTPUT_DIR + "graph_aligner_10/{gene}.gaf"
    conda:
        "envs/graphaligner.yaml"

    threads: 1

    shell:
        "/usr/bin/time -v  GraphAligner -g {input.graph} -f {input.read} -t 1 -a {output.gaf} -x vg 2> {log}"
        


