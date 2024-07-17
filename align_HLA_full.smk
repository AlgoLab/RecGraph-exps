import os

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

error_levels = ["0", "3", "5"]
rec_levels = ["0", "1", "2"]

split_read_files = []
for gene in genes:
    for rec in rec_levels:
        for err in error_levels:
            output_dir = f"output/HLA_full/genes/{gene}/{rec}/reads_{err}_split"
            if os.path.exists(output_dir):
                for read_file in os.listdir(output_dir):
                    if read_file.startswith("read_") and read_file.endswith(".fa"):
                        path = os.path.join(f"{gene}/{rec}/reads_{err}_split", read_file.split(".")[0])
                        split_read_files.append(path)

rule all:
    input:
        "bin/recgraph_a_star",
        expand("output/HLA_full/new_recgraph/{read_path}.gaf",
                gene=genes,
                rec=rec_levels,
                err=error_levels,
                read_path=split_read_files
        )
        
rule build_recgraph_a_star:
    output:
        "bin/recgraph_a_star"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        mkdir a_star
        cd a_star
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout a_star
        cargo build --release --jobs {threads}
        cd ../..
        cp a_star/RecGraph/target/release/recgraph {output}
        """

rule run_recgraph_a_star:
    input:
        fa="output/HLA_full/genes/{gene}/{rec}/reads_{err}_split/{read_file}.fa",
        gfa="output/HLA_full/genes/{gene}/graph.gfa",
        rg = "bin/recgraph_a_star"
    log:
        "output/HLA_full/new_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.log"
    output:
        "output/HLA_full/new_recgraph/{gene}/{rec}/reads_{err}_split/{read_file}.gaf"
    shell:
        "/usr/bin/time -v {input.rg} -q {input.fa} -g {input.gfa} -s 8 -r 4 -k 2 > {output} 2> {log}"
