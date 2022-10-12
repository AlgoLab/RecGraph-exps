from os.path import join as pjoin


configfile: "pandmakegraphs-config.yaml"


SEQSDIR = config["seqsdir"]
RSPOA_BIN = config["rspoa"]
MAKEPRG_BIN = config["mkprg"]

RECPAR = "1-150"

SEQS = config["seqs"]


rule run:
    input:
        expand(pjoin(SEQSDIR, "{seq}", "paperino"), seq=SEQS),


rule gunzip:
    input:
        pjoin(SEQSDIR, "{seq}.fa.gz"),
    output:
        pjoin(SEQSDIR, "{seq}.fa"),
    shell:
        """
        gunzip -k {input}
        """


rule make_graph:
    input:
        msa=pjoin(SEQSDIR, "{seq}.fa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.gfa"),
    params:
        prefix=pjoin(SEQSDIR, "{seq}", "graph"),
    threads: workflow.cores / 4
    shell:
        """
        {MAKEPRG_BIN} from_msa -i {input.msa} -o {params.prefix} -O g -t {threads}
        python3 ./scripts/clean_gfa_from_ast.py {params.prefix}.prg.gfa > {output.gfa}
        """


rule tsort_graph:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
    params:
        prefix=pjoin(SEQSDIR, "{seq}", "graph"),
    conda:
        "envs/vg.yaml"
    shell:
        """
        python3 ./scripts/increment_gfa_idx.py {input.gfa} > {params.prefix}.incr.gfa
        odgi build -g {params.prefix}.incr.gfa -o {params.prefix}.incr.gfa.og
        odgi sort -i {params.prefix}.incr.gfa.og -o {params.prefix}.incr.tsorted.gfa.og
        odgi view -i {params.prefix}.incr.tsorted.gfa.og -g > {output.gfa}
        """


rule check_acyclicity:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
    output:
        log=pjoin(SEQSDIR, "{seq}", "graph.gfa.ACYFLAG"),
    conda:
        "envs/vg.yaml"
    shell:
        """
        vg stats -A {input.gfa} > {output.log}
        """


rule remove_indels:
    input:
        msa=pjoin(SEQSDIR, "{seq}.fa"),
    output:
        fa=pjoin(SEQSDIR, "{seq}", "seqs.noindels.fa"),
    conda:
        "envs/pyt.yaml"
    shell:
        """
        python3 ./scripts/rm_indel_msa.py {input.msa} > {output.fa}
        """


rule rspoa_align_1:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "seqs.noindels.fa"),
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
    output:
        gaf=pjoin(SEQSDIR, "{seq}", "rspoa.gaf"),
    threads: workflow.cores / 2
    shell:
        """
        {RSPOA_BIN} {input.fa} {input.gfa} > {output.gaf}
        """


rule add_paths:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.gfa"),
        gaf=pjoin(SEQSDIR, "{seq}", "rspoa.gaf"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    shell:
        """
        cp {input.gfa} {output.gfa}
        cut -f 1,6 {input.gaf} | while read idx p ; do echo -e "P\\t$idx\\t$(echo $p | cut -c 2- | sed "s/>/+,/g")+\\t*" ; done >> {output.gfa}
        """


rule remove_nodes_not_in_a_path:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    shell:
        """
        python3 ./scripts/remove_nodes_not_in_a_path.py {input.gfa} > {output.gfa}
        """


rule minimal_path_cover:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    output:
        gfa=pjoin(SEQSDIR, "{seq}", RECPAR, "graph.gfa"),
        recombinations=pjoin(SEQSDIR, "{seq}", RECPAR, "recombinations.log"),
    shell:
        """
        cat {input.gfa} | python3 ./scripts/minimal_path_cover.py > {output.gfa} 2> {output.recombinations}
        """


checkpoint get_mosaics:
    input:
        log=pjoin(SEQSDIR, "{seq}", RECPAR, "recombinations.log"),
    output:
        txt=pjoin(SEQSDIR, "{seq}", RECPAR, "mosaics.txt"),
    shell:
        """
        grep "is a mosaic" {input.log} | cut -f 3 -d' ' | sort > {output.txt}
        """


###########################
###  AGGREGATE FUNCTION ###
###########################
def aggregate_input(wildcards):
    mosaics = []
    with checkpoints.get_mosaics.get(**wildcards).output[0].open() as f:
        for line in f:
            mosaic = line.strip("\n")
            mosaics.append(pjoin(SEQSDIR, wildcards.seq, RECPAR, mosaic, "sequence.fa"))
    return mosaics


rule run_after_aggregation:
    input:
        aggregate_input,
    output:
        pjoin(SEQSDIR, "{seq}", "paperino"),
    shell:
        """
        touch {output}
        """


# CHECKME this rule crashes randomly - but if you keep running the snakemake, it'll run succesfully sooner or later
# [W::fai_get_val] Reference NZ_CP008801-HQ24_RS22550-1-cytosol_aminopeptidase not found in FASTA file, returning empty sequence
# [faidx] Failed to fetch sequence in NZ_CP008801-HQ24_RS22550-1-cytosol_aminopeptidase
rule get_mosaic_sequence:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "seqs.noindels.fa"),
    output:
        fa=pjoin(SEQSDIR, "{seq}", RECPAR, "{mosaic}", "sequence.fa"),
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools faidx {input.fa} {wildcards.mosaic} > {output.fa}
        """
