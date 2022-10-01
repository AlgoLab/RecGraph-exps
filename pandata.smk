from os.path import join as pjoin

# import glob


configfile: "pandata-config.yaml"


SEQSDIR = config["seqs"]
OUTDIR = config["odir"]
GTOOLS = config["graphs"]
RSPOA_BIN = config["rspoa"]
MAKEPRG_BIN = config["mkprg"]

# SEQS = []
# for fa in glob.glob(pjoin(SEQSDIR, "*.fa")):
#     SEQS.append(fa.strip("\n").split("/")[-1].split(".")[0])

SEQS = [
    "GC00000002_131_na_aln",  # ~1000bp <- BIG
    "GC00000748_2_na_aln",  # ~300bp
    "GC00001971_2_na_aln",  # <- BIIG
    "GC00003677_na_aln",  # ~400bp
    "GC00004331_na_aln",  # ~500bp
    "GC00004341_na_aln",  # ~300bp
    "GC00004731_na_aln",  # ~300bp
    "GC00004781_na_aln",  # ~300bp
    #
    # "GC00001875_2_na_aln",  # <- BIIIG
    # "GC00004750_na_aln",  # <- BIIIIG
]


rule run:
    input:
        expand(pjoin(OUTDIR, "{seq}", "graph.gfa.ACYFLAG"), seq=SEQS),
        expand(pjoin(OUTDIR, "{seq}-15x.png"), seq=SEQS),
    output:
        png15=pjoin(OUTDIR, "all-15x.png"),
        png30=pjoin(OUTDIR, "all-30x.png"),
    shell:
        """
        python3 ./scripts/pandata_plot.py -c 15 -t rspoa-5,rspoa-9,giraffe,bwa -o {output.png15} {OUTDIR}
        python3 ./scripts/pandata_plot.py -c 30 -t rspoa-5,rspoa-9,giraffe,bwa -o {output.png30} {OUTDIR}
        """


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
        gfa=pjoin(OUTDIR, "{seq}", "graph.gfa"),
    params:
        prefix=pjoin(OUTDIR, "{seq}", "graph"),
    threads: workflow.cores / 4
    shell:
        """
        {MAKEPRG_BIN} from_msa -i {input.msa} -o {params.prefix} -O g -t {threads}
        python3 scripts/clean_gfa_from_ast.py {params.prefix}.prg.gfa > {output.gfa}
        """


rule tsort_graph:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.gfa"),
    output:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.gfa"),
    params:
        prefix=pjoin(OUTDIR, "{seq}", "graph"),
    conda:
        "envs/vg.yaml"
    shell:
        """
        python3 scripts/increment_gfa_idx.py {input.gfa} > {params.prefix}.incr.gfa
        odgi build -g {params.prefix}.incr.gfa -o {params.prefix}.incr.gfa.og
        odgi sort -i {params.prefix}.incr.gfa.og -o {params.prefix}.incr.tsorted.gfa.og
        odgi view -i {params.prefix}.incr.tsorted.gfa.og -g > {output.gfa}
        """


rule check_acyclicity:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.gfa"),
    output:
        log=pjoin(OUTDIR, "{seq}", "graph.gfa.ACYFLAG"),
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
        fa=pjoin(OUTDIR, "{seq}", "seqs.noindels.fa"),
    shell:
        """
        python3 scripts/rm_indel_msa.py {input.msa} > {output.fa}
        """


rule rspoa_align_1:
    input:
        fa=pjoin(OUTDIR, "{seq}", "seqs.noindels.fa"),
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.gfa"),
    output:
        gaf=pjoin(OUTDIR, "{seq}", "seqs.noindels.gaf"),
    threads: workflow.cores / 2
    shell:
        """
        {RSPOA_BIN} {input.fa} {input.gfa} > {output.gaf}
        """


rule add_paths:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.gfa"),
        gaf=pjoin(OUTDIR, "{seq}", "seqs.noindels.gaf"),
    output:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    shell:
        """
        cp {input.gfa} {output.gfa}
        cut -f 1,6 {input.gaf} | while read idx p ; do echo -e "P\\t$idx\\t$(echo $p | cut -c 2- | sed "s/>/+,/g")+\\t*" ; done >> {output.gfa}
        """


rule remove_nodes_not_in_a_path:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.gfa"),
    output:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    shell:
        """
        python3 ./scripts/remove_nodes_not_in_a_path.py {input.gfa} > {output.gfa}
        """


rule minimal_path_cover:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.gfa"),
    output:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
        recombinations=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.recombinations.log"
        ),
    shell:
        """
        cat {input.gfa} | python3 scripts/minimal_path_cover.py > {output.gfa} 2> {output.recombinations}
        """


checkpoint get_mosaics:
    input:
        log=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.recombinations.log"),
    output:
        txt=pjoin(OUTDIR, "{seq}", "mosaics.txt"),
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
            for coverage in ["15", "30"]:
                mosaics.append(
                    pjoin(OUTDIR, wildcards.seq, mosaic, coverage, "rspoa-5.dist.txt")
                )
                mosaics.append(
                    pjoin(OUTDIR, wildcards.seq, mosaic, coverage, "rspoa-9.dist.txt")
                )
                mosaics.append(
                    pjoin(OUTDIR, wildcards.seq, mosaic, coverage, "giraffe.dist.txt")
                )
                mosaics.append(
                    pjoin(OUTDIR, wildcards.seq, mosaic, coverage, "bwa.dist.txt")
                )
    return mosaics


rule run_after_aggregation:
    input:
        aggregate_input,
    output:
        png15=pjoin(OUTDIR, "{seq}-15x.png"),
        png30=pjoin(OUTDIR, "{seq}-30x.png"),
    shell:
        """
        python3 ./scripts/pandata_plot.py -c 15 -g {wildcards.seq} -t rspoa-5,rspoa-9,giraffe,bwa -o {output.png15} {OUTDIR}
        python3 ./scripts/pandata_plot.py -c 30 -g {wildcards.seq} -t rspoa-5,rspoa-9,giraffe,bwa -o {output.png30} {OUTDIR}
        """


rule get_mosaic_sequence:
    input:
        fa=pjoin(OUTDIR, "{seq}", "seqs.noindels.fa"),
    output:
        fa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa"),
    shell:
        """
        samtools faidx {input.fa} {wildcards.mosaic} > {output.fa}
        """


rule simulate_reads:
    input:
        fa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa"),
    output:
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
    params:
        prefix=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads"),
    conda:
        "envs/dwgsim.yaml"
    shell:
        """
        dwgsim -e 0.01 -E 0.01 -d 1 -C {wildcards.coverage} -1 150 -2 150 -r 0 -F 0 -R 0 -X 0 -y 0 -A 1 -o 1 {input.fa} {params.prefix}
        gunzip -c {params.prefix}.bwa.read1.fastq.gz > {output.fq}
        python3 scripts/rc_sample.py {params.prefix}.bwa.read2.fastq.gz >> {output.fq}
        """


rule giraffe_index:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
    output:
        gbz=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.gbz"
        ),
        dist=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.dist"
        ),
        minim=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.min"
        ),
    params:
        snarls=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.snarls"
        ),
    conda:
        "envs/vg.yaml"
    log:
        time=pjoin(OUTDIR, "{seq}", "giraffe-index.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} vg gbwt -g {output.gbz} --gbz-format -p -G {input.gfa}
        /usr/bin/time -vao {log.time} vg snarls -T {output.gbz} > {params.snarls}
        /usr/bin/time -vao {log.time} vg index -s {params.snarls} -j {output.dist} {output.gbz}
        /usr/bin/time -vao {log.time} vg minimizer -o {output.minim} -d {output.dist} {output.gbz}
        """


rule giraffe_map:
    input:
        gbz=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.gbz"
        ),
        dist=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.dist"
        ),
        minim=pjoin(
            OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.giraffe.min"
        ),
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
    output:
        gaf=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "giraffe.gaf"),
    conda:
        "envs/vg.yaml"
    log:
        time=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "giraffe.time"),
        out=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "giraffe.log"),
    shell:
        """
        /usr/bin/time -vo {log.time} vg giraffe -Z {input.gbz} -m {input.minim} -d {input.dist} -f {input.fq} -M 1 -o GAF > {output.gaf} 2> {log.out}
        """


rule bwa_index:
    input:
        fa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa"),
    output:
        sa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa.sa"),
    log:
        time=pjoin(OUTDIR, "{seq}", "{mosaic}", "bwa-index.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} bwa index {input.fa}
        """


rule bwa_mem:
    input:
        fa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa"),
        sa=pjoin(OUTDIR, "{seq}", "{mosaic}", "sequence.fa.sa"),
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
    output:
        bam=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "bwa.bam"),
    log:
        time=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "bwa.time"),
    shell:
        """
        /usr/bin/time -vo {log.time} bwa mem {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule fq2fa:
    input:
        fq="{fpath}.fq",
    output:
        fa="{fpath}.fa",
    shell:
        """
        seqtk seq -A {input.fq} > {output.fa}
        """


rule rspoa_map:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
        fa=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fa"),
    output:
        gaf=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "rspoa-{mode}.gaf"),
    log:
        time=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "rspoa-{mode}.time"),
        out=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "rspoa-{mode}.log"),
    threads: workflow.cores / 8
    shell:
        """
        /usr/bin/time -vo {log.time} {RSPOA_BIN} -m {wildcards.mode} -R 8 -r 0.25 {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule bwa_dist:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
        bam=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "bwa.bam"),
    output:
        txt=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "bwa.dist.txt"),
    shell:
        """
        python3 ./scripts/pandata_dist.py -t bwa {input.gfa} {input.fq} {input.bam} > {output.txt}
        """


rule giraffe_dist:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
        gaf=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "giraffe.gaf"),
    output:
        txt=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "giraffe.dist.txt"),
    shell:
        """
        python3 ./scripts/pandata_dist.py -t giraffe {input.gfa} {input.fq} {input.gaf} > {output.txt}
        """


rule rspoa_dist:
    input:
        gfa=pjoin(OUTDIR, "{seq}", "graph.tsorted.wpaths.clean.minimalpcover.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "reads.fq"),
        gaf=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "rspoa-{mode}.gaf"),
    output:
        txt=pjoin(OUTDIR, "{seq}", "{mosaic}", "{coverage}", "rspoa-{mode}.dist.txt"),
    shell:
        """
        python3 ./scripts/pandata_dist.py -t rspoa {input.gfa} {input.fq} {input.gaf} > {output.txt}
        """
