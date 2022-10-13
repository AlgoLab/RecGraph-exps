from os.path import join as pjoin
from glob import glob


configfile: "config.yaml"


SEQSDIR = config["seqsdir"]
OUTDIR = config["odir"]
RSPOA_BIN = config["rspoa"]

SEQS = config["seqs"]

COVS = ["15"]
Ms = [2]  # , 4, 5]
Xs = [4]  # , 6, 10]
Os = [4]  # , 6, 10]
Es = [2]  # , 3, 4]
Rs = [4]  # [8]  # 2, 3, 4, 7
rs = [0.1]  # 1, 4

rspoa_output = []
giraffe_output = []
bwa_output = []
for fpath in glob(pjoin(SEQSDIR, "*", "1-150", "*", "sequence.fa")):
    seq = fpath.split("/")[-4]
    mosaic = fpath.split("/")[-2]
    if seq in SEQS:
        prefix = pjoin(OUTDIR, seq, "1-150", mosaic)
        for c in COVS:
            bwa_output.append(pjoin(prefix, c, "bwa.dist.txt"))
            giraffe_output.append(pjoin(prefix, c, "giraffe.dist.txt"))
        for M in Ms:
            for X in Xs:
                for O in Os:
                    for E in Es:
                        for c in COVS:
                            rspoa_output.append(
                                pjoin(
                                    OUTDIR,
                                    seq,
                                    "1-150",
                                    mosaic,
                                    c,
                                    f"rspoa-5.M{M}-X{X}-O{O}-E{E}.dist.txt",
                                )
                            )
        for R in Rs:
            for r in rs:
                for c in COVS:
                    rspoa_output.append(
                        pjoin(
                            OUTDIR,
                            seq,
                            "1-150",
                            mosaic,
                            c,
                            f"rspoa-9.R{R}-r{r}.dist.txt",
                        )
                    )


rule run:
    input:
        rspoa_output,
        giraffe_output,
        bwa_output,


rule simulate_reads:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "{recpar}", "{mosaic}", "sequence.fa"),
    output:
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
    params:
        prefix=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads"),
    conda:
        "envs/dwgsim.yaml"
    shell:
        """
        dwgsim -e 0.01 -E 0.01 -d 150 -C {wildcards.coverage} -1 150 -2 150 -r 0 -F 0 -R 0 -X 0 -y 0 -A 1 -o 1 {input.fa} {params.prefix}
        gunzip -c {params.prefix}.bwa.read1.fastq.gz > {output.fq}
        python3 ./scripts/revcompl_sample.py {params.prefix}.bwa.read2.fastq.gz >> {output.fq}
        """


rule bwa_index:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "{recpar}", "{mosaic}", "sequence.fa"),
    output:
        sa=pjoin(SEQSDIR, "{seq}", "{recpar}", "{mosaic}", "sequence.fa.sa"),
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa index {input.fa}
        """


rule bwa_mem:
    input:
        fa=pjoin(SEQSDIR, "{seq}", "{recpar}", "{mosaic}", "sequence.fa"),
        sa=pjoin(SEQSDIR, "{seq}", "{recpar}", "{mosaic}", "sequence.fa.sa"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
    output:
        bam=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "bwa.bam"),
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem {input.fa} {input.fq} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule fq2fa:
    input:
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
    output:
        fa=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fa"),
    conda:
        "envs/seqtk.yaml"
    shell:
        """
        seqtk seq -A {input.fq} > {output.fa}
        """


rule rspoa5_map:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fa=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fa"),
    output:
        gaf=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-5.M{M}-X{X}-O{O}-E{E}.gaf",
        ),
    log:
        time=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-5.M{M}-X{X}-O{O}-E{E}.time",
        ),
        out=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-5.M{M}-X{X}-O{O}-E{E}.log",
        ),
    threads: 1  # workflow.cores / 8
    shell:
        """
        /usr/bin/time -vo {log.time} {RSPOA_BIN} -m 5 -M {wildcards.M} -X {wildcards.X} -O {wildcards.O} -E {wildcards.E} {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule rspoa9_map:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fa=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fa"),
    output:
        gaf=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-9.R{R}-r{r}.gaf",
        ),
    log:
        time=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-9.R{R}-r{r}.time",
        ),
        out=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-9.R{R}-r{r}.log",
        ),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} {RSPOA_BIN} -m 9 -B 0.8 -R {wildcards.R} -r {wildcards.r} {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule giraffe_index:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
    output:
        gbz=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.gbz"),
        dist=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.dist"),
        minim=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.min"),
    params:
        snarls=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.snarls"),
    log:
        time=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.time"),
    conda:
        "envs/vg.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} vg gbwt -g {output.gbz} --gbz-format -p -G {input.gfa}
        /usr/bin/time -vao {log.time} vg snarls -T {output.gbz} > {params.snarls}
        /usr/bin/time -vao {log.time} vg index -s {params.snarls} -j {output.dist} {output.gbz}
        /usr/bin/time -vao {log.time} vg minimizer -o {output.minim} -d {output.dist} {output.gbz}
        """


rule giraffe_map:
    input:
        gbz=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.gbz"),
        dist=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.dist"),
        minim=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.giraffe.min"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
    output:
        gaf=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "giraffe.gaf"),
    log:
        time=pjoin(
            OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "giraffe.time"
        ),
    conda:
        "envs/vg.yaml"
    shell:
        """
        /usr/bin/time -vo {log.time} vg giraffe -Z {input.gbz} -m {input.minim} -d {input.dist} -f {input.fq} -M 1 -o GAF > {output.gaf}
        """


rule bwa_dist:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
        bam=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "bwa.bam"),
    output:
        txt=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "bwa.dist.txt"),
    conda:
        "envs/ed.yaml"
    shell:
        """
        python3 ./scripts/compute_ed.py -f bam {input.gfa} {input.fq} {input.bam} > {output.txt}
        """


rule rspoa5_dist:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
        gaf=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-5.M{M}-X{X}-O{O}-E{E}.gaf",
        ),
    output:
        txt=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-5.M{M}-X{X}-O{O}-E{E}.dist.txt",
        ),
    conda:
        "envs/ed.yaml"
    shell:
        """
        python3 ./scripts/compute_ed.py -f gaf {input.gfa} {input.fq} {input.gaf} > {output.txt}
        """


rule rspoa9_dist:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
        gaf=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-9.R{R}-r{r}.gaf",
        ),
    output:
        txt=pjoin(
            OUTDIR,
            "{seq}",
            "{recpar}",
            "{mosaic}",
            "{coverage}",
            "rspoa-9.R{R}-r{r}.dist.txt",
        ),
    conda:
        "envs/ed.yaml"
    shell:
        """
        python3 ./scripts/compute_ed.py -f gaf {input.gfa} {input.fq} {input.gaf} > {output.txt}
        """


rule giraffe_dist:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "{recpar}", "graph.gfa"),
        fq=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "reads.fq"),
        gaf=pjoin(OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "giraffe.gaf"),
    output:
        txt=pjoin(
            OUTDIR, "{seq}", "{recpar}", "{mosaic}", "{coverage}", "giraffe.dist.txt"
        ),
    conda:
        "envs/ed.yaml"
    shell:
        """
        python3 ./scripts/compute_ed.py -f gaf {input.gfa} {input.fq} {input.gaf} > {output.txt}
        """
