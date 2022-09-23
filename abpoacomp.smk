configfile: "config.yaml"


FASTACHR = config["seq"]
WD = config["odir"]
ABPOA_BIN = config["abpoa"]  # conda version doesn't work on our server
RSPOA_BIN = config["rspoa"]
Is = range(1, config["iterations"] + 1)
NUMREADS = config["numreads"]
REGION_LENS = config["region_lens"]
SIMULATORS = config["simulators"]
TOOLS = config["tools"]


rule run:
    input:
        times=expand(
            WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/{tools}.times",
            sim=SIMULATORS,
            tools=TOOLS,
            region=Is,
            NUMREAD=NUMREADS,
            REGION_LEN=REGION_LENS,
        ),


rule exctract_rand_region:
    input:
        fa=FASTACHR,
    output:
        extracted=WD + "/results/{region}/region-{REGION_LEN}.fa",
    params:
        reglen=lambda wildcards: int(int(wildcards.REGION_LEN) * 1.1),
    run:
        shell(
            "python ./scripts/extract_random_region.py {FASTACHR} {params.reglen} > {output.extracted}"
        )


rule run_pbsim:
    input:
        region=WD + "/results/{region}/region-{REGION_LEN}.fa",
    params:
        model="/home/denti/software/pbsim2/data/P6C4.model",
        lenmax=lambda wildcards: int(int(wildcards.REGION_LEN) * 1.05),
        lenmin=lambda wildcards: int(int(wildcards.REGION_LEN) * 0.95),
        lenmean=lambda wildcards: wildcards.REGION_LEN,
        depth=lambda wildcards: wildcards.NUMREAD,  # this since region_len ~ avg_read_len
        accmax=0.9,
        accmin=0.8,
        accmean=0.85,
        prefix=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads",
    output:
        fq=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads_0001.fastq",
    conda:
        "envs/pbsim2.yaml"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        pbsim {input.region} \\
            --hmm_model {params.model} \\
            --length-max {params.lenmax} --length-min {params.lenmin} --length-mean {params.lenmean} \\
            --depth {params.depth} \\
            --accuracy-max {params.accmax} --accuracy-min {params.accmin} --accuracy-mean {params.accmean} \\
            --prefix {params.prefix}
        """


rule run_nano:
    input:
        region=WD + "/results/{region}/region-{REGION_LEN}.fa",
    params:
        c="/home/denti/software/NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy/training",
        lenmax=lambda wildcards: int(int(wildcards.REGION_LEN) * 1.05),
        lenmin=lambda wildcards: int(int(wildcards.REGION_LEN) * 0.95),
        prefix=WD + "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads",
    output:
        fa=WD
        + "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads_aligned_reads.fasta",
    threads: workflow.cores / 4
    conda:
        "envs/nanosim.yaml"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        simulator.py genome -rg {input.region} -n {wildcards.NUMREAD} \\
            -max {params.lenmax} \\
            -min {params.lenmin} \\
            -b guppy -s 0 \\
            -c {params.c} -o {params.prefix}
        """


rule fq_to_fa:
    input:
        fq=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads_0001.fastq",
    output:
        fa=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads_0001.fa",
    shell:
        """
        seqtk seq -A {input.fq} > {output.fa}
        """


rule split_reads_pbsim:
    input:
        fa=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads_0001.fa",
    output:
        fa_gfa=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads.forgfa.fa",
        fa_test=WD + "/results/{region}/pbsim/N{NUMREAD}.L{REGION_LEN}/reads.fortest.fa",
    params:
        split_for_gfa=0.3,
    run:
        totreads = sum(1 for _ in open(input.fa)) / 2
        gfa_split = int(totreads * params.split_for_gfa)
        test_split = int(totreads - gfa_split)
        gfa_split *= 2
        test_split *= 2

        shell("head -n {gfa_split} {input.fa} > {output.fa_gfa}")
        shell("tail -n {test_split} {input.fa} > {output.fa_test}")


rule split_reads_nano:
    input:
        fa=WD
        + "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads_aligned_reads.fasta",
    output:
        fa_gfa=WD + "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads.forgfa.fa",
        fa_test=WD
        + "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads.fortest.fa",
    params:
        split_for_gfa=0.3,
    run:
        totreads = sum(1 for _ in open(input.fa)) / 2
        gfa_split = int(totreads * params.split_for_gfa)
        test_split = int(totreads - gfa_split)
        gfa_split *= 2
        test_split *= 2

        shell("head -n {gfa_split} {input.fa} > {output.fa_gfa}")
        shell("tail -n {test_split} {input.fa} > {output.fa_test}")


rule build_gfa:
    input:
        fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/reads.forgfa.fa",
    output:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.gfa",
    threads: workflow.cores / 2
    shell:
        """
        /home/denti/software/abPOA-v1.4.1/bin/abpoa -r 3 {input.fa} > {output.gfa}
        """


rule sort_gfa:
    input:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.gfa",
    output:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.sorted.gfa",
    params:
        og=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.TMPODGI",
    conda:
        "envs/odgi.yaml"
    threads: workflow.cores / 2
    shell:
        """
        odgi build -g {input.gfa} -o {params.og}.og
        odgi sort -i {params.og}.og -o {params.og}.sorted.og
        odgi view -i {params.og}.sorted.og -g > {output.gfa}
        """


rule run_abpoa:
    input:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.gfa",
        fa_test=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/reads.fortest.fa",
    output:
        log=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/abpoa.times",
    params:
        wd=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/abpoa",
    threads: workflow.cores
    shell:
        """
        bash ./scripts/run_abpoa.sh {input.gfa} {input.fa_test} {params.wd} 2 {ABPOA_BIN} > {output.log}
        """


rule run_rspoa:
    input:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.sorted.gfa",
        fa_test=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/reads.fortest.fa",
    output:
        log=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.times",
    params:
        wd=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa",
    threads: workflow.cores
    shell:
        """
        bash ./scripts/run_rspoa.sh {input.gfa} {input.fa_test} {params.wd} 2 {RSPOA_BIN} > {output.log}
        """
