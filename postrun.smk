configfile: "config.yaml"

WD = config["odir"]
Is = range(1, config["iterations"] + 1)
NUMREADS = config["numreads"]
REGION_LENS = config["region_lens"]
SIMULATORS = config["simulators"]
TOOLS = config["tools"]

rule run:
    input:
        fa=expand(
            WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.run.fa",
            sim=SIMULATORS,
            region=Is,
            NUMREAD=NUMREADS,
            REGION_LEN=REGION_LENS,
        ),
        timestable=expand(
            WD + "/results/{sim}.N{NUMREAD}.L{REGION_LEN}.cmptimes.table",
            sim=SIMULATORS,
            NUMREAD=NUMREADS,
            REGION_LEN=REGION_LENS,
        ),
        levdist=expand(
            WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.levdist",
            sim=SIMULATORS,
            region=Is,
            NUMREAD=NUMREADS,
            REGION_LEN=REGION_LENS,
        ),
        diststable=expand(
            WD + "/results/{sim}.N{NUMREAD}.L{REGION_LEN}.dists.table",
            sim=SIMULATORS,
            NUMREAD=NUMREADS,
            REGION_LEN=REGION_LENS,
        ),

rule gaf2fa:
    input:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/graph.sorted.gfa",
    output:
        fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.run.fa",
    params:
        gaf=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa/*.alignment.fa",
    shell:
        """
            python ./scripts/gaf2fa.py {input.gfa} {params.gaf} > {output.fa}
        """

rule make_tables_times:
    output:
        table=WD + "/results/{sim}.N{NUMREAD}.L{REGION_LEN}.cmptimes.table",
    params:
        abpoa=WD + "/results/*/{sim}/N{NUMREAD}.L{REGION_LEN}/abpoa.times",
        rspoa=WD + "/results/*/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.times",
    shell:
        """
            python ./scripts/stats_times.py --abpoa {params.abpoa} --rspoa {params.rspoa} > {output.table}
        """

rule calc_lev:
    input:
        fortest_fa=WD+ "/results/{region}/nanosim/N{NUMREAD}.L{REGION_LEN}/reads.fortest.fa",
        rspoa_fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.run.fa",
    output:
        levdist=WD + "/results/{region}/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.levdist",
    shell:
        """
            python ./scripts/levdist.py {input.fortest_fa} {input.rspoa_fa} > {output.levdist}
        """

rule make_tables_dist:
    output:
        table=WD + "/results/{sim}.N{NUMREAD}.L{REGION_LEN}.dists.table",
    params:
        rspoa=WD + "/results/*/{sim}/N{NUMREAD}.L{REGION_LEN}/rspoa.levdist",
    shell:
        """
            python ./scripts/stats_dists.py --rspoa {params.rspoa} > {output.table}
        """
