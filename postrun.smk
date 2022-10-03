configfile: "abpoacomp-config.yaml"

WD = config["odir"]
Is = range(1, config["iterations"] + 1)
NUMREADS = config["numreads"]
READ_LENS = config["read_lens"]
SIMULATORS = config["simulators"]
TOOLS = config["tools"]

rule tables:
    input:
        diststable=expand(
            WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.dists.table",
            sim=SIMULATORS,
            NUMREAD=NUMREADS,
            READ_LEN=READ_LENS,
        ),
        timestable=expand(
            WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.cmptimes.table",
            sim=SIMULATORS,
            NUMREAD=NUMREADS,
            READ_LEN=READ_LENS,
        ),

rule calc:
    input:
        levdist=expand(
            WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/{TOOLS}.levdist",
            sim=SIMULATORS,
            region=Is,
            NUMREAD=NUMREADS,
            READ_LEN=READ_LENS,
            TOOLS=["rspoa"],
        ),
        percalns=expand(
            WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.percalns",
            sim=SIMULATORS,
            region=Is,
            NUMREAD=NUMREADS,
            READ_LEN=READ_LENS,
        ),
    output: touch(WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.calc.flag")


rule clean:
    shell:
        """
            rm /data/abpoa-comparison/results/*/*/N75.*/*.levdist
            rm /data/abpoa-comparison/results/*/*/N75.*/rspoa.run.fa
            rm /data/abpoa-comparison/results/*.table
            rm /data/abpoa-comparison/results/*/*/N75.*/rspoa.percalns
        """

rule gaf2fa:
    input:
        gfa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/graph.sorted.reversed.gfa",
    output:
        fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.run.fa",
    params:
        gaf=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa/*.alignment.fa",
    shell:
        """
            python ./scripts/gaf2fa.py {input.gfa} {params.gaf} > {output.fa}
        """

rule make_tables_times:
    output:
        table=WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.cmptimes.table",
    params:
        abpoa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/abpoa.times",
        rspoa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.times",
    shell:
        """
            python ./scripts/stats_times.py --abpoa {params.abpoa} --rspoa {params.rspoa} > {output.table}
        """

rule cat_abpoa_consensus: # TODO: remember abPOA use NOT SORTED graph
    output:
        abpoa_fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/abpoa.run.fa",
    params:
        abpoa_dir=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/abpoa",
    shell:
        """
            find {params.abpoa_dir} -name "*.consensus.fa" | sort -V | xargs cat > {output.abpoa_fa}
        """

rule calc_lev:
    input:
        fortest_fa=WD+ "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/reads.fortest.fa",
        fa=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/{TOOLS}.run.fa",
    output:
        levdist=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/{TOOLS}.levdist",
    shell:
        """
            python ./scripts/levdist.py {input.fortest_fa} {input.fa} > {output.levdist}
        """

rule calc_perc_alns:
    params:
        gaf=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa/*.alignment.fa",
    output:
        percaln=WD + "/results/{region}/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.percalns",
    shell:
        """
            python ./scripts/gaf2percaln.py {params.gaf} > {output.percaln}
        """

rule make_tables_dist:
    input:
        flag=WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.calc.flag"
    output:
        table=WD + "/results/{sim}.N{NUMREAD}.L{READ_LEN}.dists.table",
    params:
        # abpoa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/abpoa.levdist",
        rspoa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.levdist",
        rspoa_pa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.percalns",
    shell:
        """
            python ./scripts/stats_dists.py --rspoa {params.rspoa} --rspercal {params.rspoa_pa} > {output.table}
        """
