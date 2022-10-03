configfile: "abpoacomp-msa-config.yaml"

from os.path import join as pjoin
from Levenshtein import distance as lev
import glob

SAMPLES = config["samples"]
WD = config["odir"]
ABPOA_BIN = config["abpoa"]  # conda version doesn't work on our server
RSPOA_BIN = config["rspoa"]
TOOLS = config["tools"]
SAMPLES_DIR = config["sampledir"]


rule run:
    input:
        times=expand(
            WD + "/results/{sample}/{tool}.times",
            sample=SAMPLES,
            tool=TOOLS,
        ),
        out=expand(
            WD + "/results/{sample}/{tool}.levdist",
            sample=SAMPLES,
            tool=TOOLS,
        ),
        tables=expand(
            WD + "/results/{sample}/{table}.table",
            sample=SAMPLES,
            table=["cmptimes", "cmpdists"]
        ),
        # table_dist=expand(
        #     WD + "/results/{sample}/cmpdists.table",
        #     sample=SAMPLES,
        # ),
        

rule split_reads:
    input:
        fa=pjoin(SAMPLES_DIR, "{sample}") + '.fa',
    output:
        fa_gfa=WD + "/results/{sample}/reads.forgfa.fa",
        fa_test=WD + "/results/{sample}/reads.fortest.fa",
    params:
        split_for_gfa=0.33,
        prefix=WD + "/results/{sample}/reads",
    shell:
        """
        python scripts/get_rand_fas.py -i {input.fa} -n $(grep -c ">" {input.fa}) -p {params.split_for_gfa} -o {params.prefix}
        """

rule build_gfa:
    input:
        fa=WD + "/results/{sample}/reads.forgfa.fa",
    output:
        gfa=WD + "/results/{sample}/graph.gfa",
    threads: workflow.cores / 2
    shell:
        """
        /home/denti/software/abPOA-v1.4.1/bin/abpoa -m 0 -r 3 {input.fa} > {output.gfa}
        """

rule unchop_gfa:
    input:
        gfa=WD + "/results/{sample}/graph.gfa",
    output:
        gfa=WD + "/results/{sample}/graph.unchopped.gfa",
    conda:
        "envs/pggb.yaml"
    shell:
        """
        vg mod -u {input.gfa} > {output.gfa}
        """

rule sort_gfa:
    input:
        gfa=WD + "/results/{sample}/graph.unchopped.gfa",
    output:
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.gfa",
    params:
        og=WD + "/results/{sample}/graph.TMPODGI",
    conda:
        "envs/odgi.yaml"
    threads: workflow.cores / 2
    shell:
        """
        odgi build -g {input.gfa} -o {params.og}.og
        odgi sort -i {params.og}.og -o {params.og}.sorted.og
        odgi view -i {params.og}.sorted.og -g > {output.gfa}
        """

rule reverse_paths:
    input:
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.gfa",
    output:
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.reversed.gfa",
    shell:
        """
        python ./scripts/reverse_paths.py {input.gfa} > {output.gfa}
        """

rule run_abpoa:
    input:
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.reversed.gfa",
        fa_test=WD + "/results/{sample}/reads.fortest.fa",
    output:
        log=WD + "/results/{sample}/abpoa.times",
    params:
        wd=WD + "/results/{sample}/abpoa",
    threads: workflow.cores
    shell:
        """
        bash ./scripts/run_abpoa.sh {input.gfa} {input.fa_test} {params.wd} 2 {ABPOA_BIN} > {output.log}
        """


rule run_rspoa:
    input:
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.reversed.gfa",
        fa_test=WD + "/results/{sample}/reads.fortest.fa",
    output:
        log=WD + "/results/{sample}/rspoa.times",
    params:
        wd=WD + "/results/{sample}/rspoa",
    threads: workflow.cores
    shell:
        """
        bash ./scripts/run_rspoa.sh {input.gfa} {input.fa_test} {params.wd} 2 {RSPOA_BIN} > {output.log}
        """

rule gaf2fa:
    input:
        log=WD + "/results/{sample}/rspoa.times",
        gfa=WD + "/results/{sample}/graph.unchopped.sorted.reversed.gfa",
    output:
        fa=WD + "/results/{sample}/rspoa.run.fa",
    params:
        gaf=WD + "/results/{sample}/rspoa/*.alignment.fa",
    shell:
        """
            python ./scripts/gaf2fa.py {input.gfa} {params.gaf} > {output.fa}
        """

rule calc_lev_rspoa:
    input:
        fortest_fa=WD + "/results/{sample}/reads.fortest.fa",
        fa=WD + "/results/{sample}/rspoa.run.fa",
    output:
        levdist=WD + "/results/{sample}/rspoa.levdist",
    shell:
        """
            python ./scripts/levdist.py {input.fortest_fa} {input.fa} > {output.levdist}
        """

rule calc_lev_abpoa:
    input:
        log=WD + "/results/{sample}/abpoa.times",
    output:
        levdist=WD + "/results/{sample}/abpoa.levdist",
    params:
        cons=WD + "/results/{sample}/abpoa/*.consensus.msa.fa",
    run:
        out = open(output.levdist, 'w+')
        cons = glob.glob(params.cons)
        for f in cons:
            ll = []
            with open(f, 'r') as fin:
                ll = fin.readlines()
            r = ll[-3]
            c = ll[-1]
            out.write(f'{lev(r,c)}\n')
        out.close()

rule make_tables_times:
    input:
        abpoa=WD + "/results/{sample}/abpoa.times",
        rspoa=WD + "/results/{sample}/rspoa.times",
    output:
        table=WD + "/results/{sample}/cmptimes.table",
    shell:
        """
            python ./scripts/stats_times.py --abpoa {input.abpoa} --rspoa {input.rspoa} > {output.table}
        """

rule make_tables_dist:
    input:
        abpoa=WD + "/results/{sample}/abpoa.levdist",
        rspoa=WD + "/results/{sample}/rspoa.levdist",
        # rspoa_pa=WD + "/results/*/{sim}/N{NUMREAD}.L{READ_LEN}/rspoa.percalns",
    output:
        table=WD + "/results/{sample}/cmpdists.table",
    shell:
        """
            python ./scripts/stats_dists.py --abpoa {input.abpoa} --rspoa {input.rspoa} > {output.table}
        """
