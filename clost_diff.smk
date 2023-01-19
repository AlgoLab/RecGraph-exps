from os.path import join as pjoin
import re

RECGRAPH_BIN = (
    config["recgraph"]
    if "recgraph" in config
    else "/data/RecGraph/target/release/recgraph"
)
MAKEPRG_BIN = (
    config["mkprg"] if "mkprg" in config else "make_prg"
)

OUTDIR="data/cdifficile/results"

# rule all:
#     input:
#         expand(pjoin(OUTDIR,"{fa}.m{m}.B{B}.gaf"), fa=FAS, m=[9], B=[1, 0.9, 0.8])

rule all:
    input:
        expand(pjoin(OUTDIR, "simulated.m{m}.B{B}.csv"), m=[9], B=[0.95]),
        expand(pjoin(OUTDIR, "real.m{m}.B{B}.csv"), m=[9], B=[0.8])


rule make_graph:
    input:
        msa="data/cdifficile/slpa-basis.mafft.fa",
    output:
        gfa="data/cdifficile/slpa-basis.sorted.wpaths.clean.gfa",
    params:
        prefix="data/cdifficile/slpa-basis",
    shell:
        """
        {MAKEPRG_BIN} from_msa -F -i {input.msa} -o {params.prefix}.mprg -O g -t 16
        python3 ./scripts/clean_gfa_from_ast.py {params.prefix}.mprg.prg.gfa > {params.prefix}.clean.gfa
        python3 ./scripts/increment_gfa_idx.py {params.prefix}.clean.gfa > {params.prefix}.clean.incr.gfa
        odgi build -g {params.prefix}.clean.incr.gfa -o {params.prefix}.clean.incr.gfa.og
        odgi sort -i {params.prefix}.clean.incr.gfa.og -o {params.prefix}.clean.incr.tsorted.gfa.og
        odgi view -i {params.prefix}.clean.incr.tsorted.gfa.og -g > {params.prefix}.sorted.gfa
        python3 ./scripts/rm_indel_msa.py {input.msa} > {params.prefix}.noindels.fa
        {RECGRAPH_BIN} {params.prefix}.noindels.fa {params.prefix}.sorted.gfa > {params.prefix}.p.gaf
        cp {params.prefix}.sorted.gfa {params.prefix}.sorted.wpaths.gfa
        cut -f 1,6 {params.prefix}.p.gaf | while read idx p ; do echo -e "P\\t$idx\\t$(echo $p | cut -c 2- | sed "s/>/+,/g")+\\t*" ; done >> {params.prefix}.sorted.wpaths.gfa
        python3 ./scripts/remove_nodes_not_in_a_path.py {params.prefix}.sorted.wpaths.gfa > {params.prefix}.sorted.wpaths.clean.gfa
        rm {params.prefix}.clean.incr.* {params.prefix}.p.gaf {params.prefix}.mprg.* {params.prefix}.sorted.wpaths.gfa {params.prefix}.sorted.gfa {params.prefix}.noindels.fa {params.prefix}.clean.gfa
        """

rule run_recgraph:
    input:
        fa = pjoin("data/cdifficile/", "{fa}.fa"),
        gfa = "data/cdifficile/slpa-basis.sorted.wpaths.clean.gfa",
    output:
        gaf = pjoin(OUTDIR, "{fa}.m{m}.B{B}.gaf"),
    shell:  
        """
            {RECGRAPH_BIN} -m {wildcards.m} -B {wildcards.B} {input.fa} {input.gfa} > {output.gaf}
        """

rule make_csv:
    input:
        gaf = pjoin(OUTDIR, "{fa}.m{m}.B{B}.gaf"),
    output:
        csv = pjoin(OUTDIR, "{fa}.m{m}.B{B}.csv"),
    run:
        out = open(output.csv, 'w+')
        print("ReadName,RecPaths,RecPos", file=out)
        for line in open(input.gaf):
            line = line.strip()
            rid=line.split('\t')[0]
            m = re.search(r'recombination path (\d+) (\d+)', line)
            if m:
                recpath=f"{m.group(1)}>{m.group(2)}"
                bp=line.split('\t')[-1]
            else:
                recpath=bp='.'
            print(rid, recpath, bp, sep=',', file=out)
        out.close()