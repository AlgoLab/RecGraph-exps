from os.path import join as pjoin
import glob

SEQSDIR = config["seqsdir"]
RECGRAPH_BIN = (
    config["recgraph"]
    if "recgraph" in config
    else "/data/RecGraph/target/release/recgraph"
)

M = 2
X = 4
O = 4
E = 2
R = 4
B = 0.9
r = 0.1

rg0_output = []
rg4_output = []
rg8_output = []
# giraffe_output = []
ga_output = []
for fpath in glob.glob(pjoin(SEQSDIR, "*", "MPCSIM", "*", "p*", "sequence.fa")):
    seq = fpath.split("/")[-5]
    mosaic = fpath.split("/")[-3]
    p = fpath.split("/")[-2][1:]

    ga_output.append(pjoin(SEQSDIR, seq, "MPCSIM", mosaic, f"p{p}", "graphaligner.gaf"))
    # giraffe_output.append(pjoin(SEQSDIR, seq, "MPCSIM", mosaic, "giraffe.gaf"))
    rg0_output.append(
        pjoin(
            SEQSDIR,
            seq,
            "MPCSIM",
            mosaic,
            f"p{p}",
            f"recgraph-0.M{M}-X{X}-O{O}-E{E}.gaf",
        )
    )
    rg4_output.append(
        pjoin(
            SEQSDIR,
            seq,
            "MPCSIM",
            mosaic,
            f"p{p}",
            f"recgraph-4.M{M}-X{X}-O{O}-E{E}.gaf",
        )
    )
    rg8_output.append(
        pjoin(
            SEQSDIR,
            seq,
            "MPCSIM",
            mosaic,
            f"p{p}",
            f"recgraph-8.R{R}-r{r}.gaf",
        )
    )


rule run:
    input:
        rg0_output,
        rg4_output,
        rg8_output,
        ga_output,
    output:
        (SEQSDIR[:-1] if SEQSDIR.endswith("/") else SEQSDIR) + ".results.txt",
    shell:
        """
        python3 scripts/check.py {SEQSDIR} > {output}
        """


rule rg0_map:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
        fa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "p{p}", "sequence.fa"),
    output:
        gaf=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-0.M{M}-X{X}-O{O}-E{E}.gaf",
        ),
    log:
        time=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-0.M{M}-X{X}-O{O}-E{E}.time",
        ),
        out=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-0.M{M}-X{X}-O{O}-E{E}.log",
        ),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} {RECGRAPH_BIN} -m 0  -M {wildcards.M} -X {wildcards.X} -O {wildcards.O} -E {wildcards.E} {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule rg4_map:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
        fa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "p{p}", "sequence.fa"),
    output:
        gaf=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-4.M{M}-X{X}-O{O}-E{E}.gaf",
        ),
    log:
        time=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-4.M{M}-X{X}-O{O}-E{E}.time",
        ),
        out=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-4.M{M}-X{X}-O{O}-E{E}.log",
        ),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} {RECGRAPH_BIN} -m 4 -M {wildcards.M} -X {wildcards.X} -O {wildcards.O} -E {wildcards.E} {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule rg8_map:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
        fa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "p{p}", "sequence.fa"),
    output:
        gaf=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-8.R{R}-r{r}.gaf",
        ),
    log:
        time=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-8.R{R}-r{r}.time",
        ),
        out=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "recgraph-8.R{R}-r{r}.log",
        ),
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time} {RECGRAPH_BIN} -m 8 -R {wildcards.R} -r {wildcards.r} {input.fa} {input.gfa} > {output.gaf} 2> {log.out}
        """


rule graphaligner:
    input:
        gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
        fa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "p{p}", "sequence.fa"),
    output:
        gaf=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "graphaligner.gaf",
        ),
    log:
        time=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "p{p}", "graphaligner.time"),
        out=pjoin(
            SEQSDIR,
            "{seq}",
            "MPCSIM",
            "{mosaic}",
            "p{p}",
            "graphaligner.log",
        ),
    conda:
        "envs/graphaligner.yaml"
    threads: 1
    shell:
        """
        /usr/bin/time -vo {log.time}  GraphAligner -g {input.gfa} -f {input.fa} -t 1 -a {output.gaf} -x vg &> {log.out}
        """


# ###############
# ### GIRAFFE ###
# ###############
# rule giraffe_index:
#     input:
#         gfa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.gfa"),
#     output:
#         gbz=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.gbz"),
#         dist=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.dist"),
#         minim=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.min"),
#     params:
#         snarls=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.snarls"),
#     log:
#         time=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.time"),
#     conda:
#         "envs/vg.yaml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} vg gbwt -g {output.gbz} --gbz-format -p -G {input.gfa}
#         /usr/bin/time -vao {log.time} vg snarls -T {output.gbz} > {params.snarls}
#         /usr/bin/time -vao {log.time} vg index -s {params.snarls} -j {output.dist} {output.gbz}
#         /usr/bin/time -vao {log.time} vg minimizer -o {output.minim} -d {output.dist} {output.gbz}
#         """
# rule giraffe_map:
#     input:
#         gbz=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.gbz"),
#         dist=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.dist"),
#         minim=pjoin(SEQSDIR, "{seq}", "MPCSIM", "graph.giraffe.min"),
#         fa=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "sequence.fa"),
#     output:
#         gaf=pjoin(SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "giraffe.gaf"),
#     log:
#         time=pjoin(
#             SEQSDIR, "{seq}", "MPCSIM", "{mosaic}", "giraffe.time"
#         ),
#     conda:
#         "envs/vg.yaml"
#     shell:
#         """
#         /usr/bin/time -vo {log.time} vg giraffe -Z {input.gbz} -m {input.minim} -d {input.dist} -f {input.fa} -o GAF > {output.gaf}
#         """
