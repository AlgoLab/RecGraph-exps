import os.path
import re


rule all:
    input:
        expand(
            "output/cdifficile/slpa-basis.{msa}/{fa}/recgraph/{par}/full.csv",
            msa=["mafft", "gappy", "sens"],
            fa=["simulated"],
            par=[
                "mode_9-band_1-match_2-mism_4-open_4-ext_2-rec_4-recext_0.00001",
                "mode_9-band_1-match_2-mism_4-open_4-ext_2-rec_28-recext_0.00001"
            ]
        ),
        expand(
            "output/cdifficile/slpa-basis.{msa}/{fa}/jali/{mat}-{par}/full.csv",
            msa=["mafft", "gappy", "sens"],
            fa=["simulated"],
            mat=["CUSTOM_MAT"],
            par=[
                "open_4-ext_2-rec_4",
                "open_4-ext_2-rec_28",
                "open_4-ext_2-rec_48",
                "open_4-ext_0-rec_4",
                "open_4-ext_0-rec_28",
                "open_4-ext_0-rec_40",
                "open_4-ext_0-rec_48"
            ]
        )
    output:
        "output/cdifficile/full.csv"
    conda: "envs/csvkit.yaml"
    shell:
        """
        MYFI=""
        MYGR=""
        for f in {input:q}; do MYFI="$MYFI $f"; MYGR="$MYGR,${{f#output\/cdifficile\/}}"; done
        csvstack -g ${{MYGR#,}} $MYFI > {output}
        """



rule download_jali:
    output:
        "bin/jali"
    shadow: "shallow"
    shell:
        """
        wget https://bibiserv.cebitec.uni-bielefeld.de/applications/jali/resources/downloads/jali-1.3.src.tar.gz
        tar xzf jali-1.3.src.tar.gz
        cd jali1.3/
        patch < ../scripts/jali_include.patch
        make jali
        cd ..
        cp jali1.3/jali {output}
        """

rule build_recgraph:
    output:
        "bin/recgraph"
    shadow: "shallow"
    conda: "envs/rust.yaml"
    threads: 4
    shell:
        """
        git clone https://github.com/AlgoLab/RecGraph.git
        cd RecGraph/
        git checkout ef684e8fe6a9cd1218f8e418aacf195f3702e3b7
        cargo build --release --jobs {threads}
        cd ..
        cp RecGraph/target/release/recgraph {output}
        """

rule msa_auto:
    input:
        fa="data/cdifficile/slpa-basis.fa",
    output:
        msa="output/cdifficile/msa/slpa-basis.mafft.fa",
    log:
        out="output/cdifficile/msa/slpa-basis.mafft.log",
    conda: "envs/mafft.yaml"
    shell:
        """
        mafft --auto {input.fa} > {output.msa} 2> {log.out}
        """

rule msa_gappy:
    input:
        fa="data/cdifficile/slpa-basis.fa",
    output:
        msa="output/cdifficile/msa/slpa-basis.gappy.fa",
    log:
        out="output/cdifficile/msa/slpa-basis.gappy.log",
    conda: "envs/mafft.yaml"
    shell:
        """
        mafft --inputorder --anysymbol --allowshift --unalignlevel 0.8 --leavegappyregion --maxiterate 2 --retree 1 --globalpair {input.fa} > {output.msa} 2> {log.out}
        """

rule msa_sens:
    input:
        fa="data/cdifficile/slpa-basis.fa",
    output:
        msa="output/cdifficile/msa/slpa-basis.sens.fa",
    log:
        out="output/cdifficile/msa/slpa-basis.sens.log",
    conda: "envs/mafft.yaml"
    shell:
        """
        mafft --maxiterate 1000 --globalpair --op 4 --ep 2 {input.fa} > {output.msa} 2> {log.out}
        """

rule make_graph:
    input:
        recgraph="bin/recgraph",
        msa="output/cdifficile/msa/{basemsa}.fa",
    output:
        gfa="output/cdifficile/gfa/{basemsa}.sorted.wpaths.clean.gfa",
    conda: "envs/make_graph.yaml"
    threads: 4
    shadow: "shallow"
    shell:
        """
        make_prg from_msa -F -i {input.msa} -o mprg -O g -t {threads}
        python3 ./scripts/clean_gfa_from_ast.py mprg.prg.gfa > clean.gfa
        python3 ./scripts/increment_gfa_idx.py clean.gfa > clean.incr.gfa
        odgi build -g clean.incr.gfa -o clean.incr.gfa.og
        odgi sort -i clean.incr.gfa.og -o clean.incr.tsorted.gfa.og
        odgi view -i clean.incr.tsorted.gfa.og -g > sorted.gfa
        python3 ./scripts/rm_indel_msa.py {input.msa} > noindels.fa
        {input.recgraph} noindels.fa sorted.gfa > p.gaf
        cp sorted.gfa sorted.wpaths.gfa
        cut -f 1,6 p.gaf | while read idx p ; do echo -e "P\\t$idx\\t$(echo $p | cut -c 2- | sed "s/>/+,/g")+\\t*" ; done >> sorted.wpaths.gfa
        python3 ./scripts/remove_nodes_not_in_a_path.py sorted.wpaths.gfa > {output.gfa}
        """

rule run_recgraph:
    input:
        recgraph="bin/recgraph",
        fa = "output/cdifficile/split/{fa}/split.{seq}.fa",
        gfa = "output/cdifficile/gfa/{msa}.sorted.wpaths.clean.gfa",
    output:
        gaf = "output/cdifficile/{msa}/{fa}/recgraph/mode_{mod}-band_{B}-match_{match}-mism_{mism}-open_{op}-ext_{ext}-rec_{rec}-recext_{recext}/split.{seq}.gaf",
    threads: 1
    resources:
        mem="12GB"
    shell:  
        """
            {input.recgraph} \
               -m {wildcards.mod} -B {wildcards.B} \
               -O {wildcards.op} -E {wildcards.ext} \
               -R {wildcards.rec} -r {wildcards.recext} \
               -M {wildcards.match} -X {wildcards.mism} \
               {input.fa} {input.gfa} > {output.gaf}
        """

def aggregate_recgraph_input(wildcards):
    checkpoint_output = checkpoints.split_seqs.get(**wildcards).output[0]
    ret = expand(
        "output/cdifficile/{msa}/{fa}/recgraph/mode_{mod}-band_{B}-match_{match}-mism_{mism}-open_{op}-ext_{ext}-rec_{rec}-recext_{recext}/split.{seq}.gaf",
        **wildcards,
        seq=glob_wildcards(os.path.join(checkpoint_output, "split.{seq}.fa")).seq
    )
    return sorted(ret)

rule merge_recgraph:
    input:
        aggregate_recgraph_input
    output:
        "output/cdifficile/{msa}/{fa}/recgraph/mode_{mod}-band_{B}-match_{match}-mism_{mism}-open_{op}-ext_{ext}-rec_{rec}-recext_{recext}/full.gaf",
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule make_csv:
    input:
        gaf = "output/cdifficile/{msa}/{fa}/recgraph/{par}/full.gaf",
    output:
        csv = "output/cdifficile/{msa}/{fa}/recgraph/{par}/full.csv",
    run:
        out = open(output.csv, 'w+')
        print("ReadName,RecPaths,RecPos,Score,Displacement", file=out)
        for line in open(input.gaf):
            line = line.strip()
            rid=line.split('\t')[0]
            m = re.search(r'recombination path (\d+) (\d+)', line)
            if m:
                recpath=f"{int(m.group(1))+1}>{int(m.group(2))+1}"
                bp=line.split('\t')[-1]
            else:
                m = re.search(r'best path: (\d+)', line)
                if m:
                    recpath=f"{int(m.group(1))+1}"
                else:
                    recpath=''
                bp=''
            m = re.search(r'score: (\d+)', line)
            if m:
                score=f"{m.group(1)}"
            else:
                score='-'
            m = re.search(r'displacement: (\d+)', line)
            if m:
                displacement=f"{m.group(1)}"
            else:
                displacement=''

            print(rid, recpath, bp, score, displacement, sep=',', file=out)
        out.close()

checkpoint split_seqs:
    input:
        fa = "data/cdifficile/{fa}.fa",
    output:
        dir = directory("output/cdifficile/split/{fa}"),
    shell:
        """
        mkdir -p {output.dir}
        awk 'BEGIN {{ NUM = -1 }} /^>/ {{NUM = NUM+1; F=sprintf("{output.dir}/split.%.4d.fa", NUM); print > F; next; }} {{print >> F;}}' < {input.fa}
        """


rule jali:
    input:
        jali = "bin/jali",
        msa = "output/cdifficile/msa/{msa}.fa",
        fa = "output/cdifficile/split/{fa}/split.{seq}.fa",
        mat = "data/cdifficile/jali/{mat}",
    output:
        "output/cdifficile/{msa}/{fa}/jali/{mat}-open_{op}-ext_{ext}-rec_{rec}/split.{seq}.out",
    shadow: "minimal"
    threads: 1
    resources:
        mem="2GB"
    shell:
        """
        {input.jali} -p -w {input.mat} -i -{wildcards.op} -e -{wildcards.ext} -j -{wildcards.rec} {input.fa} {input.msa} > {output}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_seqs.get(**wildcards).output[0]
    ret = expand(
        "output/cdifficile/{msa}/{fa}/jali/{mat}-open_{op}-ext_{ext}-rec_{rec}/split.{seq}.out",
        msa=wildcards.msa,
        fa=wildcards.fa,
        mat=wildcards.mat,
        op=wildcards.op,
        ext=wildcards.ext,
        rec=wildcards.rec,
        seq=glob_wildcards(os.path.join(checkpoint_output, "split.{seq}.fa")).seq
    )
    return ret


rule aggregate_jali:
    input:
        aggregate_input
    output:
        txt = "output/cdifficile/{msa}/{fa}/jali/{mat}-open_{op}-ext_{ext}-rec_{rec}/full.csv",
    run:
        out = open(output.txt, 'w')
        print("ReadName,RecPaths,RecPos,Score,Displacement", file=out)
        for f in sorted(input):
            prev = None
            seqs = []
            nextline = False
            seqname = None
            for line in open(f):
                line = line.rstrip()
                if nextline and seqname is None:
                    seqname = line.split(" ")[0]
                    nextline = False
                if line.startswith(" "):
                    prev = line
                if line.startswith("SequenceNo"):
                    nums = None
                    line = line[13:]
                    if prev is None:
                        nums = [int(c) for c in line]
                    else:
                        prev = prev[13:]
                        nums = [10*int(0 if d == " " else d) + int(c) for d,c in zip(prev,line)]
                        prev = None
                    seqs += nums
                    nextline = True
                if line.startswith("The optimal alignment score is "):
                    m = re.search(r'The optimal alignment score is (\d+) using (\d+)', line)
                    if m:
                        score=f"{m.group(1)}"
                        nseq=f"{m.group(2)}"
                    else:
                        score='-'
                        nseq='-'
            pos = 0
            recs = []
            outp = [str(seqs[0])]
            for i,j in zip(seqs, seqs[1:]):
                if i != j:
                    recs.append(str(pos))
                    outp.append(str(j))
                pos += 1
            print(seqname, ">".join(outp), "-".join(recs), score, nseq, sep=",", file=out)
        out.close()
