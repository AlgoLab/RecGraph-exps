import argparse
from Bio import SeqIO
import pysam
from Levenshtein import distance as lev
from gaf2fa import load_gfa_nodes


def convert_se_giraffe(nodes, path, s, e):
    r = ""
    if ">" in path:
        path = path.split(">")[1:]
        for n in path:
            r += nodes[n]  # if n in nodes else ""
    else:
        print(path)
        raise NotImplementedError(
            "Path on reverse strand is not implemented yet.\
            Not clear how to read it."
        )
    return r[s : e + 1]


def convert_se_rspoa(nodes, path, s, e):
    r = ""
    # used_nodes = set()
    if ">" in path:
        path = path.split(">")[1:]
        r += nodes[path[0]][s:]  # if path[0] in nodes else ""
        # used_nodes.add(path[0])
        for n in path[1:-1]:
            # if n not in used_nodes:
            r += nodes[n]  # if n in nodes else ""
            # used_nodes.add(n)
        # if path[-1] not in used_nodes:
        r += nodes[path[-1]][: e + 1]  # if path[-1] in nodes else ""
    else:
        print(path)
        raise NotImplementedError(
            "Path on reverse strand is not implemented yet.\
            Not clear how to read it."
        )
    return r


def main():
    parser = argparse.ArgumentParser(description="Produce edit distance.")
    parser.add_argument(
        "graph",
        type=str,
        help="Path to graph in .gfa",
    )
    parser.add_argument(
        "sample",
        type=str,
        help="Path to sample in .fastq",
    )
    parser.add_argument(
        "alignments",
        type=str,
        help="Path to alignments in .gaf/.bam",
    )

    parser.add_argument(
        "-t",
        "--tool",
        dest="tool",
        default="rspoa",
        required=True,
        type=str,
        help="Tool (rspoa/giraffe/bwa) used to produce the alignments (default: rspoa)",
    )

    args = parser.parse_args()

    reads = {}
    for record in SeqIO.parse(args.sample, "fastq"):
        reads[record.id] = str(record.seq)
    print(len(reads))

    if args.tool == "rspoa" or args.tool == "giraffe":
        nodes = load_gfa_nodes(args.graph)
        for line in open(args.alignments):
            line = line.strip("\n").split("\t")
            idx = line[0]
            path = line[5]
            read = reads[idx]
            path_seq = ""
            if args.tool == "rspoa":
                path_seq = convert_se_rspoa(nodes, path, int(line[7]), int(line[8]))
            else:
                path_seq = convert_se_giraffe(nodes, path, int(line[7]), int(line[8]))
            print(idx, path_seq, read, lev(path_seq, read))
    elif args.tool == "bwa":
        bam = pysam.AlignmentFile(args.alignments, "rb")
        for read in bam.fetch():
            print(read.query_name, "*", "*", read.get_tag("NM"))


if __name__ == "__main__":
    main()
