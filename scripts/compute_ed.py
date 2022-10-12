import argparse
from Bio import SeqIO
import pysam
from Levenshtein import distance as lev
from gaf2fa import load_gfa_nodes


def convert(nodes, path, s, e):
    r = ""
    if ">" in path:
        path = path.split(">")[1:]
        for n in path:
            r += nodes[n]
    else:
        print(path)
        raise NotImplementedError(
            "Path on reverse strand is not implemented yet.\
            Not clear how to read it."
        )
    return r[s : e + 1]


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
        "-f",
        "--format",
        dest="fformat",
        default="gaf",
        required=True,
        type=str,
        help="Alignment (gaf/bam) format (default: gaf)",
    )

    args = parser.parse_args()

    reads = {}
    for record in SeqIO.parse(args.sample, "fastq"):
        reads[record.id] = str(record.seq)
    print(len(reads))

    if args.fformat.startswith("gaf"):
        nodes = load_gfa_nodes(args.graph)
        for line in open(args.alignments):
            line = line.strip("\n").split("\t")
            idx = line[0]
            path = line[5]
            if path == "*":
                continue  # FIXME for giraffe
            path_seq = line[-1]
            if ":" in path_seq:
                path_seq = convert(nodes, path, int(line[7]), int(line[8]))
            read = reads[idx]
            print(idx, path_seq, read, lev(path_seq, read))
    else:
        bam = pysam.AlignmentFile(args.alignments, "rb")
        for read in bam.fetch():
            print(read.query_name, "*", "*", read.get_tag("NM"))


if __name__ == "__main__":
    main()
