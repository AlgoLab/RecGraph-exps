import argparse
import sys
import os
import glob


def extract(fpath):
    data = {}
    f = open(fpath)
    _ = int(next(f)[:-1])
    for line in f:
        line = line.strip("\n").split(" ")
        r, d = line[0], int(line[3])
        if r[-2] == "/":
            r = r[:-2]
        data[r] = data[r] + d if r in data else d
    f.close()
    return data


def get_unmated(fpath):
    unmated = set()
    f = open(fpath)
    _ = int(next(f)[:-1])
    for line in f:
        line = line.strip("\n").split(" ")
        r, _ = line[0], int(line[3])
        if r[-2] == "/":
            r = r[:-2]
        if r in unmated:
            unmated.remove(r)
        else:
            unmated.add(r)
    f.close()
    return unmated


def main():
    parser = argparse.ArgumentParser(description="Dump edit distance CSV to stdout.")
    parser.add_argument(
        "WDIR",
        type=str,
        help="Path to working directory",
    )
    parser.add_argument(
        "--genes",
        dest="genes",
        required=True,
        type=str,
        help="Comma-separated list of genes or file with genes to consider (one per line)",
    )
    parser.add_argument(
        "--rspoa",
        dest="rspoa",
        default="rspoa-5.M2-X4-O4-E2,rspoa-9.R2-r0.1,rspoa-9.R4-r0.1,rspoa-9.R4-r1,rspoa-9.R8-r0.1",
        type=str,
        help="rspoa run to consider (default: rspoa-5.M2-X4-O4-E2,rspoa-9.R2-r0.1,rspoa-9.R4-r0.1,rspoa-9.R4-r1,rspoa-9.R8-r0.1)",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        dest="coverage",
        default="15",
        type=str,
        help="Coverage to consider",
    )
    args = parser.parse_args()

    genes = []
    if os.path.isfile(args.genes):
        for line in open(args.genes):
            genes.append(line.strip("\n"))
    else:
        genes = args.genes.split(",")

    reads = {}
    nfilteredpairs = 0
    for gene in genes:
        for fpath in glob.glob(
            os.path.join(args.WDIR, gene, "1-150", "*", args.coverage, "bwa.dist.txt")
        ):
            recomb = fpath.split("/")[-3]
            basedir = os.path.join(args.WDIR, gene, "1-150", recomb, args.coverage)
            unmateds = get_unmated(os.path.join(basedir, "giraffe.dist.txt"))
            nfilteredpairs += len(unmateds)

            BWA = extract(fpath)
            for r, d in BWA.items():
                if r not in unmateds:
                    reads[r] = {"bwa": d}

            giraffe = extract(os.path.join(basedir, "giraffe.dist.txt"))
            for r, d in giraffe.items():
                assert r in BWA
                if r in unmateds:
                    continue
                reads[r]["giraffe"] = d
            # RSPOA
            for rspoa in args.rspoa.split(","):
                for r, d in extract(os.path.join(basedir, f"{rspoa}.dist.txt")).items():
                    assert r in BWA
                    if r in unmateds:
                        continue
                    reads[r][rspoa] = d
    runs = args.rspoa.split(",") + ["giraffe", "bwa"]
    print("Read", "\t".join(runs), sep="\t")
    for read, Ds in reads.items():
        print(read, "\t".join([str(Ds[r]) for r in runs]), sep="\t")
    print(nfilteredpairs, "filtered pairs.", file=sys.stderr)


if __name__ == "__main__":
    main()
