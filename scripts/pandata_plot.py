import argparse
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def extract(globpath):
    nreads = 0
    data = []
    for fpath in glob.glob(globpath):
        f = open(fpath)
        nreads += int(next(f)[:-1])
        # data += [int(x.strip("\n").split(" ")[3]) for x in f]
        for line in f:
            r, _, _, d = line.strip("\n").split(" ")
            d = int(d)
            data.append(d)
            if d > 70:
                print(fpath, r)
        f.close()
    return nreads, data


def main():
    parser = argparse.ArgumentParser(description="Plot edit distance.")
    parser.add_argument(
        "ODIR",
        type=str,
        help="Path to out directory",
    )
    parser.add_argument(
        "--genes",
        dest="genes",
        required=True,
        type=str,
        help="Comma-separated list of genes or file with genes to consider - one per line",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        dest="coverage",
        required=True,
        type=str,
        help="Recombination to consider (default: *, i.e., all)",
    )
    parser.add_argument(
        "-t",
        "--tools",
        dest="tools",
        default="rspoa-9,giraffe,bwa",
        type=str,
        help="Tools to consider, comma-separated (default: rspoa-9,giraffe,bwa)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="plot.png",
        type=str,
        help="Path to output .png (default: plot.png)",
    )

    args = parser.parse_args()

    genes = []
    if os.path.isfile(args.genes):
        for line in open(args.genes):
            genes.append(line.strip("\n"))
    else:
        genes = args.genes.split(",")

    tools = args.tools.split(",")

    data = {}
    nreads = 0
    ff = True
    for tool in tools:
        toolname = tool
        # if toolname.startswith("rspoa-5"):
        #     toolname = "rspoa-5"
        # elif toolname.startswith("rspoa-9"):
        #     toolname = "rspoa-9"
        data[toolname] = []
        for gene in genes:
            n, d = extract(
                os.path.join(
                    args.ODIR, gene, "1-150", "*", args.coverage, f"{tool}.dist.txt"
                )
            )
            if ff:
                nreads += n
            data[toolname] += d
        ff = False

    for tool in data:
        data[tool] = np.array(data[tool])
        print(tool, sum(data[tool]), len(data[tool]), sum(data[tool]) / len(data[tool]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 7))

    ax1.grid(color="lightgrey", alpha=1, ls="--")
    XLIM = 10
    Xs = list(range(XLIM + 1))
    for tool in data:
        BARS = [0 for i in range(XLIM + 1)]
        for i in range(XLIM + 1):
            BARS[i] = (data[tool] <= i).sum() / nreads
        ax1.step(Xs, BARS, label=tool)
    # ax1.plot([0, 150], [0, 1], color="black", ls="-")
    ax1.set_xlabel("Edit Distance")
    ax1.set_xlim(-1, XLIM + 1)
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel("Cumulative Aligned Reads (%)")
    ax1.legend(loc=4)

    df = []
    for tool in data:
        df += [[x, tool] for x in data[tool]]
    df = pd.DataFrame(
        df,
        columns=["Edit Distance", "Tool"],
    )

    sns.boxplot(data=df, x="Edit Distance", y="Tool", ax=ax2)
    ## sns.violinplot(data=df, x="Edit Distance", y="Tool", ax=ax2)
    ## sns.swarmplot(data=df, x="Edit Distance", y="Tool", color="black", ax=ax2)
    ax2.set_xlim(-1, XLIM + 1)

    # Unwanted line plot
    # ax2.grid(color="lightgrey", alpha=1, ls="--")
    # COLORS = {"rspoa-5": "blue", "rspoa-9": "orange", "giraffe": "green"}
    # MARKERS = {"rspoa-5": "o", "rspoa-9": "o", "giraffe": "x"}
    # BWA_BARS = [0 for i in range(XLIM + 1)]
    # for i in range(XLIM + 1):
    #     BWA_BARS[i] = (data["bwa"] == i).sum()
    # print("bwa", BWA_BARS)
    # for tool in data:
    #     if tool == "bwa":
    #         continue
    #     BARS = [0 for i in range(XLIM + 1)]
    #     for i in range(XLIM + 1):
    #         BARS[i] = BWA_BARS[i] - (data[tool] == i).sum()
    #     ax2.plot(Xs, BARS, color=COLORS[tool], marker=MARKERS[tool], label=tool)
    # ax2.set_xticks(list(range(0, XLIM + 1)))
    # ax2.set_xlabel("Edit Distance")
    # ax2.set_ylabel("Difference in number of alignments wrt BWA")
    # ax2.legend(loc=4)

    plt.suptitle("Alignments Accuracy (eCDF and BoxPlot)")

    plt.tight_layout()
    plt.savefig(args.output)


if __name__ == "__main__":
    main()
