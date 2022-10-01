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
        data += [int(x.strip("\n").split(" ")[3]) for x in f]
        f.close()
    return nreads, np.array(data)


def main():
    parser = argparse.ArgumentParser(description="Plot edit distance.")
    parser.add_argument(
        "ODIR",
        type=str,
        help="Path to out directory",
    )
    parser.add_argument(
        "-g",
        "--gene",
        dest="gene",
        default="*",
        type=str,
        help="Gene to consider (default: *, i.e., all)",
    )
    parser.add_argument(
        "-r",
        "--recombination",
        dest="recombination",
        default="*",
        type=str,
        help="Recombination to consider (default: *, i.e., all)",
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
    tools = args.tools.split(",")
    data = {}
    nreads = 0
    for tool in tools:
        n, d = extract(
            os.path.join(
                args.ODIR,
                args.gene,
                args.recombination,
                args.coverage,
                f"{tool}.dist.txt",
            )
        )
        data[tool] = d
        if nreads == 0:
            nreads = n
        else:
            assert n == nreads

    for tool in data:
        print(tool, sum(data[tool]) / len(data[tool]))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 7))

    ax1.grid(color="lightgrey", alpha=1, ls="--")
    XLIM = 40
    Xs = list(range(XLIM + 1))
    for tool in data:
        BARS = [0 for i in range(XLIM + 1)]
        for i in range(XLIM + 1):
            BARS[i] = (data[tool] <= i).sum() / nreads
        ax1.step(Xs, BARS, label=tool)
    # ax.plot([0, 150], [0, 1], color="black", ls="-")
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

    # ax2.grid(color="lightgrey", alpha=1, ls="--")
    sns.violinplot(data=df, x="Edit Distance", y="Tool", ax=ax2)

    plt.suptitle("Alignments Accuracy (eCDF and ViolinPlot)")
    plt.savefig(args.output)


if __name__ == "__main__":
    main()
