import sys


def main():
    gfa_path = sys.argv[1]

    to_remove = []
    predecessors = {}
    successors = {}
    for line in open(gfa_path):
        if not line.startswith("S"):
            continue
        _, idx, seq, _ = line.split("\t")
        if seq == "*":
            to_remove.append(idx)
            predecessors[idx] = []
            successors[idx] = []
        else:
            print(line, end="")

    for line in open(gfa_path):
        if not line.startswith("L"):
            continue
        _, idx1, _, idx2, _, _ = line.split("\t")
        if idx1 in to_remove:
            successors[idx1].append(idx2)
        elif idx2 in to_remove:
            predecessors[idx2].append(idx1)
        else:
            print(line, end="")

    for idx in successors:
        while any([x in to_remove for x in successors[idx]]):
            new_successors = []
            for sidx in successors[idx]:
                if sidx in to_remove:
                    new_successors.extend(successors[sidx])
                else:
                    new_successors.append(sidx)
            successors[idx] = new_successors
    for idx in predecessors:
        while any([x in to_remove for x in predecessors[idx]]):
            new_predecessors = []
            for sidx in predecessors[idx]:
                if sidx in to_remove:
                    new_predecessors.extend(predecessors[sidx])
                else:
                    new_predecessors.append(sidx)
            predecessors[idx] = new_predecessors

    # adds new edges
    for idx in to_remove:
        if len(predecessors[idx]) == 0 or len(successors[idx]) == 0:
            continue
        for p in predecessors[idx]:
            for s in successors[idx]:
                print("L", p, "+", s, "+", "0M", sep="\t")

    # for line in gfa_path:
    #     if line.startswith("P"):
    #         _, idx1, _, idx2, _, _ = line.split("\t")
    #         if idx1 in to_remove or idx2 in to_remove:
    #             pass


if __name__ == "__main__":
    main()
