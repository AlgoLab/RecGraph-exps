def load_gfa_nodes(path):
    nodes = {}
    with open(path, "r") as fin:
        for line in fin:
            if line.startswith("S"):
                line = line.strip()
                _, nid, seq = line.split("\t")
                nodes[nid] = seq
    return nodes


def convert(nodes, path, s, e):
    r = ""
    if ">" in path:
        path = path.split(">")[1:]
        r += nodes[path[0]][s:] if path[0] in nodes else ""
        for n in path[1:-1]:
            r += nodes[n] if n in nodes else ""
        r += nodes[path[-1]][: e + 1] if path[-1] in nodes else ""
    else:
        print(path)
        raise NotImplementedError(
            "Path on reverse strand is not implemented yet. Not clear how to read it."
        )
    return r


def main():
    nodes = load_gfa_nodes(sys.argv[1])
    for read in natsort.natsorted(sys.argv[2:]):
        with open(read, "r") as fin:
            for line in fin:
                line = line.strip()
                d = line.split("\t")
                print(f">{d[0]}")
                # for n in d[5].split(">")[1:]:
                #     print(gfa[n], end='')
                print(convert(nodes, d[5], int(d[7]), int(d[8])))


if __name__ == "__main__":
    import sys
    import natsort

    main()
