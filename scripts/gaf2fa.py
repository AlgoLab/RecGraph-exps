import sys
import natsort

def load_gfa_nodes(path):
    nodes = {}
    with open(path, 'r') as fin:
        for line in fin:
            if line.startswith('S'):
                line = line.strip()
                _, nid, seq = line.split('\t')
                nodes[nid] = seq
    return nodes

def convert(nodes, path):
    r = ''
    if ">" in path:
        for n in path.split(">")[1:]:
            r += nodes[n]
    else:
        print(path)
        raise NotImplementedError("Path on reverse strand is not implemented yet. Not clear how to read it.")
    return r

def main():
    nodes = load_gfa_nodes(sys.argv[1])
    for read in natsort.natsorted(sys.argv[2:]):
        with open(read, 'r') as fin:
            for line in fin:
                line = line.strip()
                d = line.split('\t')
                print(f'>{d[0]}')
                # for n in d[5].split(">")[1:]:
                #     print(gfa[n], end='')
                print(convert(nodes, d[5]))
        print()
if __name__ == '__main__':
    main()