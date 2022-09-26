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

def main():
    gfa = load_gfa_nodes(sys.argv[1])
    for read in natsort.natsorted(sys.argv[2:]):
        with open(read, 'r') as fin:
            for line in fin:
                line = line.strip()
                d = line.split('\t')
                rid = d[0]
                print(f'>{rid}')
                for n in d[5].split(">")[1:]:
                    print(gfa[n], end='')
        print()
if __name__ == '__main__':
    main()