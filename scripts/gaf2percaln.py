import sys
import natsort

def main():
    for read in natsort.natsorted(sys.argv[1:]):
        with open(read, 'r') as fin:
            for line in fin:
                line = line.strip()
                d = line.split('\t')
                l = float(d[1])
                s = float(d[2])
                e = float(d[3])
                print(f'{(e-s)/l:.4f}')
                
if __name__ == '__main__':
    main()