from ctypes import CDLL
import sys

so = "./scripts/lev.so"
f = CDLL(so)

with open(sys.argv[1]) as fin:
    r1 = fin.readlines()

with open(sys.argv[2]) as fin:
    r2 = fin.readlines()

for ix in range(1, len(r1), 2):
    print(f.lev_dist(
        str.encode(r1[ix].rstrip()), 
        str.encode(r2[ix].rstrip()))
        )