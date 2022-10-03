import sys
from Levenshtein import distance as lev


with open(sys.argv[1]) as fin:
    r1 = fin.readlines()

with open(sys.argv[2]) as fin:
    r2 = fin.readlines()

for ix in range(1, len(r1), 2):
    print(
            lev(r1[ix].rstrip(), r2[ix].rstrip())
        )