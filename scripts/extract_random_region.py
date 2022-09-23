import sys
import random

from Bio import SeqIO


def main():
    record = next(SeqIO.parse(sys.argv[1], "fasta"))
    l = int(sys.argv[2])
    CSTART, CEND = 84000, len(record) - 200  # HARDCODED FOR CHR21

    while True:
        s = random.randint(CSTART, CEND - 2 * l)
        if str(record[s : s + l + 1].seq).count("N") < l / 500:  # HARDCODED
            SeqIO.write(record[s : s + l + 1], sys.stdout, "fasta")
            break


if __name__ == "__main__":
    main()
