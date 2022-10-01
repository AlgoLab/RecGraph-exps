import sys
import gzip
from Bio import SeqIO

with gzip.open(sys.argv[1], mode="rt") as f:
    for record in SeqIO.parse(f, "fastq"):
        new_record = record.reverse_complement()
        new_record.id = record.id
        new_record.description = record.description
        SeqIO.write(new_record, sys.stdout, "fastq")
