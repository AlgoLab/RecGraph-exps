"""
USAGE example:
python scripts/get_rand_fas.py -i /data/pandoro/ecoli_core_gene_alignments/GC00001080_3_na_aln.fa -n $(grep -c ">" /data/pandoro/ecoli_core_gene_alignments/GC00001080_3_na_aln.fa) -p 0.3 -o ecoli.test
"""

def main(args):
    randixs = np.random.choice(args.n, int(args.n * args.p), replace=False)
    ix = 0
    prefix = args.o.rstrip('.')
    forgfa_f = open(f'{prefix}.forgfa.fa', 'w+')
    fortest_f = open(f'{prefix}.fortest.fa', 'w+')
    for record in SeqIO.parse(args.i, "fasta"):
        if ix in randixs:
            forgfa_f.write(f'>{record.id}\n')
            forgfa_f.write(str(record.seq)+'\n')
        else:
            fortest_f.write(f'>{record.id}\n')
            fortest_f.write(str(record.seq)+'\n')

        ix+=1

    forgfa_f.close()
    fortest_f.close()

if __name__ == '__main__':
    from Bio import SeqIO
    import sys
    import numpy as np
    import argparse
    
    parser = argparse.ArgumentParser(description='Random selection of reads from fasta.')
    parser.add_argument('-i', type=str, required=True,
                        help='FASTA file')
    parser.add_argument('-n', type=int, required=True,
                        help='tot number of reads in FASTA file')
    parser.add_argument('-p', type=float, required=True,
                        help='percentage of reads to select for GFA (float, eg. 0.3 for 30%)')
    parser.add_argument('-o', type=str, required=True,
                        help='output prefix (will create two file *.forgfa.fa and *.fortest.fa')
    

    args = parser.parse_args()
    main(args)
