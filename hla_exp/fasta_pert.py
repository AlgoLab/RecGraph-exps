import argparse
import random

def modify_seq(seq, mod_prop):
    new_seq = ""
    for base in seq:
        if random.random() < mod_prop:
            new_base = random.choice(['A', 'C', 'G', 'T'])
            if new_base == base:
                new_base = random.choice(['', f"{base}{base}"])
            new_seq += new_base
        else:
            new_seq += base
    return new_seq

def modify_fasta_file(input_file, mod_prop):
    with open(input_file, 'r') as f_input:
        line = f_input.readline()
        while line:
            if line.startswith('>'):  
                print(line, end='')
                line = f_input.readline()
            else: 
                seq = line.strip()
                seq_mod = modify_seq(seq, mod_prop)
                print(seq_mod)
                line = f_input.readline()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modifica un file FASTA con una certa probabilità di mutazione.")
    parser.add_argument("input_file", type=str, help="Percorso del file FASTA di input.")
    parser.add_argument("--p", type=float, default=0.1, help="Probabilità di modifica della base [0.0 - 1.0]. Default: 0.1")
    args = parser.parse_args()

    modify_fasta_file(args.input_file, args.p)


