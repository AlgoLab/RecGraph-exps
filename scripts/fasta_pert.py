import argparse
import random

def modifica_sequenza(sequenza, probabilita_modifica):
    nuova_sequenza = ""
    for base in sequenza:
        if random.random() < probabilita_modifica:
            nuova_base = random.choice(['A', 'C', 'G', 'T'])
            if nuova_base == base:
                nuova_base = random.choice(['', f"{base}{base}"])
            nuova_sequenza += nuova_base
        else:
            nuova_sequenza += base
    return nuova_sequenza

def modifica_file_fasta(input_file, probabilita_modifica):
    with open(input_file, 'r') as f_input:
        linea = f_input.readline()
        while linea:
            if linea.startswith('>'):  # Intestazione della sequenza
                print(linea, end='')
                linea = f_input.readline()
            else:  # Sequenza
                sequenza = linea.strip()
                sequenza_modificata = modifica_sequenza(sequenza, probabilita_modifica)
                print(sequenza_modificata)
                linea = f_input.readline()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modifica un file FASTA con una certa probabilità di mutazione.")
    parser.add_argument("input_file", type=str, help="Percorso del file FASTA di input.")
    parser.add_argument("--p", type=float, default=0.1, help="Probabilità di modifica della base [0.0 - 1.0]. Default: 0.1")
    args = parser.parse_args()

    modifica_file_fasta(args.input_file, args.p)


