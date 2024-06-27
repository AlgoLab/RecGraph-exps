#!/bin/bash

input_file="genes_HLA.txt"
graphs_dir="https://raw.githubusercontent.com/ekg/HLA-zoo/master/graphs/seqwish/minimap2/"
seqs_dir="https://raw.githubusercontent.com/ekg/HLA-zoo/master/seqs/"

# Verifica se il file genes.txt esiste
if [ ! -f "$input_file" ]; then
    echo "Errore: Il file $input_file non esiste."
    exit 1
fi

# Leggi ogni riga del file genes.txt
while IFS= read -r gene_name
do
    gene_name_trimmed=$(echo "$gene_name" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//')
    graph_url="${graphs_dir}${gene_name_trimmed}.gfa"
    seq_url="${seqs_dir}${gene_name_trimmed}.fa"
    directory="HLA-$gene_name_trimmed"
    mkdir -p "../data/$directory"
    curl -o "../data/$directory/graph.gfa" "$graph_url"
    curl -o "../data/$directory/reads.fa" "$seq_url"

done < "$input_file"