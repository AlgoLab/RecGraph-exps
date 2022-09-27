# RSPOA-EXPS

```bash
wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

mamba create -c bioconda -n rspoaexps snakemake-minimal biopython seqtk
# install abpoa and rspoa somewhere
```

```bash
# updates config.yaml setting seq, odir, abpoa, and rspoa
snakemake [-n] [-p] -c 16 --use-conda
# requires python natsort
snakemake -s postrun.smk -c16 -R calc [-n] [-p]
snakemake -s postrun.smk -c16 -R tables [-n] [-p]

tail -n +1 /data/abpoa-comparison/results/*.N75.L10000.dists.table | most

```