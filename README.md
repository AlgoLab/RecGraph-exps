# RSPOA-EXPS

To replicate the experiments:
1. install [make_prg](https://github.com/leoisl/make_prg/tree/update) (`update` branch) and [rspoa](https://github.com/AlgoLab/rspoa)
2. download the E.Coli core genes (`.msa`) from [panX](https://pangenome.org/):
```bash
wget http://pangenome.de/dataset/Escherichia_coli/core_gene_alignments.tar.gz
tar xvfz core_gene_alignments.tar.gz
# edit strains name containing /
cd core_gene_alignments
for fa in $(ls *.fa) ; do sed -i "s/\//-/g" $fa ; done
```
3. update `config.yaml`:
   - `seqsdir`, path to `core_gene_alignments` directory
   - `odir`, desired output directory (for results)
   - `rspoa`, path to `rspoa` binary
   - `mkprg`, path to `make_prg` binary
4. build the (minimal) pangenomes graphs and get recombinant strains:
```bash
snakemake [-n] [-p] -s pandmakegraphs.smk -c 32 --use-conda
```
5. simulate and align reads:
```bash
snakemake [-n] [-p] -s pandata.smk -c 32 --use-conda
```
6. dump `.tsv` with alignments accuracy:
```bash
python3 scripts/dump_tsv.py --genes 50genes.list --rspoa rspoa-5.M2-X4-O4-E2,rspoa-9.R4-r0.1 [OUT_DIRECTORY] > alignments_accuracy.tsv
```