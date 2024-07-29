# Code for analyses in Sigeman et al. 2024
## Sigeman, H., Downing, P.A., Zhang, H. et al. The rate of W chromosome degeneration across multiple avian neo-sex chromosomes. Sci Rep 14, 16548 (2024). https://doi.org/10.1038/s41598-024-66470-7



### 1. Code to produce the zebra finch gene bed files in this script: 
`code/create_zebra_finch_bed_file.sh`


### 2. Run snakefiles for each species (but first modify the file paths in the config files)

```
snakemake -s snakemake -s snakefile_exon_intron -j 15 --cores 7 --configfile config_files/Panurus_biarmicus_config -R all --use-conda --conda-frontend conda

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Calandrella_cinerea_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Alauda_arvensis_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Cisticola_juncidis_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Sylvietta_brachyura_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Alauda_razae_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Camaroptera_brevicaudata_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Eremophila_alpestris_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp

snakemake -s snakefile_exon_intron -j 15 -R all --configfile config_files/Panurus_biarmicus_config --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -nrp
```

### 3. Summarize the results from each species (put as separate script)
`code/summarize_tables.sh`

### 4. Run R-code for plots and statistics
`Rscript1_make_tables.R`
`Rscript2_gene_loss.R`


Contact: Hanna Sigeman (hanna.sigeman@oulu.fi)
