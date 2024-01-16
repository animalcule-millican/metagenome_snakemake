#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mkdir -p /home/glbrc.org/millican/ref_db/reference_genomes/genomes/$1 /home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/$1
/home/glbrc.org/millican/repos/metagenome_snakemake/etc/get_genomes_taxonomy.py $1
#/home/glbrc.org/millican/repos/metagenome_snakemake/etc/taxify_ncbi_accessions.py $1



