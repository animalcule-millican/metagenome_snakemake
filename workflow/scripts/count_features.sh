#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate thats-so-meta
source random_directory.sh
#input=$1
#annotation=$2
#output=$3
#cpu=$4
input=/home/glbrc.org/millican/repos/metagenome_snakemake/data/map/${1}.co_assembly.sorted.bam
#annotation=/home/glbrc.org/millican/repos/metagenome_snakemake/data/annotation/${1}_functional_annotation.gff
annotation=/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/annotation/rhizo/rhizo-coassembly-annotation-feature-counting.gff
output=/home/glbrc.org/millican/repos/metagenome_snakemake/data/count/${1}.co_assembly.count
cpu=16

featureCounts -T $cpu -p --countReadPairs --tmpDir $TMPDIR -t CDS -g ID -F GFF -a "$annotation" -o "$output" "$input"

rm -r $TMPDIR 