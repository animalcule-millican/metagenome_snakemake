#!/bin/bash
source ~/.bashrc
mamba activate read-mapping
sdir=/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/mapping_reads/scripts
gff=$1
gtf=$2
bam=$3
count=$4
genelen=$5
name=$6
output=$7
/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/mapping_reads/scripts/gff2gtf.sh $gff > $gtf
htseq-count -r pos -t CDS -f bam $bam $gtf > $count
cut -f4,5,9 $gtf | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > $genelen
/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/mapping_reads/scripts/tpm_table.py -n $name -c $count -i <(echo -e "${name}\t100") -l $genelen > $output
