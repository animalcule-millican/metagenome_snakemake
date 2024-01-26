#!/bin/bash
source ~/.bashrc
mamba activate read-mapping

ann=/home/glbrc.org/millican/repos/metagenome_snakemake/data/annotation/rhizo-coassembly-annotation.gff
bam=/home/glbrc.org/millican/repos/metagenome_snakemake/data/map/${1}.co_assembly.sorted.bam
output=/home/glbrc.org/millican/repos/metagenome_snakemake/data/count/${1}.count
python -m HTSeq.scripts.count -f bam -r pos -s no -t CDS -i ID \
--additional-attr product --additional-attr ko --additional-attr=cog \
--additional-attr pfam --additional-attr ec_number --additional-attr smart \
--additional-attr tigrfam --secondary-alignments ignore --supplementary-alignments ignore \
--counts_output $output \
$bam \
$ann
