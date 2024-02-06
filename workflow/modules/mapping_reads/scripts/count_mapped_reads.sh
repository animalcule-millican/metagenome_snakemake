#!/bin/bash
source ~/.bashrc
mamba activate read-mapping

export gff=/home/glbrc.org/millican/repos/metagenome_snakemake/data/annotation/rhizo-coassembly-annotation.gff
export aln=/home/glbrc.org/millican/repos/metagenome_snakemake/data/map/${1}.co_assembly.sorted.bam
export output=/home/glbrc.org/millican/repos/metagenome_snakemake/data/count/${1}.count.txt

htseq-count -f bam -r pos -s no -t CDS -i ID --additional-attr product_name --additional-attr ko --additional-attr cog --additional-attr pfam --additional-attr ec_number --additional-attr tigrfam $aln $gff > $output



