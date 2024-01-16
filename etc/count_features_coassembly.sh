#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
gzcofile="/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/bam/rhizo/rhizo-coassembly.sam.bam.gz"
cofile="/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/bam/rhizo/rhizo-coassembly.sam.bam"
gff="/home/glbrc.org/millican/projects/TMP_transfer_data/ptmp/adina/millican/oil/data/annotation/rhizo/rhizo-coassembly-annotation-feature-counting.gff"
gunzip -c $gzcofile > $cofile
samtools index -@ 32 $cofile
featureCounts -T 32 -p --countReadPairs --tmpDir $TMPDIR -t CDS -g ID -F GFF -a $gff -o /home/glbrc.org/millican/repos/metagenome_snakemake/co_assembly.count $cofile

rm -rf $TMPDIR