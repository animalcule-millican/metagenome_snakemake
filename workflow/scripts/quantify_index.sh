#!/bin/bash
source ~/.bashrc
mamba activate salmon-env
source random_directory.sh

salmon index -t $1 -i $2 -k 31 --tmpdir $TMPDIR --threads $3 

rm -r $TMPDIR


~/metagenome/oil/WT_ome_39_FD/IMG_Data/Ga0501086_contigs.fna
~/metagenome/oil/WT_ome_39_FD/IMG_Data/Ga0501086_functional_annotation.gff

tar xzf ~/metagenome/oil/WT_ome_39_FD/IMG_Data/3300053751.tar.gz 3300053751/3300053751.a.fna -C $TMPDIR
salmon index -t 3300053751/3300053751.a.fna -i Ga0501086.quant -k 31 tmpdir $TMPDIR --threads 16
rename.sh in=~/metagenome/oil/WT_ome_39_FD/Filtered_Raw_Data/52683.2.418583.GTCGGATT-CTTACGAG.filter-METAGENOME.fastq.gz out1=$TMPDIR/Ga0501086_R1.fastq.gz out2=$TMPDIR/Ga0501086_R2.fastq.gz int=t
mamba deactivate thats-so-meta; mamba activate salmon-env
salmon quant -i Ga0501086.quant \
-1 $TMPDIR/Ga0501086_R1.fastq.gz \
-2 $TMPDIR/Ga0501086_R2.fastq.gz \
-g ~/metagenome/oil/WT_ome_39_FD/IMG_Data/Ga0501086_functional_annotation.gff \
-p 16 -o quants/${samp}_quant \
--meta 