#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate thats-so-meta
source random_directory.sh
input=$1
annotation=$2
output=$3
cpu=$4

featureCounts -T $cpu -p --countReadPairs --tmpDir $TMPDIR -t CDS -g ID -F GFF -a "$annotation" -o "$output" "$input"

rm -r $TMPDIR 