#!/bin/bash
source ~/.bashrc
mamba activate thats-so-meta
source random_directory.sh
#export BB_rename=/home/glbrc.org/millican/miniconda/miniconda/envs/thats-so-meta/bin/rename.sh
#name=$(basename -s $4 $1)

rename.sh in=$1 out=$2 out2=$3 int=t #  prefix=$name 

rm -r $TMPDIR