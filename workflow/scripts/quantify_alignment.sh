#!/bin/bash
source ~/.bashrc
mamba activate salmon-env
source random_directory.sh

salmon quant --alignments $1 --targets $2 --output $3 --threads $4 --mappingCacheMemoryLimit 100000000 --meta

rm -r $TMPDIR