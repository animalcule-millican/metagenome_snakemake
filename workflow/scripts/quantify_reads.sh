#!/bin/bash
source ~/.bashrc
mamba activate salmon-env
source random_directory.sh

salmon quant -i $1 -1 $2 -2 $3 -o $4 --meta -g $5 --threads $6

rm -r $TMPDIR