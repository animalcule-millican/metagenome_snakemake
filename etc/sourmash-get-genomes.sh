#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate branchwater
sourmash scripts get-genomes --output-directory $1 $2