#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate branchwater
source random_directory.sh
echo $PYTHONPATH
csvfile=$TMPDIR/mag_list.csv

/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/scripts/sketch_mags.py $1 $2 $3 $4
#/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/scripts/build_mag_list.py $1 $csvfile

#sourmash scripts manysketch -o $2 -c $3 -p k=21,k=31,k=51,scaled=1000,abund $csvfile

trap "rm -r $TMPDIR" SIGKILL
trap "rm -r $TMPDIR" EXIT
trap "rm -r $TMPDIR" ERR

rm -r $TMPDIR