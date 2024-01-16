#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate read-mapping
source random_directory.sh
oldPERL5LIB=/opt/bifxapps/perl516/share/perl5:/opt/bifxapps/perl516/lib64/perl5:/opt/bifxapps/perl516/lib/perl5
export PERL5LIB=/home/glbrc.org/millican/mambaforge/envs/read-mapping/lib/perl5
contig=$1
index=$2

bowtie2-build -f $fasta $index

rm -r $TMPDIR
export PERL5LIB=$oldPERL5LIB