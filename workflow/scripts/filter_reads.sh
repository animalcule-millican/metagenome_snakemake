#!/bin/bash
source ~/repos/Slime_Py/workflow/env/slimepy.env
mamba activate slime-py
source random_directory.sh

cd $TMPDIR

rqcfilter2.sh -Xmx30g in=$1 out=$2 \
rqcfilterdata=$RQCFilterData \
barcodefilter=f \
qtrim=r \
trimq=5 \
silvalocal=f \
usetmpdir=t \
tmpdir="$TMPDIR" \
merge=f

cd $HOME

trap handle_signal EXIT
trap handle_signal SIGTERM
trap handle_signal ERR
trap handle_signal SIGINT
trap handle_signal SIGQUIT
trap handle_signal SIGKILL
