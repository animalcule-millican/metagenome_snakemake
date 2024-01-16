#!/bin/bash
source ~/repos/Slime_Py/workflow/env/slimepy.env
mamba activate slime-py
bbcms.sh in=$1 out=$2 bits=4 hashes=3 k=31

if [ -f $2 ]; then
    rm $1
fi