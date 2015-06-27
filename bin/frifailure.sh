#!/bin/bash

# Author: Tyler Camp
# Date: 2015/6/9
# frifailure: a tool to control the FRI-Failure analysis pipeline


# Execute seqproc.py
# This script creates the FASTA files that will be used by MAFFT

python ../python/seqproc.py


# All files that need to be aligned with MAFFT are in the templates folder

cd ../
PROJECT_DIR=$(readlink -f ${0%/*}) # neat trick to get abolute path from bash
ALIGN_DIR=$PROJECT_DIR/alignments/
TEMPLATE_DIR=$PROJECT_DIR/templates/
cd $TEMPLATE_DIR


mkdir $ALIGN_DIR # Create directory to store alignment output from MAFFT

for file in $PROJECT_DIR/templates/*
do
    TEMP=$(basename $file)
    ALIGN="alignment_$TEMP"
    file_dest=$ALIGN_DIR/$ALIGN
    mafft $file > $file_dest
done

# Execute analysis.py
# This script performs the actual analyses of the alignment files

cd ../
python ./python/analysis.py


