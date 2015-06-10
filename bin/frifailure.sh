#!/bin/bash

# Author: Tyler Camp
# Date: 2015/6/9
# frifailure: a tool to control the FRI-Failure analysis pipeline


# Execute seqproc.py
# This script creates the FASTA files that will be used by MAFFT
python ../python/seqproc.py


#Go through all directories and create relevant alignments with MAFFT
cd ../sequences/
SCRIPT_DIR=$(readlink -f ${0%/*})

for dir in $SCRIPT_DIR/*
do
    echo "Making directory in $dir/alignments"
    mkdir $dir/alignments
    for subdir in $dir/*
    do
	if [ -d $subdir ] && [ $subdir = "$dir/templates" ]
	    then 
       	    for file in $subdir/*
       	    do
		TEMP=$(basename $file)
		ALIGN="alignment_$TEMP"
		file_loc=$dir/templates/$TEMP
		file_dest=$dir/alignments/$ALIGN
		echo "Running mafft on $file_loc and $file_dest"
		mafft $file_loc > $file_dest
	    done
	    
	fi
    done
    
done

# Execute analysis.py
# This script performs the actual analyses of the alignment files

cd ../
python ./python/analysis.py


