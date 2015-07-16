# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:05:04 2015

@author: tyler
"""

from Bio import SeqIO
import os
import subprocess
import time
import Analysis2



read_path = "./seqs/"
template_path = "./templates/"
mob_path = "./elements/"
plasmid_path = "./plasmids/"
alignment_path = "./alignments/"
curr_output_path = "./output/"
log_path = "./output/"
print "Create fewer alignments."
start = time.time()

err_check = Analysis2.seqproc_2(read_path, plasmid_path, template_path, mob_path)
# Run files through MAFFT
for dirName, subdirList, fileList in os.walk(template_path):
    for f in fileList:
        outfile_path = alignment_path+f+"_alignment"
        with open(outfile_path, "w") as out_file:  
           print "MAFFT:\nInput:\t"+f+"\nOutput:\t"+outfile_path
           print "Creating alignment..."
           check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet", template_path+f], stdout=out_file)
           
    # Analyze the alignments
err_check = Analysis2.analysis_control(alignment_path, curr_output_path, log_path)
    
   
stop = time.time()

print "TIME:", stop-start      