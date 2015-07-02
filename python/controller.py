# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:23:44 2015

@author: tyler
"""

import os
import subprocess
import datetime


import Analysis
   
def controller():
    """This is the controller for the FRI-Failure analysis pipeline."""
    

    """Setup: Create necessary directories"""

    mob_path = "../mob_elements/"
    read_path = "../reads/"
    plasmid_path = "../plasmids/"
    template_path = "../templates/" # location of .fasta files that will be run through MAFFT
    alignment_path = "../alignments/"
    output_path = "../output/"
    date_time_string = str(datetime.datetime.today()).replace(" ", "-")
    curr_output_path = "../output/genomediff_" + date_time_string 
    
    if not(os.access(template_path, os.F_OK)):
        os.mkdir(template_path)    
    if not(os.access(alignment_path, os.F_OK)):
        os.mkdir(alignment_path)
    if not(os.access(output_path, os.F_OK)):
        os.mkdir(output_path)
    if not(os.access(curr_output_path, os.F_OK)):
        os.mkdir(curr_output_path)    
    
    log_path = curr_output_path + "/logfile.log"
    with open(log_path, "a") as log:
        log.write("Log file created: " + date_time_string+"\n")
       
    # Create templates for alignment
#    err_check = Analysis.seqproc(read_path, plasmid_path, template_path, mob_path)
#    if (err_check):
#        with open(log_path, "a") as log:
#            log.write("Error in Analysis.seqproc\n")
#        return 1
   
    
    # Run files through MAFFT
#    for dirName, subdirList, fileList in os.walk(template_path):
#        for f in fileList:
#            outfile_path = alignment_path+f+"_alignment"
#            with open(outfile_path, "w") as out_file:  
#                print "MAFFT:\nInput:\t"+f+"\nOutput:\t"+outfile_path
#               print "Creating alignment..."
#              check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet", template_path+f], stdout=out_file)
#             if (check):
#                with open(log_path, "a") as log:
#                   log.write("Mafft error (input=" + f + ", output=" + outfile_path + ")\n")
    
    # Analyze the alignments
    err_check = Analysis.analysis_control(alignment_path, curr_output_path, log_path)
    if (err_check):
        with open(log_path, "a") as log:
            log.write("Error in Analysis.analysis_control\n")
    return 1

def main():
    controller()
if (__name__ == "__main__"):
    main()
