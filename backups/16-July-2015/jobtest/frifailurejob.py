# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler
"""

    
from Bio import AlignIO
from jobmanager import jobmanager
import os
import subprocess
import sys
from decimal import *


def controller(job):
    """Controller for FRI-Failure Analysis pipeline."""
    read_path = job.input_dirs["reads"]
    plasmid_path = job.input_dirs["plasmids"]
    template_path = job.output_dirs["templates"]
    mob_path = job.input_dirs["mobile_elements"]
    # Using module Analysis
    Analysis = job.process_module    
    
    # Create templates to be aligned by MAFFT
    Analysis.seqproc(read_path, plasmid_path, template_path, mob_path)
    # Add templates to input dictionary
    job.input_dirs["templates"] = job.output_dirs["templates"]
    
    #Run files through MAFFT    
    alignment_path = job.output_dirs["alignments"]
    print "Creating alignments with MAFFT:"
    for dirName, subdirList, fileList in os.walk(template_path):
        file_num = len(fileList)
        total_complete = 0                    
        for f in fileList:
            percent_complete = Decimal(total_complete)/Decimal(file_num)
            percent_complete = int(percent_complete * 100)
            print str(percent_complete) + "% complete...\r",            
            sys.stdout.flush()
            outfile_path = alignment_path+f+"_alignment"
            with open(outfile_path, "w") as out_file:  
               check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet", template_path+f], stdout=out_file)
               sys.stdout.flush()             
               if (check):
                   with open(job.logfile_name, "a") as log:
                       log.write("Mafft error (input=" + f + ", output=" + outfile_path + ")\n")
                       
            
            sys.stdout.flush()
            total_complete += 1
        print "100% complete."
  # Add alignments to input dictionary
    job.input_dirs["alignments"] = job.output_dirs["alignments"]
    
    # Analyze alignments


    mob_evidence_dict = dict()
    if not(os.access(alignment_path, os.F_OK)):
        err_str =  "Error: cannot access", alignment_path
        with open(job.logfile_name, "a") as log:
            log.write(err_str + "\n")
        return 1
    
    curr_output_path = job.master_output_dir
    print "Gathering mobile element homology:"
    for dirName, subdirList, fileList in os.walk(alignment_path):
        file_num = len(fileList)
        total_complete = 0                    
        for f in fileList:           
            if "rm" in f:
                percent_complete = Decimal(total_complete)/Decimal(file_num)
                percent_complete = int(percent_complete * 100)
                print str(percent_complete) + "% complete...\r",
                # read/mobile element alignment
                # need to determine if any homology to do nt analysis
                rm_alignments = AlignIO.read(alignment_path+f, "fasta")                
                evidence = Analysis.mob_info(rm_alignments)
                """Evidence contains
                [0]: id of Sanger read
                [1]: id of mobile element
                [2]: list of mutations identified by analyzing mobile element
                     against Sanger read and validity score
                [3]: offset from beginning of Sanger read at which to truncate analysis"""  
                    
                if not(evidence is None):
                    key = evidence[0]
                    value = [evidence[1], evidence[2], evidence[3]]
                    mob_evidence_dict[key]= value
                total_complete += 1
        print "100% complete."
    print "Analyzing alignments of plasmids/reads:"
    for dirName, subdirList, fileList in os.walk(alignment_path):
        file_num = len(fileList)
        total_complete = 0                    
        for f in fileList:            
            if "pr" in f:  # plasmid/read alignment 
                percent_complete = Decimal(total_complete)/Decimal(file_num)
                percent_complete = int(percent_complete * 100)
                print str(percent_complete) + "% complete...\r",
                mutation_list = list()
                pr_alignments = AlignIO.read(alignment_path+f, "fasta")
                print pr_alignments[0].id + "\r", 
                sys.stdout.flush()
                plasmid_indices = Analysis.get_indices(pr_alignments[0])
                #print "plasmid_indices:", plasmid_indices
                for x in range(1, len(pr_alignments)): # only get Sange reads
                    read_indices = Analysis.get_indices(pr_alignments[x])
                    #print "read_indices:", read_indices
                    start = max(plasmid_indices[0], read_indices[0])
                    stop = min(plasmid_indices[1], read_indices[1])
                    #print "initially, stop was", stop
                    key = pr_alignments[x].id
                    if (key in mob_evidence_dict):
                        stop = start + mob_evidence_dict[key][2] # offset should be from beginning
                        #print "...but then, stop was", stop, "for key", key
                        mutation_string = ("MOB\t.\t.\t"+key+"\t"+str(stop) + "\t" + 
                                            mob_evidence_dict[key][0])
                        mutation_list.append(mutation_string)                                                
                        if (mob_evidence_dict[key][1]): # append any mutations within the mobile element
                            """mob_evidence_dict[key] is the evidence list,
                            [1] is the list of mutations, validity,
                            [0] gets the list of mutations"""                            
                            for mutation in mob_evidence_dict[key][1][0]:
                                mutation_list.append(mutation)
                    
                    #print "analyze:", pr_alignments[0].id, ",", pr_alignments[x].id, ",", start, ",", stop
                                        
                    mutation_list_2 = Analysis.analyze_seq(pr_alignments[0], pr_alignments[x], start, stop)
                    for mutation in mutation_list_2[0]: # mutation list at index 0
                        mutation_list.append(mutation)
                total_complete += 1
                
                name = pr_alignments[0].id + ".gd"
                with open(curr_output_path+name, "w") as out_file:
                    out_file.write("#=GENOME_DIFF 1.0\n")
                    for mutation in mutation_list:
                        out_file.write(mutation + "\n")
                del(mutation_list)     
        print "100% complete."
    
def main():
    input_dirs = dict()    
    output_dirs = dict()
    input_dirs["mobile_elements"] = "/home/tyler/Documents/research/FRI-failure/mob_elements/" # location of mobile elements
    input_dirs["reads"] = "/home/tyler/Documents/research/FRI-failure/reads/" # location of Sanger reads
    input_dirs["plasmids"] = "/home/tyler/Documents/research/FRI-failure/plasmids/"
    output_dirs["templates"] = "templates/" # location of .fasta files that will be run through MAFFT
    output_dirs["alignments"] = "alignments/" # path to store and access alignments produced by MAFFT
    myJob = jobmanager(input_dirs, output_dirs, "logfile.log", "Analysis")
    controller(myJob)
    
main()