# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler

This file is part of fri-failure-analysis.

    fri-failure-analysis is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    fri-failure-analysis is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with fri-failure-analysis.  If not, see <http://www.gnu.org/licenses/>.
"""

    
from Bio import AlignIO
#Legacy import
#from Bio import SeqIO
from jobmanager import jobmanager
import os
import subprocess
import sys
import re
from decimal import *

def phase_1(job):
    """FRI-Failure Analysis: Phase 1
    
    Phase 1 involves processing sequencing files. Each file written includes two sequences: a template 
    (first sequence) and a target (second sequence). Filenames containing "pr" denote plasmid/Sanger read
    files. Filenames containing "rm" denote read/mobile element files."""
    
    # Define necessary names   
    read_path = job.input_dirs["reads"]
    plasmid_path = job.input_dirs["plasmids"]
    template_path = job.output_dirs["templates"]
    mob_path = job.input_dirs["mobile_elements"]
    Analysis = job.process_module
    
    Analysis.seqproc(read_path, plasmid_path, template_path, mob_path)
    # Add templates to input dictionary
    job.input_dirs["templates"] = job.output_dirs["templates"]   
    return
def phase_2(job):
    """FRI-Failure Analysis: Phase 2
    
    Phase 2 creates the alignments with the MAFFT alignment program. File input comes from files written
    during phase 1."""

    # Define necessary directories    
    alignment_path = job.output_dirs["alignments"]
    template_path = job.input_dirs['templates']
    
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
    
    return
def phase_3(job):
    """FRI-Failure Analysis: Phase 3
    
    Phase 3 analyzes the alignments created during phase 2 using module Analysis. Read
    homology to mobile elements is gathered and used during the analysis. All mutations are
    written out in GenomeDiff format currently."""
    
    # Define necessary names
    alignment_path = job.output_dirs['alignments']
    curr_output_path = job.master_output_dir
    Analysis = job.process_module
    
    """As of 2015/8/26, we are no longer tracking mutations in the insertion region.

    # Define dictionary to map reads to plasmids
    pr_map = dict()
    working_dir = job.output_dirs['templates']
    for dirName, subdirList, fileList in os.walk(working_dir):
        for f in fileList:
            if ("pr" in f):
                name_lst= list()
                for seq in SeqIO.parse(working_dir+f, "fasta"):
                    name_lst.append(seq.id)
                pr_map[name_lst[1]] = name_lst[0] # map reads onto plasmids
                del(name_lst)
    """            
    mob_evidence_dict = dict()
    if not(os.access(alignment_path, os.F_OK)):
        err_str =  "Error: cannot access", alignment_path
        with open(job.logfile_name, "a") as log:
            log.write(err_str + "\n")
        return 1
    
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
               
                """Legacy function call.               
               evidence = Analysis.mob_info(rm_alignments, pr_map)
               """
                evidence = Analysis.mob_info(rm_alignments)
                """Evidence contains
                [0]: id of Sanger read
                [1]: id of mobile element
                [2]: list of mutations identified by analyzing mobile element
                     against Sanger read and validity score
                [3]: fraction of correctly matched residues
                [4]: offset from beginning of Sanger read at which to truncate analysis"""  
                    
                if not(evidence is None):
                    key = evidence[0]
                    value = [evidence[1], evidence[2], evidence[4]]
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
                pr_alignment = AlignIO.read(alignment_path+f, "fasta")
                print pr_alignment[0].id + "\r", 
                sys.stdout.flush()
                plasmid_indices = Analysis.get_indices(pr_alignment[0])
                #print "plasmid_indices:", plasmid_indices
                read_indices = Analysis.get_indices(pr_alignment[1])
                #print "read_indices:", read_indices
                start = max(plasmid_indices[0], read_indices[0])
                stop = min(plasmid_indices[1], read_indices[1])
                #print "initially, stop was", stop
                key = pr_alignment[1].id
                if (key in mob_evidence_dict):
                    # index in reference sequence for mob element calculated by:
                        # adding start_index of subsequent analysis
                        # subtracting total inserted residues
                    strand = 1                        
                    stop = start + mob_evidence_dict[key][2] # offset should be from beginning
                    mob_ele_id = mob_evidence_dict[key][0]
                    if ("reverse_complement" in mob_ele_id):
                        temp = re.split("-", mob_ele_id)
                        mob_ele_id = temp[0]
                        strand = -1
                    # compensate for insertions before mobile element insertion to properly map to reference seq
                    validity = Analysis.analyze_seq(pr_alignment[0], pr_alignment[1], start, stop, mutation_list)
                    ins_count = 0
                    for mutation in mutation_list:
                        tmp = re.split("\t", mutation)
                        if ("INS" in tmp[0]):
                            ins_count += len(tmp[5]) # add length of any insertions
                            
                    mutation_string = ("MOB\t.\t.\t"+pr_alignment[0].id+"\t"+str(stop-ins_count) + "\t" + 
                                        mob_ele_id + "\t" + str(strand))
                    mutation_list.append(mutation_string)          

                    """As of 2015/8/26, we are no longer tracking mutations in the insertion region.                                      
                    
                    if (mob_evidence_dict[key][1]): # append any mutations within the mobile element
                        #mob_evidence_dict[key] is the list of mutations
                        for mutation in mob_evidence_dict[key][1]:
                            mutation_list.append(mutation)
                            
                    """
                
                mutation_list_2 = list() # will contain mutation in Sanger read
                validity = Analysis.analyze_seq(pr_alignment[0], pr_alignment[1], start, stop, mutation_list_2)
                for mutation in mutation_list_2: # mutation list at index 0
                    mutation_list.append(mutation)
                
                read_name = pr_alignment[1].id
                output_path = job.output_dirs["genomediff"]
                with open(output_path+read_name + ".gd", "w") as out_file:
                    out_file.write("#=GENOME_DIFF 1.0\n")
                    for mutation in mutation_list:
                        out_file.write(mutation + "\n")
                del(mutation_list) 
                total_complete += 1
        
        print "100% complete."
    
    return
def controller(job):
    """Controller for FRI-Failure Analysis pipeline."""
    # Using module Analysis  
    phase_1(job)
    phase_2(job)
    phase_3(job)
    return
    
def main():
    input_dirs = dict()    
    output_dirs = dict()
    input_dirs["mobile_elements"] = "mob_elements/" # location of mobile elements
    input_dirs["reads"] = "reads/" # location of Sanger reads
    input_dirs["plasmids"] = "plasmids/"
    output_dirs["templates"] = "templates/" # location of .fasta files that will be run through MAFFT
    output_dirs["alignments"] = "alignments/" # path to store and access alignments produced by MAFFT
    output_dirs["genomediff"] = "genomediff/" # path to store genomediff files
    myJob = jobmanager(input_dirs, output_dirs, "logfile.log", "Analysis")

    # Start pipeline
    controller(myJob)
    
main()
