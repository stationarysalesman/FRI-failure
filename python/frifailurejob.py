# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler
"""

    
from Bio import AlignIO
#Legacy import
#from Bio import SeqIO
from jobmanager import jobmanager
import os
import subprocess
import sys
import re
import argparse
from decimal import *


def phase_1(job):
    """FRI-Failure Analysis: Phase 1
    
    Phase 1 involves processing sequencing files. Each file written includes two sequences: a template 
    (first sequence) and a target (second sequence). Filenames containing "pr" denote plasmid/Sanger read
    files. File names containing "rm" denote read/mobile element files."""
    
    # Define necessary names   
    read_path = job.input_dirs["reads"]
    plasmid_path = job.input_dirs["plasmids"]
    template_path = job.output_dirs["templates"]
    mob_path = job.input_dirs["mobile_elements"]
    Analysis = job.process_module
    phred_cutoff = 10 #d efault
    if 'cutoff' in job.arg_map.keys() and job.arg_map['cutoff']:
        phred_cutoff = int(job.arg_map['cutoff'])
    Analysis.seqproc(read_path, plasmid_path, template_path, mob_path, phred_cutoff)
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
            """Calculate how much progress has been made and display for the user."""
            percent_complete = Decimal(total_complete)/Decimal(file_num)
            percent_complete = int(percent_complete * 100)
            print str(percent_complete) + "% complete...\r",            
            sys.stdout.flush()
            """Run MAFFT on all template files."""
            outfile_path = alignment_path+f+"_alignment"
            with open(outfile_path, "w") as out_file:  
               check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet", "--op", "3", template_path+f], stdout=out_file)
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
    
    """Create global mappings for this file."""
    alignment_path = job.output_dirs['alignments']
    curr_output_path = job.master_output_dir
    Analysis = job.process_module
    mob_evidence_dict = dict()
    
    """Check that we can access the directory containing alignments."""
    if not(os.access(alignment_path, os.F_OK)):
        err_str =  "Error: cannot access", alignment_path
        with open(job.logfile_name, "a") as log:
            log.write(err_str + "\n")
        return 1

    """Gathers mobile element homology by stepping through the alignment
    path and analyzing each alignment for mutations. Inside the loop,
    map a key (which is a sample ID) onto a named tuple that contains
    information about the mobile element with which the sample shares
    the highest degree of homology."""
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
               
              
                mobile_evidence_tuple = Analysis.mob_info(rm_alignments)
 
                    
                if not(evidence is None):
                    key = mobile_evidence_tuple['sangerID']
                    value = mobile_evidence_tuple
                    mob_evidence_dict[key]= value
                total_complete += 1
        print "100% complete."

    """This loop performs the actual analysis of samples against their
    respective template plasmids (i.e., this is where we generate most
    of the most important information to our future analyses).
    
    For each alignment file, we first check if we had previously identified
    homology between the sample and some mobile element. If so, we 
    truncate analysis at the start of the insertion."""
    print "Analyzing alignments of plasmids/reads:"
    for dirName, subdirList, fileList in os.walk(alignment_path):
        
        """Variables used to display progress to the user"""
        file_num = len(fileList)
        total_complete = 0
        
        for f in fileList:            
            if "pr" in f:  # plasmid/read alignment
                """Display current progress."""
                percent_complete = Decimal(total_complete)/Decimal(file_num)
                percent_complete = int(percent_complete * 100) # Estimate percentage complete
                print str(percent_complete) + "% complete...\r", # Move cursor to beginning

                """Begin analysis of sample and template alignment."""
                mutation_list = list()
                pr_alignment = AlignIO.read(alignment_path+f, "fasta") # Read sample file
                print pr_alignment[0].id + "\r",
                # TODO: test effects of removing flush() calls on output
                sys.stdout.flush()
                
                """Get indices at which to begin and terminate analysis.
                To determine the start, we use the maximum of the two
                sequences; for the end, we use the minimum of the two.
                This prevents mistakenly identifying mutations in the 
                portions of the alignment where the sequences do not
                overlap."""
                plasmid_indices = Analysis.get_indices(pr_alignment[0])
                read_indices = Analysis.get_indices(pr_alignment[1])
                start = max(plasmid_indices[0], read_indices[0])
                stop = min(plasmid_indices[1], read_indices[1])

                """Determine whether or not the sample exhibits substantial
                homology to mobile elements. If so, we can use the start
                of the insertion to truncate analysis, thereby preventing
                erroneous mutations from appearing in our final analysis."""
                key = pr_alignment[1].id
                if (key in mob_evidence_dict):
                    # index in reference sequence for mob element calculated by:
                        # adding start_index of subsequent analysis
                        # subtracting total inserted residues
                    strand = 1
                    """Compute the index at which to terminate analysis. This
                    should be based on two things:
                    1) The index (on the sample) of the insertion
                    2) The offset at which the sample first appears in this
                    (current) alignment file
                    Taking these two factors into account yields the index
                    at which the mobile element insertion occurs in the current
                    alignment.

                    However, when documenting any mutations, we are
                    interested in where those mutations occur in relation to 
                    the original template sequence only. To this end, we must
                    subtract from the index obtained above the count of insertions
                    that occur in our sample prior to said mobile element insertion."""
                    stop = start + mob_evidence_dict[key]['start']
                    mob_ele_id = mob_evidence_dict[key]['elementID']
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
                else:
                    validity = Analysis.analyze_seq(pr_alignment[0], pr_alignment[1], start, stop, mutation_list)


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
    # Specify argument storage variables
    phred_cutoff = 10 # default
    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutoff", type=int, help="average PHRED score (over 3 nts) to cutoff Sanger reads")
    args = parser.parse_args()
    # User has specified a custom cutoff
    if args.cutoff:
        phred_cutoff = args.cutoff
    input_dirs = dict()    
    output_dirs = dict()
    input_dirs["mobile_elements"] = "../mob_elements/" # location of mobile elements
    input_dirs["reads"] = "../reads/" # location of Sanger reads
    input_dirs["plasmids"] = "../plasmids/"
    output_dirs["templates"] = "templates/" # location of .fasta files that will be run through MAFFT
    output_dirs["alignments"] = "alignments/" # path to store and access alignments produced by MAFFT
    output_dirs["genomediff"] = "genomediff/" # path to store genomediff files
    myJob = jobmanager(input_dirs, output_dirs, "logfile.log", "Analysis", vars(args))

    # Start pipeline
    controller(myJob)
    
main()
