# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler
"""

    
from Bio import AlignIO
import EvolutionaryAnalysis
from jobmanager import jobmanager
import os
import subprocess
import sys
import re
import argparse
from multiprocessing import Process

from decimal import *

"""Global flags."""

"""Number of logical cores. Note: on modern architectures, each physical 
core may have multiple logical cores. These logical cores can fetch 
and execute multiple instruction streams in parallel by combining multiple 
ALUs or other execution circuitry onto a shared core."""
NUM_CORES = 8

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
    phred_cutoff = 10 #default
    if 'cutoff' in job.arg_map.keys() and job.arg_map['cutoff']:
        phred_cutoff = int(job.arg_map['cutoff'])
    Analysis.seqproc(read_path, plasmid_path, template_path, mob_path, phred_cutoff)
    # Add templates to input dictionary
    job.input_dirs["templates"] = job.output_dirs["templates"]   
    return


def create_alignment(file_list, infile_path, alignment_path, log_path):
    """Create several alignments file_list, and output them to outfile_path.

    This function handles the generation of alignments with MAFFT, and options 
    should be changed here. This function runs concurrently on different cores
    to achieve higher throughput."""
    for f in file_list:
        with open(alignment_path+f, "w") as out_file:
           check = subprocess.call(["/usr/lib/mafft/bin/mafft", "--quiet",
                                    infile_path+f], stdout=out_file)
           sys.stdout.flush()             
           if (check):
               with open(log_path, "a") as log:
                   log.write("Mafft error (input=" + f + ", output=" + outfile_path + ")\n")


def phase_2(job):
    """FRI-Failure Analysis: Phase 2
    
    Phase 2 creates the alignments with the MAFFT alignment program. File input comes from files written
    during phase 1."""

    """Define necessary directories"""    
    alignment_path = job.output_dirs["alignments"]
    template_path = job.input_dirs['templates']
    log_path = job.logfile_name
    
    print "Creating alignments with MAFFT."
    for dirName, subdirList, fileList in os.walk(template_path):
        file_list = fileList

    """Partition the list of files. Each core should get roughly
    the same number of files to process."""
    partition_size = len(file_list)/NUM_CORES
    """Create a list of file lists."""
    partition_list = list()
    index = 0

    """Partition the list by multiplying an index by the partition size to 
    get the starting index, and add the partition_size to get the end 
    index. If the number of files does not divide evenly by NUM_CORES, 
    we must add the tail end of the file list outside of the loop."""
    while (index < NUM_CORES):
        start_index = index*partition_size
        partition_list.append(file_list[start_index:start_index+partition_size])
        index += 1
    if (len(file_list) % NUM_CORES):
        """Allocate the leftover files over all the processors."""
        start_index = index*partition_size
        leftover_list = file_list[start_index:len(file_list)]
        mod_index = 0
        for f in leftover_list:
            """Modulo the index to wrap around the file_list."""""
            partition_list[mod_index % len(file_list)].append(f)
            mod_index += 1

    assert not(len(partition_list) % NUM_CORES)

    """Create list of processes so we can join them."""
    proc_list = list()
    
    """Create a new process for each partition of the file list."""
    for template_list in partition_list:        
        proc = Process(target=create_alignment,
                       args=(template_list, template_path, alignment_path,
                             log_path))
        proc_list.append(proc)
        proc.start()

    """Join the processes."""
    for proc in proc_list:
        proc.join()

    """Add alignments directory to the input directories map."""
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
    EvolutionaryAnalysis = job.process_module
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
        """Build necessary file lists for analysis."""
        plasmid_read_filelist = list()
        read_mobile_filelist = list()
        for f in fileList:
            if "pr" in f:
                plasmid_read_filelist.append(f)
            elif "rm" in f:
                read_mobile_filelist.append(f)
                
    analyze_mobile_homology(read_mobile_filelist, curr_output_path, alignment_path, mob_evidence_dict)
    
    print "Analyzing alignments of plasmids/reads:"
    mutation_lst = parse_samples(plasmid_read_filelist, curr_output_path, alignment_path, mob_evidence_dict)

    return


def analyze_mobile_homology(fileList, output_path, alignment_path, mob_evidence_dict):
    """Gather data relating to homology between samples and mobile elements.

    This function uses a namedtuple to store information about homology. After calling
    helper function mob_info(), update the mob_evidence_dict (a dictionary) with keys 
    (sample names) and values (namedtuple containing homology data)."""
    file_num = len(fileList)
    total_complete = 0  
    for f in fileList:           
        """Display current progress."""
        percent_complete = Decimal(total_complete)/Decimal(file_num)
        percent_complete = int(percent_complete * 100)
        print str(percent_complete) + "% complete...\r",
        """Here, we determine mobile element homology to inform our
        analysis later in the work flow."""
        rm_alignments = AlignIO.read(alignment_path+f, "fasta")                

        """Generate a named tuple that contains information about
        the strongest (if any) mobile element insertion."""
        mobile_evidence_tuple = EvolutionaryAnalysis.mob_info(rm_alignments)

        if not(mobile_evidence_tuple is None):
            key = mobile_evidence_tuple.sangerID
            value = mobile_evidence_tuple
            mob_evidence_dict[key]= value
            
        total_complete += 1
    print "100% complete."
    
    return


def parse_samples(fileList, output_path, alignment_path, mob_evidence_dict):
    """Traverse the list of alignments containing plasmids and reads and 
    analyze them for mutations.

    This function utilizes a helper function analyze_sample to perform the 
    analysis. Then, parse_samples outputs the results to a file named for 
    the current sample."""
    
    """Variables used to display progress to the user"""
    file_num = len(fileList)
    total_complete = 0
    """Build a list of files that contain alignments between plasmids
    and reads."""

    for f in fileList:
        """Display current progress."""
        percent_complete = Decimal(total_complete)/Decimal(file_num)
        percent_complete = int(percent_complete * 100) # Estimate percentage complete
        print str(percent_complete) + "% complete...\r", # Move cursor to beginning
        """Define names for file I/O"""
        mutation_list = list()
        
        outfile_name = analyze_sample(alignment_path+f, mob_evidence_dict, mutation_list)
        
        """At this point, we have the list of mutations. Now we write them to disk."""
        if not outfile_name:
            print "ERROR: no outfile name"
            continue

        with open(output_path+outfile_name+ ".gd", "w") as outfile:
            outfile.write("#=GENOME_DIFF 1.0\n")
            for mutation in mutation_list:
                outfile.write(mutation + "\n")
        total_complete += 1
    return



def analyze_sample(infile_path, mob_evidence_dict, mutation_list):
    """This function performs the actual analysis of samples against their
    respective template plasmids (i.e., this is where we generate most
    of the most important information to our future analyses).

    For each alignment file, we first check if we had previously identified
    homology between the sample and some mobile element. If so, we 
    truncate analysis at the start of the insertion.

    This function assigns two of its arguments, outfile_name and mutation_list,
    to new values. The caller expects them to have a certain state. Modify these 
    values with care."""
    
    """Begin analysis of sample and template alignment."""
    pr_alignment = AlignIO.read(infile_path, "fasta") # Read sample file
    outfile_name = pr_alignment[1].id

    """Get indices at which to begin and terminate analysis.
    To determine the start, we use the maximum of the two
    sequences; for the end, we use the minimum of the two.
    This prevents mistakenly identifying mutations in the 
    portions of the alignment where the sequences do not
    overlap."""
    plasmid_indices = EvolutionaryAnalysis.get_indices(pr_alignment[0])
    read_indices = EvolutionaryAnalysis.get_indices(pr_alignment[1])
    start = max(plasmid_indices[0], read_indices[0])
    stop = min(plasmid_indices[1], read_indices[1])

    """Determine whether or not the sample exhibits substantial
    homology to mobile elements. If so, we can use the start
    of the insertion to truncate analysis, thereby preventing
    erroneous mutations from appearing in our final analysis."""
    key = pr_alignment[1].id
    if (key in mob_evidence_dict):
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
        stop = start + mob_evidence_dict[key].startIndex
        mob_ele_id = mob_evidence_dict[key].elementID
        analyze_mobile_ins(mob_ele_id, pr_alignment[0], pr_alignment[1], start, stop, mutation_list) 
    else: 
        validity = EvolutionaryAnalysis.analyze_seq(pr_alignment[0], pr_alignment[1], start, stop, mutation_list)

    return outfile_name
    
    
    


def analyze_mobile_ins(mob_ele_id, template, sample, start, stop, mutation_list):
    """Analyze a sequence, compensating for mobile element insertions.
    This function finalizes parameters for analysis and proceeds with the analyze function. 
    This function performs no I/O; it is left to the caller."""
    strand = 1    
    if ("reverse_complement" in mob_ele_id):
        temp = re.split("-", mob_ele_id)
        mob_ele_id = temp[0]
        strand = -1
    """Compensate for insertions before mobile element insertion to properly map 
    to reference seq"""
    validity = EvolutionaryAnalysis.analyze_seq(template, sample,
                                                start, stop, mutation_list)
    ins_count = 0
    for mutation in mutation_list:
        tmp = re.split("\t", mutation)
        if ("INS" in tmp[0]):
            ins_count += len(tmp[5]) # add length of any insertions

    mutation_string = ("MOB\t.\t.\t"+template.id+"\t"+str(stop-ins_count) + "\t" + 
                       mob_ele_id + "\t" + str(strand))
    mutation_list.append(mutation_string)       

    return
        

    


def controller(job):
    """Controller for FRI-Failure Analysis pipeline."""
    
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
    
if __name__ == "main":
    main()
