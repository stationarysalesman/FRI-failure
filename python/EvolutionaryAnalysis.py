# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:21:00 2015

@author: tyler
FRI-Failure Analysis Pipeline
Copyright Tyler Camp 2015
module: EvolutionaryAnalysis
"""
from Bio import SeqIO
from Bio import AlignIO
import os
import re
import sys
from decimal import *
from collections import namedtuple
from Mutation import Mutation
from Analysis import Analysis
from datetime import datetime
import time
import errno
import random


class EvolutionaryAnalysis(Analysis):
    """EvolutionaryAnalysis: A subclass of Analysis used to perform evolution
    and stability based analyses.

    Unique instance variables:

    ancestor_seq_filepath: location on disk of ancestor sequence file
    sample_lst: list of sample sequence file paths"""

    def __init__(self, input_dir, output_dir):
        self.output_dir = input_dir
        self.input_dir = output_dir

        # Default values for instance variables
        self.ancestor_seq_path = None
        self.sample_lst = None
        self.tmp_dir = self.init_temp_dir(type(self).__name__)

    def start_analysis(self):
        """Start the analysis for an Evolutionary study"""

        """Identify sequences to be analyzed"""
        file_list = list()
        for d, sd, fl in os.walk(self.input_dir):
            for f in fl:
                if ".gb" in f:
                    self.ancestor_seq_path = self.input_dir + "/" + f
                elif ".fasta" in f or ".ab1" in f:
                    self.sample_lst.append(f)
                else:
                    print "Ignoring file: ", f

        """Create miscellaneous temporary directories"""
        os.mkdir(self.tmp_dir+"templates/")
        os.mkdir(self.tmp_dir+"alignments/")
        os.mkdir(self.tmp_dir+"output/")

        """Begin work flow"""




        print "Success."
        return





    """input_dir: director containing input"""

    _input_dir = None

    def get_input_dir(self):
        return self._input_dir

    def set_inputdirs(self, val):
        self._input_dir = val

    input_dir = property(get_input_dir, set_inputdirs)

    """output_dir: directory containing output"""

    _output_dir = None

    def get_output_dir(self):
        return self._outputDir

    def set_output_dir(self, val):
        self._output_dir = val

    output_dir = property(get_output_dir, set_output_dir)

    """log file"""
    @property
    def log_file(self):
        return None

    def write_log(self, data):
        return None

    """Analysis metadata"""

    @property
    def start_time(self):
        return None

    @property
    def end_time(self):
        return None


    """Ancestor sequence (Biopython Sequence object)"""
    _ancestor_seq_obj = None
    @property
    def ancestor_seq_obj(self):
        return self._ancestor_seq_obj



"""Define the collection of data that carries information about mobile
element insertions."""

MobileElementInfo = namedtuple('MobileElementInfo',
                               ['sangerID', #id of Sanger read 
                                'elementID', #id of mobile element
                                'startIndex', #index the mobile element inserted 
                                'validity', #proportion of correctly matched n.t.
                                'mutations']) #list of differences bt. read & mob
"""Declare flags."""
NO_CALL = 1
NO_MUT = 2

def trim_n(seq, phred_cutoff):
    """Trim the 5' and 3' ends of Sanger reads."""
    # Trim 5' end
    i = 0
    skip = False
    done = False

    while (i < len(seq.seq)) and not(done):
        if (((seq.letter_annotations.values()[0][i]+
              seq.letter_annotations.values()[0][i+1]+
              seq.letter_annotations.values()[0][i+2])/3)
              > 30):
                 seq = seq[i:]
                 done = True
                 break
        i += 1

                
    # Trim 3' end
    j = i
    done = False
    while (j < len(seq.seq)) and not(done):
        if (((seq.letter_annotations.values()[0][j]+
              seq.letter_annotations.values()[0][j+1]+
              seq.letter_annotations.values()[0][j+2])/3) 
              < phred_cutoff):
                 seq = seq[i:j]
                 done = True
                 break
        j += 1
    return seq

def seqproc(read_path, plasmid_path, template_path, mob_path, phred_cutoff):
    """Create templates that will be later aligned with MAFFT"""
    
    # Create list containing mobile elements
    mob_list = list()
    for dirName, subdirList, fileList in os.walk(mob_path):
        for f in fileList:
            mob_seq = SeqIO.read(dirName+f, "genbank")
            mob_rc = mob_seq.reverse_complement()
            mob_rc.id = mob_seq.id + "-reverse_complement"
            mob_list.append(mob_seq)
            mob_list.append(mob_rc)
            
    for dirName, subdirList, fileList in os.walk(plasmid_path):
        for plasmid_file in fileList:
            read_list = list()
            initials = ""
            plasmid = SeqIO.read(plasmid_path+plasmid_file, "genbank") # object containing plasmid info     
            temp = plasmid_file   
            dest = []
            dest = re.split("_", temp) #get initials from plasmid file name         
            initials = dest[0].lower()
            for dirName2, subdirList2, fileList2 in os.walk(read_path):
                for read_file in fileList2:
                    if (initials in read_file.lower()):
                        out_file = SeqIO.read(read_path+read_file, "abi")
                        trimmed_file = trim_n(out_file, phred_cutoff) # trim n's
                        read_list.append(trimmed_file)
                   
    # At this point we have a list of all Sanger reads corresponding to our current plasmid file.
    # Ends have been trimmed of N's.            
            
            for i, sequence in enumerate(read_list):
                for j, mob_seq in enumerate(mob_list):
                    temp_lst = [sequence, mob_seq]
                    SeqIO.write(temp_lst, template_path+"rm"+str(i)+str(j)+"_"+initials, "fasta")
                    del(temp_lst)
            for i, sequence in enumerate(read_list):
                temp_lst = [plasmid, sequence]
                SeqIO.write(temp_lst, template_path+"pr"+str(i)+"_"+initials, "fasta")
                del(temp_lst)
    return
    
def get_seq_list(seq_path):
    """Return a list of Seq objects from GenBank files located in seq_path."""
    seq_list = []
    for dirName, subdirList, fileList in os.walk(seq_path):
        for f in fileList:
            seq_list.append(SeqIO.read(dirName+f, "genbank"))
    return seq_list
def get_indices(seq):
    """Determine start and stop indices for analysis.
    
    Obtain start and stop indices for analyzing a single alignment
    within a multiple sequence alignment file. This sequence should be
    smaller than your template."""
    start = 0
    stop = len(seq.seq)-1
    while ((seq.seq[stop] == "-") and (stop > 0)):
        stop -= 1
    if (stop == 0):
        print "ERROR: stop_index = 0 for sequence", seq.id
        return
    while (seq.seq[start] == "-") and (start < len(seq.seq)):
        start += 1
    if (start >= len(seq.seq)):
        print "ERROR: start >= length for sequence", seq.id
        return
    if (start > stop):
        print "ERROR: start > stop for sequence", seq.id
        return
    if (stop > len(seq.seq)):
        print "ERROR: stop > len(seq.seq) for sequence", seq.id
    #print "returning", start, stop, "for seq", seq.id
    return [start, stop]
def get_sub_rep(seq_frag):
    """Return the length of the smallest repeated sequence in the
       sequence fragment, and how many times the sequence is repeated."""
    repeat_num = 1
    for x in range(1, len(seq_frag)):
        i = x
        j= 0
        while (i<len(seq_frag) and str(seq_frag[j:i]) == str(seq_frag[i:i+x])):
            repeat_num += 1
            i += x
            j += x
        if (repeat_num * x == len(seq_frag)):
            return [seq_frag[:x], repeat_num]
        else:
            repeat_num = 1
            
    return [seq_frag, 1]
    
 # determine number of repeats, if any
def build_repeat_string(seq, x, mut_length, ins, tandem_seq):
    """Determine number of tandem repeats in seq, starting at x."""
    
    """Starting from seq[x], determine the number of repeated units
    of length mut_length directly downstream or upstream of
    seq[x:x+mut_length].
    arg: ins is boolean (1 -> insertion, 0 -> deletion"""
    mut_str = ""
    seq_str = str.replace("-", "", str(seq))           
    sub_rep = get_sub_rep(seq_str[x:x+mut_length])
    repeat_seq = sub_rep[0]
    delta_units = sub_rep[1] # number of repeated units deleted/inserted
    # determine if read sequence contains repeats
    # now we must find the bounds of the read repeat sequence
    repeat_len = len(repeat_seq)    
    start_index = x
    stop_index = x+repeat_len
    skip = mut_length-1
    
    """This is an error condition and will cause the process to exit."""
    if (repeat_len <= 0): 
        sys.exit("ERROR: repeat_len <= 0")
    
    """First, we move toward the 5' end by repeat_len each iteration"""
    while (str(seq_str[start_index:start_index+repeat_len])
    == str(seq_str[start_index-repeat_len:start_index]) and start_index>0):
        start_index -= repeat_len
    
    """Now, we move toward the 3' end in the same manner"""
    while (str(seq_str[stop_index-repeat_len:stop_index])
    == str(seq_str[stop_index:stop_index+repeat_len]) and
           (stop_index+repeat_len)<len(seq_str)):
        stop_index += repeat_len
        
    repeat_ref_num = (stop_index-start_index)/repeat_len
    if (repeat_ref_num == 1 and delta_units == 1):
        """Not repeat mediated."""
        return
    elif (repeat_len * repeat_ref_num >= 5):
        """Meets GenomeDiff criteria for RMD/RMI"""
        repeat_new_copies = 0
        if (ins):
            repeat_new_copies = repeat_ref_num + delta_units
        else:
            repeat_new_copies = repeat_ref_num - delta_units
            
        tandem_seq = repeat_seq
        mut_str = ("\trepeat_seq="+ repeat_seq + "\trepeat_len="+str(repeat_len)+
                  "\t"+"repeat_ref_num="+str(repeat_ref_num)+
                  "\trepeat_new_copies="+str(repeat_new_copies))
        return mut_str
    else:
        return


def identify_disjoint_rmd(template_seq, x, del_length):
    """Determine whether or not a disjoint repeat mediated deletion 
    occurred at this location in the template.

    There are two ways for the sequences to align if this occurs. 
    The alignment program can either place the deleted region first, 
    or it can align the complement of the repeat sequence first.
    The ASCII art below gives examples of both cases. Repeat sequences
    are bracketed.

    Case 1: Deletion first
    Template: [ctagag]-----[ctagag]tag
    Sample:    ------ ----- ctagag tag

    Case 2: One copy of repeat sequence first
    Template: [ctagag]-----[ctagag]tag
    Sample:   [ctagag]----- ------ tag

    This function will account for both possibilities: if the deletion 
    occurs first, we build the repeat sequence starting at the 5' end 
    and compare to the region just downstream of the deletion region, 
    building the length of the repeat sequence one at a time up to a 
    reasonable limit. If a copy of the repeat sequence is aligned first,
    we build the sequence starting at the 3' end of the deletion region,
    and we compare to the region just upstream of the deletion region."""

    repeat_length = 3
    repeat_seq = None

    """If the deletion length is less than our minimum repeat length, 
    this algorithm will fail with false positives, so we terminate if 
    the deletion length is too small."""
    if (del_length < 3):
        return
    """Case 1: deletion aligned first. Start building from the 5' end 
    (inside the deletion region) and compare to the region just downstream 
    of the deletion region."""

    while (repeat_length < 10):
        upstr_index = x
        downstr_index = x + repeat_length
        candidate_seq = template_seq[upstr_index:downstr_index]
        compare_upstr_index = x + del_length
        compare_downstr_index = compare_upstr_index + repeat_length
        compare_seq = template_seq[compare_upstr_index:compare_downstr_index]
        if str(candidate_seq) == str(compare_seq):
            repeat_seq = candidate_seq
            repeat_length += 1
        else:
            break

    if not(repeat_seq is None):
        return repeat_seq


    """Case 2: a copy of the repeat sequence aligned first. Start building 
    from the 3' end (inside the deletion region) and compare to the region 
    just upstream of the deletion region."""
    
    
    """The repeat sequence should be, at minimum, 3 nucleotides in length."""
    repeat_length = 3
    while (repeat_length < 10):
        downstr_index = x+del_length
        upstr_index = downstr_index-repeat_length
        candidate_seq = template_seq[upstr_index:downstr_index]
        compare_upstr_index = x-repeat_length
        compare_downstr_index = x
        compare_seq = template_seq[compare_upstr_index:compare_downstr_index]
        if str(candidate_seq) == str(compare_seq):
            repeat_seq = candidate_seq
        else:
            break

    
    if not(repeat_seq is None):
        return repeat_seq
    
    return
    

def eval_del(mutObject, template, target, x, ref_index):
    """Evaluate a previously identified deletion.

    At entry, the caller has identified that there is a deletion at 
    index x in the alignment. Assign a Mutation object with 
    information based on what kind of deletion it is."""

    mutObject.set_mutation_type("DEL")

    template_seq = template.seq
    template_id = template.id
    target_seq = target.seq

    """Determine the length of the deletion and store it in the
    Mutation object."""
    z = x
    while (z < len(target_seq) and target_seq[z] == "-"):
        z += 1
    del_length = z - x    
    mutObject.set_length(del_length)

    """Build the initial output string."""
    mutation_str = ("DEL\t.\t.\t"+ template_id + "\t" +
                    str(ref_index) + "\t" + str(del_length))

    """Determine if this is a RMD, which can be one of two cases:
    1) tandem repeats
    2) between disjoint repeats
    These cases are mutually exclusive."""

    """Determine if this is an RMD caused by tandem repeats."""
    tandem_seq = None

    """tandem_rmd_str holds a GenomeDiff mutation string if there is
    a tandem RMD here."""
    tandem_rmd_str = build_repeat_string(template_seq, x, del_length, 0,
                                        tandem_seq)
    """If tandem repeats, we append the new information."""
    if (tandem_rmd_str):
        mutation_str += tandem_rmd_str
    """Regardless of status, set the mutation string of the Mutation
    object."""
    mutObject.set_string(mutation_str)
        
    """Determine if this is a disjoint repeat mediated deletion.
    Mutually exclusive with tandem repeats."""
    if not(tandem_rmd_str):
        disj_rmd_seq = identify_disjoint_rmd(template_seq, x, del_length)
        if (disj_rmd_seq):
            errata = "\tbetween="+str(disj_rmd_seq)
            mutObject.set_string(mutation_str+errata)
    else:
        mutObject.set_string(mutation_str)                    

    return 0

def eval_ins(mutObject, template, target, x, ref_index):
    """Evaluate a previously identified insertion.

    At entry, the caller has identified that there is an insertion at 
    index x in the alignment. Assign a  Mutation object with 
    information based on what kind of insertion it is."""


    mutObject.set_mutation_type("INS")
    
    template_seq = template.seq
    template_id = template.id
    target_seq = target.seq
    z = x
    ignore_count = 0
    while (z < len(template_seq) and template_seq[z] == "-"):
        z += 1
    ins_length = z - x

    mutation_str = ("INS\t.\t.\t" + template_id + "\t" +
                       str(ref_index) + "\t" +
                       str(target_seq[x:x+ins_length]))
    tandem_seq = None
    errata = None

    """tandem_rmi_str holds a GenomeDiff mutation string if there is
    a tandem RMI here."""
    tandem_rmi_str = build_repeat_string(target_seq, x, ins_length, 1,
                                 tandem_seq)
    """If tandem repeats, we append the new information."""
    if (tandem_rmi_str):
        mutation_str += tandem_rmi_str
    """Regardless of status, set the mutation string of the Mutation
    object."""
    mutObject.set_string(mutation_str)

    return 0


def eval_snp(mutObject, template, target, x, ref_index):
    """Evaluate a previously identified SNP.

    At entry, the caller has identified that there is a SNP at 
    index x in the alignment. Assign a  Mutation object with 
    information based on the SNP's location."""

    mutObject.set_mutation_type("SNP")
    template_id = template.id
    target_nt = target.seq[x]
    mutObject.set_string("SNP\t.\t.\t"+template_id+"\t"+str(ref_index)+"\t"+target_nt)

    return 0

def eval_nt(template, target, x, ref_index):
    """Evaluate whether or not a mutation occurred here,
    and call helpers to identify them.

    This function returns a Mutation object that contains 
    information relevant to the identified mutation. All 
    information needed by the caller should be returned 
    in this object."""

    template_nt = template.seq[x]
    target_nt = target.seq[x]
    """Some kind of discrepancy here."""
    if (target_nt != template_nt):

        """Define Mutation object that will be passed to evaluation 
        functions."""
        mutObject = Mutation()
        """Sequencing could not reach a consensus for any nucleotide here,
        so we don't count these as mutations."""
        if (target_nt == "n"):
            del(mutObject)
            return NO_CALL
        elif (target_nt == "-" and template_nt != "-"):
            """Deletion"""
            check = eval_del(mutObject, template, target, x, ref_index)
            return mutObject
        elif (template_nt == "-" and target_nt != "-"):
            """Insertion"""
            check = eval_ins(mutObject, template, target, x, ref_index)
            return mutObject
        else:
            check = eval_snp(mutObject, template, target, x, ref_index)
            return mutObject
            
    else:
        return NO_MUT
        
def analyze_seq(template, target, start_index, stop_index, mutation_list):
    """Analyze an alignment between template and target
    between the given indices."""
    template_seq = template.seq
    target_seq = target.seq
    template_id = template.id
    if ("reverse_complement" in template_id):
        template_id = (re.split("-", template_id))[0]
 
    
    ins_count = 0 # total number of nts inserted  
    del_count = 0 # total number of nts deleted    
    skip = 0
    
    true_count = 0
    for x in range(start_index, stop_index):   
        # calculate "true index" to map mutations onto reference sequence
        # only things that affect index into reference are start index and total number of inserted bases
        ref_index = x - ins_count
        if (skip):
            skip -= 1
            continue
        template_nt = template.seq[x]
        target_nt = target.seq[x]
            
        """Evaluate the alignment at index x, identifying any mutations 
        and returning them."""
        mutObject = eval_nt(template, target, x, ref_index)

        """Did a mutation occur?"""
        if (mutObject == NO_MUT):
            true_count += 1
            continue
        """Was there a base call for the target sequence at this location?"""
        if (mutObject == NO_CALL):
            continue    
        
        mut_type = mutObject.get_mutation_type()
        """Do we need to skip nucleotides following this mutation?
        We can calculate this for in/dels by subtracting 1 from the 
        length of the in/del."""
        if (mut_type == "DEL" or mut_type == "INS"):
            skip = mutObject.get_length() - 1

        """Update del_count and ins_count"""
        if (mut_type == "DEL"):
            del_count += 1
        elif (mut_type == "INS"):
            ins_count += 1
            
        """Add the Mutation object's string field to the mutation list."""
        mutation_list.append(mutObject.get_string())
        
        """Delete the object."""
        del(mutObject)

    validity = Decimal(true_count)/Decimal(stop_index-start_index) # percent nts that match in analysis frame
    return validity

def mob_info(alignments):
    """Returns a list of info about an alignment analysis based on
    sequence homology to mobile elements."""
    
    read = alignments[0] # This is the Sanger read
    
    mob_evidence_list = list()
    for r in range(1, len(alignments)):
        mob = alignments[r] # this is the current mobile element we want to analyze    
        # Analyze the alignment only in the region of the mobile element
        mob_indices = get_indices(mob)
        read_indices = get_indices(read)
        # We want indices that capture only portions of the alignment where both sequences
        # are aligned. This may overlook edge insertions/deletions.
        start = max(mob_indices[0], read_indices[0])
        stop = min(mob_indices[1], read_indices[1])
        if ((stop - start) < 50): # not long enough for any meaningful analysis
            return
        # call analyze_seq on our read/mobile element to get info about sequence homology
        current_element_mutation_list = list()
        validity = analyze_seq(mob, read, start, stop, current_element_mutation_list)
        current_mobile_element_info = MobileElementInfo(read.id, mob.id, start, validity,
                                                        current_element_mutation_list)
        mob_evidence_list.append(current_mobile_element_info)
        del current_element_mutation_list 
    curr_info_element = ""
    curr_len = 9999999 # do not use this value for "big data" alignments/analyses

    """Shouldn't attempt to access items in empty list"""
    if not(mob_evidence_list):
        return
    """The goal is to return the strongest possible evidence of a mobile
    element insertion. To this end, we must do two things:
    1) Identify the collection of information which has the fewest differences
       between its sample and mobile element
    2) Verify that this pair actually exhibits a significant degree of homology"""
    for info_element in mob_evidence_list: # find the mobile element with the least number of mutations      
        if (len(info_element.mutations) < curr_len):
            curr_len = len(info_element.mutations)
            curr_info_element = info_element      
    if (curr_info_element and curr_info_element.validity >= .45): # is the mobile element actually homologous?
        del(mob_evidence_list)
        return curr_info_element
    
    return                                     

"""Build a dictionary of information relating to mobile element
homology, such that we can determine if a given sample exhibits
substantial homology to any mobile elements (note: only 
informatino about the mobile element which is most similar to 
the sample is stored in the helper function mob_info). 

The dictionary that this function returns can be indexed by 
sample ID, i.e. each sample ID is a key in the dictinoary,
which stores a MobileElementInfo namedtuple."""
def build_mob_evidence_dict(alignment_path):
    mob_evidence_dict = dict()
    for dirName, subdirList, fileList in os.walk(alignment_path):
        for f in fileList:
            if "rm" in f:
                # read/mobile element alignment
                # need to determine if any homology to do nt analysis
                rm_alignments = AlignIO.read(alignment_path+f, "fasta")
                print ("Gathering mobile element homology data for read " +
                rm_alignments[0].id + " and " + rm_alignments[1].id + "..."),
                mobile_element_evidence = mob_info(rm_alignments)
                if not(mobile_element_evidence is None):
                    key = mobile_element_evidence['sangerID']
                    value = mobile_element_evidence
                    mob_evidence_dict[key]= value
                print "done."
    return mob_evidence_dict
