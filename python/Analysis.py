# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:21:00 2015

@author: tyler
FRI-Failure Analysis Pipeline
Copyright Tyler Camp 2015
module: Analysis
"""
from Bio import SeqIO
from Bio import AlignIO
import os
import re
import sys
from decimal import *
from collections import namedtuple
from Mutation import Mutation

"""Define the collection of data that carries information about mobile
element insertions."""

MobileElementInfo = namedtuple('MobileElementInfo',
                               ['sangerID', #id of Sanger read 
                                'elementID', #id of mobile element
                                'startIndex', #index the mobile element inserted 
                                'validity', #proportion of correctly matched n.t.
                                'mutations']) #list of differences bt. read & mob

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
def build_repeat_string(seq, x, mut_length, ins, tandem_seq, mut_str):
    """Determine number of tandem repeats in seq, starting at x."""
    
    """Starting from seq[x], determine the number of repeated units
    of length mut_length directly downstream or upstream of
    seq[x:x+mut_length].
    arg: ins is boolean (1 -> insertion, 0 -> deletion"""

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
        return -1
    
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
        return 0    
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
        return 1
    else:
        return 0


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
    Sample:    ------ ----- gatctc atc

    Case 2: One copy of repeat sequence first
    Template: [ctagag]-----[ctagag]tag
    Sample:   [gatctc]----- ------ atc

    This function will account for both possibilities: if the deletion 
    occurs first, we build the repeat sequence starting at the 5' end 
    and compare to the region just downstream of the deletion region, 
    building the length of the repeat sequence one at a time up to a 
    reasonable limit. If a copy of the repeat sequence is aligned first,
    we build the sequence starting at the 3' end of the deletion region,
    and we compare to the region just upstream of the deletion region."""

    repeat_length = 3
    repeat_seq = None
    """Case 1: deletion aligned first. Start building from the 5' end 
    (inside the deletion region) and compare to the region just downstream 
    of the deletion region."""

    while (repeat_length < 10):
        upstr_index = x
        downstr_index = x + repeat_length
        candidate_seq = template_seq[upstr_index:downstr_index]
        compare_upstr_index = downstr_index
        compare_downstr_index = downstr_index + repeat_length
        compare_seq = template_seq[compare_upstr_index:compare_downstr_index]
        if candidate_seq == compare_seq:
            repeat_seq = candidate_seq
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
        if candidate_seq == compare_seq:
            repeat_seq = candidate_seq
        else:
            break

    
    if not(repeat_seq is None):
        return repeat_seq
    
    return
    

def eval_del(template, target, x, ref_index, STOP):
    """Evaluate a previously identified deletion.

    At entry, the caller has identified that there is a deletion at 
    index x in the alignment. Create a Mutation object with 
    information based on what kind of deletion it is."""

    """Create Mutation object to store mutation data."""
    mutObject = Mutation()
    mutObject.set_type("DEL")

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
    errata = None
    is_tandem_rmd = build_repeat_string(template_seq, x, del_length, 0,
                                        tandem_seq, errata)
    
    """Error condition. For now, we want to terminate the process."""
    if (is_tandem_rmd < 0):
        sys.exit("Error: repeat_len <= 0")
    elif (is_tandem_rmd):
        mutation_str += errata
        mutObject.set_string(mutation_str)
        
    """Determine if this is a disjoint repeat mediated deletion"""
    disj_rmd_seq = identify_disjoint_rmd(template_seq, x, del_length)
    if (disj_rmd_seq):
        errata = "\tbetween="+disj_rmd_seq
        mutObject.set_string(mutation_str+errata)
    else:
        mutObject.set_string(mutation_string)                    
    
    """Check if analysis must stop in caller"""
    if ((skip + x) >= stop_index):
        STOP = 1
           

    return mutObject

def eval_ins(template, target, x, ref_index, STOP):
    """Evaluate a previously identified insertion.

    At entry, the caller has identified that there is an insertion at 
    index x in the alignment. Create a  Mutation object with 
    information based on what kind of insertion it is."""

    """Create Mutation object to store mutation data."""
    mutObject = Mutation()
    mutObject.set_type("INS")
    
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
    is_tandem_rmi = build_repeat_string(target_seq, x, ins_length, 1,
                                 tandem_seq, errata)
    
    """Error condition. For now, we want to terminate the process."""
    if (is_tandem_rmi < 0):
        sys.exit("Error: repeat_len <= 0")
    elif (is_tandem_rmi):
        mutation_str += errata
        mutObject.set_string(mutation_str)
    
    if (errata):
        mutation_list.append(mutation_string+errata)
    else:
        mutation_list.append(mutation_string)
    skip = errata[1]
    if ((skip + x) >= stop_index):
       return Decimal(true_count)/Decimal(stop_index-start_index)
def eval_nt(template, target, x, ref_index, STOP):
    """Evaluate whether or not a mutation occurred here,
    and call helpers to identify them.

    This function returns a Mutation object that contains 
    information relevant to the identified mutation. All 
    information needed by the caller should be returned 
    in this object."""

    template_nt = template_seq[x]
    target_nt = target_seq[x]
    """Some kind of mutation occurred here."""
    if (target_nt != template_nt):
            
        """Sequencing could not reach a consensus for any nucleotide here,
        so we don't count these as mutations."""
        if (target_nt == "n"):
            del(mutObject)
            return
        elif (target_nt == "-" and template_nt != "-"):
            """Deletion"""
            mutObject = eval_del(template, target, x, ref_index, STOP)
            return mutObject
        elif (template_nt == "-" and target_nt != "-"):
            """Insertion"""
            mutObject = eval_ins(template, target, x, ref_index, STOP)
            
        else:
            mutation_string = "SNP\t.\t.\t"+template_id+"\t"+str(ref_index)+"\t"+target_nt
            mutation_list.append(mutation_string)
    else:
        true_count += 1
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
    
    """Flag used to track whether or not to continue analysis. Based on 
    in/dels occuring on the trailing ends of the reads (unlikely edge 
    case)."""
    STOP = 0
    
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
        mutObject = eval_nt(template, target, x, ref_index, STOP)

        """Did a mutation occur?"""
        if not(mutObject):
            true_count += 1
            continue
        

        mut_type = mutObject.get_type()
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

        """Do we need to stop analysis (due to an in/del occuring at the 
        end of an alignment)? To do this, we can prematurely compute the 
        percent of nucleotides that were successfully matched and return 
        this value."""
        if STOP:
            return Decimal(true_count)/Decimal(stop_index-start_index)

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
                if not(evidence is None):
                    key = mobile_element_evidence['sangerID']
                    value = mobile_element_evidence
                    mob_evidence_dict[key]= value
                print "done."
    return mob_evidence_dict
