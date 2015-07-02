# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:21:00 2015

@author: tyler
"""
from Bio import SeqIO
from Bio import AlignIO
import os
import re

def trim_n(seq):
    """Trims the 5' and 3' ends of Sanger reads."""
    # Trim 5' end
    i = 0
    skip = False
    done = False
    while (i < len(seq.seq)) and not(done):
        if (seq.seq[i] != "N"):
            for j in range(10):
                if (seq.seq[i+j+1] == "N"):
                    skip = True
                    i += 1
                    break
            if (skip):
                i += 1
                skip = False
                continue
            elif not(skip):
                seq = seq[i:]
                done = True
        else:
            i += 1
                
                
    # Trim 3' end
    j = i
    done = False
    while (j < len(seq.seq)) and not(done):
        if (((seq.letter_annotations.values()[0][j]+
              seq.letter_annotations.values()[0][j+1]+
              seq.letter_annotations.values()[0][j+2])/3) 
              < 10):
                 seq = seq[i:j]
                 done = True
                 break
        j += 1
    return seq

def seqproc(read_path, plasmid_path, template_path, mob_path):
    """Create templates that will be later aligned with MAFFT"""
    
    # Create list containing mobile elements
    mob_list = list()
    for dirName, subdirList, fileList in os.walk(mob_path):
        for f in fileList:
            mob_list.append(SeqIO.read(dirName+f, "genbank"))
            
    for dirName, subdirList, fileList in os.walk(plasmid_path):
        for plasmid_file in fileList:
            read_list = list()
            initials = ""
            plasmid_read_list = [] # will contain plasmid sequence and all Sanger reads
            plasmid = SeqIO.read(plasmid_path+plasmid_file, "genbank") # object containing plasmid info     
            temp = plasmid_file   
            dest = []
            dest = re.split("_", temp) #get initials from plasmid file name         
            initials = dest[0].lower()
            for dirName2, subdirList2, fileList2 in os.walk(read_path):
                for read_file in fileList2:
                    if (initials in read_file.lower()):
                        out_file = SeqIO.read(read_path+read_file, "abi")
                        trimmed_file = trim_n(out_file) # trim n's
                        read_list.append(trimmed_file)
                   
    # At this point we have a list of all Sanger reads corresponding to our current plasmid file.
    # Ends have been trimmed of N's. 
            plasmid_read_list.append(plasmid)
            plasmid_read_list += read_list
       
            for i, sequence in enumerate(read_list):
                newList = list()
                newList.append(sequence)
                newList += mob_list
                SeqIO.write(newList, template_path+"/rm"+str(i)+"_"+initials, "fasta")
                del(newList)                
            SeqIO.write(plasmid_read_list, template_path+"/pr_"+initials, "fasta")
            del(read_list)
            del(plasmid_read_list)
    return

    
def get_seq_list(seq_path):
    """Returns a list of Seq objects from GenBank files located
    in seq_path."""
    seq_list = []
    for dirName, subdirList, fileList in os.walk(seq_path):
        for f in fileList:
            seq_list.append(SeqIO.read(dirName+f, "genbank"))
    return seq_list
def get_indices(seq):
    """Obtain start and stop indices for analyzing a single alignment
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
    """Returns the length of the smallest repeated sequence in the
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
def build_repeat_string(seq, x, mut_length):
    seq_str = str.replace("-", "", str(seq))
            
    sub_rep = get_sub_rep(seq_str[x:x+mut_length])
    repeat_seq = sub_rep[0]
    delta_units = sub_rep[1] # number of repeated units deleted/inserted
    # determine if read sequence contains repeats
    # now we must find the bounds of the read repeat sequence
    repeat_len = len(repeat_seq)
    
    start_index = x
    stop_index = x+repeat_len
   # print "start_index:", start_index
   # print "stop_index:", stop_index
   # print "repeat_len:", repeat_len
    skip = mut_length-1
    if (repeat_len <= 0):
        print "Error: repeat_len <= 0"
        return ["", skip]
    # first, we move toward the 5' end by repeat_len each iteration
    while (str(seq_str[start_index:start_index+repeat_len])
    == str(seq_str[start_index-repeat_len:start_index]) and start_index>0):
        #print seq_str[start_index:start_index+repeat_len], "equals",seq_str[start_index-repeat_len:start_index]
        start_index -= repeat_len
    #print "new start_index:", start_index
    
    # now, we move toward the 3' end in the same manner
    while (str(seq_str[stop_index-repeat_len:stop_index])
    == str(seq_str[stop_index:stop_index+repeat_len]) and (stop_index+repeat_len)<len(seq_str)):
        #print seq_str[stop_index-repeat_len],"equals",seq_str[stop_index:stop_index+repeat_len]
        
        stop_index += repeat_len
        
    #print "new stop_index:", stop_index
    #print "entire repeated sequence:", seq_str[start_index:stop_index]
    repeat_ref_num = (stop_index-start_index)/repeat_len
    #print "repeat_ref_num:", repeat_ref_num
    if (repeat_ref_num == 1 and delta_units == 1):
            # not repeat mediated
        return ["", skip]       
    elif (repeat_len * repeat_ref_num >= 5):
        # meets GenomeDiff criteria for RMD/RMI
        repeat_new_copies = repeat_ref_num + delta_units
        errata = ("\trepeat_seq="+ repeat_seq + "\trepeat_len="+str(repeat_len)+"\t"+
        "repeat_ref_num="+str(repeat_ref_num)+"\trepeat_new_copies="+
        str(repeat_new_copies))
        return [errata, skip]
    else:
        return ["", skip]
        
def analyze_seq(template, target, start_index, stop_index):
    """Analyze an alignment between template and target
    between the given indices."""
    mutation_dict = dict()
    template_seq = template.seq
    target_seq = target.seq
    target_id = target.id    

    #print "Analyzing", target_id,"against element", template.id
    #print "start:", start_index
    #print "stop:", stop_index
    
    skip = 0   
    read_start = 0
    #while (target_seq[read_start] == "-"):
     #   read_start += 1
        
    for x in range(start_index, stop_index):   
        if (skip):
            skip -= 1
            continue
        

        template_nt = template.seq[x]
        target_nt = target.seq[x]
        if (target_nt != template_nt):
            if (target_nt == "n"):
                # useless!
               # print "Encountered 'n' in read", target_seq_id, "at nt", x
                continue
            elif (target_nt == "-" and template_nt != "-"):
                # deletion
                #print "Deletion at nt", x, ": plasmid:", template_nt, ", read:", target_nt
                z = x
                ignore_count = 0
                while (z < len(target_seq) and target_seq[z] == "-"):
                    if (template_seq[z] == "-"): # interrupted by virtue of being a multiple seq alignment
                        ignore_count += 1                      
                    z += 1
                del_length = z - x - ignore_count
               # print "del_length:", del_length
                mutation_string = ("DEL\t.\t.\t"+
                str(target.id) + "\t" +
                str(x-read_start) + "\t" + str(del_length))
                errata = build_repeat_string(template_seq, x, del_length)
                if (errata[0] != ""):
                    key = hash(mutation_string+errata[0])
                    if not(key in mutation_dict):
                        mutation_dict[key] = mutation_string+errata[0]
                else:
                    key = hash(mutation_string)
                    if not(key in mutation_dict):
                        mutation_dict[key] = mutation_string                    
                skip = errata[1] + ignore_count # number of nt to skip  
                if ((skip + x) >= stop_index):
                   return mutation_dict
                    
            elif (template_nt == "-" and target_nt != "-"):
                #print "Insertion in", target.id, "at nt", x, ": plasmid:", template_nt, ", read:", target_nt
                z = x
                ignore_count = 0
                while (z < len(template_seq) and template_seq[z] == "-"):
                    if (target_seq[z] == "-"):
                        ignore_count += 1
                    z += 1
                ins_length = z - x - ignore_count
               # print "ins_length:", ins_length
                mutation_string = ("INS\t.\t.\t"+
                str(target.id) + "\t" +
                str(x-read_start) + "\t" + str(ins_length))
                errata = build_repeat_string(target_seq, x, ins_length)
                if (errata[0] != ""):
                    key = hash(mutation_string+errata[0])
                    if not(key in mutation_dict):
                        mutation_dict[key] = mutation_string+errata[0]
                else:
                    key = hash(mutation_string)
                    if not(key in mutation_dict):
                        mutation_dict[key] = mutation_string
                skip = errata[1] + ignore_count
                if ((skip + x) >= stop_index):
                   return mutation_dict
            else:
                #print "SNP at nt", x, ": plasmid:", template_nt, ", read:", target_nt
                
                mutation_string = "SNP\t.\t.\t"+target_id+"\t"+str(x-read_start)+"\t"+target_nt
                key = hash(mutation_string)
                if not(key in mutation_dict):
                    mutation_dict[key] = mutation_string
   
    return mutation_dict

def mob_info(alignments):
    """Returns a list of info about an alignment analysis based on
    sequence homology to mobile elements."""

    read = alignments[0] # This is the Sanger read
    mob_mut_list = list()
    for r in range(1, len(alignments)):
        mob = alignments[r] # this is the current mobile element we want to analyze    
        # Analyze the alignment only in the region of the mobile element
        mob_indices = get_indices(mob)
        #print "mob_indices:", mob_indices
        read_indices = get_indices(read)
        #print "read_indices:", read_indices
        # We want indices that capture only portions of the alignment where both sequences
        # are aligned. This may overlook edge insertions/deletions.
        start = max(mob_indices[0], read_indices[0])
        stop = min(mob_indices[1], read_indices[1])
        if ((stop - start) < 200): # not long enough for any meaningful analysis
            return
        #print "stop decided:", stop
        # call analyze_seq on our read/mobile element to get info about sequence homology
        temp = analyze_seq(mob, read, start, stop)
        mob_seq_rev = mob[::-1]
    
        
        temp_rev = analyze_seq(mob_seq_rev, read, start, stop)
        if ("is5" in mob.id.lower()):
            print len(temp_rev.keys())
            for rec in temp_rev:
                print temp_rev[rec]
        mob_mut_list.append([read.id, mob.id, temp, start])
        mob_mut_list.append([read.id, mob.id, temp_rev, start])

    curr_lst = list()
    curr_len = 9999999 # do not use this value for "big data" alignments/analyses
    if not(mob_mut_list):
        return
    for lst in mob_mut_list: # find the mobile element with the least number of mutations      
        if (len(lst[2].keys()) < curr_len):
            #print "new list"
            curr_len = len(lst[2].keys())
            curr_lst = lst
    if (curr_lst and curr_len< 50): # is the mobile element actually homologous?
        del(mob_mut_list)
        print curr_lst
        return curr_lst
    
   
    return                                     
    
def analysis_control(alignment_path, curr_output_path, log_path):

    mob_evidence_dict = dict()
    if not(os.access(alignment_path, os.F_OK)):
        err_str =  "Error: cannot access", alignment_path
        with open(log_path, "a") as log:
            log.write(err_str + "\n")
        return 1
    
    for dirName, subdirList, fileList in os.walk(alignment_path):
        for f in fileList:
            if "rm" in f and "md" in f:
                # read/mobile element alignment
                # need to determine if any homology to do nt analysis
                rm_alignments = AlignIO.read(alignment_path+f, "fasta")
                print "Gathering mobile element homology data for read", rm_alignments[0].id, "...",
                evidence = mob_info(rm_alignments)
                if not(evidence is None):
                    key = evidence[0]
                    value = [evidence[1], evidence[2], evidence[3]]
                    mob_evidence_dict[key]= value
                print "done."

    for dirName, subdirList, fileList in os.walk(alignment_path):
        for f in fileList:
            if "pr" in f and "md" in f:  # plasmid/read alignment 
                mutation_dict = dict()
                pr_alignments = AlignIO.read(alignment_path+f, "fasta")
                print "Analyzing alignments for plasmid:", pr_alignments[0].id, "...", 
                plasmid_indices = get_indices(pr_alignments[0])
                for x in range(1, len(pr_alignments)): # only get Sange reads
                    read_indices = get_indices(pr_alignments[x])
                    start = max(plasmid_indices[0], read_indices[0])
                    stop = min(plasmid_indices[1], read_indices[1])
                    #print "initially, stop was", stop
                    key = pr_alignments[x].id
                    if (key in mob_evidence_dict):
                        stop = mob_evidence_dict[key][2]
                        #print "...but then, stop was", stop, "for key", key
                        mutation_string = "MOB\t.\t.\t"+key+"\t"+mob_evidence_dict[key][0]
                        key2 = hash(mutation_string)
                        if not(key2 in mutation_dict):
                            mutation_dict[key2] = mutation_string
                        if (mob_evidence_dict[key][1]): # append any mutations within the mobile element
                            for key3 in mob_evidence_dict[key][1].keys():
                                key4 = hash(key3)
                                if not(key4 in mutation_dict):
                                    mutation_dict[key4] = mob_evidence_dict[key][1][key3]
                    
                    #print "analyze:", pr_alignments[0].id, ",", pr_alignments[x].id, ",", start, ",", stop
                    mutation_dict_2 = analyze_seq(pr_alignments[0], pr_alignments[x], start, stop)
                    if (mutation_dict_2):
                        mutation_dict.update(mutation_dict_2)
                name = pr_alignments[0].id + ".gd"
                with open(curr_output_path+"/"+name, "w") as out_file:
                    out_file.write("#=GENOME_DIFF 1.0\n")
                    for v in mutation_dict.values():
                        out_file.write(v + "\n")

                del(mutation_dict)
                print "done."
    return
