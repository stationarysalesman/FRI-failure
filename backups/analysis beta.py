# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 13:45:07 2015

@author: tyler
"""

from Bio import SeqIO
from Bio import AlignIO
import os
import re

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
def homology_trunc(alignments):
    """Return an index at which to truncate nt analysis based on
    sequence homology to mobile elements."""
    for alignment in alignments:
        seq_id = alignment[0].id
        seq_length = len(alignment[0].seq)
        y = 0
        z = 0
        for x in range(len(alignment[0])):
            z = 0
            while (alignment[0][x] != "-"):
                z += 1
                x += 1
            if (z > y): y = z
    if (y > 100): return {seq_id: y}
    else: return {seq_id: seq_length}
def analyze(alignments, trunc_index_dict):
    """Analyze alignment file and output a genomediff file."""
    mutation_list = list()
    mutation_index = 1
    evidence_list = list()
    evidence_index = 1
    purines = ["a", "g"]
    pyrimidines = ["t", "c"]
    for alignment in alignments:
        alignment_length = len(alignment[0])
        plasmid_seq = alignment[0].seq              
        for y in range(1, len(alignment)):
            read_seq = alignment[y].seq
            # get start/stop indices for current sequence 
            stop_index = len(read_seq)-1
            while ((read_seq[stop_index] == "-") and (stop_index > 0)):
                stop_index -= 1
            if (stop_index == 0):
                print "ERROR: stop_index = 0 for read sequence", alignment[y].id
                return
            start_index = 0
            while (read_seq[start_index] == "-") and (start_index < len(read_seq)):
                start_index += 1
            if (start_index == len(read_seq)-1):
                print "ERROR: start_index = length-1 for read sequence", alignment[y].id
                return
            if (start_index > stop_index):
                print "ERROR: start_index > stop_index for read sequence", alignment[y].id
                return
            if (stop_index > alignment_length):
                print "ERROR: stop_index > alignment_length for read sequence", alignment[y].id
                return
            print "start_index:", start_index, "\nstop_index:", stop_index  
            skip = 0
            for x in range(start_index, stop_index+1):   
                if (skip):
                    skip -= 1
                    continue
                
                plasmid_nt = plasmid_seq[x]
                plasmid_nt_grp = plasmid_nt in purines
                read_nt = alignment[y][x]
                read_nt_grp = read_nt in purines
                if (read_nt != plasmid_nt):
                    if (read_nt == "n"):
                        # useless!
                        print "Encountered 'n' in read", alignment[y].id, "at nt", x
                        continue
                    elif (read_nt == "-" and plasmid_nt != "-"):
                        # deletion
                        print "Deletion at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                        z = x
                        while (read_seq[z] == "-"):
                            z += 1
                        del_length = z-x
                        print "del_length:", del_length
                        mutation_string = ("DEL\t" +str(mutation_index) + "\t" +
                        str(alignment[y].id) + "\t" +
                        str(x-start_index) + "\t" + str(del_length))
                        # determine number of repeats, if any
                        sub_rep = get_sub_rep(plasmid_seq[x:x+del_length])
                        repeat_seq = sub_rep[0]
                        deleted_units = sub_rep[1] # number of repeated units deleted
                        # determine if ancestor sequence contains repeats
                        # now we must find the bounds of the ancestor repeat sequence
                        repeat_len = len(repeat_seq)
                        start_index = x
                        read_start_index = start_index
                        stop_index = x+repeat_len
                        print "start_index:", start_index
                        print "stop_index:", stop_index
                        print "repeat_len:", repeat_len
                        skip = del_length-1
                            
                        # first, we move toward the 5' end by repeat_len each iteration
                        while (str(plasmid_seq[start_index:start_index+repeat_len])
                        == str(plasmid_seq[start_index-repeat_len:start_index]) and start_index>0):
                            print plasmid_seq[start_index:start_index+repeat_len], "equals",plasmid_seq[start_index-repeat_len:start_index]
                            start_index -= repeat_len
                        
                        print "new start_index:", start_index
                        
                        # now, we move toward the 3' end in the same manner
                        while (str(plasmid_seq[stop_index-repeat_len:stop_index])
                        == str(plasmid_seq[stop_index:stop_index+repeat_len]) and (stop_index+repeat_len)<len(plasmid_seq)):
                            print plasmid_seq[read_start_index:stop_index],"equals",plasmid_seq[stop_index:stop_index+repeat_len]
                            
                            stop_index += repeat_len
                            
                        print "new stop_index:", stop_index
                        print "entire repeated sequence:", plasmid_seq[start_index:stop_index]
                        repeat_ref_num = (stop_index-start_index)/repeat_len
                        print "repeat_ref_num:", repeat_ref_num
                        if (repeat_ref_num == 1 and deleted_units == 1):
                                # not repeat mediated
                            mutation_list.append(mutation_string)
                                
                            mutation_index +=1
                            break
                            
                        elif (repeat_len * repeat_ref_num >= 5):
                            # meets GenomeDiff criteria for RMD
                            repeat_new_copies = repeat_ref_num - deleted_units
                            errata = ("\trepeat_seq="+ repeat_seq + "\trepeat_len="+str(repeat_len)+"\t"+
                            "repeat_ref_num="+str(repeat_ref_num)+"\trepeat_new_copies="+
                            str(repeat_new_copies))
                            mutation_list.append(mutation_string+errata)
                            mutation_index += 1
                            
                            
                            
                            
                    elif (plasmid_nt == "-" and read_nt != "-"):
                        print "Insertion at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                        z = x
                        while (plasmid_seq[z] == "-"):
                            z += 1
                        ins_length = z-x
                        print "ins_length:", ins_length
                        mutation_string = ("INS\t" +str(mutation_index) + "\t" +
                        str(alignment[y].id) + "\t" +
                        str(x-start_index) + "\t" + str(ins_length))
                        # determine number of repeats, if any
                        sub_rep = get_sub_rep(read_seq[x:x+ins_length])
                        repeat_seq = sub_rep[0]
                        inserted_units = sub_rep[1] # number of repeated units inserted
                        # determine if read sequence contains repeats
                        # now we must find the bounds of the read repeat sequence
                        repeat_len = len(repeat_seq)
                        start_index = x
                        read_start_index = start_index
                        stop_index = x+repeat_len
                        print "start_index:", start_index
                        print "stop_index:", stop_index
                        print "repeat_len:", repeat_len
                        skip = ins_length-1
                            
                        # first, we move toward the 5' end by repeat_len each iteration
                        while (str(read_seq[start_index:start_index+repeat_len])
                        == str(read_seq[start_index-repeat_len:start_index]) and start_index>0):
                            print read_seq[start_index:start_index+repeat_len], "equals",read_seq[start_index-repeat_len:start_index]
                            start_index -= repeat_len
                        
                        print "new start_index:", start_index
                        
                        # now, we move toward the 3' end in the same manner
                        while (str(read_seq[stop_index-repeat_len:stop_index])
                        == str(read_seq[stop_index:stop_index+repeat_len]) and (stop_index+repeat_len)<len(plasmid_seq)):
                            print read_seq[read_start_index:stop_index],"equals",read_seq[stop_index:stop_index+repeat_len]
                            
                            stop_index += repeat_len
                            
                        print "new stop_index:", stop_index
                        print "entire repeated sequence:", read_seq[start_index:stop_index]
                        repeat_ref_num = (stop_index-start_index)/repeat_len
                        print "repeat_ref_num:", repeat_ref_num
                        if (repeat_ref_num == 1 and inserted_units == 1):
                                # not repeat mediated
                            mutation_list.append(mutation_string)
                                
                            mutation_index +=1
                            break
                            
                        elif (repeat_len * repeat_ref_num >= 5):
                            # meets GenomeDiff criteria for RMD
                            repeat_new_copies = repeat_ref_num + inserted_units
                            errata = ("\trepeat_seq="+ repeat_seq + "\trepeat_len="+str(repeat_len)+"\t"+
                            "repeat_ref_num="+str(repeat_ref_num)+"\trepeat_new_copies="+
                            str(repeat_new_copies))
                            mutation_list.append(mutation_string+errata)
                            mutation_index += 1
                            
                            
                            
                            
                    elif (read_nt_grp != plasmid_nt_grp):
                        print "Transversion at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                    
                    
    print plasmid_seq            
    for mut in mutation_list:
        print mut        
    return
    
path = "../sequences/"
trunc_index_dict = dict()
for dirName, subdirList, fileList in os.walk(path):
    for subdir in subdirList:
        for dirName2, subdirList2, filelist2 in os.walk(path + subdir):
            for subdir2 in subdirList2:
                if subdir2 == "alignments":
                    for dirName3, subdirList3, filelist3 in os.walk(path+subdir+"/"+subdir2):
                        for f in filelist3:
                            if "rm" in f:
                                # read/mobile element alignment
                                # need to determine if any homology to do nt analysis
                                rm_alignments = AlignIO.parse(path+subdir+"/"+subdir2+"/"+f, "fasta")
                                trunc_index_dict.update(homology_trunc(rm_alignments))
                                print trunc_index_dict
                            else:
                                continue
                            
for dirName, subdirList, fileList in os.walk(path):
   for subdir in subdirList:
        for dirName2, subdirList2, filelist2 in os.walk(path + subdir):
            for subdir2 in subdirList2:
                if subdir2 == "alignments":
                    for dirName3, subdirList3, filelist3 in os.walk(path+subdir+"/"+subdir2):
                        for f in filelist3:  
                            if "pr" in f:
                                alignment_seqs = AlignIO.parse(path+subdir+"/"+subdir2+"/"+f, "fasta")
                                analyze(alignment_seqs, trunc_index_dict)
                            else:
                                # temporary - remove in release
                                alignment_seqs = AlignIO.parse(path+subdir+"/"+subdir2+"/"+f, "fasta")
                                analyze(alignment_seqs, trunc_index_dict)
        