# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 13:45:07 2015

@author: tyler
"""

from Bio import AlignIO
import os
import datetime

 # determine number of repeats, if any
def build_repeat_string(seq, x, mut_length):
    sub_rep = get_sub_rep(seq[x:x+mut_length])
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
        
    # first, we move toward the 5' end by repeat_len each iteration
    while (str(seq[start_index:start_index+repeat_len])
    == str(seq[start_index-repeat_len:start_index]) and start_index>0):
        #print seq[start_index:start_index+repeat_len], "equals",seq[start_index-repeat_len:start_index]
        start_index -= repeat_len
    
    #print "new start_index:", start_index
    
    # now, we move toward the 3' end in the same manner
    while (str(seq[stop_index-repeat_len:stop_index])
    == str(seq[stop_index:stop_index+repeat_len]) and (stop_index+repeat_len)<len(seq)):
        #print seq[stop_index-repeat_len],"equals",seq[stop_index:stop_index+repeat_len]
        
        stop_index += repeat_len
        
    #print "new stop_index:", stop_index
    #print "entire repeated sequence:", seq[start_index:stop_index]
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

    homolog_id = ""
    start_index = 0
    stop_index = 0
    read_id = alignments[0].id
    read_length = len(alignments[0].seq)
    y = 0 # will hold stop index
    z = 0 # reset every time a gap is encountered
    for r in range(1, len(alignments)):
        for x in range(len(alignments[r])):
            current_read = alignments[r].seq
            z = 0
            temp_start_index = x
            while (current_read[x] != "-"):
                # count how many nt match between Sanger read and any of the mob eles.
                z += 1
                x += 1
            temp_stop_index = x
            if (z > y): 
                y = z
                start_index = temp_start_index
                stop_index = temp_stop_index
     # Continue above loop until all mobile element alignments have been evaluated,
     # then select the best one           
    if (y > 100):
        # determine which mobile element id matches
        for a in alignments[a]:
            if (alignments[a][start_index:stop_index] == alignments[0][start_index:stop_index]):
               homolog_id = alignments[a].id 
               return ({seq_id: y}, homolog_id) 
    else: return ({seq_id: seq_length}, "")
     
     
def analyze(alignments, trunc_index_list):
    """Analyze alignment file and output a genomediff file."""
    mutation_list = list()
    alignment_length = len(alignments[0])
    plasmid_seq = alignments[0].seq     
    for y in range(1, len(alignments)):
        read_seq = alignments[y].seq
        read_seq_id = alignments[y].id
        # get start/stop indices for current sequence 
        stop_index = len(read_seq)-1
        q = 0
        for q in range(0, len(trunc_index_list)):
            for w in range(0, len(trunc_index_list[q])):
                if (trunc_index_list[q][1] != ''):      
                    # mobile element
                    stop_index = trunc_index_list[q][w][read_seq_id]
                    homolog_id = trunc_index_list[1]
                    mutation_string = "MOB\t.\t.\t"+read_seq_id+"\t"+homolog_id
                    mutation_list.append(mutation_string)
        while ((read_seq[stop_index] == "-") and (stop_index > 0)):
            stop_index -= 1
        if (stop_index == 0):
            print "ERROR: stop_index = 0 for read sequence", read_seq_id
            return
        start_index = 0
        while (read_seq[start_index] == "-") and (start_index < len(read_seq)):
            start_index += 1
        if (start_index == len(read_seq)-1):
            print "ERROR: start_index = length-1 for read sequence", read_seq_id
            return
        if (start_index > stop_index):
            print "ERROR: styearart_index > stop_index for read sequence", read_seq_id
            return
        if (stop_index > alignment_length):
            print "ERROR: stop_index > alignment_length for read sequence", read_seq_id
            return
        #print "start_index:", start_index, "\nstop_index:", stop_index  
        skip = 0
        
        for x in range(start_index, stop_index):   
            if (skip):
                skip -= 1
                continue
            
            plasmid_nt = plasmid_seq[x]
            read_nt = alignments[y][x]
            if (read_nt != plasmid_nt):
                if (read_nt == "n"):
                    # useless!
                   # print "Encountered 'n' in read", read_seq_id, "at nt", x
                    continue
                elif (read_nt == "-" and plasmid_nt != "-"):
                    # deletion
                    #print "Deletion at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                    z = x
                    while (read_seq[z] == "-"):
                        z += 1
                    del_length = z-x
                   # print "del_length:", del_length
                    mutation_string = ("DEL\t.\t.\t"+
                    str(alignments[y].id) + "\t" +
                    str(x-start_index) + "\t" + str(del_length))
                    errata = build_repeat_string(plasmid_seq, x, del_length)
                    if (errata[0] != ""):
                        mutation_list.append(mutation_string+errata[0])
                    else:
                        mutation_list.append(mutation_string)
                    skip = errata[1] # number of nt to skip   
                        
                elif (plasmid_nt == "-" and read_nt != "-"):
                   # print "Insertion at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                    z = x
                    while (plasmid_seq[z] == "-"):
                        z += 1
                    ins_length = z-x
                   # print "ins_length:", ins_length
                    mutation_string = ("INS\t.\t.\t"+
                    str(alignments[y].id) + "\t" +
                    str(x-start_index) + "\t" + str(ins_length))
                    errata = build_repeat_string(read_seq, x, ins_length)
                    if (errata[0] != ""):
                        mutation_list.append(mutation_string+errata[0])
                    else:
                        mutation_list.append(mutation_string)
                    skip = errata[1]
                      
                else:
                    #print "SNP at nt", x, ": plasmid:", plasmid_nt, ", read:", read_nt
                    
                    mutation_string = "SNP\t.\t.\t"+read_seq_id+"\t"+str(x-start_index)+"\t"+read_nt
                    mutation_list.append(mutation_string)
    return mutation_list

# Create directory that will contain genome diff files
def main():
    date_time_string = str(datetime.datetime.today())
    output_dir = "../genomediff_" + date_time_string 
    os.mkdir(output_dir)
    path = "../alignments/"
    if not(os.access(path, os.F_OK)):
        print "Error: cannot access ../alignments/"
        return
    trunc_index_list = list()
    for dirName, subdirList, fileList in os.walk(path):
        for f in fileList:
            if "rm" in f:
                # read/mobile element alignment
                # need to determine if any homology to do nt analysis
                rm_alignments = AlignIO.read(path+f, "fasta")
                trunc_index_list.append(homology_trunc(rm_alignments))
                #print trunc_index_list
            else:
                continue
       
    for dirName, subdirList, fileList in os.walk(path):
        for f in fileList: 
            if "pr" in f:
                alignment_seqs = AlignIO.read(path+f, "fasta")
                file_input = analyze(alignment_seqs, trunc_index_list)
                # Create GenomeDiff file   
                name = alignment_seqs[0].id + ".gd"
                with open(output_dir+"/"+name, "w") as out_file:
                    out_file.write("#=GENOME_DIFF 1.0\n")
                    for mutation in file_input:
                        out_file.write(mutation + "\n")
    return  
                          
main()                            

        