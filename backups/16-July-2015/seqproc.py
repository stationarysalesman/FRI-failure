# -*- coding: utf-8 -*-
"""
Created on Fri May 29 17:30:16 2015

@author: tyler
"""

from Bio import SeqIO
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
            
            
            
#Generate FASTA file containing sequences to be aligned


#Get list of mobile elements (only do this once)
mob_list = []
mob_path = "../mob_elements/"
for dirName, subdirList, fileList in os.walk(mob_path):
    for f in fileList:
        mob_list.append(SeqIO.read(dirName+f, "genbank"))


read_path = "../reads/"
plasmid_path = "../plasmids/"
template_path = "../templates/" # location of .fasta files that will be run through MAFFT


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