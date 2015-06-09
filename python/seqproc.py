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
    i = 0
    done = False
    while (i < len(seq.seq)) and not(done):
        if (((seq.letter_annotations.values()[0][i]+
              seq.letter_annotations.values()[0][i+1]+
              seq.letter_annotations.values()[0][i+2])/3) 
              < 10):
                 seq = seq[:i]
                 done = True
        i += 1
    return seq
            
            
            
#Generate FASTA file containing sequences to be aligned


#Get list of mobile elements (only do this once)
mob_list = []
mob_path = "../mob_elements/"
for dirName, subdirList, fileList in os.walk(mob_path):
    for f in fileList:
         mob_list.append(next(SeqIO.parse(dirName+"/"+f, "genbank")))

# Get list of plasmid/sanger reads
# We will do this for every dataset
read_list = [] # will contain Sanger reads
plasmid = [] # will contain plasmid
path = "../sequences/"
for dirName, subdirList, fileList in os.walk(path):
    if (not subdirList):
        template_dir = dirName+"/templates"
        os.mkdir(template_dir) # this is where .fasta files to be aligned using MAFFT are saved
        a = re.split("/", dirName)
        name= a[-1] + ".fasta" #we will name all analysis files based on dir name
        for f in fileList:
            if ".gb" in f: # this should be plasmid/template
                plasmid.append(next(SeqIO.parse(dirName+"/"+f, "genbank"))) #get GenBank file from sequences directory
            elif ("VF" in f) or ("VR" in f): # these should be Sanger reads
                read_file = next(SeqIO.parse(dirName+"/"+f, "abi"))
                trimmed_file = trim_n(read_file)
                read_list.append(trimmed_file) #add trimmed Sanger read files
        plasmid_read_list = plasmid + read_list
        for i, sequence in enumerate(read_list):
            newList = list()
            newList.append(sequence)
            newList += mob_list
            SeqIO.write(newList, template_dir+"/rm"+str(i)+"_"+name, "fasta")
            del(newList)
        SeqIO.write(plasmid_read_list, template_dir+"/pr_"+name, "fasta")