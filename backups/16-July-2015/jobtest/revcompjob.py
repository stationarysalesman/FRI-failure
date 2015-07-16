# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler
"""


from Bio import SeqIO    
from jobmanager import jobmanager
import os

def controller(job):
    
    input_dir = job.input_dirs["sequences"]
    output_dir = job.master_output_dir
    revcomp = job.process_module
    outfile_name = "outfile.txt"
    seq_list = list()
    output_list = list()
    for dirName, subdirList, fileList in os.walk(input_dir):
        for f in fileList:
            seq_list.append(SeqIO.read(input_dir+f, "fasta"))
    for seq in seq_list:
        output_list.append(revcomp.revcomp(seq))
    with open(output_dir+outfile_name, "w") as out:
        SeqIO.write(output_list, out, "fasta")
    return
    
def main():
    input_dirs = dict()    
    input_dirs["sequences"] = "/home/tyler/Documents/research/sequences/"
    input_dirs["plasmids"] = "/home/tyler/Documents/research/plasmids/"
    myJob = jobmanager(input_dirs, [], "logfile.log", "revcomp")
    controller(myJob)
    
main()