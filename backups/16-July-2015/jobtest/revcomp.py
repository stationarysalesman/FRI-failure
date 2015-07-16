# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:44:45 2015

@author: tyler
"""

from Bio.SeqRecord import SeqRecord

def revcomp(seq):
    if isinstance(seq, SeqRecord):
        return seq.reverse_complement()
    else:
        print "shit"
        return