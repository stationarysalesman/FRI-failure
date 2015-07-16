# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 11:15:19 2015

@author: tyler
"""

from multiprocessing import Pool
import time

def f(x):
    return (x * 73) / 7 + (x*117)

if __name__ == '__main__':
    start = time.time()
    for y in range(100):
        for x in range(9999999):
            f(x)
        stop= time.time()
    total = stop-start
    print "Standard:", total
    
    start = time.time()
    p = Pool()
    for y in range(100):    
        p.map(f, range(9999999))
        stop = time.time()
    total = stop-start
    print "Multithreaded:", total