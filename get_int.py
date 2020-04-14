#!/root/anaconda3/bin/python
from __future__ import print_function
import sys, os
import argparse
from itertools import islice
from subprocess import call
import numpy as np
import math
# access timings
import time
# import external file with functions
from func import get_dist, get_jval
from func import get_q1_val, get_q2_val, get_q3_val, get_q4_val

#conversion factors
ang2bohr = 1.8897161646320724
bohr2ang = 0.52918

# read input from command line
parser = argparse.ArgumentParser()
parser.add_argument('--atoms', type=str)
parser.add_argument('--fix', type=str)
args = parser.parse_args() 
var = args.atoms.split(',')
idx = [int(x) for x in var]
atom1 = idx[0]
atom2 = idx[1]
# check if a fixatom has been defined
if args.fix != "F":
    fixatom = args.fix
    idxfix = int(args.fix)
    print("fixatom", fixatom)
else:
    print("no fixatom set")

dist = get_dist(idx)
# step = float(0.1)
# step = float(0.05)
step = float(0.2)
start = time.time()
q1 = get_q1_val(step,idx,idxfix,dist) 
q2 = get_q2_val(step,idx,idxfix,dist) 
q3 = get_q3_val(step,idx,idxfix,dist) 
q4 = get_q4_val(step,idx,idxfix,dist) 

end = time.time()
walltime = end -start  # in seconds
j = float(q1 + q2 + q3 + q4) 

print(" ")
print("valq1 final =", "%8.5f" %(q1) ) 
print("valq2 final =", "%8.5f" %(q2) ) 
print("valq3 final =", "%8.5f" %(q3) ) 
print("valq4 final =", "%8.5f" %(q4) ) 
print("SUM =        ", "%8.5f" %(j) ) 

print(" ")
print("wall time in seconds", walltime)
print("wall time in minutes", float(walltime/60.0))

# store result in output file
fout="result.txt"
f1 = open(fout,"w") 
f1.write("bond atom 1 = "+ str(atom1) + "\n") 
f1.write("bond atom 2 ="+ str(atom2) + "\n") 
f1.write("fixpoint = "+ fixatom + "\n")
# f1.write("Q1 = "+ str(q1) + "\n")
f1.write("Q1 = "+ f"{q1:.5f}" + "\n")
f1.write("Q2 = "+ f"{q2:.5f}" + "\n")
f1.write("Q3 = "+ f"{q3:.5f}" + "\n")
f1.write("Q4 = "+ f"{q4:.5f}" + "\n")
f1.write("J  = "+ f"{j:.5f}" + "\n")
f1.close()




