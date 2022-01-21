# -*- coding: utf-8 -*-
"""
script to calculate lead scores from tab files

script input:
    infile_control          *.tab file as from sam2tab.pl
    infile_treatment        *.tab file as from sam2tab.pl
    read_depth_control      as number or short cut 
    read_depth_treatment    as number or short cut 
    name_outfile_scores 
"""
print("input: infile_control, infile_treatment, read_depth_control, read_depth_treatment, name_outfile_scores\n")


import sys
from itertools import chain
from tools import scores as s 
from tools.box import write_results 

# read input variables
infile_1, infile_2 = sys.argv[1], sys.argv[2]
depth1, depth2 = sys.argv[3], sys.argv[4]
outfile_scores = sys.argv[5]

# read default read depths
dephs_short = {'rd1':9123314, 'rd2':8204588, 'rd3':15240566, 'rd4':9511133}
if depth1.isalpha: depth1 = dephs_short[depth1]
if depth2.isalpha: depth2 = dephs_short[depth2]
        
     
# read in raw counts
factor=100000000.0 #(1)
rc_1, rc_2 = s.read_tabfile(infile_1),  s.read_tabfile(infile_2) 
(infile_2)  #(3)

# normalize counts
normed1 = {key:s.norm_counts(rc_1[key], int(depth1), factor) for key in rc_1.keys() }
normed2 = {key:s.norm_counts(rc_2[key], int(depth2), factor) for key in rc_2.keys() } 
all_norm = list(chain(*normed1.values())) + list(chain(*normed2.values())) 
m_all = s.find_iteration(all_norm) #(5)

# calculate scores
scores={key:s.make_scores(normed1[key], normed2[key], m_all) for key in normed1.keys()}  

# write scores to output
write_results(scores, outfile_scores)
 
