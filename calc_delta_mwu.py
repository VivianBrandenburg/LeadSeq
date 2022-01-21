# -*- coding: utf-8 -*-

import sys
from tools.box import read_scores, write_scores
from tools.mwu import filter_scores, compare_with_mwu

# data
INF_l25, INF_l37 =  sys.argv[1], sys.argv[2]
OUTF_delta, OUTF_mwu = sys.argv[3], sys.argv[4]

# read in data
l25 = read_scores(INF_l25)
l37 = read_scores(INF_l37)


# delta scores
filtered = [key for key in l25.keys() if filter_scores(l25[key]) and filter_scores(l37[key])]
lead_delta = {key:[x-y for x,y in zip(l25[key], l37[key])] for key in filtered}


# size of window which is to compared to the whole transcript.
winSize = 23 # has to be odd number > 21
ends = int(winSize/2) 


# mwu
counter = 0
corr_dic, pval_dic = {}, {}
for key, value in lead_delta.items():
    counter += 1
    print(counter, key)
    corr_dic[key], pval_dic[key] =  compare_with_mwu(value, winSize, ends)

# write results
write_scores(pval_dic, OUTF_delta)
write_scores(lead_delta, OUTF_mwu)
  



