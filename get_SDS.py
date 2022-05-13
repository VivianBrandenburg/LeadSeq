# -*- coding: utf-8 -*-

import sys
import pandas as pd
from tools.box import read_scores, read_UTR_file, get_SDS
from tools.plots import plot_SDS


# INF_delta, INF_mwu = sys.argv[1], sys.argv[2]
# INF_5UTR =sys.argv[3]
# OUTF = sys.argv[4]

INF_delta, INF_mwu = 'results/delta.txt', 'results/mwu.txt'
INF_5UTR = 'data/5UTR_annotation.txt'
OUTF = 'results/SDS.csv'


# read in data
delta = read_scores(INF_delta)
mwu = read_scores(INF_mwu)
UTRs = read_UTR_file(INF_5UTR)


# find delta 
SDS_delta, SDS_mwu = get_SDS(UTRs, delta, mwu)

# sort results
res = pd.DataFrame({'tag':SDS_delta.keys(), 'delta': SDS_delta.values(), 'mwu':SDS_mwu.values()})
res = res.sort_values(by=['mwu', 'delta'])

# write results 
res.to_csv(OUTF, index=False)

# plot results
res['label'] = (res.delta<0) & (res.mwu<-0.5)
res['label'] =['candidate' if x==True else 'RBS' for x in res['label'] ]
plot_SDS(res, 'plots/SDS.png')


