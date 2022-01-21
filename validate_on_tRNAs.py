import pandas as pd
from tools.box import read_scores, eval_trnas
from tools.plots import plot_validation, plot_score_distribution
import sys

# =============================================================================
# prepare data
# =============================================================================


INF_l25, INF_l37, INF_tRNAs = sys.argv[1], sys.argv[2], sys.argv[3]


# read in tRNAs and scores
tRNA = pd.read_table(INF_tRNAs, names=['AS', 'ID', 'tRNA_num', 'start',  'stop', 'seq', 'struc'], index_col=False)
l25 = read_scores(INF_l25)            
l37 = read_scores(INF_l37)            

# combine scores and tRNAs
tRNA['l25'] = [l25[x] for x in tRNA.ID]
tRNA['l37'] = [l37[x] for x in tRNA.ID]


# =============================================================================
# selectivity and sensitivity
# =============================================================================


# calculate selectivity/sensitivity
valid_all = pd.DataFrame()   
cut = 0.5
for i in tRNA.index:
    res_l25 = eval_trnas(tRNA.loc[i].ID, tRNA.loc[i].struc, tRNA.loc[i].l25, cut)
    res_l37 = eval_trnas(tRNA.loc[i].ID, tRNA.loc[i].struc, tRNA.loc[i].l37, cut)
    valid_all = valid_all.append(dict(res_l25, **{'sample':'25 °C'}), ignore_index = True)
    valid_all = valid_all.append(dict(res_l37, **{'sample':'37 °C'}), ignore_index = True)

# plot 
valid_all = pd.melt(valid_all, id_vars = ['ID', 'cutoff', 'sample'])
plot_validation(valid_all[valid_all.cutoff==0.5], 'plots/performance.png')


# =============================================================================
# distribution of scores at paired/unpaired nts
# =============================================================================

# extract scores and pairing state 
scores_l25 = pd.DataFrame({'scores':[i for s in list(tRNA.l25) for i in s],
                           'paired':['paired' if z==1 else 'unpaired' for z in [i for s in [[int(x) for x in list(y)] for y in tRNA.struc] for i in s]],
                           'sample':'25 °C',})
scores_l37 = pd.DataFrame({'scores':[i for s in list(tRNA.l37) for i in s],
                           'paired':['paired' if z==1 else 'unpaired' for z in [i for s in [[int(x) for x in list(y)] for y in tRNA.struc] for i in s]],
                           'sample':'37 °C',})

# plot distribution 
plot_score_distribution(scores_l25, 'l25', 'plots/l25_distribution.png')
plot_score_distribution(scores_l37, 'l37', 'plots/l37_distribution.png')



