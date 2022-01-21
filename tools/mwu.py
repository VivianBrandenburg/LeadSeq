# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import mannwhitneyu


# filter out transcripts which have not more then 10 informative lead-scores
def filter_scores(scores):
    scores = [x for x in scores if x>0.5]
    return len(scores)>10
    
 
    
# sliding window, subroutine of compare_with_mwu
from collections import deque
def window(seq, n=3):
    it = iter(seq)
    win = deque((next(it, None) for _ in range(1)), maxlen=n)
    yield win
    append = win.append
    for e in it:
        append(e)
        yield win
    for _ in range(len(win)-1):
        win.popleft()
        yield win
    

# compare  windows against rest with mwu
def compare_with_mwu(deltas, winSize, ends):
    filt_up=np.percentile([x for x in deltas if x > 0 ], 30) if max(deltas)>0 else  0 #find top quantile
    filt_down=np.percentile([x for x in deltas if x < 0 ], 70) if min(deltas)<0 else  0 #find bottom quantile
    
    deltas_no_center = [x for x in deltas if x > filt_up or x < filt_down] # filter out center deltas for whole transcript
    deltas_no_center_mean = np.mean(deltas_no_center)
    
    corr_s, pval_s = [], []
    for win_delta in window(deltas, winSize): #slide window acros transcript
        if len(win_delta)>=21 and min(win_delta) != max(win_delta):
                win_no_center = [x for x in win_delta if x > filt_up or x < filt_down] #filter out center deltas in window
                if win_no_center:
                    corr, pval = mannwhitneyu(win_no_center, deltas_no_center) #compare filtered window with filtered transcript
                else:
                    corr, pval = 0, 1.0
                if pval > 0.05: pval = 1.0 # p-values above 0.05 will be zero after log10
                if pval != 0: pval = np.log10(pval) * -1 # high p-vals are high confidence level
                if np.mean(win_no_center) < deltas_no_center_mean: pval = pval *-1 # negative pval -> l37 scores are higher -> structure melts with higher temp   
                 
                corr_s.append(corr)
                pval_s.append(pval)
        else:
            corr_s.append(0)
            pval_s.append(0)
    return corr_s[ends:-ends], pval_s[ends:-ends] 
    
