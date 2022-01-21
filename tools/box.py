# -*- coding: utf-8 -*-

import numpy as np



def read_scores(infile):
    with open(infile, 'r') as inf:
        data = {}
        for line in inf:
            ID, scores = line.strip().split('\t')
            data[ID] = [ float(x) for x in scores.split(';')]
    return data



#write output in file
def write_scores(inDic, outFile):
    with open (outFile, "w") as outf:
        for key, value in inDic.items():
            valueS = ';'.join([str(x) for x in value])
            outf.write('%s\t%s\n' % (key,valueS))



# =============================================================================
# evaluate on tRNAs 
# =============================================================================

def eval_trnas(ID, y_true, y_pred, cut):
    y_true = [int(x) for x in list(y_true)]
    
    TP = sum([true == 0 and pred > cut for true, pred in zip(y_true, y_pred)])
    FP = sum([true == 1 and pred > cut for true, pred in zip(y_true, y_pred)])
    TN = sum([true == 1 and pred <= cut for true, pred in zip(y_true, y_pred)])
    FN = sum([true == 0 and pred <= cut for true, pred in zip(y_true, y_pred)])
    
    selectivity = round(TN/(TN+FP),3)
    sensitivity = round(TP/(TP+FN),3)
    return{ 'ID': ID,  'cutoff': cut,  'selectivity':selectivity, 'sensitivity':sensitivity}



# =============================================================================
# find shine-dalgarno-sequence
# =============================================================================

def read_UTR_file(UTR_file):
    out= []
    with open(UTR_file, 'r') as inf:
        for line in inf:
            tab = line.replace('\n', '').split('\t')
            name, UTR_start, UTR_stop= tab[0].replace('YPK_TSS','TSS'), int(tab[4]), int(tab[5])
            dic = {'name': name,  'UTR_start':UTR_start, 'UTR_stop':UTR_stop, 'strand':tab[1]}
            length_UTR = UTR_stop - UTR_start +1
            if length_UTR > 14:    
                dic.update({'SDS_start':length_UTR - 14, 'SDS_stop':length_UTR - 6, 'UTR_start':1, 'UTR_stop':length_UTR})
                out.append(dic)
            else:
                continue
    return out
 


def get_SDS(UTR_listIn, deltaIn, mwuIn):
    SDS_delta, SDS_mwu = {}, {}
    for i in UTR_listIn:
        if i['name'] in deltaIn.keys():
            name = i['name']
            SDS_delta[name]  = np.mean(deltaIn[name][(i['SDS_start']-1):i['SDS_stop']])
            SDS_mwu[name]  = np.mean(mwuIn[name][(i['SDS_start']-1):i['SDS_stop']])
    return     SDS_delta, SDS_mwu
   
