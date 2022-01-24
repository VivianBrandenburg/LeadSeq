# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import textwrap



# =============================================================================
# sequence handling
# =============================================================================


reverse_tab = str.maketrans("ACUG", "UGAC")
def reverse_complement(seq):
    return seq.translate(reverse_tab)[::-1]


def TU_converter(seq):
    if seq[:1000].find('U') == -1:
        return seq.replace('T', 'U')
    if seq[:1000].find('T') == -1:
        return seq.replace('U', 'T')

# =============================================================================
# general data reads
# =============================================================================

def read_gff(INF):
    gff_fields = ['seqname' ,'source', 'feature', 'start', 'end', 'score', 'strand', 'frame' , 'attribute']
    gff = pd.read_table(INF, comment='#', names = gff_fields)
    return gff


def read_fna(INF, convert_TU=False):
    with open(INF) as inf:
        fna = inf.read().strip()
        fna = [x.split('\n')  for x in fna.split('>') if x]
                
        fna = {x[0].split(' ')[0]:''.join(x[1:]) for x in fna}
        if convert_TU:
            return {key:TU_converter(value) for key, value in fna.items()}
    return fna

    
# =============================================================================
# convert to transcriptome
# =============================================================================

def get_attrtibutes(att):
    att = [x.split('=') for x in att.split(';')]
    return {x[0]:x[1] for x in att}

def genome2transcript(genome, gff, feature='CDS' ,gene_ids='locus_tag'):
    gff = gff[gff.feature==feature]
    res = {}
    for i in gff.index:
        gene = gff.loc[i]
        
        seqname = gene.seqname
        if seqname not in genome.keys(): raise KeyError('sequence names in gff are not present in fasta-file')
        seq = genome[seqname]
        
        transcript = seq[gene.start:gene.end+1]
        if not transcript:
            print('\nABORTING: gff and length of genome/scores do not match\n')
            return res
                        
        if type(transcript[0]) == int or float  and  gene.strand == '-':
                transcript = transcript[::-1]
            
        elif type(transcript[0]) == str and  gene.strand == '-':
            transcript = reverse_complement(transcript)
       
        attributes = get_attrtibutes(gene.attribute)
        name = attributes[gene_ids]
        
        res[name]=transcript
    return res


def write_transcriptome(inDic, outFile):
    with open (outFile, "w") as outf:
        for key, value in inDic.items():
            value_wrap = [value[i:i+80] for i in range(0, len(value), 80)]
            outf.write('\n'.join(['>'+key] + value_wrap )+ '\n')


# =============================================================================
# lead-scores read/write   
# =============================================================================

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
            name, UTR_start, UTR_stop= tab[0], int(tab[4]), int(tab[5])
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
   
