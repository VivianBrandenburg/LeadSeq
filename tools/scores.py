import numpy as np 
  
#(3): read in raw counts from infile1 and infile2        
def read_tabfile(infile):
    out_dic={}
    with open(infile, "r") as inf:
        for line in inf:
            tabs = line.replace("\n","").split("\t")
            out_nums= tabs[1].split(";")
            out_nums=[float(x) for x in out_nums]
            out_dic[tabs[0]]=out_nums
    return out_dic
                
#(4): normalize counts
def norm_counts(in_rc, tiefe, factor):
    return [round(((x+5)*factor)/(tiefe),5) for x in in_rc] 
     

#(5): find m, which is needed for calculating the scores 
def find_iteration (in_list): 
    for _ in range(3):
        m = np.median(in_list)      
        in_list=[x for x in in_list if x>m] 
    return m                       


#(6): calculate lead-scores     
def make_scores(in1, in2, m_here):
    out = list()
    for i,j in zip(in1,in2):   
        score = round(np.log2((j+m_here)/(i+m_here)),5)    # score = log2((NormCountTreat + m) / (NormCountControl + m))   
        out.append(score)
    out = [x if x >0  else 0 for x in out]   # cut scores at minimum of zero (negative scores are meaningless because of experimental procedure)
    if max(out) == min(out):
        out = [0 for x in out]
    return out                             


#write output in file
def write_scores(inDic, outFile):
    with open (outFile, "w") as outf:
        for key, value in inDic.items():
            valueS = ';'.join([str(x) for x in value])
            outf.write('%s\t%s\n' % (key,valueS))
