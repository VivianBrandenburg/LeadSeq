# -*- coding: utf-8 -*-

from tools.box import read_fna, read_gff, read_scores, genome2transcript
from tools.box import write_scores, write_transcriptome
import argparse


parser = argparse.ArgumentParser(usage='%(prog)s genome_file gff_file outfile [-h] [-f] [-n]')
parser.add_argument("genome_file", type=str, help="infile of genome or scores")
parser.add_argument("gff_file", type=str, help="gff infile")
parser.add_argument("outfile", type=str, help="outfile")
parser.add_argument('-f', '--feature' ,type=str, help="gff feature to use. Default is 'CDS'", default='CDS', metavar='' ,)
parser.add_argument("-n", '--name', type=str, help="identifier to use from gff. Default is 'locus_tag'", default='locus_tag', metavar='')
args = parser.parse_args()



# check which input we got
with open(args.genome_file, 'r') as f:
    mode = 'seq' if ''.join(f.readlines(2000)).find('>') != -1 else 'scores'
    
# read input
genome = read_fna(args.genome_file, True) if mode =='seq' else read_scores(args.genome_file)
gff = read_gff(args.gff_file)

# split into transcriptome
transcriptome = genome2transcript(genome, gff, args.feature, args.name)


# write outfile
if mode == 'seq': write_transcriptome(transcriptome, args.outfile)
else: write_scores(transcriptome, args.outfile)  

