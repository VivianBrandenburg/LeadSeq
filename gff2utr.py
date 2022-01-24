# -*- coding: utf-8 -*-

from tools.box import read_gff, get_attrtibutes
import argparse
import pandas as pd

parser = argparse.ArgumentParser(usage='%(prog)s genome_file gff_file outfile [-h] [-f] [-n]')
parser.add_argument("gff_file", type=str, help="gff infile")
parser.add_argument("outfile", type=str, help="outfile")
parser.add_argument('-f', '--feature' ,type=str, help="gff feature to use. Default is '5UTR'", default='5UTR', metavar='' ,)
parser.add_argument("-n", '--name', type=str, help="name to use from gff attributes. Default is 'locus_tag'", default='locus_tag', metavar='')
args = parser.parse_args()



# read input
gff = read_gff(args.gff_file)

# make utr
gff = gff[gff.feature == args.feature]
attributes = [get_attrtibutes(att) for att in gff.attribute]
ID = [att[args.name] for att in attributes]

# add parent gene if possible
genes = []
for att in attributes:
    try: genes.append(att['gene'])
    except KeyError:  genes.append('.')

lengths = [abs(x-y) for x,y in zip(gff.start, gff.end)]


utr = pd.DataFrame({ 
    'ID':ID, 'strand': gff.strand, 'parent': genes,
    'genes_involved':'.', 'start': gff.start, 'end':gff.end,
    'length':lengths
    })

utr.to_csv(args.outfile, sep='\t', index=False, header=False)
