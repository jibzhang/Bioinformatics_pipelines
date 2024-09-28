#!/usr/bin/env python

'''
For a directory containing myCpG/methCall files, this script will find the average percent methylation for DMRs (inputed as a seperate bedfile), by doing the following:
1) Preparing an empty table indexed by the DMRs (which are inputed as a bed-3 file without a header)
2) Bedtools intersect of myCpG/methCall and DMR file to obtain a table of DMC that overlap DMR
3) Find percent methylation for each CpG that overlaps a DMR
3) Find average percentage methylation for each DMR
4) Create a dataframe that contains the average percent meth for each sample. Each row is a DMR, and each column is a sample.
if there is no overlap, a NA value will be printed.
'''

from __future__ import division, print_function
import sys
import pybedtools
import glob
import os
import pandas as pd
import argparse
from collections import Counter
import numpy as np
from datetime import datetime

parser = argparse.ArgumentParser(description='This program finds the average percentage of methylation of myCpG files for DMRs')
parser.add_argument('-e', '--region', help='region file containing DEGs and DMRs', required=True)
parser.add_argument('-i', '--inputDirectory', help='Directory that contains myCpG/methCall files', required=True)
parser.add_argument('-o', '--output', help='file with the methylation matrix in the region')

args = parser.parse_args()
DMRs = pd.read_table(args.region, index_col=None, sep="\t")
outdir = os.path.dirname(args.region)
DMRs.columns = ['Chrom','Start','End', "gene"]
DMRs.index = DMRs.Chrom + '.' + DMRs.Start.astype(str) + '.' + DMRs.End.astype(str)
a = pybedtools.BedTool.from_dataframe(DMRs)

### Find average percent methylation for each DMR
path = args.inputDirectory
listing =os.listdir(path)
#group = args.output.split('_')[0]
contrast = args.output.split('_')[0]
group = contrast.split("vs")
#listing = [item for item in listing if group in item]
listing = [item for item in listing if any(g in item for g in group)]

first = pd.read_table(path + "/" + listing[0], header=None, index_col=None, sep="\t")
first.columns = ['chrom','start','end', 'PercMeth', 'freqC', 'freqT']
first.index = first.chrom + '.' + first.start.astype(str) + '.' + first.end.astype(str)
methylation = pd.DataFrame(index=first.index)

for methCall in listing:
    fullpath= os.path.join(path, methCall)
    substrings=methCall.split('_')
    outputName="_".join(substrings[1:4])
    with open(fullpath,'r') as b:
        myDF = pd.read_table(b, header=None, index_col=None, sep="\t")
    myDF.columns=['chrom', 'start', 'end', 'PercMeth', 'freqC', 'freqT']
    myDF = myDF[(myDF.iloc[:, 4] + myDF.iloc[:, 5]) >= 5]
    myDF.index = myDF.chrom + '.' + myDF.start.astype(str) + '.' + myDF.end.astype(str)
    myDF.rename(columns={'PercMeth': outputName}, inplace=True)
    methylation = methylation.join(myDF.iloc[:,3], how='left')

methylation = methylation.dropna()
methylation['ID']=methylation.index
methylation[['chr', 'start', 'end']]=methylation['ID'].str.split('.', expand=True)
cols=methylation.columns.tolist()
cols=cols[-3:] + cols[:-4]
methylation=methylation.loc[:,cols]

d=pybedtools.BedTool.from_dataframe(methylation) ##this is adding decimal places to freqC and freqT
x=a.intersect(d, wa=True, wb=True)
ovlpDF=pybedtools.BedTool.to_dataframe(x, names=['chrom','start','end','gene','chr', 'Start', 'base' ] + list(methylation.columns[3:])) # freqC and freqT-columns 8/9 are float64, need to convert to int
ovlpDF=ovlpDF.drop(ovlpDF.columns[[4, 5, 6]], axis=1)
#Take mean percent meth for each DMR:
avgMeth=ovlpDF.groupby(['chrom','start','end','gene']).mean().reset_index()
avgMeth.to_csv(os.path.join(outdir, args.output), sep='\t', na_rep='NA', index=False)