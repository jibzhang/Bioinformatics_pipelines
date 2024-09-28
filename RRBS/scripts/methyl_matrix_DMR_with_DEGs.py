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
parser.add_argument('-e', '--region', help='file containing DEGs and DMRs', required=True)
parser.add_argument('-i', '--inputDirectory', help='Directory that contains myCpG/methCall files', required=True)
parser.add_argument('-o', '--output', help='Bedfile that contains chr, start, and end for DMR')

args = parser.parse_args()
DMRs = pd.read_table(args.region, index_col=None, sep="\t")
outdir = os.path.dirname(args.region)
DMRs.columns = ['chrom','start','end', "gene"]
DMRs.index = DMRs.chrom + '.' + DMRs.start.astype(str) + '.' + DMRs.end.astype(str) + '.' + DMRs.gene
FinalDF = pd.DataFrame(index=DMRs.index)
a = pybedtools.BedTool.from_dataframe(DMRs)

### Find average percent methylation for each DMR
path = args.inputDirectory
listing =os.listdir(path)

for methCall in listing:
    fullpath= os.path.join(path, methCall)
    substrings=methCall.split('_')
    outputName="_".join(substrings[1:4])
    with open(fullpath,'r') as b:
        myDF = pd.read_table(b, header=None, index_col=None, sep="\t")
    myDF.columns=['chr', 'base', 'base', 'PercMeth', 'freqC', 'freqT']
    d=pybedtools.BedTool.from_dataframe(myDF) ##this is adding decimal places to freqC and freqT
    x=a.intersect(d, wa=True, wb=True)
    ovlpDF=pybedtools.BedTool.to_dataframe(x) # freqC and freqT-columns 8/9 are float64, need to convert to int
    ovlpDF.columns = ['chrom','start','end','gene','chr', 'base', 'base', 'PercMeth', 'freqC', 'freqT']
    
    #Take mean percent meth for each DMR:
    avgMeth=ovlpDF.groupby(['chrom','start','end','gene'])['PercMeth'].mean().reset_index()
    avgMeth.index=avgMeth.chrom + '.' + avgMeth.start.astype(str) + '.' + avgMeth.end.astype(str) + '.' + avgMeth.gene
    avgMeth.rename(columns={'PercMeth':outputName}, inplace=True)
    FinalDF=FinalDF.join(avgMeth.iloc[:,4], how='left')

### Export final matrix that contains percent methylation
FinalDF['ID']=FinalDF.index
FinalDF = FinalDF.dropna()
FinalDF[['chr', 'start', 'end', 'gene']]=FinalDF['ID'].str.split('.', expand=True)
cols=FinalDF.columns.tolist()
cols=cols[-4:] + cols[:-4]
FinalDF=FinalDF.loc[:,cols]
FinalDF.drop(['ID'], axis=1, inplace=True)
FinalDF.to_csv(os.path.join(outdir, args.output), sep='\t', na_rep='NA', index=False)
