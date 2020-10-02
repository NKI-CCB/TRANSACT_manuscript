# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui
with code from Tycho Bismeijer

2020/01/03

READ DOWNLOADED DATA

This script reads files downloaded from TCGA using GDC tool.
"""

import os
import gzip
import pandas as pd


def read_metadata(pth):
    return pd.read_csv(pth, sep='\t')


def read_one_methyl_file(s, raw_folder):
    methyl_columns = ['Chromosome', 'Start', 'End', 'Gene_Symbol', 'Gene_Type']
    if '-' not in s:
        return False
    
    files = os.listdir(raw_folder + s)
    files = [f for f in files if '.txt' in f]
    
    if len(files) > 1:
        print('PROBLEME: MORE THAN ONE FILE FOR %s'%(s))
        
    methyl_df = pd.read_csv(raw_folder + s + '/' + files[0],
                            sep='\t')
    
    data_df = methyl_df[methyl_df.columns[:2]]
    data_df.columns = ['REF', 'beta_values']
    data_df = data_df.set_index('REF').T
    
    methyl_carac_df = methyl_df[methyl_df.columns[2:]]
    
    for c in methyl_columns:
        if c not in methyl_carac_df.columns:
            return s, data_df, methyl_carac_df
    return s, data_df, methyl_carac_df[methyl_columns]


def read_one_rnaseq_file(s, raw_folder):
    if '-' not in s:
        return False
    
    files = os.listdir(raw_folder + s)
    files = [f for f in files if 'htseq.counts.gz' in f or 'FPKM' in f or 'htseq_counts.txt.gz' in f]
    
    if len(files) > 1:
        print('PROBLEME: MORE THAN ONE FILE FOR %s'%(s))
    
    counts = gzip.open(raw_folder + s + '/' + files[0], 'rb')
    counts = pd.read_csv(raw_folder + s + '/' + files[0], compression='gzip', sep='\t', header=None)
    counts.columns = ['gene', 'count']
    counts = counts[counts.gene.str.contains('ENSG')]
    
    return s, counts, pd.DataFrame()

def read_one_miRNA_file(s, raw_folder):
    if '-' not in s:
        return s, pd.DataFrame(), pd.DataFrame()
        
    f = os.listdir(raw_folder + s)
    f = [e for e in f if 'quantification' in e][0]
    return s, pd.read_csv('%s%s/%s'%(raw_folder, s, f), sep='\t'), pd.DataFrame()