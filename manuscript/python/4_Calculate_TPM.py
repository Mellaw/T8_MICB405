# -*- coding: utf-8 -*-

#Sunset Enhancers: Tracing H3K27 Acetylation on Closed Chromatin in Myeloid Lineage Differentiation
#Authors: Melanie Law, Helena Sokolovska, Andrew Murtha, Kitoosepe Martens, Annice Li, and Kalen Dofher

"""
Created on Mon Nov 29 15:55:02 2021

@author: andym
"""

import os
import pandas as pd
import numpy as np

# =============================================================================
# 
# =============================================================================

folder = 'C:/Users/Helena/405_Linux_outputs/'
gene_length = pd.read_csv('C:/Users/andym/OneDrive/Documents/MICB405/gene_size.tsv', sep = '\t')
gene_length.columns= ['ensemble_id','gene_size']

gene_length['gene_size_kb'] = gene_length['gene_size'] / 1000

cts = pd.Series([f.split('_')[1] for f in os.listdir(folder) if '.tsv' not in f]).unique().tolist()

for ct in cts:
    dfs = [pd.read_csv(os.path.join(folder,f), header = None, sep = '\t', names = ['ensemble_id','gene_id','read_count']) for f in os.listdir(folder) if '_%s_'%ct in f]
    for i,df in enumerate(dfs):
        df = df.merge(gene_length[['ensemble_id','gene_size_kb']], on = 'ensemble_id')
        df['RPK'] = df['read_count'] / df['gene_size_kb']
        rpk_sf = df['RPK'].sum() / 1000000
        df['TPM'] = df['RPK']/rpk_sf
        df = df[['ensemble_id','gene_id','TPM']]
        dfs[i] = df.copy()
    df = dfs[0]
    for df1 in dfs[1:]:
        df = df.merge(df1, on = ['ensemble_id','gene_id'])
    df = df.set_index(['ensemble_id','gene_id'])
    df['TPM'] = df.median(axis = 1)
    df = df[['TPM']]
    df.to_csv(os.path.join(folder,'%s_geneExpression.tsv'%ct), sep = '\t', header = None)
