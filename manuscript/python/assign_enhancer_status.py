# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:42:09 2022

@author: andym
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import seaborn as sns
import os
import math

data_dir = os.path.join(getcwd(), "data/macs2_comb")

def get_atac(ct):
    s = 'atac'
    df = pd.read_csv(os.path.join(data_dir, 'ATAC-Seq_%s_minq.bdg' % ct), sep = '\t')
    df = df.rename(columns = {'significant':'%s_%s_significant' % (ct,s)})
    df = df[['chr','start','end','enhancer','%s_%s_significant' % (ct,s)]]
    return df

def get_ac(ct):
    s = 'ac'
    df = pd.read_csv(os.path.join(data_dir, 'H3K27Ac_%s_minq.bdg' % ct), sep = '\t')
    df = df.rename(columns = {'significant':'%s_%s_significant' % (ct,s)})
    df = df[['chr','start','end','enhancer','%s_%s_significant' % (ct,s)]]
    return df

def get_me(ct):
    s = 'me'
    df = pd.read_csv(os.path.join(data_dir, 'H3K4me1_%s_minq.bdg' % ct), sep = '\t')
    df = df.rename(columns = {'significant':'%s_%s_significant' % (ct,s)})
    df = df[['chr','start','end','enhancer','%s_%s_significant' % (ct,s)]]
    return df

atac = pd.read_csv(os.path.join(data_dir, 'ATAC-Seq_LSKBMC_minq.bdg'), sep = '\t')
me = pd.read_csv(os.path.join(data_dir, 'H3K4me1_MPP_minq.bdg'), sep = '\t')
ac = pd.read_csv(os.path.join(data_dir, 'H3K27ac_MPP_minq.bdg'), sep = '\t')
atac.head()

### Get only the closed DNA> This means FE < 0 and not in mpp_ac_peaks
enhs = pd.read_csv(os.path.join(getcwd(), 'data/references','mm10_enhancerAll.bed'), sep= '\t', header = None)
enhs.columns = ['chr','start','end','enhancer']
for s,df in zip(['atac','ac','me'],[atac,ac,me]):
    df = df.rename(columns = {'significant':'MPP_%s_significant' % s})
    enhs = enhs.merge(df[['chr','start','end','enhancer','MPP_%s_significant' % s]], on = ['chr','start','end','enhancer'], how = 'left')
    
    
### Get the acetylated enhancers
for ct in ['CMP','GMP','Bcell','Mono','Gran']:
    atac = get_atac(ct)
    me = get_me(ct)
    ac = get_ac(ct)
    for df in [atac,me,ac]:
        enhs = enhs.merge(df, on = ['chr','start','end','enhancer'], how = 'left')
enhs = enhs.fillna(0)

enhs.to_csv('data/enhancers_status.tsv', sep = '\t', index = None)
