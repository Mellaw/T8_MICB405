# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 16:38:21 2021

@author: andym
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Patch
import scipy.stats as stats

# =============================================================================
# Helpers
# =============================================================================

labels = ['Closed, H3K27ac-',
          'Open, H3K27ac-',
          'Closed, H3K27ac+',
          'Open, H3K27ac+']

# Return 0 if closed and unacetylated. 1 if open and unacetylated.
# 2 if closed and acetylated. 3 if open and acetylated
def get_ct_status(ct, row):
    atac = row['%s_atac_significant' % ct]
    ac =  row['%s_ac_significant' % ct]
    if atac == 1 and ac == 1:
        return 3
    elif atac == 0 and ac == 1:
        return 2
    elif atac == 1 and ac == 0:
        return 1
    else:
        return 0
    

# =============================================================================
# Import enhancer status
# =============================================================================

enhs = pd.read_csv('C:/Users/andym/OneDrive/Documents/MICB405/enhancers_status.tsv', sep = '\t')

cell_types = ['MPP','CMP','GMP','Bcell','Mono','Gran']

for ct in cell_types:
    enhs['%s_status' % ct] = enhs.apply(lambda row: get_ct_status(ct, row), axis = 1)
    
enhs_status = enhs.copy().set_index(['chr','start','end','enhancer'])
enhs_status = enhs_status[[col for col in enhs_status if 'status' in col]].reset_index()

# =============================================================================
# 
# =============================================================================

g_e = pd.read_csv('C:/Users/andym/OneDrive/Documents/MICB405/enhancer_gene_map.txt', sep = '\t')
e_g_map = dict(zip(g_e['enhancer'],g_e['gene']))

fig,ax = plt.subplots()

enhs_status['gene_id'] = enhs_status['enhancer'].map(e_g_map)
xbase = 0
xticks = []
colors = ['grey','blue','green','k']

for ct in cell_types:
    print(ct)
    xs = np.arange(xbase,xbase+4,1)
    xticks.append(xs.mean())
    
    rna = pd.read_csv('C:/Users/andym/OneDrive/Documents/MICB405/htseq_output/%s_geneExpression.tsv'%ct,sep = '\t', header = None, names = ['ensbl_ID','gene_id','TPM'])
    ct_enh_status = enhs_status[['enhancer','%s_status'%ct,'gene_id']].merge(rna, on = 'gene_id')
    boxes = [ct_enh_status[ct_enh_status['%s_status'%ct]==status]['TPM'] for status in [0,1,2,3]]
    ax.boxplot(boxes, positions = xs, showfliers = False)
    for i,(l,x) in enumerate(zip(boxes,xs)):
        sc_x = pd.Series([x]*len(l)).apply(lambda x: x+np.random.uniform(-0.1,0.1))
        ax.scatter(sc_x,l,s = 2, lw = 0,alpha = 0.1,color = colors[i])
    for i in [0,1,3]:
        print(2,i,stats.mannwhitneyu(boxes[2],boxes[i],alternative='two-sided'))
    xbase = xbase + 5
    
    
ax.set_yscale('log')
ax.set_xticks(xticks)
ax.set_xticklabels(cell_types)
ax.tick_params(labelsize = 6)
ax.set_ylabel('TPM', fontsize = 6)

handles = [Patch(color = c) for c in colors]

ax.legend(handles, labels, fontsize = 6)

ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_linewidth(0.8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/nearest_geneExpression_ByEnhancerStatus.pdf')
fig.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/nearest_geneExpression_ByEnhancerStatus.png')