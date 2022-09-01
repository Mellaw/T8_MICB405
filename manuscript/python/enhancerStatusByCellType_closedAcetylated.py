# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 14:02:53 2021

@author: andym
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Patch

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
enhs_status = enhs_status[[col for col in enhs_status if 'status' in col]]

enhs_status = enhs_status[enhs_status['MPP_status'] == 2]

enh_valueCounts = enhs_status.apply(pd.value_counts) / len(enhs_status)
enh_valueCounts = enh_valueCounts.sort_index().fillna(0)

# =============================================================================
# Plot frequencies
# =============================================================================

fig,[ax1,ax2] = plt.subplots(ncols = 2,figsize = (7,2.5),sharey=True)

xbase = 0
xticks = []

colors = ['grey','blue','green','k']

for ct in cell_types:
    xs = np.arange(xbase,xbase+4,1)
    xticks.append(xs.mean())
    ax1.bar(xs, enh_valueCounts['%s_status'%ct], color = colors)
    for x in range(len(xs)):
        frac = enh_valueCounts.at[x,"%s_status" % ct]
        if ct != 'MPP':
            ax1.text(xs[x],frac+.02,"%.1f" % (frac*100),rotation = 90, ha = 'center', fontsize = 6)
    xbase = xbase + 5

ax1.set_xticks(xticks)
ax1.set_xticklabels(cell_types)
ax1.tick_params(labelsize = 6)
ax1.set_ylabel('Enhancer fraction', fontsize = 6)
ax1.set_ylim(0,1)

handles = [Patch(color = c) for c in colors]

# ax1.legend(handles, labels, fontsize = 6)

ax1.spines['left'].set_linewidth(0.8)
ax1.spines['bottom'].set_linewidth(0.8)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# =============================================================================
# Helpers
# =============================================================================

labels = ['Closed, H3K4me1-',
          'Open, H3K4me1-',
          'Closed, H3K4me1+',
          'Open, H3K4me1+']

# Return 0 if closed and unmethylated. 1 if open and unmethylated.
# 2 if closed and methylated. 3 if open and methylated
def get_ct_status(ct, row):
    atac = row['%s_atac_significant' % ct]
    me = row['%s_me_significant' % ct]
    if atac == 1 and me == 1:
        return 3
    elif atac == 0 and me == 1:
        return 2
    elif atac == 1 and me == 0:
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
enhs_status = enhs_status[[col for col in enhs_status if 'status' in col]]

enhs_status = enhs_status[enhs_status['MPP_status'] == 2]

enh_valueCounts = enhs_status.apply(pd.value_counts) / len(enhs_status)
enh_valueCounts = enh_valueCounts.sort_index().fillna(0)

# =============================================================================
# Plot frequencies
# =============================================================================


xbase = 0
xticks = []

colors = ['grey','blue','green','k']

for ct in cell_types:
    xs = np.arange(xbase,xbase+4,1)
    xticks.append(xs.mean())
    ax2.bar(xs, enh_valueCounts['%s_status'%ct], color = colors )
    for x in range(len(xs)):
        frac = enh_valueCounts.at[x,"%s_status" % ct]
        if ct != 'MPP':
            ax2.text(xs[x],frac+.02,"%.1f" % (frac*100),rotation = 90, ha = 'center', fontsize = 6)
    xbase = xbase + 5

ax2.set_xticks(xticks)
ax2.set_xticklabels(cell_types)
ax2.tick_params(labelsize = 6)

handles = [Patch(color = c) for c in colors]

ax2.spines['left'].set_linewidth(0.8)
ax2.spines['bottom'].set_linewidth(0.8)

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax2.legend(handles, labels, fontsize = 6)

plt.tight_layout()
plt.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_MPPclosedAcetylated_byCellType.pdf')
plt.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_MPPclosedAcetylated_byCellType.png')