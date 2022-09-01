# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:12:41 2021

@author: andym
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
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
# 
# =============================================================================

enhs = pd.read_csv('C:/Users/andym/OneDrive/Documents/MICB405/enhancers_status.tsv', sep = '\t')

cell_types = ['MPP','CMP','GMP','Bcell','Mono','Gran']

for ct in cell_types:
    enhs['%s_status' % ct] = enhs.apply(lambda row: get_ct_status(ct, row), axis = 1)
    
enhs_status = enhs.copy().set_index(['chr','start','end','enhancer'])
enhs_status = enhs_status[[col for col in enhs_status if 'status' in col]].reset_index()

# enhs_status = enhs_status[enhs_status['MPP_status'] == 2]

enhs = enhs[['chr','start','end','enhancer']+[col for col in enhs.columns if '_me_' in col]].merge(enhs_status, on = ['chr','start','end','enhancer'])

# =============================================================================
# 
# =============================================================================

fig,ax = plt.subplots()

xbase = 0
xticks = []

colors = ['grey','blue','green','k']

for ct in cell_types:
    xs = np.arange(xbase,xbase+4,1)
    xticks.append(xs.mean())
    enhs_ct = enhs[[col for col in enhs.columns if ct in col]]
    enhs_ct_me = enhs_ct.groupby('%s_status'%ct).sum()
    enhs_ct_status = enhs_ct.groupby('%s_status'%ct).count()
    enhs_ct_me['me_percent'] = enhs_ct_me['%s_me_significant'%ct] / enhs_ct_status['%s_me_significant'%ct]*100
    ax.bar(xs, enhs_ct_me['me_percent'], color = colors )
    for i in [0,1,2,3]:
        frac = enhs_ct_me.at[i,"me_percent"]
        ax.text(xs[i],frac+2,"%.1f" % (frac),rotation = 90, ha = 'center', fontsize = 6)
    # frac = enhs_ct_me.at[1,"me_percent"]
    # ax.text(xs[1],frac+.02,"%.1f" % (frac),rotation = 90, ha = 'center', fontsize = 6)
    xbase = xbase + 5


ax.set_xticks(xticks)
ax.set_xticklabels(cell_types)
ax.tick_params(labelsize = 8)

handles = [Patch(color = c) for c in colors]

ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_linewidth(0.8)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_ylabel('Percent H3K4me1')

ax.legend(handles, labels, fontsize = 6)

plt.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_methylation_percent_byAC_ATAC.pdf')
plt.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_methylation_percent_byAC_ATAC.png')