# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 16:41:53 2021

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

enhs_status = enhs_status[(enhs_status == 2).any(axis=1)]

enh_melt = enhs_status.reset_index().melt(id_vars = ['chr','start','end','enhancer'], 
                                           var_name = 'cell_type', value_name = 'status')

enh_melt['cell_type'] = enh_melt['cell_type'].str.split('_').str[0]

# =============================================================================
# Set order
# =============================================================================

enh_order = enhs_status.copy()

for i,col in enumerate(enh_order.columns):
    enh_order[col] = 4**(len(enh_order.columns)-i-1) * (enh_order[col]+1)
    
enh_order['sum'] = enh_order.sum(axis = 1)
enh_order = enh_order.sort_values('sum').reset_index()
y_order = dict(zip(enh_order['enhancer'],np.arange(len(enh_order))))
x_order = dict(zip(cell_types,np.arange(len(cell_types))))

colors = ['grey','blue','green','k']
color_dict = dict(zip([0,1,2,3],colors))


enh_melt['y'] = enh_melt['enhancer'].map(y_order)
enh_melt['x'] = enh_melt['cell_type'].map(x_order)
enh_melt['color'] = enh_melt['status'].map(color_dict)

# =============================================================================
# 
# =============================================================================

fig,ax = plt.subplots(figsize = (3,5))

ax.bar(enh_melt['x'],1.,bottom=enh_melt['y'],color=enh_melt['color'])

ax.set_ylim(0,len(enhs_status))

ax.set_ylabel("Enhancer")
ax.tick_params(left = False, labelleft = False, bottom = False)
ax.set_xticks(list(x_order.values()))
ax.set_xticklabels(list(x_order.keys()))

ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)

fig.tight_layout()

fig.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_MPPclosedAcetylated_byCellType_heatmap_v2.png', dpi = 500)
fig.savefig('C:/Users/andym/OneDrive/Documents/MICB405/figures/enhancer_MPPclosedAcetylated_byCellType_heatmap_v2.pdf')