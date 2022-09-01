import os
import numpy as np
import pandas as pd

folder = "/projects/micb405/analysis/GROUP8/comb_bdg"

for f in os.listdir(folder):
	if not '_q.bdg' in f: continue;
	cell_type = f.split("_q")[0]
	print(cell_type)
	## Import compiled bedgraph file
	df = pd.read_csv(os.path.join(folder,f), sep = '\t', header = None)
	df.columns = ['chr','start','end','enhancer','chr_dup','start_dup','end_dup','-log10q']
	df['overlap_percent'] = (df[['end','end_dup']].min(axis = 1) - df[['start','start_dup']].max(axis = 
1))/(df['end_dup']-df['start_dup'])
	df = df[df['overlap_percent'] >= 0.5]

	## Groupby enhancer and control or treatment. Take median values
	df = df[['chr','start','end','enhancer','-log10q']]
	df_grp = df.groupby(['enhancer']).max().reset_index()

	## Sort both files

	df_grp['chr_num'] = df_grp['chr'].str.extract('(\d+)').astype(int)

	df_grp = df_grp.sort_values(['chr_num','start'], ascending = True).drop_duplicates().drop('chr_num', axis = 1)
	df_grp['q'] = df_grp['-log10q'].apply(lambda x: 10**(-1*x))
	df_grp['significant'] = 0
	df_grp.loc[df_grp['q'] <= 0.05, 'significant'] = 1
	df_grp = df_grp[['chr','start','end','enhancer','-log10q','q','significant']]

	df_grp.to_csv(os.path.join(folder,"%s_minq.bdg" % cell_type), sep = '\t', index = None)



