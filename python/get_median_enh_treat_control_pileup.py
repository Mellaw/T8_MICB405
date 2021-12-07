import os
import numpy as np
import pandas as pd

folder = "/projects/micb405/analysis/GROUP8/comb_bdg"

for f in os.listdir(folder):
	if not 'FE' in f: continue;
	cell_type = f.split("_FE")[0]
	print(cell_type)
	## Import compiled bedgraph file
	df = pd.read_csv(os.path.join(folder,f), sep = '\t', header = None)
	df.columns = ['chr','start','end','enhancer','chr_dup','start_dup','end_dup','FE']
	df['overlap_percent'] = (df[['end','end_dup']].min(axis = 1) - df[['start','start_dup']].max(axis = 
1))/(df['end_dup']-df['start_dup'])
	df = df[df['overlap_percent'] >= 0.5]

	## Groupby enhancer and control or treatment. Take median values
	df = df[['chr','start','end','enhancer','FE']]
	df_grp = df.groupby(['enhancer']).median().reset_index()

	## Merge chromosomes back onto dataframes
	df_enh =  pd.read_csv('/projects/micb405/analysis/GROUP8/mm10_enhancerAll.bed', sep = '\t', header = None)
	df_enh.columns = ['chr','start','end','enhancer']
	df_enh['start'] = df_enh['start'].astype(int)
	df_enh['end'] = df_enh['end'].astype(int)

	df_grp = df_grp.merge(df_enh, on = ['start','end','enhancer'], how = 'left')

	## Sort both files

	df_grp['chr_num'] = df_grp['chr'].str.extract('(\d+)').astype(int)

	df_grp = df_grp.sort_values(['chr_num','start'], ascending = True).drop_duplicates().drop('chr_num', axis = 1)
	df_grp['log_FE'] = df_grp['FE'].apply(lambda x: np.log2(x))

	df_grp = df_grp[['chr','start','end','enhancer','FE','log_FE']]

	df_grp.to_csv(os.path.join(folder,"%s_FE.bdg" % cell_type), sep = '\t', index = None)



