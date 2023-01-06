#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import os

ct_names = {'Granulocyte Monocyte Progenitor (GMP)': 'GMP',
            'Monocytes':'Mono',
            'Granulocytes':'Gran',
            'Common Myeloid Progenitor (CMP)':'CMP',
            'Multipotent Progenitor (MPP)':'MPP', 
            'Bone Marrow Macrophages':'Macro',
	    'Lin-Sca1+Kit- bone marrow cells':'LSKBMC',
	    'B cells':'Bcell'}

cwd = getcwd()
srr_runs = pd.read_csv(os.path.join(cwd, "srr_metadata/srr_to_download.tsv"), sep = '\t')
srr_runs = srr_runs.set_index('Run')
srr_runs['ct_file'] = srr_runs['cell_type'].map(ct_names)

fq_dir = os.path.join(cwd, 'fastq')

for f in os.listdir(fq_dir):
    row = srr_runs.loc[f.split('.')[0]]
    ct_name = row['ct_file']
    exp_type = row['Assay Type']
    if exp_type == 'ChIP-Seq':
        exp_type = row['chip_antibody']
    if exp_type == 'OTHER':
        exp_type = 'ATAC-Seq'
    i=1
    new_file = '%s_%s_%i.fastq.gz' % (exp_type,ct_name,i)
    while os.path.isfile(os.path.join(fq_dir,new_file)):
        i+=1
        new_file = '%s_%s_%i.fastq.gz' % (exp_type,ct_name,i)
    os.rename(os.path.join(fq_dir,f),os.path.join(fq_dir,new_file))
    print(f,new_file)
