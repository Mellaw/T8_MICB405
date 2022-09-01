import pandas as pd

cell_types = ['Granulocyte Monocyte Progenitor (GMP)', 'Common Myeloid Progenitor (CMP)','Multipotent Progenitor (MPP)',
              'Monocytes','Granulocytes', 'Bone Marrow Macrophages','Lin-Sca1+Kit- bone marrow cells','B cells']

# Read data
srr_runs = pd.read_csv('/projects/micb405/analysis/GROUP8/srr_metadata/SraRunTable.txt')
# Keep cell types
srr_runs = srr_runs[srr_runs['cell_type'].isin(cell_types)]
# Keep wanted experiments
srr_runs = srr_runs[~srr_runs['chip_antibody'].isin(['H3K4me2', 'H3K4me3'])]

srr_runs.to_csv('/projects/micb405/analysis/GROUP8/srr_metadata/srr_to_download.tsv', index = None, sep = '\t')
rna = srr_runs[srr_runs['Assay Type'] == 'RNA-Seq'].copy()
rna.to_csv('/projects/micb405/analysis/GROUP8/srr_metadata/rna_to_download.tsv', index = None, sep = '\t')

