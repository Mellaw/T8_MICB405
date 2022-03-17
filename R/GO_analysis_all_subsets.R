#Purpose:
#Filtering gene subsets by nearest enhancer status (open/closed chromatin, H3K27ac-/+)
#Writing gene subset catalogues
#GO term analysis of gene subsets

library(DESeq2)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)
library(topGO)
library(xlsx)

setwd("/Users/Helena/405_Linux_outputs/")

#import enhancer + nearest gene data tables
enhancer_status <- read.table(file = 'enhancers_status.tsv', sep = '\t', header = TRUE)
enhancer_gene_map <- read.table(file = 'enhancer_gene_map.txt', sep = '\t', header = TRUE)

#import expression data
MPP_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/MPP_geneExpression.tsv', sep = '\t', header = TRUE)
CMP_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/CMP_geneExpression.tsv', sep = '\t', header = TRUE)
GMP_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/GMP_geneExpression.tsv', sep = '\t', header = TRUE)
B_cell_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/Bcell_geneExpression.tsv', sep = '\t', header = TRUE)
Mono_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/Mono_geneExpression.tsv', sep = '\t', header = TRUE)
Gran_TPM <- read.table(file = 'drive-download-20220310T081839Z-001/Gran_geneExpression.tsv', sep = '\t', header = TRUE)

#merge enhancers + nearest genes
enhancer_status_with_genes <- merge(enhancer_gene_map, enhancer_status, by="enhancer")

#merge with ENSEMBL IDs + expression data
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes, MPP_TPM, by="gene", all.x = TRUE)
enhancer_status_with_genes_ensembl <- subset(enhancer_status_with_genes_ensembl, select = -c(ENSEMBL))
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes_ensembl, CMP_TPM, by="gene", all.x = TRUE)
enhancer_status_with_genes_ensembl <- subset(enhancer_status_with_genes_ensembl, select = -c(ENSEMBL))
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes_ensembl, GMP_TPM, by="gene", all.x = TRUE)
enhancer_status_with_genes_ensembl <- subset(enhancer_status_with_genes_ensembl, select = -c(ENSEMBL))
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes_ensembl, B_cell_TPM, by="gene", all.x = TRUE)
enhancer_status_with_genes_ensembl <- subset(enhancer_status_with_genes_ensembl, select = -c(ENSEMBL))
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes_ensembl, Mono_TPM, by="gene", all.x = TRUE)
enhancer_status_with_genes_ensembl <- subset(enhancer_status_with_genes_ensembl, select = -c(ENSEMBL))
enhancer_status_with_genes_ensembl <- merge(enhancer_status_with_genes_ensembl, Gran_TPM, by="gene", all.x = TRUE)

enhancers_and_genes <- relocate(enhancer_status_with_genes_ensembl, ENSEMBL, .before = "gene")

write.table(enhancers_and_genes, file = "enhancer_status_with_genes_and_expression.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)


### filtering gene subsets by nearest enhancer status (open/closed chromatin, H3K27ac-/+) ###

#open/ac-
MPP_closed_unac <- filter(enhancers_and_genes, MPP_atac_significant == 0, MPP_ac_significant == 0)
CMP_closed_unac <- filter(enhancers_and_genes, CMP_atac_significant == 0, CMP_ac_significant == 0)
GMP_closed_unac <- filter(enhancers_and_genes, GMP_atac_significant == 0, GMP_ac_significant == 0)
B_cell_closed_unac <- filter(enhancers_and_genes, Bcell_atac_significant == 0, Bcell_ac_significant == 0)
Mono_closed_unac <- filter(enhancers_and_genes, Mono_atac_significant == 0, Mono_ac_significant == 0)
Gran_closed_unac <- filter(enhancers_and_genes, Gran_atac_significant == 0, Gran_ac_significant == 0)
#closed/ac+
MPP_closed_ac <- filter(enhancers_and_genes, MPP_atac_significant == 0, MPP_ac_significant == 1)
CMP_closed_ac <- filter(enhancers_and_genes, CMP_atac_significant == 0, CMP_ac_significant == 1)
GMP_closed_ac <- filter(enhancers_and_genes, GMP_atac_significant == 0, GMP_ac_significant == 1)
B_cell_closed_ac <- filter(enhancers_and_genes, Bcell_atac_significant == 0, Bcell_ac_significant == 1)
Mono_closed_ac <- filter(enhancers_and_genes, Mono_atac_significant == 0, Mono_ac_significant == 1)
Gran_closed_ac <- filter(enhancers_and_genes, Gran_atac_significant == 0, Gran_ac_significant == 1)
#open/ac-
MPP_open_unac <- filter(enhancers_and_genes, MPP_atac_significant == 1, MPP_ac_significant == 0)
CMP_open_unac <- filter(enhancers_and_genes, CMP_atac_significant == 1, CMP_ac_significant == 0)
GMP_open_unac <- filter(enhancers_and_genes, GMP_atac_significant == 1, GMP_ac_significant == 0)
B_cell_open_unac <- filter(enhancers_and_genes, Bcell_atac_significant == 1, Bcell_ac_significant == 0)
Mono_open_unac <- filter(enhancers_and_genes, Mono_atac_significant == 1, Mono_ac_significant == 0)
Gran_open_unac <- filter(enhancers_and_genes, Gran_atac_significant == 1, Gran_ac_significant == 0)
#open/ac+
MPP_open_ac <- filter(enhancers_and_genes, MPP_atac_significant == 1, MPP_ac_significant == 1)
CMP_open_ac <- filter(enhancers_and_genes, CMP_atac_significant == 1, CMP_ac_significant == 1)
GMP_open_ac <- filter(enhancers_and_genes, GMP_atac_significant == 1, GMP_ac_significant == 1)
B_cell_open_ac <- filter(enhancers_and_genes, Bcell_atac_significant == 1, Bcell_ac_significant == 1)
Mono_open_ac <- filter(enhancers_and_genes, Mono_atac_significant == 1, Mono_ac_significant == 1)
Gran_open_ac <- filter(enhancers_and_genes, Gran_atac_significant == 1, Gran_ac_significant == 1)

### write gene subset catalogues ###

#closed/ac-
write.table(MPP_closed_unac, file='MPP_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(CMP_closed_unac, file='CMP_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(GMP_closed_unac, file='GMP_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(B_cell_closed_unac, file='B_cell_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Mono_closed_unac, file='Mono_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Gran_closed_unac, file='Gran_closed_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
#closed/ac+
write.table(MPP_closed_ac, file='MPP_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(CMP_closed_ac, file='CMP_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(GMP_closed_ac, file='GMP_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(B_cell_closed_ac, file='B_cell_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Mono_closed_ac, file='Mono_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Gran_closed_ac, file='Gran_closed_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
#open/ac-
write.table(MPP_open_unac, file='MPP_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(CMP_open_unac, file='CMP_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(GMP_open_unac, file='GMP_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(B_cell_open_unac, file='B_cell_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Mono_open_unac, file='Mono_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Gran_open_unac, file='Gran_open_unac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
#open/ac+
write.table(MPP_open_ac, file='MPP_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(CMP_open_ac, file='CMP_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(GMP_open_ac, file='GMP_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(B_cell_open_ac, file='B_cell_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Mono_open_ac, file='Mono_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)
write.table(Gran_open_ac, file='Gran_open_ac_gene_catalogue.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE)

### GO term analysis of gene subsets (nearest enhancer open/closed chromatin, H3K27ac-/+) ###

#add ENTREZ IDs to all genes (gene universe for GO analysis)
enhancers_and_genes$ENTREZ <- mapIds(
  org.Mm.eg.db,
  keys = enhancers_and_genes$ENSEMBL,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#closed/ac-
GO_MPP_closed_unac <- {
  
  MPP_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
     
  gene_subset <- MPP_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
   
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                                      geneIds = gene_subset,
                                      universeGeneIds = all_genes,
                                      annotation = "org.Mm.eg.db",
                                      ontology = "BP",
                                      pvalueCutoff = 0.05,
                                      conditional = FALSE,
                                      testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_unac <- {
  
  CMP_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_unac <- {
  
  GMP_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_unac <- {
  
  B_cell_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_unac <- {
  
  Mono_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_unac <- {
  
  Gran_closed_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#closed/ac+
GO_MPP_closed_ac <- {
  
  MPP_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_ac <- {
  
  CMP_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_ac <- {
  
  GMP_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_ac <- {
  
  B_cell_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_ac <- {
  
  Mono_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_ac <- {
  
  Gran_closed_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac-
GO_MPP_open_unac <- {
  
  MPP_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_unac <- {
  
  CMP_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_unac <- {
  
  GMP_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_unac <- {
  
  B_cell_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_unac <- {
  
  Mono_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_unac <- {
  
  Gran_open_unac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_unac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_unac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_unac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac+
GO_MPP_open_ac <- {
  
  MPP_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_ac <- {
  
  CMP_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_ac <- {
  
  GMP_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_ac <- {
  
  B_cell_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_ac <- {
  
  Mono_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_ac <- {
  
  Gran_open_ac$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_ac$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_ac %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_ac_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}

### GO term analysis of genes with highest and lowest expression in each subset, defined by top and top quartile of TPM ###

#extract top quantile of gene expression in each subset

#closed/ac-
MPP_closed_unac_TPM_quantiles <- quantile(MPP_closed_unac$MPP_TPM, na.rm = TRUE)
MPP_closed_unac_top_quantile <- subset(MPP_closed_unac, MPP_TPM >= MPP_closed_unac_TPM_quantiles[4])
CMP_closed_unac_TPM_quantiles <- quantile(CMP_closed_unac$CMP_TPM, na.rm = TRUE)
CMP_closed_unac_top_quantile <- subset(CMP_closed_unac, CMP_TPM >= CMP_closed_unac_TPM_quantiles[4])
GMP_closed_unac_TPM_quantiles <- quantile(GMP_closed_unac$GMP_TPM, na.rm = TRUE)
GMP_closed_unac_top_quantile <- subset(GMP_closed_unac, GMP_TPM >= GMP_closed_unac_TPM_quantiles[4])
B_cell_closed_unac_TPM_quantiles <- quantile(B_cell_closed_unac$B_cell_TPM, na.rm = TRUE)
B_cell_closed_unac_top_quantile <- subset(B_cell_closed_unac, B_cell_TPM >= B_cell_closed_unac_TPM_quantiles[4])
Mono_closed_unac_TPM_quantiles <- quantile(Mono_closed_unac$Mono_TPM, na.rm = TRUE)
Mono_closed_unac_top_quantile <- subset(Mono_closed_unac, Mono_TPM >= Mono_closed_unac_TPM_quantiles[4])
Gran_closed_unac_TPM_quantiles <- quantile(Gran_closed_unac$Gran_TPM, na.rm = TRUE)
Gran_closed_unac_top_quantile <- subset(Gran_closed_unac, Gran_TPM >= Gran_closed_unac_TPM_quantiles[4])
#closed/ac+
MPP_closed_ac_TPM_quantiles <- quantile(MPP_closed_ac$MPP_TPM, na.rm = TRUE)
MPP_closed_ac_top_quantile <- subset(MPP_closed_ac, MPP_TPM >= MPP_closed_ac_TPM_quantiles[4])
CMP_closed_ac_TPM_quantiles <- quantile(CMP_closed_ac$CMP_TPM, na.rm = TRUE)
CMP_closed_ac_top_quantile <- subset(CMP_closed_ac, CMP_TPM >= CMP_closed_ac_TPM_quantiles[4])
GMP_closed_ac_TPM_quantiles <- quantile(GMP_closed_ac$GMP_TPM, na.rm = TRUE)
GMP_closed_ac_top_quantile <- subset(GMP_closed_ac, GMP_TPM >= GMP_closed_ac_TPM_quantiles[4])
B_cell_closed_ac_TPM_quantiles <- quantile(B_cell_closed_ac$B_cell_TPM, na.rm = TRUE)
B_cell_closed_ac_top_quantile <- subset(B_cell_closed_ac, B_cell_TPM >= B_cell_closed_ac_TPM_quantiles[4])
Mono_closed_ac_TPM_quantiles <- quantile(Mono_closed_ac$Mono_TPM, na.rm = TRUE)
Mono_closed_ac_top_quantile <- subset(Mono_closed_ac, Mono_TPM >= Mono_closed_ac_TPM_quantiles[4])
Gran_closed_ac_TPM_quantiles <- quantile(Gran_closed_ac$Gran_TPM, na.rm = TRUE)
Gran_closed_ac_top_quantile <- subset(Gran_closed_ac, Gran_TPM >= Gran_closed_ac_TPM_quantiles[4])
#open/ac-
MPP_open_unac_TPM_quantiles <- quantile(MPP_open_unac$MPP_TPM, na.rm = TRUE)
MPP_open_unac_top_quantile <- subset(MPP_open_unac, MPP_TPM >= MPP_open_unac_TPM_quantiles[4])
CMP_open_unac_TPM_quantiles <- quantile(CMP_open_unac$CMP_TPM, na.rm = TRUE)
CMP_open_unac_top_quantile <- subset(CMP_open_unac, CMP_TPM >= CMP_open_unac_TPM_quantiles[4])
GMP_open_unac_TPM_quantiles <- quantile(GMP_open_unac$GMP_TPM, na.rm = TRUE)
GMP_open_unac_top_quantile <- subset(GMP_open_unac, GMP_TPM >= GMP_open_unac_TPM_quantiles[4])
B_cell_open_unac_TPM_quantiles <- quantile(B_cell_open_unac$B_cell_TPM, na.rm = TRUE)
B_cell_open_unac_top_quantile <- subset(B_cell_open_unac, B_cell_TPM >= B_cell_open_unac_TPM_quantiles[4])
Mono_open_unac_TPM_quantiles <- quantile(Mono_open_unac$Mono_TPM, na.rm = TRUE)
Mono_open_unac_top_quantile <- subset(Mono_open_unac, Mono_TPM >= Mono_open_unac_TPM_quantiles[4])
Gran_open_unac_TPM_quantiles <- quantile(Gran_open_unac$Gran_TPM, na.rm = TRUE)
Gran_open_unac_top_quantile <- subset(Gran_open_unac, Gran_TPM >= Gran_open_unac_TPM_quantiles[4])
#open/ac+
MPP_open_ac_TPM_quantiles <- quantile(MPP_open_ac$MPP_TPM, na.rm = TRUE)
MPP_open_ac_top_quantile <- subset(MPP_open_ac, MPP_TPM >= MPP_open_ac_TPM_quantiles[4])
CMP_open_ac_TPM_quantiles <- quantile(CMP_open_ac$CMP_TPM, na.rm = TRUE)
CMP_open_ac_top_quantile <- subset(CMP_open_ac, CMP_TPM >= CMP_open_ac_TPM_quantiles[4])
GMP_open_ac_TPM_quantiles <- quantile(GMP_open_ac$GMP_TPM, na.rm = TRUE)
GMP_open_ac_top_quantile <- subset(GMP_open_ac, GMP_TPM >= GMP_open_ac_TPM_quantiles[4])
B_cell_open_ac_TPM_quantiles <- quantile(B_cell_open_ac$B_cell_TPM, na.rm = TRUE)
B_cell_open_ac_top_quantile <- subset(B_cell_open_ac, B_cell_TPM >= B_cell_open_ac_TPM_quantiles[4])
Mono_open_ac_TPM_quantiles <- quantile(Mono_open_ac$Mono_TPM, na.rm = TRUE)
Mono_open_ac_top_quantile <- subset(Mono_open_ac, Mono_TPM >= Mono_open_ac_TPM_quantiles[4])
Gran_open_ac_TPM_quantiles <- quantile(Gran_open_ac$Gran_TPM, na.rm = TRUE)
Gran_open_ac_top_quantile <- subset(Gran_open_ac, Gran_TPM >= Gran_open_ac_TPM_quantiles[4])

#closed/ac-
GO_MPP_closed_unac_top_quantile <- {
  
  MPP_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_unac_top_quantile <- {
  
  CMP_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_unac_top_quantile <- {
  
  GMP_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_unac_top_quantile <- {
  
  B_cell_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_unac_top_quantile <- {
  
  Mono_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_unac_top_quantile <- {
  
  Gran_closed_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#closed/ac+
GO_MPP_closed_ac_top_quantile <- {
  
  MPP_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_ac_top_quantile <- {
  
  CMP_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_ac_top_quantile <- {
  
  GMP_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_ac_top_quantile <- {
  
  B_cell_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_ac_top_quantile <- {
  
  Mono_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_ac_top_quantile <- {
  
  Gran_closed_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac-
GO_MPP_open_unac_top_quantile <- {
  
  MPP_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_unac_top_quantile <- {
  
  CMP_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_unac_top_quantile <- {
  
  GMP_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_unac_top_quantile <- {
  
  B_cell_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_unac_top_quantile <- {
  
  Mono_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_unac_top_quantile <- {
  
  Gran_open_unac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_unac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_unac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_unac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac+
GO_MPP_open_ac_top_quantile <- {
  
  MPP_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_ac_top_quantile <- {
  
  CMP_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_ac_top_quantile <- {
  
  GMP_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_ac_top_quantile <- {
  
  B_cell_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_ac_top_quantile <- {
  
  Mono_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_ac_top_quantile <- {
  
  Gran_open_ac_top_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_ac_top_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_ac_top_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_ac_top_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}

#extract bottom quantile of gene expression in each subset

#closed/ac-
MPP_closed_unac_bottom_quantile <- subset(MPP_closed_unac, MPP_TPM <= MPP_closed_unac_TPM_quantiles[2])
CMP_closed_unac_bottom_quantile <- subset(CMP_closed_unac, CMP_TPM <= CMP_closed_unac_TPM_quantiles[2])
GMP_closed_unac_bottom_quantile <- subset(GMP_closed_unac, GMP_TPM <= GMP_closed_unac_TPM_quantiles[2])
B_cell_closed_unac_bottom_quantile <- subset(B_cell_closed_unac, B_cell_TPM <= B_cell_closed_unac_TPM_quantiles[2])
Mono_closed_unac_bottom_quantile <- subset(Mono_closed_unac, Mono_TPM <= Mono_closed_unac_TPM_quantiles[2])
Gran_closed_unac_bottom_quantile <- subset(Gran_closed_unac, Gran_TPM <= Gran_closed_unac_TPM_quantiles[2])
#closed/ac+
MPP_closed_ac_bottom_quantile <- subset(MPP_closed_ac, MPP_TPM <= MPP_closed_ac_TPM_quantiles[2])
CMP_closed_ac_bottom_quantile <- subset(CMP_closed_ac, CMP_TPM <= CMP_closed_ac_TPM_quantiles[2])
GMP_closed_ac_bottom_quantile <- subset(GMP_closed_ac, GMP_TPM <= GMP_closed_ac_TPM_quantiles[2])
B_cell_closed_ac_bottom_quantile <- subset(B_cell_closed_ac, B_cell_TPM <= B_cell_closed_ac_TPM_quantiles[2])
Mono_closed_ac_bottom_quantile <- subset(Mono_closed_ac, Mono_TPM <= Mono_closed_ac_TPM_quantiles[2])
Gran_closed_ac_bottom_quantile <- subset(Gran_closed_ac, Gran_TPM <= Gran_closed_ac_TPM_quantiles[2])
#open/ac-
MPP_open_unac_bottom_quantile <- subset(MPP_open_unac, MPP_TPM <= MPP_open_unac_TPM_quantiles[2])
CMP_open_unac_bottom_quantile <- subset(CMP_open_unac, CMP_TPM <= CMP_open_unac_TPM_quantiles[2])
GMP_open_unac_bottom_quantile <- subset(GMP_open_unac, GMP_TPM <= GMP_open_unac_TPM_quantiles[2])
B_cell_open_unac_bottom_quantile <- subset(B_cell_open_unac, B_cell_TPM <= B_cell_open_unac_TPM_quantiles[2])
Mono_open_unac_bottom_quantile <- subset(Mono_open_unac, Mono_TPM <= Mono_open_unac_TPM_quantiles[2])
Gran_open_unac_bottom_quantile <- subset(Gran_open_unac, Gran_TPM <= Gran_open_unac_TPM_quantiles[2])
#open/ac+
MPP_open_ac_bottom_quantile <- subset(MPP_open_ac, MPP_TPM <= MPP_open_ac_TPM_quantiles[2])
CMP_open_ac_bottom_quantile <- subset(CMP_open_ac, CMP_TPM <= CMP_open_ac_TPM_quantiles[2])
GMP_open_ac_bottom_quantile <- subset(GMP_open_ac, GMP_TPM <= GMP_open_ac_TPM_quantiles[2])
B_cell_open_ac_bottom_quantile <- subset(B_cell_open_ac, B_cell_TPM <= B_cell_open_ac_TPM_quantiles[2])
Mono_open_ac_bottom_quantile <- subset(Mono_open_ac, Mono_TPM <= Mono_open_ac_TPM_quantiles[2])
Gran_open_ac_bottom_quantile <- subset(Gran_open_ac, Gran_TPM <= Gran_open_ac_TPM_quantiles[2])

#closed/ac-
GO_MPP_closed_unac_bottom_quantile <- {
  
  MPP_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_unac_bottom_quantile <- {
  
  CMP_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_unac_bottom_quantile <- {
  
  GMP_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_unac_bottom_quantile <- {
  
  B_cell_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_unac_bottom_quantile <- {
  
  Mono_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_unac_bottom_quantile <- {
  
  Gran_closed_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#closed/ac+
GO_MPP_closed_ac_bottom_quantile <- {
  
  MPP_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_closed_ac_bottom_quantile <- {
  
  CMP_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_closed_ac_bottom_quantile <- {
  
  GMP_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_closed_ac_bottom_quantile <- {
  
  B_cell_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_closed_ac_bottom_quantile <- {
  
  Mono_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_closed_ac_bottom_quantile <- {
  
  Gran_closed_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_closed_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_closed_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_closed_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac-
GO_MPP_open_unac_bottom_quantile <- {
  
  MPP_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_unac_bottom_quantile <- {
  
  CMP_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_unac_bottom_quantile <- {
  
  GMP_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_unac_bottom_quantile <- {
  
  B_cell_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_unac_bottom_quantile <- {
  
  Mono_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_unac_bottom_quantile <- {
  
  Gran_open_unac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_unac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_unac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_unac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
#open/ac+
GO_MPP_open_ac_bottom_quantile <- {
  
  MPP_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = MPP_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- MPP_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='MPP_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_CMP_open_ac_bottom_quantile <- {
  
  CMP_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = CMP_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- CMP_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='CMP_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_GMP_open_ac_bottom_quantile <- {
  
  GMP_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = GMP_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- GMP_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='GMP_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_B_cell_open_ac_bottom_quantile <- {
  
  B_cell_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = B_cell_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- B_cell_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='B_cell_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Mono_open_ac_bottom_quantile <- {
  
  Mono_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Mono_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Mono_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Mono_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}
GO_Gran_open_ac_bottom_quantile <- {
  
  Gran_open_ac_bottom_quantile$ENTREZ <- mapIds(
    org.Mm.eg.db,
    keys = Gran_open_ac_bottom_quantile$ENSEMBL,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  gene_subset <- Gran_open_ac_bottom_quantile %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  all_genes <- enhancers_and_genes %>% 
    as.data.frame() %>% 
    pull(ENTREZ) %>% 
    unique()
  
  GO_run <- hyperGTest(new("GOHyperGParams",
                           geneIds = gene_subset,
                           universeGeneIds = all_genes,
                           annotation = "org.Mm.eg.db",
                           ontology = "BP",
                           pvalueCutoff = 0.05,
                           conditional = FALSE,
                           testDirection = "over"))
  
  write.table(summary(GO_run), file='Gran_open_ac_bottom_quantile_GO_analysis.tsv', quote=FALSE, sep='\t', row.names = FALSE, col.names = TRUE) 
  
}