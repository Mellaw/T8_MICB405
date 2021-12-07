# Kitty Martens
# GO analysis of tpm data

# 5 Dec 2021

#load necessary packages
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)
library(dplyr)

dir <- "TPM"

#import gene-enhancer list
Enhancers <- read_tsv("all_closest_gene_expression.tsv")
Enhancers 

#Import enhancer status
Enhancer_Status <- read_tsv("enhancer_status_byCellType.tsv")
Enhancer_Status

# assign file path to a variable
sample_metadata <- read_csv(file.path(dir, "TPMsamples.csv"))
sample_metadata

# assign file paths to TPM data to a variable
files <- file.path(dir, sample_metadata$Filepath)
files
files[1]

# extract tsv's from the folder, extract top quantile of TPM
#for each cell type

#MPP
MPP_TPM <- read_tsv(files[1], 
                 col_names = c("EnsemblID", "Gene", "TPM"))
#add ENTREZ ids
MPP_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = MPP_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
MPP_TPM_quantiles <- quantile(MPP_TPM$TPM)

MPP_TPM <- subset(MPP_TPM, TPM > MPP_TPM_quantiles[4])

#extract only genes found in the list of enhancers
MPP_TPM <-  subset(MPP_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
MPP_TPM_enhancers <- subset(Enhancers, Gene %in% c(MPP_TPM$Gene))
MPP_TPM_enhancers <- merge(MPP_TPM_enhancers, MPP_TPM, by = "Gene")

#extract only enhancers that are closed/unacetylated (status = 0)
MPP_TPM_Closed_Inactive <- subset(Enhancer_Status, MPP_status == 0)
MPP_TPM_CI <- MPP_TPM_enhancers %>% 
  filter(Enhancer %in% MPP_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
MPP_TPM_Closed_Active <- subset(Enhancer_Status, MPP_status == 1)
MPP_TPM_CA <- MPP_TPM_enhancers %>% 
  filter(Enhancer %in% MPP_TPM_Closed_Active$enhancer)

#GO
GOTPM_MPP <- hyperGTest(new("GOHyperGParams",
                                    geneIds = MPP_TPM_CA$entrez,
                                    universeGeneIds = MPP_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_MPP  %>% summary() %>% head(10)

write.csv(summary(GOTPM_MPP), file = "GOTPM_MPP.csv")



#CMP
CMP_TPM <- read_tsv(files[2],
                    col_names = c("EnsemblID", "Gene", "TPM"))

#add ENTREZ ids
CMP_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = CMP_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
CMP_TPM_quantiles <- quantile(CMP_TPM$TPM)

CMP_TPM <- subset(CMP_TPM, TPM > CMP_TPM_quantiles[4])

#extract only genes found in the list of enhancers
CMP_TPM_enhancers <- subset(Enhancers, Gene %in% c(CMP_TPM$Gene))
CMP_TPM <-  subset(CMP_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
CMP_TPM_enhancers <- merge(CMP_TPM_enhancers, CMP_TPM, by = "Gene", all.y = TRUE)

#extract only enhancers that are closed/unacetylated (status = 0)
CMP_TPM_Closed_Inactive <- subset(Enhancer_Status, CMP_status == 0)
CMP_TPM_CI <- CMP_TPM_enhancers %>% 
  filter(Enhancer %in% CMP_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
CMP_TPM_Closed_Active <- subset(Enhancer_Status, CMP_status == 1)
CMP_TPM_CA <- CMP_TPM_enhancers %>% 
  filter(Enhancer %in% CMP_TPM_Closed_Active$enhancer)

#GO
GOTPM_CMP <- hyperGTest(new("GOHyperGParams",
                                    geneIds = CMP_TPM_CA$entrez,
                                    universeGeneIds = CMP_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_CMP  %>% summary() %>% head(10)

write.csv(summary(GOTPM_CMP), file = "GOTPM_CMP.csv")

#GMP
GMP_TPM <- read_tsv(files[3],
                   col_names = c("EnsemblID", "Gene", "TPM"))

#add ENTREZ ids
GMP_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = GMP_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
GMP_TPM_quantiles <- quantile(GMP_TPM$TPM)

GMP_TPM <- subset(GMP_TPM, TPM > GMP_TPM_quantiles[4])

#extract only genes found in the list of enhancers
GMP_TPM_enhancers <- subset(Enhancers, Gene %in% c(GMP_TPM$Gene))
GMP_TPM <-  subset(GMP_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
GMP_TPM_enhancers <- merge(GMP_TPM_enhancers, GMP_TPM, by = "Gene", all.y = TRUE)

#extract only enhancers that are closed/unacetylated (status = 0)
GMP_TPM_Closed_Inactive <- subset(Enhancer_Status, GMP_status == 0)
GMP_TPM_CI <- GMP_TPM_enhancers %>% 
  filter(Enhancer %in% GMP_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
GMP_TPM_Closed_Active <- subset(Enhancer_Status, GMP_status == 1)
GMP_TPM_CA <- GMP_TPM_enhancers %>% 
  filter(Enhancer %in% GMP_TPM_Closed_Active$enhancer)

#GO
GOTPM_GMP <- hyperGTest(new("GOHyperGParams",
                                    geneIds = GMP_TPM_CA$entrez,
                                    universeGeneIds = GMP_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_GMP %>% summary() %>% head(10)

write.csv(summary(GOTPM_GMP), file = "GOTMP_GMP.csv")

#Gran

Gran_TPM <- read_tsv(files[4],
                     col_names = c("EnsemblID", "Gene", "TPM"))

#add ENTREZ ids
Gran_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = Gran_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
Gran_TPM_quantiles <- quantile(Gran_TPM$TPM)

Gran_TPM <- subset(Gran_TPM, TPM > Gran_TPM_quantiles[4])

#extract only genes found in the list of enhancers
Gran_TPM <-  subset(Gran_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
Gran_TPM_enhancers <- subset(Enhancers, Gene %in% c(Gran_TPM$Gene))
Gran_TPM_enhancers <- merge(Gran_TPM_enhancers, Gran_TPM, by = "Gene", all.y = TRUE)

#extract only enhancers that are closed/unacetylated (status = 0)
Gran_TPM_Closed_Inactive <- subset(Enhancer_Status, Gran_status == 0)
Gran_TPM_CI <- Gran_TPM_enhancers %>% 
  filter(Enhancer %in% Gran_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
Gran_TPM_Closed_Active <- subset(Enhancer_Status, Gran_status == 1)
Gran_TPM_CA <- Gran_TPM_enhancers %>% 
  filter(Enhancer %in% Gran_TPM_Closed_Active$enhancer)

#GO
GOTPM_Gran <- hyperGTest(new("GOHyperGParams",
                                    geneIds = Gran_TPM_CA$entrez,
                                    universeGeneIds = Gran_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_Gran %>% summary() %>% head(10)

write.csv(summary(GOTPM_Gran), file = "GOTPM_Gran.csv")

#Mono

Mono_TPM <- read_tsv(files[6],
                     col_names = c("EnsemblID", "Gene", "TPM"))

#add ENTREZ ids
Mono_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = Mono_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
Mono_TPM_quantiles <- quantile(Mono_TPM$TPM)

Mono_TPM <- subset(Mono_TPM, TPM > Mono_TPM_quantiles[4])

#extract only genes found in the list of enhancers
Mono_TPM <-  subset(Mono_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
Mono_TPM_enhancers <- subset(Enhancers, Gene %in% c(Mono_TPM$Gene))
Mono_TPM_enhancers <- merge(Mono_TPM_enhancers, Mono_TPM, by = "Gene", all.y = TRUE)

#extract only enhancers that are closed/unacetylated (status = 0)
Mono_TPM_Closed_Inactive <- subset(Enhancer_Status, Mono_status == 0)
Mono_TPM_CI <- Mono_TPM_enhancers %>% 
  filter(Enhancer %in% Mono_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
Mono_TPM_Closed_Active <- subset(Enhancer_Status, Mono_status == 1)
Mono_TPM_CA <- Mono_TPM_enhancers %>% 
  filter(Enhancer %in% Mono_TPM_Closed_Active$enhancer)

#GO
GOTPM_Mono <- hyperGTest(new("GOHyperGParams",
                                    geneIds = Mono_TPM_CA$entrez,
                                    universeGeneIds = Mono_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_Mono  %>% summary() %>% head(10)

write.csv(summary(GOTPM_Mono), file = "GOTPM_Mono.csv")

#Bcell

Bcell_TPM <- read_tsv(files[7],
                      col_names = c("EnsemblID", "Gene", "TPM"))

#add ENTREZ ids
Bcell_TPM$entrez <- mapIds(
  org.Mm.eg.db,
  keys = Bcell_TPM$EnsemblID,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#extract top quartile
Bcell_TPM_quantiles <- quantile(Bcell_TPM$TPM)

Bcell_TPM <- subset(Bcell_TPM, TPM > Bcell_TPM_quantiles[4])

#extract only genes found in the list of enhancers
Bcell_TPM <-  subset(Bcell_TPM, Gene %in% c(Enhancers$Gene))

#extract only enhancers that are associated with the named genes
Bcell_TPM_enhancers <- subset(Enhancers, Gene %in% c(Bcell_TPM$Gene))
Bcell_TPM_enhancers <- merge(Bcell_TPM_enhancers, Bcell_TPM, by = "Gene", all.y = TRUE)

#extract only enhancers that are closed/unacetylated (status = 0)
Bcell_TPM_Closed_Inactive <- subset(Enhancer_Status, Bcell_status == 0)
Bcell_TPM_CI <- Bcell_TPM_enhancers %>% 
  filter(Enhancer %in% Bcell_TPM_Closed_Inactive$enhancer)


#extract only enhancers that are closed/acetylated (status = 1)
Bcell_TPM_Closed_Active <- subset(Enhancer_Status, Bcell_status == 1)
Bcell_TPM_CA <- Bcell_TPM_enhancers %>% 
  filter(Enhancer %in% Bcell_TPM_Closed_Active$enhancer)

#GO
GOTPM_Bcell <- hyperGTest(new("GOHyperGParams",
                                    geneIds = Bcell_TPM_CA$entrez,
                                    universeGeneIds = Bcell_TPM_CI$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

GOTPM_Bcell %>% summary() %>% head(10)

write.csv(summary(GOTPM_Bcell), file = "GOTPM_Bcell.csv")




