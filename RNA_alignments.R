# Kitty Martens
# RNA alignments in deseq2 - following STAR
# 22 Nov 2021

#load necessary packages
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(GO.db)
library(GOstats)

dir <- "HTSeq"

# assign file path to a variable
sample_metadata <- read_csv(file.path(dir, "samples.csv"))
sample_metadata

# assign HTSeq files to a variable
files <- file.path(dir, sample_metadata$sample)
files

# create a new data table with the variables sampleName, fileName and cell_type
sample_df <- data.frame(sampleName = sample_metadata$column1,
                        fileName = files,
                        condition = sample_metadata$cell_type,
                        run = sample_metadata$run)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sample_df,
                                       design = ~ condition)

# set control condition using the relevel function
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "MPP")

# run DESeq on ddsHTSeq
dds <- DESeq(ddsHTSeq)

# visualize clustering with a PCA plot
rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")

# visualize  using a distance matrix
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)

rld@colData$condition
rld@colData$run


colnames(sample_dist_matrix) <- paste(rld@colData$condition, rld@colData$run, sep = " _ " )
rownames(sample_dist_matrix) <- paste(rld@colData$condition, rld@colData$run, sep = " _ " )
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

# pulling out the names of genes identified as more or less abundant
all.names <-  results(dds)

res_Bcell_vs_MPP <- results(dds, name = "condition_Bcell_vs_MPP")

res_CMP_vs_MPP <- results(dds, name = "condition_CMP_vs_MPP")

res_GMP_vs_MPP <- results(dds, name = "condition_GMP_vs_MPP")

res_Gran_vs_MPP <- results(dds, name = "condition_Gran_vs_MPP")

res_Mono_vs_MPP <- results(dds, name = "condition_Mono_vs_MPP")



# subset out only differentially abundant genes
significant_res_Bcell <- subset(res_Bcell_vs_MPP, padj < 0.05)
write.csv(as.data.frame(significant_res_Bcell), file = "Bcell_vs_MPP_05.csv")

significant_res_CMP <- subset(res_CMP_vs_MPP, padj < 0.05)
write.csv(as.data.frame(significant_res_CMP), file = "CMP_vs_MPP_05.csv")

significant_res_GMP <- subset(res_GMP_vs_MPP, padj < 0.05)
write.csv(as.data.frame(significant_res_GMP), file = "GMP_vs_MPP_05.csv")

significant_res_Gran <- subset(res_Gran_vs_MPP, padj < 0.05)
write.csv(as.data.frame(significant_res_Gran), file = "Gran_vs_MPP_05.csv")

significant_res_Mono <- subset(res_Mono_vs_MPP, padj < 0.05)
write.csv(as.data.frame(significant_res_Mono), file = "Mono_vs_MPP_05.csv")


all.names$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(all.names),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

all.names$entrez <- mapIds(
  org.Mm.eg.db,
  keys = rownames(all.names),
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

all_genes <- all.names %>% 
  as.data.frame() %>% 
  pull(entrez) %>% 
  unique()

# GO term analysis
GOanalysis <- function (val) {

# add Ensembl gene ids as a new row
val$symbol <- mapIds(
  org.Mm.eg.db,
  keys = rownames(val),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

val$entrez <- mapIds(
  org.Mm.eg.db,
  keys = rownames(val),
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

all_genes <- val %>% 
  as.data.frame() %>% 
  pull(entrez) %>% 
  unique()

genes_upregulated <- val %>% 
  as.data.frame() %>% 
  filter(log2FoldChange > 4) %>% 
  pull(entrez) %>% 
  unique()

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated,
                                    universeGeneIds = all_genes,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.01,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = paste("GOup/", val, ".csv"))

}

GOanalysis(significant_res_Bcell)
write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Bcell.csv")

GOanalysis(significant_res_CMP)
write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_CMP.csv")
          
GOanalysis(significant_res_GMP)
write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_GMP.csv")
          
GOanalysis(significant_res_Gran)
write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Gran.csv")

GOanalysis(significant_res_Mono)
write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Mono.csv")


