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

# assign closed acetylated genes identified in GREAT to a variable
file <- read_tsv("closest_gene_expression_narm.tsv")
file

# create a new data frame to contain the enhancer files
sample_df <- data.frame(file)

#add ENSEMBL id's
sample_df$ensembl <- mapIds(
  org.Mm.eg.db,
  keys = sample_df$gene,
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
)

#add ENTREZ ids
sample_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = sample_df$ensembl,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

#create a list of all unique genes
all_genes_enhancers <- sample_df %>% 
  as.data.frame() %>% 
  pull(entrez) %>% 
  unique()



# subset out only differentially abundant genes
significant_res_Bcell <- subset(sample_df, bcell_padj < 0.05)
write_csv(as.data.frame(significant_res_Bcell), file = "Bcell_vs_MPP_05_Enhancer.csv")

significant_res_CMP <- subset(sample_df, cmp_padj< 0.05)
write_csv(as.data.frame(significant_res_CMP), file = "CMP_vs_MPP_05_Enhancer.csv")

significant_res_GMP <- subset(sample_df, gmp_padj< 0.05)
write_csv(as.data.frame(significant_res_GMP), file = "GMP_vs_MPP_05_Enhancer.csv")

significant_res_Gran <- subset(sample_df, gran_padj < 0.05)
write_csv(as.data.frame(significant_res_Gran), file = "Gran_vs_MPP_05_Enhancer.csv")

significant_res_Mono <- subset(sample_df, mono_padj < 0.05)
write_csv(as.data.frame(significant_res_Mono), file = "Mono_vs_MPP_05_Enhancer.csv")


#subset by log2fold change for each cell type
genes_upregulated_bcell <- significant_res_Bcell %>% 
  as.data.frame() %>% 
  filter(bcell_log2fold > 1) %>%  
  pull(entrez) %>% 
  unique()

write_csv(as.data.frame(genes_upregulated_bcell), file = "Bcell_vs_MPP_05_1_Enhancer.csv")

genes_upregulated_cmp <- significant_res_CMP %>% 
  as.data.frame() %>% 
  filter(cmp_log2fold > 1) %>% 
  pull(entrez) %>% 
  unique()

write_csv(as.data.frame(genes_upregulated_cmp), file = "CMP_vs_MPP_05_1_Enhancer.csv")

genes_upregulated_gmp <- significant_res_GMP %>% 
  as.data.frame() %>% 
  filter(gmp_log2fold > 1) %>% 
  pull(entrez) %>% 
  unique()

write_csv(as.data.frame(genes_upregulated_gmp), file = "GMP_vs_MPP_05_1_Enhancer.csv")

genes_upregulated_gran <- significant_res_Gran %>% 
  as.data.frame() %>% 
  filter(gran_log2fold > 1) %>% 
  pull(entrez) %>% 
  unique()

write_csv(as.data.frame(genes_upregulated_gran), file = "Gran_vs_MPP_05_1_Enhancer.csv")

genes_upregulated_mono <- significant_res_Mono %>% 
  as.data.frame() %>% 
  filter(mono_log2fold > 1) %>% 
  pull(entrez) %>% 
  unique()

write_csv(as.data.frame(genes_upregulated_mono), file = "Mono_vs_MPP_05_4_Enhancer.csv")

#create hypergtest object for each cell type and print to file
go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated_bcell,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)


write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Bcell.csv")

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated_cmp,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)


write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_CMP.csv")
          
go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated_gmp,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_GMP.csv")

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated_gran,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Gran.csv")

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = genes_upregulated_mono,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = "GOupsignificant_res_Mono.csv")

go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = sample_df$entrez,
                                    universeGeneIds = all_genes_enhancers,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = "GOupsignificant_all_Enhancers.csv")


#I did this in a very bad way - the all.names variable comes from the script RNA_alignments.R
go_bp_upregulated <- hyperGTest(new("GOHyperGParams",
                                    geneIds = all_genes_enhancers_enhancers[1-228],
                                    universeGeneIds = all.names$entrez,
                                    annotation = "org.Mm.eg.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    conditional = FALSE,
                                    testDirection = "over"))

go_bp_upregulated %>% summary() %>% head(10)

write.csv(summary(go_bp_upregulated), file = "GOupsignificant_all_Enhancers.csv")

