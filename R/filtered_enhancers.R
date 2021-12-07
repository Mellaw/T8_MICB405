#Filtered enhancers 

# libraries  --------------------------------------------------------------
library(tidyverse)
library(dplyr)

enhancers <-  read.table(file = 'rna_align/enhancers_status.tsv',
                         sep = '\t', header = TRUE)
# Subset enhancers --------------------------------------------------------

mpp_acet_closed <- enhancers %>% 
  filter(MPP_atac_significant == 0 &
           MPP_ac_significant == 1)

mpp_acet_closed_great <-  mpp_acet_closed %>% 
  select(chr, start, end, enhancer)

#Export enhancers as .bed for visualization in GREAT 
write_tsv(mpp_acet_closed_great, "./rna_align/mpp_acet_closed_great.bed", col_names = FALSE)


# Filter GREAT ---------------------------------------------------------

great <- read.table(file = 'rna_align/GREAT/single_closest_gene.txt',
                    sep = '\t', header = FALSE)
colnames(great)<- c("enhancer", "gene_dist")

gene_dist <- great %>% 
  mutate(gene_ls = strsplit(gene_dist, split = " ",fixed =TRUE))

gene_dist_ls <- gene_dist$gene_ls

index <- 1:388
gene_nm <- c()
bp_dist <- c()

for (i in index) {
  gene_nm <- append(gene_nm, unlist(gene_dist_ls[i])[1])
}

for (i in index) {
  bp_dist <- append(bp_dist, unlist(gene_dist_ls[i])[2])
}

gene_dist$gene = gene_nm
gene_dist$bp_dist = bp_dist
gene_dist <- gene_dist %>% 
  select(-c("gene_ls"))

# Filter lineages ---------------------------------------------------------
bcell <- read.table(file = 'rna_align/RNA_gene_counts_all/Bcell_vs_MPP_all.txt'
                    , sep = '\t', header = TRUE)
cmp <- read.table(file = 'rna_align/RNA_gene_counts_all/CMP_vs_MPP_all.txt'
                  , sep = '\t', header = TRUE)
gmp <- read.table(file = 'rna_align/RNA_gene_counts_all/GMP_vs_MPP_all.txt'
                  , sep = '\t', header = TRUE)
gran <- read.table(file = 'rna_align/RNA_gene_counts_all/Gran_vs_MPP_all.txt'
                   , sep = '\t', header = TRUE)


colnames(bcell) <- c("gene", "baseMean", "bcell_log2fold", "lfcSE", "stat",
                    "pvalue", "bcell_padj")
colnames(cmp) <- c("gene", "baseMean", "cmp_log2fold", "lfcSE", "stat",
                   "pvalue", "cmp_padj")
colnames(gmp) <- c("gene", "baseMean", "gmp_log2fold", "lfcSE", "stat",
                   "pvalue", "gmp_padj")
colnames(gran) <- c("gene", "baseMean", "gran_log2fold", "lfcSE", "stat",
                    "pvalue", "gran_padj")

#Filter cell lineages for desired columns

bcell_sub <- bcell %>% 
  filter(gene %in% gene_nm) %>% 
  select(gene, bcell_log2fold, bcell_padj)

cmp_sub <- cmp %>% 
  filter(gene %in% gene_nm) %>% 
  select(gene, cmp_log2fold, cmp_padj)

gmp_sub <- gmp %>% 
  filter(gene %in% gene_nm) %>% 
  select(gene, gmp_log2fold, gmp_padj)

gran_sub <- gran %>% 
  filter(gene %in% gene_nm) %>% 
  select(gene, gran_log2fold, gran_padj)

#Merge cell linage information with gene & enhancer table 
gene_dist = merge(gene_dist, bcell_sub, by = "gene", all.x = TRUE)
gene_dist = merge(gene_dist, cmp_sub, by = "gene", all.x = TRUE)
gene_dist = merge(gene_dist, gmp_sub, by = "gene", all.x = TRUE)
gene_dist = merge(gene_dist, gran_sub, by = "gene", all.x = TRUE)

gene_dist <- gene_dist %>% 
  select(-c(gene_dist))

write_tsv(gene_dist, "./rna_align/closest_gene_expression.tsv", col_names = TRUE)

summary_narm <- na.omit(gene_dist)

write_tsv(summary_narm, "./rna_align/closest_gene_expression_narm.tsv", col_names = TRUE)

# enhancers that become open acet -----------------------------------------

cmp_acet_open <- enhancers %>% 
  filter(CMP_atac_significant == 1 &
           CMP_ac_significant == 1)

gmp_acet_open <- enhancers %>% 
  filter(GMP_atac_significant == 1 &
           GMP_ac_significant == 1)

gran_acet_open <- enhancers %>% 
  filter(Gran_atac_significant == 1 &
           Gran_ac_significant == 1)

Bcell_acet_open <- enhancers %>% 
  filter(Bcell_atac_significant == 1 &
           Bcell_ac_significant == 1)

#Add T/F for each enhancer in summary_narm if enhancer becomes open/acet for
#each cell type

open_enhancer <- summary_narm %>% 
  mutate(cmp_open = ifelse(enhancer %in% cmp_acet_open$enhancer, TRUE, FALSE)) %>% 
  mutate(gmp_open = ifelse(enhancer %in% gmp_acet_open$enhancer, TRUE, FALSE)) %>% 
  mutate(gran_open = ifelse(enhancer %in% gran_acet_open$enhancer, TRUE, FALSE)) %>% 
  mutate(bcell_open = ifelse(enhancer %in% Bcell_acet_open$enhancer, TRUE, FALSE))

write_tsv(open_enhancer, "./rna_align/closest_gene_open_ac_narm.tsv", col_names = TRUE)

#count how many True values 
open_enhancer %>% summarize(n = sum(cmp_open  == 1))
open_enhancer %>% summarize(n = sum(gmp_open  == 1))
open_enhancer %>% summarize(n = sum(gran_open  == 1))
open_enhancer %>% summarize(n = sum(bcell_open  == 1))
                          