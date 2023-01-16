# LAB REPORT 

### EXERCISE 8.1 - Download packages needed ############################################

# R version 4.2.2

library(BiocManager)
# Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.2
# (2022-10-31 ucrt)

library(tximport)
# Version 1.26.0

library(DESeq2)
# Version 1.38.1

library(biomaRt)
# Version 2.54.0

library(pheatmap)
# Version 1.0.12

library(tidyverse)
# tidyverse 1.3.2
# ??? ggplot2 3.4.0      ??? purrr   0.3.5 
# ??? tibble  3.1.8      ??? dplyr   1.0.10
# ??? tidyr   1.2.1      ??? stringr 1.5.0 
# ??? readr   2.1.3      ??? forcats 0.5.2 

library(magrittr)
# Version 2.0.3

library(RColorBrewer)
# Version 1.1-3

#install.packages("ggrepel")
library(ggrepel)
# Version 0.9.2

#install.packages("ggpubr")
library(ggpubr)
# Version 0.5.0

### EXERCISE 8.2 - Read data in for DESeq2 processing ############################################

#Set correct working directory
setwd("C:/Users/fongf/OneDrive - Newcastle University/Newcastle MRes/MMB8052 - Bioinformatics/Assignment - Lab report")

# Download sample table with (3 conditions x 4 replicates) data
sample_table <- read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')

# Extract names of sample
files <- pull(sample_table, Run)

# Edit samples with file path to extract data files from "count" folder
# "count" folder - result of RNA quantification using Salmon
files <- paste0('counts/', files, '/quant.sf')

# names () - Functions to get or set the names of an object.
names(files) <- pull(sample_table, Run)

# Gene map - to map genes using Ensembl IDs
# Contains 2 columns: ensmust,ensmusg

gene_map <- read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')

# Import & summarize transcript-level estimates for transcript & gene-level analysis
txi <- tximport(files, 
                type='salmon',
                tx2gene=gene_map,
                ignoreTxVersion=TRUE)
# List of 4
# Include "abundance", "counts", "lengths", "countsFromAbundance" (no)

### EXERCISE 8.3 - DESeq2 processing ############################################

# DESeq2 - normalisation & stabilisation for processing Salmon count output

# Prepare data - convert tximport into DESeq2 input 
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)

# OPTIONAL - PRE-FILTERING: Remove rows with low count numbers
# keep <- rowSums(counts(dds)) >= 3
# dds <- dds[keep,]

# SIZE FACTOR NORMALISATION
# Normalise experiments by size factors - "median of ratios" method
# log eliminates all genes transcribed in only one sample type
# Geometric mean smooth over outlier reads
# Median downplays genes that soak up lots of reads
# Puts emphasis on more moderately expressed genes
dds <- estimateSizeFactors(dds)
# combines size factor normalisation with average transcript length normalisation

# DISPERSION
# genewise measure of the variance of the count data - model biological noise
# Variance typically decreases with increasing read count
# Maximum likelihood estimate (MLE)
# Shrink gene-wise dispersion estimates towards value predicted by curve
# Limits detection of false positives
dds <- estimateDispersions(dds)

# WALD TEST FOR NEGATIVE BINOMIAL DISTRIBUTION
# RNA-seq raw count data 'naturally' follows a negative binomial distribution (Poisson-like)
# Wald test (a.k.a. Wald Chi-Squared Test) 
# find out if explanatory variables in a model are significant
# Null hypothesis for the test is: some parameter = some value
dds <- nbinomWaldTest(dds)

# EXPORT DISPERSION PLOT (AS TIFF IMAGE)

# tiff("Dispersion_estimate_plot_01.tiff", height = 750, width = 1000)

options(scipen=999)
plotDispEsts(dds,
             xlab = "Mean of normalised counts",
             ylab = "Dispersion",
             legend = T,
             log = "xy",
             cex = 0.7)
options(scipen=0)

# dev.off() 

# ACCESSOR METHODS FOR DESEQ2

# View raw counts & normalised counts
counts(dds)
counts(dds, normalized = T)

# View sample table
col_dds <- colData(dds)

# View matrix of normalisation factors
normalizationFactors(dds)


### EXERCISE 8.4 - Data quality ####################################

# REGULARISED LOGARITHM TRANSFORMATION
rld <- rlog(dds)

# PCA PLOT
# Base plot with DESeq2 function
pca_data <- plotPCA(rld, intgroup ='Group', returnData = TRUE)

# Modify and export PCA plot

#tiff("PCA_plot_03.tiff", height = 750, width = 1000, res = 100)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size =3) +
  xlab("PC1: 60% variance") +
  ylab("PC2: 18% variance") +
  ggtitle("PCA plot") +
  coord_fixed() +
  scale_color_manual(values = c(
    "Allo24h" = "tomato1", "Allo2h" = "yellowgreen", Naive = "steelblue3")) +
  theme_bw()

#dev.off()

# SAMPLE DISTANCE ANALYSIS
# Measure Euclidean distance 

# To set alternate row names for distance plot from example code
# Made copy to retain original data
dds_alt <- dds

# Reset rownames from SRR number to sample names
rownames(colData(dds_alt)) <- colData(dds_alt)$Sample_Name

rld_alt <- rlog(dds_alt)

# Convert to sample distance matrix
sample_distance_alt <- dist(t(assay(rld_alt)), method='euclidian')

sample_distance_matrix_alt <- as.matrix(sample_distance_alt)

# Annotations for distance heatmap
heatmap_annotation_alt <- data.frame(group=colData(dds_alt)[,c('Group')],
                                     row.names=rownames(colData(dds_alt))) %>% 
  rename("Group" = group)

# Specify colours for annotation to match PCA plot
annot_colors <- list(Group=c(Allo24h="tomato1",Allo2h="yellowgreen", Naive="steelblue3"))


# CONSTRUCT AND EXPORT DISTANCE HEATMAP

#tiff("Distance_heatmap_02.tiff", height = 1000, width = 1000, res = 100)

pheatmap(sample_distance_matrix_alt,
         clustering_distance_rows=sample_distance_alt,
         clustering_distance_cols=sample_distance_alt,
         annotation_col = heatmap_annotation_alt,
         annotation_row = heatmap_annotation_alt,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="YlGnBu")))(100),
         border_color = NA,
         annotation_colors = annot_colors,
         annotation_names_row = F,
         annotation_names_col = F,
         fontsize = 10) 

#dev.off()


### EXERCISE 8.5 - Detect DEGs #################################################

# COMPARISON BTWN ALLO24H vs NAIVE

# Results table
results_table <- results(dds, contrast= c('Group', 'Allo24h','Naive'))

summary(results_table)
# out of 33071 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3682, 11%
# LFC < 0 (down)     : 3030, 9.2%
# outliers [1]       : 48, 0.15%
# low counts [2]     : 13589, 41%
# (mean count < 3)

# LFC = shrunken log2 fold changes

# Convert to tibble & filter out rows with NA anywhere 
results_tibble <- as_tibble(results_table, rownames='ensembl_gene_id')
# n = 55858

filtered_results <- filter(results_tibble,
                           complete.cases(results_tibble))
# n = 19434

# Add columns for -log10(padj), absolute log2 foldchange, and significance
filtered_results <- mutate(filtered_results, logPVal = -log10(padj))

filtered_results_01 <- filtered_results %>% 
  mutate("Abs_log2foldchange" = abs(log2FoldChange)) %>% 
  mutate("Significant" = case_when(padj<0.05 & Abs_log2foldchange>1 ~ "Both",
                                   padj<0.05 ~ "p-value only",
                                   padj>=0.05 ~ "Neither")) 

# COMPARISON BTWN ALLO2H vs NAIVE

# Results table
results_table_2h <- results(dds, contrast= c('Group', 'Allo2h','Naive'))

summary(results_table_2h)
# out of 33071 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2794, 8.4%
# LFC < 0 (down)     : 2785, 8.4%
# outliers [1]       : 48, 0.15%
# low counts [2]     : 12358, 37%

# Convert to tibble & filter out rows with NA anywhere 
results_tibble_2h <- as_tibble(results_table_2h, rownames='ensembl_gene_id')
# n = 55858

filtered_results_2h <- filter(results_tibble_2h,
                              complete.cases(results_tibble_2h))
# n = 20665

# Add columns for -log10(padj), absolute log2 foldchange, and significance
filtered_results_2h_2 <- filtered_results_2h %>% 
  mutate(logPVal = -log10(padj)) %>% 
  mutate("Abs_log2foldchange" = abs(log2FoldChange)) %>% 
  mutate("Significant" = case_when(padj<0.05 & Abs_log2foldchange>1 ~ "Both",
                                   padj<0.05 ~ "p-value only",
                                   padj>=0.05 ~ "Neither")) 


### EXERCISE 8.6 - Annotating DEG results #########################

# Biomart - system for extracting tables of data from Ensembl

ensembl108 <- useEnsembl(biomart="ensembl", version=108)

ensembl108 <- useDataset("mmusculus_gene_ensembl", mart=ensembl108)

annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                 'start_position', 'end_position',
                                 'strand', 'gene_biotype',
                                 'external_gene_name',
                                 'description'),
                    filters = 'ensembl_gene_id', values =
                      filtered_results$ensembl_gene_id,
                    mart = ensembl108)


# Add annotations for genes differentially expressed in Allo24h vs Naive

annot_results_01 <- left_join(filtered_results_01, annotation)

annot_results_01 <- arrange(annot_results_01, padj) 
# reorder genes by smallest adjusted p-value to largest 

# Export annotated table
# write.csv(annot_results_01, "Allo24h_vs_naive_filtered_gene_list.csv", row.names = F)

# Filter for those fulfilling criteria for foldchange and padj
degs_01 <- annot_results_01 %>% 
  filter(Significant == "Both")
# n = 1510


# Add annotations for genes differentially expressed in Allo2h vs Naive

annot_results_02 <- left_join(filtered_results_2h_2, annotation)

annot_results_02 <- arrange(annot_results_02, padj)
# reorder genes by smallest adjusted p-value to largest 

# Export annotated table
# write.csv(annot_results_02, "Allo2h_vs_naive_filtered_gene_list.csv", row.names = F)

# Filter for those fulfilling criteria for foldchange and padj
degs_02 <- annot_results_02 %>% 
  filter(Significant == "Both")
# n = 1268


# VOLCANO PLOT WITH ANNOTATIONS

# Allo24h vs naive

plot_24h <- 
  ggplot(annot_results_01, aes(x=log2FoldChange, y=logPVal, color = Significant)) +
  geom_point(size =0.5, show.legend = F) + 
  geom_text_repel(aes(label=case_when(
    ((log2FoldChange < -7 | 
        log2FoldChange > 7) &
       Significant == "Both")|
      -log10(padj) >100 ~ as.character(external_gene_name))),
    hjust=0.5,vjust=0, 
    show.legend = F,
    color = "maroon",
    max.overlaps = 20) +
  theme_bw()+
  scale_color_manual(values = c("Both" = "red", 
                                "p-value only" = "grey48", 
                                "Neither" = "black")) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "steelblue") +
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", 
             color = "steelblue") +
  ggtitle("Allogeneic 24 hours vs Naive") +
  xlab("log2 Fold Change") +
  ylab("-log10 (padj)") +
  theme(plot.title = 
          element_text(color="black", size=10, face="bold",hjust = 0.5),
        text = element_text(size = 10, color = "black")) +
  xlim(-21, 25) +
  ylim(0, 200)

# Allo2h vs naive

plot_2h <- 
  ggplot(annot_results_02, aes(x=log2FoldChange, y=logPVal, color = Significant)) +
  geom_point(size =0.5, show.legend = F) + 
  geom_text_repel(aes(label=case_when(
    ((log2FoldChange < -7.5 | 
        log2FoldChange > 10) &
       Significant == "Both") |
      -log10(padj) >50 ~ as.character(external_gene_name))),
    hjust=0.5,vjust=0, 
    #position=position_jitter(width=2,height=1),
    show.legend = F,
    color = "maroon",
    max.overlaps = 20) +
  theme_bw()+
  scale_color_manual(values = c("Both" = "red", 
                                "p-value only" = "grey48", 
                                "Neither" = "black")) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "steelblue") +
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", 
             color = "steelblue") +
  ggtitle("Allogeneic 2 hours vs Naive") +
  xlab("log2 Fold Change") +
  ylab("-log10 (padj)") +
  theme(plot.title = 
          element_text(color="black", size=10, face="bold",hjust = 0.5),
        text = element_text(size = 10, color = "black")) +
  xlim(-21, 25) +
  ylim(0, 200)

# Export figure with both volcano plots by ggarrange()
# tiff("Volcano_plots_03.tiff",height = 750, width = 1500, res = 100)

ggarrange(plot_24h, plot_2h,
          ncol = 2, nrow = 1)


# dev.off()

### EXERCISE 9.1 - Heatmap by foldchange for DEG #############################

# ALLO24H vs NAIVE

#Use filtered results - produce heatmap for top 250 DEG for Allo24h vs naive
degs_01_1 <- degs_01 %>% 
  arrange(desc(Abs_log2foldchange))

# Export final DEG list
# write.csv(degs_01_1, "Allo24h_vs_naive_DEG_list.csv", row.names = F)

# Filter only top 250 DEGs for heatmap 
# - gene names annotated for rows
# - sample names annotated for columns

degs_01_2 <- degs_01_1 %>% 
  slice_head(n = 250)

heatmap_genes <- degs_01_2$ensembl_gene_id

heatmap_data_alt <- assay(rld_alt)[heatmap_genes,]

rownames(heatmap_data_alt) <- degs_01_2$external_gene_name[1:250]

column_annotation_alt <- data.frame(group = colData(dds_alt)[,c('Group')],
                                    row.names = rownames(colData(dds_alt)))
# Use dds_alt with rownames changed to sample names instead of SRR number

annot_colors_2 <- list(group=c(Allo24h="tomato1",Allo2h="yellowgreen", Naive="steelblue3"))


# Export top 250 DEG heatmap for Allo24h vs naive

# tiff("DEG_heatmap_03.tiff", width = 1700, height = 2500, res = 180)

pheatmap(heatmap_data_alt,
         annotation_col = column_annotation_alt,
         annotation_colors = annot_colors_2,
         scale = "row", 
         fontsize = 5,
         cutree_rows = 3,
         annotation_names_col = F)

# dev.off()



# ALLO2H vs NAIVE

# *repeated as above

#Use filtered results - produce heatmap for top 250 DEG for Allo2h vs naive
degs_02_1 <- degs_02 %>% 
  arrange(desc(Abs_log2foldchange))

# Export final DEG list
#write.table(degs_02_1, "Allo2h_vs_naive_DEG_list", row.names = F)

# Filter only top 250 DEGs for heatmap 
degs_02_2 <- degs_02_1 %>% 
  slice_head(n = 250)

heatmap_genes <- degs_02_2$ensembl_gene_id

heatmap_data_alt <- assay(rld_alt)[heatmap_genes,]
rownames(heatmap_data_alt) <- degs_02_2$external_gene_name[1:250]
column_annotation_alt <- data.frame(group = colData(dds_alt)[,c('Group')],
                                    row.names = rownames(colData(dds_alt)))

# Export top 250 DEG heatmap for Allo2h vs naive

# tiff("DEG_heatmap_04.tiff", width = 1700, height = 2500, res = 180)

pheatmap(heatmap_data_alt,
         annotation_col = column_annotation_alt,
         annotation_colors = annot_colors_2,
         scale = "row", 
         fontsize = 5,
         cutree_rows = 3,
         annotation_names_col = F)

# dev.off()
