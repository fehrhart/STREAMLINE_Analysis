#Part 1 DESeq2 - from raw counts to differentially expressed genes
# This script is to analyse RNAseq data from the STREAMLINE using DESeq2 package. 

#clean workspace
rm(list=ls())

#Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("tidyr")
BiocManager::install("dplyr")
BiocManager::install("stringr")
BiocManager::install("ggfortify")

#load required libraries
library(DESeq2)
library(ggplot2)
library(ggfortify)
library(tidyr)
library(dplyr)
library(stringr)
library(dbplyr)
library(biomaRt)

#Set the working directory by copy-pasting your path
setwd("C:/Users/friederike.ehrhart/Documents/STREAMLINE/DESeq analysis")

#load raw count data and create some overview plots.
data = read.csv(file = '8vs11.csv', row.names = 1, header = TRUE)

#Some raw count data visualisation with boxplot and PCA
boxplot(log(data))
PCA <- prcomp(t(data))
autoplot(PCA, label = TRUE, label.size = 3)

#Load the metadata file. 
metadata <- read.csv(file = 'm8vs11.csv', header = TRUE)

#Create the "dds" object from count data and metadata. The experimental design compares via diagnosis (COS vs control). This is a design if you have only one variable to compare experimental groups.
dds <- DESeqDataSetFromMatrix(countData = data,
                             colData = metadata,
                             design = ~ Disease, tidy = FALSE)

#run the DESeq2 function. This will take a while (up to 10 min are normal).
dds <- DESeq(dds)

#result output
resultsNames(dds)
DEG <- results(dds)
write.csv(DEG, file = "DEG_N8vs11_2.csv")

#Add HGNC symbols
#Select a database and a dataset to define your "ensembl mart" object
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
Ensembl_hsa_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
colnames(Ensembl_hsa_genes) <-c('X','HGNC')

DEG_N8 = read.csv("DEG_N8vs11_2.csv", header = TRUE)

#Merge with the DEG.csv dataset . Note that not all Ensembl IDs will get a HGNC ID.
DEG_N8_HGNC <- merge(DEG_N8, Ensembl_hsa_genes, by="X")

#export the complete result file to CSV for later analysis
write.csv(DEG_N8_HGNC, file = "DEG_N8vs11_2.csv")

#Volcano plot
par(mfrow=c(1,1))
with(DEG, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
with(subset(DEG, padj<0.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DEG, padj<0.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#compare normalised counts vs. raw counts 
dds_n <- estimateSizeFactors(dds); 
dds_n <- counts(dds_n, normalized=TRUE)
boxplot(log(dds_n))

PCA <- prcomp(t(dds_n))
autoplot(PCA, label = TRUE, label.size = 3)
