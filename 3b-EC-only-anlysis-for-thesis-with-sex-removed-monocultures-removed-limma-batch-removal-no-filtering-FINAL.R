# Load in the data
library(tidyverse)

load("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Preprocessing_of_all_nosex.RData")

# EC only
eccountdata <- allcountdata[,c(4,5,6,10,11,12,14,20,21,23,25,30,31,33,35,40,41,43,45,50,51,53,55,60,61,63,65,70,71,
                               73,75)]
ecsampleinfo <- allsampleinfo[c(4,5,6,10,11,12,14,20,21,23,25,30,31,33,35,40,41,43,45,50,51,53,55,60,61,63,65,70,71,
                                73,75),]

# remove EC2016_BX064_EC_Res # 3741
eccountdata <- eccountdata[,c(-6)]
ecsampleinfo <- ecsampleinfo[-c(6),]

# remove the monocultures
eccountdata <- eccountdata[,-c(4, 5, 7, 11, 15,19, 23, 27)]
ecsampleinfo <- ecsampleinfo[-c(4, 5, 7, 11, 15,19, 23, 27),]

# load packages 
library(DESeq2)
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)

# Vectors for colours
DiagnosisCol <- match(ecsampleinfo$Diagnosis, c("Res", "VeRA", "JRep")) + 1

# Data transformation 
vstcounts <- vst(eccountdata)
vst_pcDat <- prcomp(t(vstcounts))

# PCA to confirm data has been loaded in correctly
autoplot(vst_pcDat, 
         data = ecsampleinfo, 
         colour = 'EC_batch',
         shape = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) +
 geom_text_repel(aes(x=PC1, y=PC2, label=Sample_name), box.padding = 0.8)

autoplot(vst_pcDat, 
         data = ecsampleinfo, 
         colour = 'EC_Sex',
         shape = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3) 

# Run DESeq2 

# design suggesting diagnosis causes differences in groups
design <- as.formula(~ EC_batch + Diagnosis)

# what is the y intercept as standard? 
modelMatrix <- model.matrix(design, data = ecsampleinfo)
modelMatrix

# change y intercept to resolving
ecsampleinfo$Diagnosis <- factor(ecsampleinfo$Diagnosis, levels = c("Res", "VeRA", "JRep"))
modelMatrix <- model.matrix(design, data = ecsampleinfo)
modelMatrix

# To access the p values

# Running DESeq2 with filtering

dds <- DESeqDataSetFromMatrix(countData = eccountdata,
                              colData = ecsampleinfo,
                              design = design)

dds <- estimateSizeFactors(dds)

# Filter out the lower expressed genes 
# e.g. idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
# would filter out genes where there are less than 3 samples with normalized counts greater than or equal to 5
idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

# The results and histograms
JRep_vs_Res <- results(dds,contrast=c("Diagnosis", "JRep", "Res"), alpha=0.05)
hist(JRep_vs_Res$pvalue, breaks=20, col="grey" )

VeRA_vs_Res <- results(dds,contrast=c("Diagnosis", "VeRA", "Res"), alpha=0.05)
hist(VeRA_vs_Res$pvalue, breaks=20, col="grey" )

VeRA_vs_JRep <- results(dds,contrast=c("Diagnosis", "VeRA", "JRep"), alpha=0.05)
hist(VeRA_vs_JRep$pvalue, breaks=20, col="grey" )

# DESeq with no filtering
dds <- DESeqDataSetFromMatrix(countData = eccountdata,
                              colData = ecsampleinfo,
                              design = design)

dds <- DESeq(dds)

#plot of the gene dispersion estimates
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/Gene_dispersion_estimate_plot.tiff",
     width = 500, height = 300)
plotDispEsts(dds)
dev.off()

# To save at this point
save(dds, eccountdata, ecsampleinfo, modelMatrix, vst, vst_pcDat, vstcounts, design, DiagnosisCol,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data sex genes and monoculture removed.RData")

load("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EEC data sex genes and monoculture removed.RData")

# More PCAs, with biplots, loadings etc. 
vst <- assay(vst(dds))

# To colour the graphs accoridng to Diagnosis, Joint etc. 
my_metadata <- data.frame(row.names = colnames(dds))
my_metadata$Diagnosis <- dds$Diagnosis
my_metadata$Joint <- dds$Joint
my_metadata$Fibroblast_Sex <- dds$Fibroblast_Sex
my_metadata$EC_Sex <- dds$EC_Sex
my_metadata$EC_donor <- dds$EC_batch
my_metadata$Age <- dds$age
my_metadata$Culture <- dds$Co_vs_mono
my_metadata$Fibroblast_ID <- dds$Fibroblast_ID

# To remove the lower 10% of variables based on variance
p <- pca(vst, metadata = my_metadata, removeVar = 0.1)

# To get gene symbol rather than ensembl gene names
library('org.Hs.eg.db')
newnames <- mapIds(org.Hs.eg.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="ENSEMBL")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by donor sex removed.tiff",
     width = 700, height = 400)
PCAtools::pairsplot(p, colby = 'EC_donor', pointSize=3)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC PCA coloured by EC sex sex removed.tiff",
     width = 700, height = 400)
PCAtools::pairsplot(p, colby = 'EC_Sex', pointSize=3)
dev.off()

# Biplots according to donor and sex of EC
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC biplot coloured by EC donor sex removed.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'EC_donor', pointSize=3.5, showLoadings = TRUE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'EC_donor', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC biplot coloured by EC sex sex removed.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, "PC3", "PC5", colby = 'EC_Sex', pointSize=3.5, showLoadings = TRUE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'EC sex', legendPosition = 'right')
dev.off()

# Loadings plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC original loadings plot.tiff",
     width = 500, height = 600)
PCAtools::plotloadings(p, labSize = 3)
dev.off()

# Effect of fibroblast in culture pairsplots
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by diagnosis sex removed.tiff",
     width = 750, height = 400)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=2.5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by Joint sex removed.tiff",
     width = 750, height = 400)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=2.5)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by fibroblast sex, sex removed.tiff",
     width = 750, height = 400)
PCAtools::pairsplot(p, colby = 'Fibroblast_Sex', pointSize=2.5)
dev.off()

# Biplots of the pairsplots, just to get legends
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC biplot coloured by diagnosis sex removed.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC biplot coloured by fibro sex, sex removed.tiff",
     width = 500, height = 400)
PCAtools::biplot(p, colby = 'Fibroblast_Sex', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast_Sex', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC biplot coloured by fibro joint, sex removed.tiff",
     width = 600, height = 400)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = TRUE, lab = NULL, 
                  sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
library(org.Hs.eg.db)

# To see what they're called
columns(org.Hs.eg.db)

topVarGenes <- head( order( rowVars(vst), decreasing=TRUE ), 50 )
subset <- vst[ topVarGenes, ]
newnames_subset <- mapIds(org.Hs.eg.db,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="ENSEMBL")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
fibroblast_sex_col <- as.data.frame(dds$Fibroblast_Sex)
culture_conditions <- as.data.frame(dds$Co_vs_mono)
EC_sex_col <- as.data.frame(dds$EC_Sex)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(5, 20, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC 50 variable genes heatmap no sex no monoculture.tiff",
     width = 600, height = 800)
Heatmap(subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex,
                                           Joint = dds$Joint, EC_donor = dds$EC_batch, EC_Sex = dds$EC_Sex, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"))))
dev.off()

save(dds, eccountdata, ecsampleinfo, modelMatrix, vst, vst_pcDat, vstcounts, design, DiagnosisCol, diagnosis_col, 
     EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, my_metadata, p, subset, design, newnames, newnames_subset,
     topVarGenes, col_fun,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data no filter sex removed monoculture removed.RData")

load("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data no filter sex removed monoculture removed.RData")

# Look at the results 
JRep_vs_Res <- results(dds,contrast=c("Diagnosis", "JRep", "Res"), alpha=0.05)
JRep_vs_Res

hist(JRep_vs_Res$pvalue, breaks=20, col="grey" )
sum(JRep_vs_Res$padj < 0.05, na.rm = TRUE)

VeRA_vs_Res <- results(dds,contrast=c("Diagnosis", "VeRA", "Res"), alpha=0.05)
VeRA_vs_Res

hist(VeRA_vs_Res$pvalue, breaks=20, col="grey" )
sum(VeRA_vs_Res$padj < 0.05, na.rm = TRUE)

VeRA_vs_JRep <- results(dds,contrast=c("Diagnosis", "VeRA", "JRep"), alpha=0.05)
VeRA_vs_JRep

hist(VeRA_vs_JRep$pvalue, breaks=20, col="grey" )
sum(VeRA_vs_JRep$padj < 0.05, na.rm = TRUE)

# Export results as csv

# Get ensembl names as rownames
# Although this says VeRA_vs_Res, it's the same gene list for all
ensembl <- as.vector(rownames(VeRA_vs_Res))

entrez <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL'))
symbol <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'SYMBOL', 'ENSEMBL'))

# JRep_vs_Res ----
# Get DEGs as a data frame
annotJRep_vs_Res <- as.data.frame(JRep_vs_Res) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotJRep_vs_Res <- cbind (annotJRep_vs_Res, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotJRep_vs_Res), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC_JRep_vs_Res_all.csv")

# export the significant results
JRep_vs_Res_sig <- as.data.frame(annotJRep_vs_Res)
JRep_vs_Res_sig <- JRep_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
JRep_vs_Res_sig <- JRep_vs_Res_sig[order(JRep_vs_Res_sig$log2FoldChange),]
write.csv(as.data.frame(JRep_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC_JRep_vs_Res_sig.csv")

# VeRA_vs_Res ----
# Get DEGs as a data frame
annotVeRA_vs_Res <- as.data.frame(VeRA_vs_Res) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Res <- cbind (annotVeRA_vs_Res, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_Res), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC_VeRA_vs_Res_all.csv")

# export the significant results
VeRA_vs_Res_sig <- as.data.frame(annotVeRA_vs_Res)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig[order(VeRA_vs_Res_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC_VeRA_vs_Res_sig.csv")

# VeRA vs Jrep ----
# Get DEGs as a data frame
annotVeRA_vs_JRep <- as.data.frame(VeRA_vs_JRep) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_JRep <- cbind (annotVeRA_vs_JRep, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_JRep), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC_VeRA_vs_JRep_all.csv")
# export the significant results
VeRA_vs_JRep_sig <- as.data.frame(annotVeRA_vs_JRep)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig[order(VeRA_vs_JRep_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_JRep_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC_VeRA_vs_JRep_sig.csv")

# To get gene names for heatmaps and volcano plots
library("EnsDb.Hsapiens.v79")

# To get volcano plots
library(EnhancedVolcano)

# VeRA_vs_Res
# Shrink LF change for volcano plot
VeRA_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Res'), 
                                res=VeRA_vs_Res, type = 'ashr')
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_Res_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_Res_shrunk) <- newnames
rownames(VeRA_vs_Res_shrunk)[rownames(VeRA_vs_Res_shrunk) == "ENSG00000206989"] <- "SNORD63"

# JRep_vs_Res
# Shrink LF change for volcano plot
JRep_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','JRep','Res'), 
                                res=JRep_vs_Res, type = 'ashr')
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(JRep_vs_Res_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(JRep_vs_Res_shrunk) <- newnames

rowname(JRep_vs_Res_shrunk) <- newnames
# Not all gene names are genes are mapped, so manually mapped this one 
rownames(JRep_vs_Res_shrunk)[rownames(JRep_vs_Res_shrunk) == "ENSG00000206989"] <- "SNORD63"

# VeRA_vs_JRep
# Shrink LF change for volcano plot
VeRA_vs_JRep_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','JRep'), 
                                 res=VeRA_vs_JRep, type = 'ashr')
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_JRep_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_JRep_shrunk) <- newnames
rownames(VeRA_vs_JRep_shrunk)[rownames(VeRA_vs_JRep_shrunk) == "ENSG00000206989"] <- "SNORD63"

# Volcano plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC volcanos padj.tiff",
     width = 1400, height = 400)
a <- EnhancedVolcano(JRep_vs_Res_shrunk,
                     lab = rownames(JRep_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'JRep vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
b <- EnhancedVolcano(VeRA_vs_Res_shrunk,
                     lab = rownames(VeRA_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
c <- EnhancedVolcano(VeRA_vs_JRep_shrunk,
                     lab = rownames(VeRA_vs_JRep_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs JRep',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
cowplot::plot_grid(a, b, c,
                   ncol = 3, nrow = 1)
dev.off()

# Heatmap of the Jrep vs Res genes (those between Jrep and VeRA are included in this)
Jrep_vs_Res_df <- as.data.frame(JRep_vs_Res)

Jrep_vs_Res_adj_sig <- Jrep_vs_Res_df %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 )

Jrep_vs_Res_adj_sig_ls <- rownames(Jrep_vs_Res_adj_sig)

# adding the DEG between VeRA and JRep
append(Jrep_vs_Res_adj_sig_ls, "ENSG00000206989")

Jrep_of_Res_subset_adj <- vst[Jrep_vs_Res_adj_sig_ls, ]

# New names for the vst object
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(vst),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)

rownames(vst) <- newnames

# new names for the subsetted genes
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(Jrep_vs_Res_adj_sig),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(Jrep_vs_Res_adj_sig) <- newnames

Jrep_vs_Res_adj_sig_ls <- rownames(Jrep_vs_Res_adj_sig)
Jrep_vs_Res_adj_sig_ls <- append(Jrep_vs_Res_adj_sig_ls, "SNORD63")
Jrep_of_Res_subset_adj <- vst[Jrep_vs_Res_adj_sig_ls, ]

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(breaks = seq(4, 10, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))  

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC sig DEGs heatmap.tiff",
     width = 500, height = 400)
Heatmap(Jrep_of_Res_subset_adj,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun, 
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis, EC_donor = dds$EC_batch, EC_sex = dds$EC_sex, col = list
                                                                                                (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "chocolate1", "JRep" = "darkslateblue"),
                                                                                                  EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                                                                               "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                                                                               "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                                                                                  EC_sex = c("F" = "pink", "M" = "Blue"))))
dev.off()

# genes from the PCA
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/Genes from heatmap.tiff", width = 1200, height = 900)
a <- plotCounts(dds, gene='ENSG00000180785', intgroup="Diagnosis", returnData = TRUE) # OR51E1
b <- plotCounts(dds, gene='ENSG00000236427', intgroup="Diagnosis", returnData = TRUE) # AL589986.2
c <- plotCounts(dds, gene='ENSG00000133055', intgroup="Diagnosis", returnData = TRUE) # MYBPH
d <- plotCounts(dds, gene='ENSG00000203435', intgroup="Diagnosis", returnData = TRUE) # E2F3P2
e <- plotCounts(dds, gene='ENSG00000165694', intgroup="Diagnosis", returnData = TRUE) # FRMD7
f <- plotCounts(dds, gene='ENSG00000105664', intgroup="Diagnosis", returnData = TRUE) # COMP
g <- plotCounts(dds, gene='ENSG00000181195', intgroup="Diagnosis", returnData = TRUE) # PENK
h <- plotCounts(dds, gene='ENSG00000248491', intgroup="Diagnosis", returnData = TRUE) # AC093772.1
i <- plotCounts(dds, gene='ENSG00000227158', intgroup="Diagnosis", returnData = TRUE) # AC073621.1
j <- plotCounts(dds, gene='ENSG00000206989', intgroup="Diagnosis", returnData = TRUE) # SNORD63
k <- plotCounts(dds, gene='ENSG00000223727', intgroup="Diagnosis", returnData = TRUE) # AC034195.1
l <- plotCounts(dds, gene='ENSG00000180875', intgroup="Diagnosis", returnData = TRUE) # GREM2
A <- ggplot(a, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("OR51E1")
B <- ggplot(b, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("AL589986.2")
c <- ggplot(c, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("MYBPH")
D <- ggplot(d, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("E2F3P2")
E <- ggplot(e, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("FRMD7")
FF <- ggplot(f, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("COMP")
G <- ggplot(g, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("PENK")
H <- ggplot(h, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("AC093772.1")
I <- ggplot(i, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("AC073621.1")
J <- ggplot(i, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("SNORD63")
K <- ggplot(i, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("AC034195.1")
L <- ggplot(i, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("GREM2")
cowplot::plot_grid(A, B, c, D, E, FF, G, H, I, J,K,L,
                   ncol = 3, nrow = 4)
dev.off()

# Make some files for GSEA
# To get volcano plots
library(EnhancedVolcano)

# VeRA_vs_Res
# Shrink LF change for volcano plot
VeRA_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','Res'), 
                                res=VeRA_vs_Res, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_Res_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_Res_shrunk) <- newnames
rownames(VeRA_vs_Res_shrunk)[rownames(VeRA_vs_Res_shrunk) == "ENSG00000206989"] <- "SNORD63"

# JRep_vs_Res
# Shrink LF change for volcano plot
JRep_vs_Res_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','JRep','Res'), 
                                res=JRep_vs_Res, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(JRep_vs_Res_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(JRep_vs_Res_shrunk) <- newnames

rowname(JRep_vs_Res_shrunk) <- newnames
# Not all gene names are genes are mapped, so manually mapped this one 
rownames(JRep_vs_Res_shrunk)[rownames(JRep_vs_Res_shrunk) == "ENSG00000206989"] <- "SNORD63"

# VeRA_vs_JRep
# Shrink LF change for volcano plot
VeRA_vs_JRep_shrunk <- lfcShrink(dds, contrast = c('Diagnosis','VeRA','JRep'), 
                                 res=VeRA_vs_JRep, type = 'ashr')

library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_JRep_shrunk),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_JRep_shrunk) <- newnames
rownames(VeRA_vs_JRep_shrunk)[rownames(VeRA_vs_JRep_shrunk) == "ENSG00000206989"] <- "SNORD63"

# Volcano plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC volcanos padj.tiff",
     width = 1400, height = 400)
a <- EnhancedVolcano(JRep_vs_Res_shrunk,
                     lab = rownames(JRep_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'JRep vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
b <- EnhancedVolcano(VeRA_vs_Res_shrunk,
                     lab = rownames(VeRA_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
c <- EnhancedVolcano(VeRA_vs_JRep_shrunk,
                     lab = rownames(VeRA_vs_JRep_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs JRep',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 5)
cowplot::plot_grid(a, b, c,
                   ncol = 3, nrow = 1)
dev.off()

# Files for GSEA
counts <- counts(dds)
write.csv(counts,file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC_dds_counts_table.csv")
write.csv(ecsampleinfo$Sample_name,file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC_dds_sample_names.csv")

# Save this before batch correction
save(dds, diagnosis_col, EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, modelMatrix, my_metadata, p, subset, vst, vst_pcDat,
     vstcounts, design, DiagnosisCol, newnames, newnames_subset, topVarGenes, col_fun, counts, annotJRep_vs_Res, annotVeRA_vs_Res, annotVeRA_vs_JRep,
     Jrep_of_Res_subset_adj, JRep_vs_Res, Jrep_vs_Res_adj_sig, Jrep_vs_Res_df, JRep_vs_Res_sig, VeRA_vs_JRep, VeRA_vs_JRep_sig, VeRA_vs_Res, 
     VeRA_vs_Res_sig, ensembl, entrez, Jrep_vs_Res_adj_sig_ls, newnames, newnames_subset, symbol, topVarGenes, col_fun,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data no filter sex removed before limma correction.RData")

# Do some batch removal with limma ----

# PCA of before batch correction
vsd <- vst(dds)
plotPCA(vsd, "Diagnosis")

design0 <- model.matrix(~ 0  + ecsampleinfo$Diagnosis, data=data.frame(assay(vsd)))

# Removing batch 
assay(vsd) <- limma::removeBatchEffect(assay(vsd),batch=as.vector(ecsampleinfo$EC_batch), design = design0)
plotPCA(vsd, "Diagnosis") 

# To remove the lower 10% of variables based on variance
p <- pca(assay(vsd), metadata = my_metadata, removeVar = 0.1)

# To get gene symbol rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

# Limma corrected biplot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC Biplot coloured by EC_DONOR sex removed batch corrected limma LABELLED .tiff",
     width = 600, height = 400)
PCAtools::biplot(p, colby = 'EC_donor', pointSize=3.5, showLoadings = FALSE, lab = NULL, encircle = FALSE,
                 sizeLoadingsNames = 5, colLegendTitle = 'EC_donor', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC Biplot coloured by diagnosis sex removed batch corrected limma LABELLED .tiff",
     width = 600, height = 400)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, encircle = TRUE,
                 sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')
dev.off()

# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
VST <- assay(vsd)
topVarGenes <- head( order( rowVars(VST), decreasing=TRUE ), 100 )
subset <- VST[ topVarGenes, ]

newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames)
rownames(subset) <- newnames_subset

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(0, 20, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC 100 variable genes heatmap no sex genes limma batch corrected with filter.tiff",
     width = 600, height = 1000)
Heatmap(subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)),col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor  = dds$EC_batch, EC_Sex = dds$EC_Sex,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "Unknown" = "gold2"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("F" = "hotpink", "M" = "navyblue"))))
dev.off()

save(dds, diagnosis_col, EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, modelMatrix, my_metadata, p, subset, vst, vst_pcDat, vstcounts, 
     design, DiagnosisCol, newnames, newnames_subset, topVarGenes, col_fun, dds, eccountdata, ecsampleinfo, modelMatrix, vst, vst_pcDat, vstcounts, design, DiagnosisCol, diagnosis_col, 
     EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, my_metadata, p, subset, design, newnames, newnames_subset,
     topVarGenes, col_fun,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data no filter sex removed limma correction.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC data no filter sex removed limma correction.RData")

