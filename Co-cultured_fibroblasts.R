# Loading data in ----

library (tidyverse)

# load in the preproccessing section, with sex removed
load("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Preprocessing_of_all_nosex.RData")

allsampleinfo <- read_tsv("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Sample_info/EC_fibroblast_sample_info_correct_sex.txt")

# Select just the fibroblasts
fibrocountdata <- allcountdata[,c (-4,-5,-6,-10,-11,-12,-14,-20,-21,-23,-25,-30,-31,-33,-35,-40,-41,
                                   -43,-45,-50,-51,-53,-55,-60,-61,-63,-65,-70,-71,-73,-75)]
fibrosampleinfo <- allsampleinfo[-c(4,5,6,10,11,12,14,20,21,23,25,30,31,33,35,40,41,43,45,50,51,53,55,60,61,63,65,70,71,73,75),]

# select only co-cultured fibrblasts 
fibrocountdata <- fibrocountdata[,c (-4:-6, -10:-12, -16:-18, -22:-24, -28:-30, -34:-36, -40:-42, -46:-48)]
fibrosampleinfo <- fibrosampleinfo[-c(4:6, 10:12, 16:18, 22:24, 28:30, 34:36, 40:42, 46:48),]

# Now remove BX076, aka 3761 and 3766
fibrocountdata <- fibrocountdata[,c(-10)]
fibrosampleinfo <- fibrosampleinfo[-c(10),]

# load packages 
library(DESeq2)
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)
library (cowplot) # to organise plots of figures

# Quick look at the PCA with outliers removed
vstcounts <- vst(fibrocountdata)
vst_pcDat <- prcomp(t(vstcounts))

# Quick sanity check everything is working
autoplot(vst_pcDat, 
         data = fibrosampleinfo, 
         colour = 'Diagnosis',
         main = 'Fibroblasts VST transform PCA',
         max.overlaps = 10,
         size = 3)

autoplot(vst_pcDat, 
         data = fibrosampleinfo, 
         colour = 'Co_or_mono',
         shape = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3)  +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample_name), box.padding = 0.8)


# Then run DESeq2 ----

# design suggesting diagnosis causes differences in groups
design <- as.formula(~ Diagnosis)

# add in Fibro donor too?? 

# what is the y intercept as standard? 
modelMatrix <- model.matrix(design, data = fibrosampleinfo)
modelMatrix

# change y intercept to resolving
fibrosampleinfo$Diagnosis <- factor(fibrosampleinfo$Diagnosis, levels = c("Res", "VeRA", "JRep"))
modelMatrix <- model.matrix(design, data = fibrosampleinfo)
modelMatrix

dds <- DESeqDataSetFromMatrix(countData = fibrocountdata,
                              colData = fibrosampleinfo,
                              design = design)

dds <- DESeq(dds)

save(fibrocountdata, fibrosampleinfo, modelMatrix, vst_pcDat, vstcounts, design, dds,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibroblasts_sex_genes_removed_dds_co-culture.RData")

# plot the dispersion estimates
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Dispersion_plot.tiff",
     width = 500, height = 300)
plotDispEsts(dds)
dev.off()

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


# To get gene names rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

# Biplot and pairplot of fibroblats 

# Diagnosis
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts PCA coloured by diagnosis no sex biplot.tiff",
     width = 750, height = 500)
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right',
                 encircle = TRUE)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts pairsplot coloured by diagnosis no sex biplot.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=3)
dev.off()

# Joint 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts PCA coloured by Joint no sex biplot.tiff",
     width = 750, height = 500)
PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts pairsplot coloured by Joint no sex biplot.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=3)
dev.off()

# Fibroblast sex
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts PCA coloured by Fibroblast sex no sex biplot.tiff",
     width = 750, height = 500)
PCAtools::biplot(p, colby = 'Fibroblast_Sex', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast_Sex', legendPosition = 'right',
                 ellipse = TRUE, ellipseLevel = 0.9)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts pairsplot coloured by Fiborblast_Sex no sex pairsplot.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Fibroblast_Sex', pointSize=3)
dev.off()

# EC sex
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts PCA coloured by EC sex no sex biplot.tiff",
     width = 750, height = 500)
PCAtools::biplot(p, "PC3","PC4", colby = 'EC_Sex', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right',
                 ellipse = TRUE, ellipseLevel = 0.9)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured fibroblasts pairsplot coloured by EC_Sex no sex pairsplot.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'EC_Sex', pointSize=3)
dev.off()

PCAtools::pairsplot(p, colby = 'EC_donor', pointSize=3)
PCAtools::biplot(p, colby = 'EC_donor', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'EC_donor', legendPosition = 'right',
                 ellipse = TRUE, ellipseLevel = 0.9)

#PCAtools::biplot(p,x = "PC3", y = "PC4", colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right',
 #                encircle = TRUE,   encircleAlpha = 1/6)
# Or can draw ellipse in - these are statistically sig, the encirle ones aren't 
#PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right',
                # ellipse = TRUE,   ellipseAlpha = 1/6, ellipseLevel = 0.9)


# PCA plots coloured by other variables
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured Fibroblasts PCA coloured by different factors no sex.tiff",
     width = 1000, height = 1000)
a <- PCAtools::biplot(p, colby = 'EC_donor', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'EC donor', legendPosition = 'right')
b <- PCAtools::biplot(p, colby = 'Joint', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')
c <- PCAtools::biplot(p, colby = 'Fibroblast_ID', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast donor', legendPosition = 'right')
d <- PCAtools::biplot(p, colby = 'Culture', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Culture', legendPosition = 'right')
e <- PCAtools::biplot(p, colby = 'Fibroblast_Sex', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast Sex', legendPosition = 'right')
f <- PCAtools::biplot(p, colby = 'EC_Sex', pointSize=3, showLoadings = FALSE, lab = NULL, sizeLoadingsNames = 5, colLegendTitle = 'EC Sex', legendPosition = 'right')
plot_grid(a,b,c,d,e,f, 
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)
dev.off()

# Loadings plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblasts PCA loadings plot no sex.tiff",
     width = 600, height = 600)
PCAtools::plotloadings(p, labSize = 3)
dev.off()

# pairplots - to see if effects of these variables are seen in other PCs
# doesn't seem to be, therefore not shown in results
#PCAtools::pairsplot(p, colby = 'EC_donor', pointSize=3)
#PCAtools::pairsplot(p, colby = 'Joint', pointSize=3)

PCAtools::pairsplot(p, colby = 'Fibroblast_Sex', pointSize=3)
PCAtools::pairsplot(p, colby = 'EC_Sex', pointSize=3)

# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
topVarGenes <- head( order( rowVars(vst), decreasing=TRUE ), 100 )
subset <- vst[ topVarGenes, ]
newnames_subset <- mapIds(EnsDb.Hsapiens.v79,
                          keys = rownames(subset),
                          column = c('SYMBOL'),
                          keytype="GENEID")
newnames_subset <- ifelse(is.na(newnames_subset) | duplicated(newnames_subset),
                          names(newnames_subset), newnames)
rownames(subset) <- newnames_subset

diagnosis_col <- as.data.frame(dds$Diagnosis)
fibroblast_sex_col <- as.data.frame(dds$Fibroblast_Sex)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(0, 20, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Co-cultured Fibroblasts 100 variable genes heatmap no sex.tiff",
     width = 800, height = 1200)
Heatmap(subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor= dds$EC_batch,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                                                             "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                                                             "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "Unknown" = "gold2"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"))))
dev.off()

# Also to look at the z-scored values
subset - rowMeans((subset)) -> subset_meanSubtract
subset_meanSubtract/rowSds(as.matrix(subset)) ->
  heatmap_data_zscores

Heatmap(heatmap_data_zscores,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), #col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor= dds$EC_batch,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "Unknown" = "gold2"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"))))



save(dds, fibrocountdata, fibrosampleinfo, modelMatrix, vst_pcDat, vstcounts, design,  
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibroblasts_sex_removed.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibroblasts_sex_removed.RData")

# Then look at the results 

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

# To get the names
# load
library('org.Hs.eg.db')

# To see what they're called
columns(org.Hs.eg.db)

# Get ensembl names as rownames, although this says VeRA_vs_Res, it's the same gene list for all
ensembl <- as.vector(rownames(VeRA_vs_Res))

entrez <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'ENTREZID', 'ENSEMBL'))
symbol <- as.vector(mapIds(org.Hs.eg.db, ensembl, 'SYMBOL', 'ENSEMBL'))

# JRep_vs_Res
# Get DEGs as a data frame
annotJRep_vs_Res <- as.data.frame(JRep_vs_Res) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotJRep_vs_Res <- cbind(annotJRep_vs_Res, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotJRep_vs_Res), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_JRep_vs_Res_all.csv")

# export the significant results
JRep_vs_Res_sig <- as.data.frame(annotJRep_vs_Res)
JRep_vs_Res_sig <- JRep_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
JRep_vs_Res_sig <- JRep_vs_Res_sig[order(JRep_vs_Res_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "JRep", "Norm"), alpha = 0.05)
write.csv(as.data.frame(JRep_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_JRep_vs_Res_sig.csv")

# VeRA_vs_Res
# Get DEGs as a data frame
annotVeRA_vs_Res <- as.data.frame(VeRA_vs_Res) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_Res <- cbind (annotVeRA_vs_Res, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_Res), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_VeRA_vs_Res_all.csv")

# export the significant results
VeRA_vs_Res_sig <- as.data.frame(annotVeRA_vs_Res)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig[order(VeRA_vs_Res_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_VeRA_vs_Res_sig.csv")

# VeRA vs Jrep ----
# Get DEGs as a data frame
annotVeRA_vs_JRep <- as.data.frame(VeRA_vs_JRep) %>% 
  rownames_to_column("GeneID") 

# Add entrez IDs and gene symbols to the end
annotVeRA_vs_JRep <- cbind (annotVeRA_vs_JRep, entrez, symbol)

# export the results as csv
write.csv(as.data.frame(annotVeRA_vs_JRep), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_VeRA_vs_JRep_all.csv")
# export the significant results
VeRA_vs_JRep_sig <- as.data.frame(annotVeRA_vs_JRep)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig %>% dplyr::filter (padj < 0.05 & log2FoldChange > 1 | padj < 0.05 & log2FoldChange < -1)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig[order(VeRA_vs_JRep_sig$log2FoldChange),]
#<- results(macro_dds_filter, contrast=c("Diagnosis", "VeRA", "Norm"), alpha = 0.05)
write.csv(as.data.frame(VeRA_vs_JRep_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/Fibro_co_VeRA_vs_JRep_sig.csv")


# THen get volcano plots
library(EnhancedVolcano)

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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblast co volcanos 0.1.tiff",
     width = 1500, height = 500)
a <- EnhancedVolcano(JRep_vs_Res_shrunk,
                     lab = rownames(JRep_vs_Res_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'JRep vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 4)
b <- EnhancedVolcano(VeRA_vs_Res_shrunk,
                     lab = rownames(VeRA_vs_Res_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs Res',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 4)
c <- EnhancedVolcano(VeRA_vs_JRep_shrunk,
                     lab = rownames(VeRA_vs_JRep_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs JRep',
                     pCutoff = 0.1,
                     FCcutoff = 1,
                     labSize = 4)
cowplot::plot_grid(a, b, c, 
                   labels = c("", "", "") ,
                   ncol = 3, nrow = 1)
dev.off()


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblast co volcanos 0.05.tiff",
     width = 1500, height = 500)
a <- EnhancedVolcano(JRep_vs_Res_shrunk,
                     lab = rownames(JRep_vs_Res_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'JRep vs Res',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 5)
b <- EnhancedVolcano(VeRA_vs_Res_shrunk,
                     lab = rownames(VeRA_vs_Res_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs Res',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 5)
c <- EnhancedVolcano(VeRA_vs_JRep_shrunk,
                     lab = rownames(VeRA_vs_JRep_shrunk),
                     xlim = c(-5,5),
                     ylim = c(0,3.5),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs JRep',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 5)
cowplot::plot_grid(a, b, c, 
                   labels = c("", "", "") ,
                   ncol = 3, nrow = 1)
dev.off()

plotCounts(dds, gene = "ENSG00000120708", intgroup = "Diagnosis")


COMP <- plotCounts(dds, gene = "ENSG00000105664", intgroup = "Diagnosis", returnData = TRUE)
IL36RN <- plotCounts(dds, gene = "ENSG00000136695", intgroup = "Diagnosis", returnData = TRUE)
IL36B <- plotCounts(dds, gene = "ENSG00000136696", intgroup = "Diagnosis", returnData = TRUE)

IL6 <- plotCounts(dds, gene = "ENSG00000136244", intgroup = "Diagnosis", returnData = TRUE)
SOCS3 <- plotCounts(dds, gene = "ENSG00000136696", intgroup = "Diagnosis", returnData = TRUE)
SOCS1 <- plotCounts(dds, gene = "ENSG00000185338", intgroup = "Diagnosis", returnData = TRUE)
TGFB1 <- plotCounts(dds, gene = "ENSG00000105329", intgroup = "Diagnosis", returnData = TRUE)

TNF <- plotCounts(dds, gene = "ENSG00000232810", intgroup = "Diagnosis", returnData = TRUE)
ACVRL1 <- plotCounts(dds, gene = "ENSG00000139567", intgroup = "Diagnosis", returnData = TRUE)

CXCL10 <- plotCounts(dds, gene = "ENSG00000169245", intgroup = "Diagnosis", returnData = TRUE)
CXCL9 <- plotCounts(dds, gene = "ENSG00000138755", intgroup = "Diagnosis", returnData = TRUE)


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblast co interest genes.tiff",
     width = 800, height = 500)
a <- ggplot(COMP, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("COMP")+
  ylab("Normalised count")
b <- ggplot(IL36RN, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL36RN")+
  ylab("Normalised count")
c <-  ggplot(IL36B, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL36B")+
  ylab("Normalised count")
d <- ggplot(IL6, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("IL6")+
  ylab("Normalised count")
e <- ggplot(SOCS3, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SOCS3")+
  ylab("Normalised count")
f <- ggplot(SOCS1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("SOCS1")+
  ylab("Normalised count")
g <- ggplot(TGFB1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TGFB1")+
  ylab("Normalised count")
h <- ggplot(TNF, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("TNF")+
  ylab("Normalised count")
i <- ggplot(ACVRL1, aes(x=Diagnosis, y=count, colour = fibrosampleinfo$Diagnosis)) + 
  theme(legend.position = "none") +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ggtitle("ACVRL1")+
  ylab("Normalised count")
cowplot::plot_grid(a, b, c, d, e, f, g, h, i, 
                   labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I") ,
                   ncol = 3, nrow = 3)
dev.off()

# To get heat maps with the sig genes ----

vst <- assay(vst(dds))

library(EnsDb.Hsapiens.v79)

# New names for the vst object
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(vst),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)

rownames(vst) <- newnames

#VeRA vs JRep

VeRA_vs_JRep_df <- as.data.frame(VeRA_vs_JRep)

VeRA_vs_JRep_adj_sig <- VeRA_vs_JRep_df %>% dplyr::filter (padj < 0.05 & abs(log2FoldChange > 1 | log2FoldChange < -1))

# new names for the adjusted sig subsetted genes
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(VeRA_vs_JRep_adj_sig),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(VeRA_vs_JRep_adj_sig) <- newnames

VeRA_vs_JRep_adj_sig_ls <- rownames(VeRA_vs_JRep_adj_sig)

VeRA_vs_JRep_subset_adj <- vst[VeRA_vs_JRep_adj_sig_ls, ]

# JRep vs Res

JRep_vs_Res_df <- as.data.frame(JRep_vs_Res)

JRep_vs_Res_adj_sig <- JRep_vs_Res_df %>% dplyr::filter (padj < 0.05 & abs(log2FoldChange > 1 | log2FoldChange < -1))

newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(JRep_vs_Res_adj_sig),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(JRep_vs_Res_adj_sig) <- newnames

JRep_vs_Res_adj_sig_ls <- rownames(JRep_vs_Res_adj_sig)

JRep_vs_Res_subset_adj <- vst[JRep_vs_Res_adj_sig_ls, ]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

col_fun <- colorRamp2(breaks = seq(9, 12, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101)) 


col_fun <- colorRamp2(breaks = seq(0, 15, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblast co heatmaps of sig.tiff",
     width = 500, height = 500)
Heatmap(VeRA_vs_JRep_subset_adj,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)),col = col_fun, column_title = "VeRA vs JRep",
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor= dds$EC_batch,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "Unknown" = "gold2"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"))))
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/Fibroblast co heatmaps of sig JRep vs Res.tiff",
     width = 500, height = 600)
Heatmap(JRep_vs_Res_subset_adj,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)),col = col_fun,column_title = "JRep vs Res",
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor= dds$EC_batch,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "Unknown" = "gold2"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"))))
dev.off()

# Now try with cluster profiler
library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code')

# For this you need entrez IDs
JRep_vs_Res$ensembl <- sapply( strsplit( rownames(JRep_vs_Res), split="\\+" ), "[", 1 )

library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = JRep_vs_Res$ensembl,
                  mart = ensembl )

idx <- match( filtered_Jrep_vs_Res$ensembl, genemap$ensembl_gene_id )
filtered_Jrep_vs_Res$entrez <- genemap$entrezgene[ idx ]
filtered_Jrep_vs_Res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

# So that any pathview stuff goes to the right location 
setwd("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only")

# JRep vs Res 

# Select only the significant genes, for pathway analysis using loose p values
sigGenesJRep_vs_Res <- annotJRep_vs_Res %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

# VeRA vs JRep

# Select only the significant genes, for pathway analysis using loose p values
sigGenesVeRA_vs_JRep <- annotVeRA_vs_JRep %>% 
  drop_na(entrez, padj) %>% 
  dplyr::filter(pvalue < 0.1 & abs(log2FoldChange) > 1) %>% 
  pull(entrez)

# Just to look at 
sigGenesJRep_vs_Res_df <- as.data.frame(sigGenesJRep_vs_Res)

library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code')

# KEGG stuff

kk_Jrep_vs_Res <- clusterProfiler::enrichKEGG(gene = sigGenesJRep_vs_Res, organism = 'hsa')
kk_VeRA_vs_JRep <- enrichKEGG(gene = sigGenesVeRA_vs_JRep, organism = 'hsa')

# Print the top head
kk_printed <- print(head(kk, n=20))

# Make this a dot plot
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Fibroblast only/Final/KEGG dotplot.tiff",
     width = 750, height = 300)
a <- dotplot(kk_VeRA_vs_JRep)
b <- dotplot(kk_Jrep_vs_Res)
cowplot::plot_grid(a, b, 
                   labels = c("A", "B") ,
                   ncol = 2, nrow = 1)
dev.off()

# Make files for GSEA

counts <- counts(dds)
write.csv(counts,file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/fibroblast_dds_counts_table.csv")

write.csv(fibrosampleinfo$Sample_name,file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/Fibroblasts/fibroblast_dds_sample_names.csv")
