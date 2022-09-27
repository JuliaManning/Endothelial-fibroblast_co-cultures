# ComBAT Seq to remove batch variation ----

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

# Using Combat-seq to remove batch effects
library(sva)

# For running the batch variation, EC2016 has to be removed as it's by itself
eccountdata <- eccountdata[,c(-4)]
ecsampleinfo <- ecsampleinfo[-c(4),]

# Quick PCA to check there's still batch variation and that labels are correct
vstcounts <- vst(eccountdata)
pcDat <- prcomp(t(vstcounts))
DiagnosisCol <- match(ecsampleinfo$Diagnosis, c("Res", "VeRA", "JRep")) + 1

autoplot(pcDat, 
         data = ecsampleinfo, 
         colour = 'EC_batch',
         shape = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 10,
         size = 3)

my_batch <- as.vector(ecsampleinfo$EC_batch)
my_group <- ecsampleinfo$Diagnosis
adjusted <- ComBat_seq(counts = eccountdata, batch=my_batch, group=my_group)

dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                       colData = ecsampleinfo,
                                       design = as.formula(~ Diagnosis))

dds <- DESeq(dds)

#plot of the gene dispersion estimates
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/Gene_dispersion_estimate_plot_combat.tiff",
     width = 500, height = 300)
plotDispEsts(dds)
dev.off()

# PCA of this to see if batch correction has worked
vsd <- vst(dds)
plotPCA(vsd, "Diagnosis")
plotPCA(vsd, "EC_batch")
plotPCA(vsd, "Joint")
plotPCA(vsd, "Fibroblast_sex")
plotPCA(vsd, "EC_Sex")

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

vst <- assay(vst(dds))

# To remove the lower 10% of variables based on variance
p <- pca(vst, metadata = my_metadata, removeVar = 0.1)

# To get gene symbol rather than ensembl gene names
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC PCA coloured by diagnosis and batch sex removed filtered combat.tiff",
     width = 600, height = 800)
PCAtools::biplot(p, colby = 'EC_donor', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                      sizeLoadingsNames = 5, colLegendTitle = 'EC_donor', legendPosition = 'right')
PCAtools::biplot(p, colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, 
                      sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right', encircle = TRUE)
cowplot::plot_grid(a, b,
                   labels = c("", "") ,
                   ncol = 1, nrow = 2)
dev.off()


tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by EC donor sex removed filtered combat.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'EC_donor', pointSize=3)
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC pairsplot coloured by Diagnosis sex removed filtered combat.tiff",
     width = 1000, height = 500)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=3)
dev.off()

PCAtools::pairsplot(p, colby = 'EC_Sex', pointSize=3)
PCAtools::pairsplot(p, colby = 'Fibroblast_Sex', pointSize=3)
PCAtools::pairsplot(p, colby = 'Joint', pointSize=3)

PCAtools::biplot(p, colby = 'EC_Sex', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'EC_Sex', legendPosition = 'right')

PCAtools::biplot(p, colby = 'Fibroblast_Sex', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast_Sex', legendPosition = 'right')

PCAtools::biplot(p, colby = 'Joint', pointSize=3.5, showLoadings = FALSE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'Joint', legendPosition = 'right')

# genes from the PCA, 
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/Genes from biplot combat.tiff", width = 1200, height = 750)
a <- plotCounts(dds, gene='ENSG00000120708', intgroup="Diagnosis", returnData = TRUE) # TGFBI
b <- plotCounts(dds, gene='ENSG00000115461', intgroup="Diagnosis", returnData = TRUE) # IGFBP5
c <- plotCounts(dds, gene='ENSG00000123610', intgroup="Diagnosis", returnData = TRUE) # TNFAIP6
d <- plotCounts(dds, gene='ENSG00000163359', intgroup="Diagnosis", returnData = TRUE) # COL6A3
e <- plotCounts(dds, gene='ENSG00000187608', intgroup="Diagnosis", returnData = TRUE) # ISG15
f <- plotCounts(dds, gene='ENSG00000198888', intgroup="Diagnosis", returnData = TRUE) # MT-ND1
g <- plotCounts(dds, gene='ENSG00000198899', intgroup="Diagnosis", returnData = TRUE) # MT-ATP6
h <- plotCounts(dds, gene='ENSG00000228253', intgroup="Diagnosis", returnData = TRUE) # MT-ATP8
i <- plotCounts(dds, gene='ENSG00000198695', intgroup="Diagnosis", returnData = TRUE) # MT-ND6
A <- ggplot(a, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("TGFBI")
B <- ggplot(b, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("IGFBP5")
C <- ggplot(c, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("TNFAIP6")
D <- ggplot(d, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("COL6A3")
E <- ggplot(e, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("ISG15")
FF <- ggplot(f, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("MT-ND1")
G <- ggplot(g, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("MT-ATP6")
H <- ggplot(h, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("MT-ATP8")
I <- ggplot(i, aes(x=Diagnosis, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  ggtitle("MT-ND6")
cowplot::plot_grid(A, B, C, D, E, FF, G, H, I, K,
                   labels = c("", "", "", "", "", "", "", "", "", "", "") ,
                   ncol = 3, nrow = 3)
dev.off()


# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
vst <- assay(vst(dds))
library(EnsDb.Hsapiens.v79)
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
EC_sex_col <- as.data.frame(dds$EC_Sex)

library("RColorBrewer")
library(circlize)
library(ComplexHeatmap)

# Colour palette function for heatmap
col_fun <- colorRamp2(breaks = seq(5, 20, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC 100 variable genes heatmap no sex no monoculture combat.tiff",
     width = 600, height = 900)
Heatmap(subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex, EC_donor = dds$EC_batch,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, EC_Sex = dds$EC_Sex, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick"),
                                             EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                          "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                          "EC0242" = "magenta4", "EC0221"= "cornflowerblue"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4"))))
dev.off()


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

JRep_vs_Res_sig <- as.data.frame(JRep_vs_Res)
JRep_vs_Res_sig <- JRep_vs_Res_sig %>% dplyr::filter (padj < 0.05)
JRep_vs_Res_sig <- JRep_vs_Res_sig[order(JRep_vs_Res_sig$log2FoldChange),]

write.csv(as.data.frame(JRep_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/JRep_vs_Res_EC_sig_combat.csv")

VeRA_vs_Res_sig <- as.data.frame(VeRA_vs_Res)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig %>% dplyr::filter (padj < 0.05)
VeRA_vs_Res_sig <- VeRA_vs_Res_sig[order(VeRA_vs_Res_sig$log2FoldChange),]

write.csv(as.data.frame(VeRA_vs_Res_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/VeRA_vs_Res_EC_sig_combat.csv")

VeRA_vs_JRep_sig <- as.data.frame(VeRA_vs_JRep)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig %>% dplyr::filter (padj < 0.05)
VeRA_vs_JRep_sig <- VeRA_vs_JRep_sig[order(VeRA_vs_JRep_sig$log2FoldChange),]

write.csv(as.data.frame(VeRA_vs_JRep_sig), 
          file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/VeRA_vs_JRep_EC_sig_FC_combat.csv")

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

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/Final/EC combat volcanos.tiff",
     width = 1500, height = 500)
a <- EnhancedVolcano(VeRA_vs_Res_shrunk,
                     lab = rownames(VeRA_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs Res',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 3.5)
b <- EnhancedVolcano(JRep_vs_Res_shrunk,
                     lab = rownames(JRep_vs_Res_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'JRep vs Res',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 3.5)
c <- EnhancedVolcano(VeRA_vs_JRep_shrunk,
                     lab = rownames(VeRA_vs_JRep_shrunk),
                     x = 'log2FoldChange',
                     y = 'padj',
                     title = 'VeRA vs JRep',
                     pCutoff = 0.05,
                     FCcutoff = 1,
                     labSize = 3.5)
cowplot::plot_grid(b, a, c, 
                   labels = c("", "", "") ,
                   ncol = 3, nrow = 1)
dev.off()


save(adjusted, diagnosis_col, EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, my_metadata, pcDat, subset, 
     VeRA_vs_JRep, VeRA_vs_JRep_shrunk, VeRA_vs_JRep_sig, VeRA_vs_Res, VeRA_vs_Res_sig, VeRA_vs_Res_shrunk, vsd, vst
     ,vstcounts, DiagnosisCol, my_batch, my_group, newnames, newnames_subset, topVarGenes, col_fun,
     file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC-fibroblast analysis combat-seq.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC-fibroblast analysis combat-seq.RData")


