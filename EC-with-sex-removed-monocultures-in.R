# Loading data in ----
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

# load packages 
library(DESeq2)
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)

# Vectors for colours
DiagnosisCol <- match(ecsampleinfo$Diagnosis, c("EC_monoculture", "Res", "VeRA", "JRep")) + 1

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
ecsampleinfo$Diagnosis <- factor(ecsampleinfo$Diagnosis, levels = c("EC_monoculture", "Res", "VeRA", "JRep"))
modelMatrix <- model.matrix(design, data = ecsampleinfo)
modelMatrix

# DESeq
dds <- DESeqDataSetFromMatrix(countData = eccountdata,
                              colData = ecsampleinfo,
                              design = design)

dds <- DESeq(dds)

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
library(EnsDb.Hsapiens.v79)
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames

PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=3)
PCAtools::pairsplot(p, colby = 'EC_Sex', pointSize=3)
PCAtools::pairsplot(p, colby = 'Fibroblast_Sex', pointSize=3)

# Plots put here, and made pretty in powerpoint
PCAtools::pairsplot(p, colby = 'EC_donor', pointSize=3)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=3)

PCAtools::plotloadings(p)

# To get the legends
PCAtools::biplot(p, "PC2", "PC5", colby = 'EC_donor', pointSize=3.5, showLoadings = TRUE, lab = NULL, encircle = TRUE,
                 sizeLoadingsNames = 5, colLegendTitle = 'EC donor', legendPosition = 'right')
PCAtools::biplot(p,colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, encircle = TRUE,
                 sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC biplot coloured by Diagnosis.tiff",
     width = 600, height = 400)
PCAtools::biplot(p, "PC2", "PC5", colby = 'Diagnosis', pointSize=3.5, showLoadings = TRUE, lab = NULL, encircle = TRUE,
                 sizeLoadingsNames = 5, colLegendTitle = 'Diagnosis', legendPosition = 'right')
dev.off()

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC pairsplot coloured by sex.tiff",
     width = 1000, height = 600)
PCAtools::pairsplot(p, colby = 'Diagnosis', pointSize=3)
dev.off()

PCAtools::bioplot(p, colby = 'Fibroblast_Sex', pointSize=3.5, showLoadings = TRUE, lab = NULL, 
                 sizeLoadingsNames = 5, colLegendTitle = 'Fibroblast_Sex', legendPosition = 'right')


# Heatmap of the most variable genes
# Selects the top 100 most variable genes across the samples
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

# Colour palette function for hetamap
col_fun <- colorRamp2(breaks = seq(5, 20, length.out=101),
                      colorRampPalette(rev(brewer.pal(n=11, "RdYlBu")))(101))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC only output/EC 100 variable genes heatmap no sex.tiff",
     width = 1000, height = 1200)
Heatmap(subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), col = col_fun,
        top_annotation = HeatmapAnnotation(Diagnosis=dds$Diagnosis,Fibroblast_Sex = dds$Fibroblast_Sex,
                                           Joint = dds$Joint, Culture = dds$Co_vs_mono, EC_Sex = dds$EC_Sex, col = list
                                           (Diagnosis = c("Res" = "chartreuse4", "VeRA" = "dodgerblue3", 
                                                          "JRep" = "firebrick", "EC_monoculture" = "grey"),
                                             EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                             Fibroblast_Sex = c("M" = "darkorange2", "F"= "purple4", "EC_Mono" = "grey"),
                                             Joint = c("mcp" = "seagreen2", "ankle" = "lightpink2", "knee" = "lightgoldenrod2", "EC_mono" = "grey"),
                                             Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4", "EC_Mono" = "grey"))))
dev.off()

save(dds, diagnosis_col, EC_sex_col, eccountdata, ecsampleinfo, fibroblast_sex_col, modelMatrix, my_metadata, p, subset, vst, vst_pcDat, vstcounts, design, 
     DiagnosisCol, newnames, newnames_subset, topVarGenes, col_fun, file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC_preprocessing_sex_removed.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Outputs/EC/EC_preprocessing_sex_removed.RData")

# For Poppy ----
plotCounts(dds, gene='ENSG00000196611', intgroup = "EC_batch") # MMP1
plotCounts(dds, gene='ENSG00000087245', intgroup = "EC_batch") # MMP2
plotCounts(dds, gene='ENSG00000100985', intgroup = "EC_batch") # MMP9
plotCounts(dds, gene='ENSG00000137801', intgroup = "EC_batch") # THBS1

# New names for the vst object
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(vst),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)

rownames(vst) <- newnames

# subset from raw, transformed
vstcounts <- vst(eccountdata)

Poppy_genes <- c("ENSG00000196611", "ENSG00000087245", "ENSG00000100985", "ENSG00000137801")
Poppy_subset <- vstcounts[ Poppy_genes, ]

# new names for the subsetted genes
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(Poppy_subset),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(Poppy_subset) <- newnames

diagnosis_col <- as.data.frame(dds$Diagnosis)
fibroblast_sex_col <- as.data.frame(dds$Fibroblast_Sex)
culture_conditions <- as.data.frame(dds$Co_vs_mono)

Heatmap(Poppy_subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), cluster_columns = FALSE,
        top_annotation = HeatmapAnnotation( EC_donor = dds$EC_batch, EC_Sex = dds$EC_Sex, Culture = dds$Co_or_mono, col = list
                                            ( EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                              EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                           "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                           "EC0242" = "magenta4", "EC0221"= "cornflowerblue")),
                                            Culture = c("Co-culture"="deeppink2", "Mono-culture" = "royalblue4")))

# For Imy

# subset from raw, transformed
vstcounts <- vst(eccountdata)

# For CD55, THY1
Imy_genes <- c("ENSG00000196352", "ENSG00000154096", "ENSG00000161638")

Imy_subset <- vstcounts[ Imy_genes, ]

# new names for the subsetted genes
newnames <- mapIds(EnsDb.Hsapiens.v79,
                   keys = rownames(Imy_subset),
                   column = c('SYMBOL'),
                   keytype="GENEID")

newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(Imy_subset) <- newnames

Heatmap(Imy_subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), 
        top_annotation = HeatmapAnnotation( EC_donor = dds$EC_batch, EC_Sex = dds$EC_Sex, Culture = dds$Co_or_mono, col = list
                                            ( EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                              EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                           "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                           "EC0242" = "magenta4", "EC0221"= "cornflowerblue")),
                                            Culture = c("Co"="deeppink2", "Mono" = "royalblue4")))

Heatmap(Imy_subset,  column_labels = dds$Sample_description, column_names_gp = gpar(fontsize = c(8), fontface = "bold"), 
        row_names_gp = gpar(fontsize = c(8)), cluster_columns = FALSE, cluster_rows = FALSE,
        top_annotation = HeatmapAnnotation( EC_donor = dds$EC_batch, EC_Sex = dds$EC_Sex, Culture = dds$Co_or_mono, col = list
                                            ( EC_Sex = c("M" = "navyblue", "F" = "hotpink"),
                                              EC_donor = c("EC2015" = "firebrick",  "EC2016" = "palegreen3", "EC0229" = "goldenrod", 
                                                           "EC0223" = "darkcyan", "EC0271" = "hotpink3", "EC0270" = "seagreen", 
                                                           "EC0242" = "magenta4", "EC0221"= "cornflowerblue")),
                                            Culture = c("Co"="deeppink2", "Mono" = "royalblue4")))

