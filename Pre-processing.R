
# Loading data in ----

library (tidyverse)

# loading in original smaple info - all data
allsampleinfo <- read_tsv("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Sample_info/EC_fibroblast_sample_info_correct_sex.txt")

# loading in the counts table
seqdata <- read_tsv("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Featurecounts/Featurecounts.txt", comment = "#")

# Transform the data to matrix of counts
allcountdata <- seqdata %>%
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  rename_all(str_remove, "/rds/projects/m/mcgetthm-ra-icase/Fibroblast-endothelial/Analysis/Aligned/") %>%# remove this string  from the column names
  select(allsampleinfo$File_names) %>% # keep sample columns using sampleinfo
  as.matrix()

# load packages
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)

# Plot the RIN for each sample
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/RINs.tiff", 
     width = 1000, height = 450)
par(mar=c(11,5,2,2))
barplot(allsampleinfo$RIN, names.arg = allsampleinfo$Sample_description, las = 2, col = "grey", main = "RIN")
dev.off()

# Look at the raw data
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Raw counts distribution.tiff", 
     width = 1000, height = 450)
par(mar=c(11,5,1,1))
boxplot(allcountdata, main='Raw counts distribution', las=2, names = allsampleinfo$Sample_description)
dev.off()

par(mar=c(1,1,1,1))

tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Mean vs Sd.tiff")
plot(rowMeans(allcountdata), rowSds(allcountdata), 
     main='Raw Counts : Mean Vs SD', 
     xlim=c(0,10000), 
     ylim=c(0,5000))
dev.off()

#Vectors

CellCol <- match(allsampleinfo$EC_or_Fibroblast, c("EC", "Fibroblast")) + 1
DiagnosisCol <- match(allsampleinfo$Diagnosis, c("EC_monoculture", "Res", "VeRA", "JRep")) + 1

# Transform data with variance stabalisation (VST)
vstcounts <- vst(allcountdata)

# Plot the counts distribution
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/VST counts distribution fibroblasts.tiff", 
     width = 2000, height = 450)
par(mar=c(11,5,1,1))
boxplot(vstcounts, 
        xlab="", 
        ylab="VST counts",
        las=2,
        col=CellCol,
        subset = allsampleinfo
        names = allsampleinfo$Sample_description) 
dev.off()

# Plot the mean vs Sd
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/VST mean vs SD.tiff")
par(mar=c(2,2,2,2))
plot(rowMeans(vstcounts), rowSds(vstcounts), 
     main='VST Counts : Mean Vs SD')
dev.off()

# PCA
vst_pcDat <- prcomp(t(vstcounts))

# plot PCA
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/VST PCA alldata.tiff",
     width = 700, height = 450)
autoplot(vst_pcDat, 
         data = allsampleinfo, 
         shape = 'EC_or_Fibroblast',
         colour = "Diagnosis",
         main = 'VST transform PCA',
         max.overlaps = 30,
         size = 3) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample_name), box.padding = 0.8)
dev.off()

save(allcountdata, allsampleinfo, file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Preprocessing_of_all.RData")

load(file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Preprocessing_of_all.RData")


```
