# Remove the Y chromosome and XIST gene
library(tidyverse)

# loading in original smaple info - all data
allsampleinfo <- read_tsv("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Sample_info/EC_fibroblast_sample_info_correct_sex.txt")

# loading in the counts table
seqdata <- read_tsv("~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/Featurecounts/Featurecounts.txt", comment = "#")

typeof(seqdata)

# to remove the Y chromosome genes
seqdata <- seqdata[!grepl("Y", seqdata$Chr),]

# And to remove the XIST and TSIX genes
seqdata <-seqdata[!grepl("ENSG00000229807", seqdata$Geneid),]
seqdata <-seqdata[!grepl("ENSG00000270641", seqdata$Geneid),]

allcountdata <- seqdata %>%
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  rename_all(str_remove, "/rds/projects/m/mcgetthm-ra-icase/Fibroblast-endothelial/Analysis/Aligned/") %>% 
  # remove this string  from the column names
  select(allsampleinfo$File_names) %>% # keep sample columns using sampleinfo
  as.matrix()

# load packages
library(ggfortify)
library (ggrepel)
library (DESeq2)
library (PCAtools)
library(limma)

# Transform data with variance stabalisation (VST)
vstcounts <- vst(allcountdata)

# PCA
vst_pcDat <- prcomp(t(vstcounts))

# plot PCA
tiff(filename = "~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/VST PCA alldata no sex.tiff",
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

save(allcountdata, allsampleinfo, file="~/OneDrive/Documents/PhD/Projects/Endothelial RNA seq/EC + fibro output/Preprocessing_of_all_nosex.RData")



