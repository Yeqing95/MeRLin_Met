setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/original/")
source("~/Documents/Wistar/Haiyin/barcode extractor/utils_sc.R")
library(ggplot2)

df_primary <- data.frame(len = 10:30, timepoint = "Primary",
                         unique_bc = c(length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 10)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 11)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 12)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 13)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 14)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 15)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 16)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 17)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 18)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 19)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 20)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 21)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 22)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 23)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 24)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 25)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 26)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 27)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 28)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 29)$barcode)),
                                       length(unique(extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 30)$barcode))))

df_lung <- data.frame(len = 10:30, timepoint = "Lung",
                         unique_bc = c(length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 10)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 11)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 12)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 13)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 14)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 15)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 16)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 17)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 18)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 19)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 20)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 21)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 22)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 23)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 24)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 25)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 26)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 27)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 28)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 29)$barcode)),
                                       length(unique(extract_barcode(bam = "Lung/Lung_barcode.bam", len = 30)$barcode))))

df_liver <- data.frame(len = 10:30, timepoint = "Liver",
                      unique_bc = c(length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 10)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 11)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 12)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 13)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 14)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 15)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 16)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 17)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 18)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 19)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 20)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 21)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 22)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 23)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 24)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 25)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 26)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 27)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 28)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 29)$barcode)),
                                    length(unique(extract_barcode(bam = "Liver/Liver_barcode.bam", len = 30)$barcode))))


df <- rbind(df_primary, df_lung, df_liver)  

ggplot(data = df, mapping = aes(x = len, y = unique_bc, group = timepoint, color = timepoint)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = ~timepoint, ncol = 3) +
  theme_bw()

#######################################################################################################
primary <- extract_barcode(bam = "Primary/Tumor_barcode.bam", len = 10)
lung <- extract_barcode(bam = "Lung/Lung_barcode.bam", len = 10)
liver <- extract_barcode(bam = "Liver/Liver_barcode.bam", len = 10)

write.table(x = primary, file = "Primary/Primary_bc.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(x = lung, file = "Lung/Lung_bc.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(x = liver, file = "Liver/Liver_bc.txt", sep = "\t", quote = F, row.names = F, col.names = T)


