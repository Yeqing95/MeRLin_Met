setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/APA/")
library(MAAPER)
options(future.globals.maxSize = 24000 * 1024^2)

pas_annotation = readRDS(file = "human.PAS.hg38.rds")
gtf = "/Volumes/herlynm/linux/ychen/refSeq/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
bam_c1 = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/scRNAseq/Tumor/outs/possorted_genome_bam.bam"
bam_c2 = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/scRNAseq/Liver/outs/possorted_genome_bam.bam"
#bam_c2 = "/Volumes/herlynm/linux/ychen/Haiyin/metastasis/scRNAseq/Lung/outs/possorted_genome_bam.bam"

maaper(gtf, # full path of the GTF file
       pas_annotation, # PAS annotation
       output_dir = "./Liver_vs_Primary", # output directory
       bam_c1, bam_c2, # full path of the BAM files
       read_len = 280, # read length
       ncores = 1  # number of cores used for parallel computation 
)

# plot
setwd("~/Documents/Wistar/Haiyin/Metastasis/data/batch1/processed/APA/")
library(ggplot2)

APA1 <- read.table("/Volumes/herlynm/linux/ychen/Haiyin/metastasis/APA/Liver_vs_Primary/gene.txt", header =T)
APA2 <- read.table("/Volumes/herlynm/linux/ychen/Haiyin/metastasis/APA/Lung_vs_Primary/gene.txt", header =T)

APA1$REDu.pval.adj <- p.adjust(APA1$REDu.pval, method = "BH")
APA1$REDi.pval.adj <- p.adjust(APA1$REDi.pval, method = "BH")
APA2$REDu.pval.adj <- p.adjust(APA2$REDu.pval, method = "BH")
APA2$REDi.pval.adj <- p.adjust(APA2$REDi.pval, method = "BH")

APA1$REDu.regu <- "No"
APA1$REDu.regu[which(APA1$REDu.pval.adj < 0.05 & APA1$REDu > log(1.2))] <- "Lengthened"
APA1$REDu.regu[which(APA1$REDu.pval.adj < 0.05 & APA1$REDu < -log(1.2))] <- "Shortened"
APA2$REDu.regu <- "No"
APA2$REDu.regu[which(APA2$REDu.pval.adj < 0.05 & APA2$REDu > log(1.2))] <- "Lengthened"
APA2$REDu.regu[which(APA2$REDu.pval.adj < 0.05 & APA2$REDu < -log(1.2))] <- "Shortened"

APA1$REDi.regu <- "No"
APA1$REDi.regu[which(APA1$REDi.pval.adj < 0.05 & APA1$REDi > log(1.2))] <- "Lengthened"
APA1$REDi.regu[which(APA1$REDi.pval.adj < 0.05 & APA1$REDi < -log(1.2))] <- "Shortened"
APA2$REDi.regu <- "No"
APA2$REDi.regu[which(APA2$REDi.pval.adj < 0.05 & APA2$REDi > log(1.2))] <- "Lengthened"
APA2$REDi.regu[which(APA2$REDi.pval.adj < 0.05 & APA2$REDi < -log(1.2))] <- "Shortened"

APA1$Comparation <- "Liver_vs_Primary"
APA2$Comparation <- "Lung_vs_Primary"

common_genes <- Reduce(intersect, list(APA1$gene, APA2$gene))

df <- rbind(subset(APA1, gene %in% common_genes), 
            subset(APA2, gene %in% common_genes))

df.REDu <- as.data.frame(table(df$Comparation, df$REDu.regu))
colnames(df.REDu) <- c("Comparation", "Type", "Counts")
df.REDu <- subset(df.REDu, Type != "No")
df.REDu$Counts <- ifelse(df.REDu$Type == "Shortened", -df.REDu$Counts, df.REDu$Counts)
ggplot(data = df.REDu, mapping = aes(x = Counts, y = Comparation, fill = Type)) +
  geom_bar(stat="identity", width = 0.7) +
  theme_classic()

df2 <- subset(df, REDu.regu %in% c("Lengthened", "Shortened"))
tab <- table(df2$Comparation, df2$REDu.regu)
fisher.test(tab)



df.REDi <- as.data.frame(table(df$Comparation, df$REDi.regu))
colnames(df.REDi) <- c("Comparation", "Type", "Counts")
df.REDi <- subset(df.REDi, Type != "No")
df.REDi$Counts <- ifelse(df.REDi$Type == "Shortened", -df.REDi$Counts, df.REDi$Counts)
ggplot(data = df.REDi, mapping = aes(x = Counts, y = Comparation, fill = Type)) +
  geom_bar(stat="identity", width = 0.7) +
  theme_classic()

