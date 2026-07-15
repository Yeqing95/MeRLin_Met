setwd("~/Documents/Wistar/TCGA/")
library(readxl)
library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(GSVA)
library(tibble)
library(ggpubr)

#load clinical information
Clinical <- read_excel(path = "SKCM/SKCM_Clinical_Index.xlsx")
Clinical <- Clinical %>%
  subset(vital_status != "Not Reported") %>%
  select(c("submitter_id", "ajcc_pathologic_stage", "age_at_index", "gender",
           "vital_status", "days_to_death", "days_to_last_follow_up",
           "treatments_pharmaceutical_treatment_or_therapy")) %>%
  rename(Patient = submitter_id, ajcc_stage = ajcc_pathologic_stage, age = age_at_index,
         treatment = treatments_pharmaceutical_treatment_or_therapy)

Clinical$ajcc_stage <- str_remove(string = Clinical$ajcc_stage, pattern = "[A|B|C]")
Clinical$ajcc_stage <- str_replace_na(string = Clinical$ajcc_stage, replacement = "Not Reported")

Clinical$days <- ifelse(test = Clinical$vital_status == "Dead", 
                        yes = Clinical$days_to_death, no = Clinical$days_to_last_follow_up)

#load gene expression profile
RNAseq <- read.csv(file = "SKCM/SKCM_RNAseq_TPM_Counts.csv", check.names = F)
gene_list <- read_excel("SKCM/gene_list.xlsx")

Signature <- read_excel(path = "../Haiyin/Metastasis/2025-09-11 MET_DEG/2025-09-04 Group Signature.xlsx", sheet = "Group 2_Lipid")
Lipid <- as.character(na.omit(Signature$Gene))
Signature <- read_excel(path = "../Haiyin/Metastasis/2025-09-11 MET_DEG/2025-09-04 Group Signature.xlsx", sheet = "Group 1_NC")
NC <- as.character(na.omit(Signature$Gene))

Signature <- read_excel(path = "../Haiyin/MeRLin/Data/Source data/Marine_signatures.xlsx")
NC <- Signature$`Neural Crest-like`

RNAseq$gene_id <- gene_list$gene_name[match(x = RNAseq$gene_id, table = gene_list$gene_id)]
RNAseq <- RNAseq[-which(duplicated(RNAseq$gene_id)),]
rownames(RNAseq) <- NULL
RNAseq <- column_to_rownames(RNAseq, var = "gene_id")
RNAseq <- RNAseq[-which(rowSums(RNAseq) == 0),]

param <- ssgseaParam(exprData = as.matrix(RNAseq), geneSets = list(Lipid = Lipid, NC = NC))
gsva_res <- gsva(param = param)

df <- as.data.frame(t(gsva_res))
df <- rownames_to_column(df, var = "Patient")
df$Patient <- str_split_fixed(string = df$Patient, pattern = "-[0-9]{2}[A|B]", n = 2)[,1]

#load mutation profile
wes <- read_excel("SKCM/SKCM_Masked_Somatic_Mutation.xlsx")
wes <- subset(wes, Hugo_Symbol == "BRAF" & HGVSp_Short %in% c("p.V640E", "p.V640M")) #p.V640G
braf_pat <- str_split_fixed(string = wes$Tumor_Sample_Barcode, pattern = "-[0-9]{2}[A|B]", n = 2)[,1]

df$BRAF <- "WT"
df$BRAF[which(df$Patient %in% braf_pat)] <- "Mut"

#plot
df <- merge(x = Clinical, y = df, by = "Patient")
df2 <- subset(df, BRAF == "Mut")

df_surv <- select(.data = df2, c("NC", "vital_status", "days"))
df_surv <- subset(df_surv, !is.na(days) & days > 0)
df_surv$vital_status <- ifelse(test = df_surv$vital_status == "Dead", yes = 1, no = 0)

res.cut <- surv_cutpoint(df_surv, time = "days", event = "vital_status", variables = "NC", minprop = 0.1)
plot(res.cut, "NC", palette = "npg")

res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(days, vital_status) ~ NC, data = res.cat)

ggsurvplot(fit,
           pval = TRUE, 
           font.x = c(12, "bold", "Black"),
           font.y = c(12, "bold", "Black"),
           font.tickslab = c(9, "plain", "Black"),
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#F8766D","#619CFF"),
           data = res.cat,
           fun = "pct"
)

res.cat2 <- df_surv
mid <- median(res.cat2$NC)
res.cat2$NC[which(res.cat2$NC >= mid)] <- "high"
res.cat2$NC[which(res.cat2$NC < mid)] <- "low"

fit2 <- survfit(Surv(days, vital_status) ~ NC, data = res.cat2)

ggsurvplot(fit2,
           pval = TRUE, 
           font.x = c(12, "bold", "Black"),
           font.y = c(12, "bold", "Black"),
           font.tickslab = c(9, "plain", "Black"),
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#F8766D","#619CFF"),
           data = res.cat2,
           fun = "pct"
)


ggplot(data = df, mapping = aes(x = ajcc_stage, y = NC, color = BRAF)) +
  geom_boxplot()


wilcox.test(x = df$NC[which(df$BRAF == "Mut" & df$ajcc_stage == "Stage III")], 
            y = df$NC[which(df$BRAF == "Mut" & df$ajcc_stage == "Stage IV")])

sig_results <- data.frame(
  group1 = c("Stage II", "Stage II", "Stage III"),
  group2 = c("Stage III", "Stage IV", "Stage IV"),
  p.adj = c(0.6412, 0.7208, 0.4804),
  label = c("ns", "ns", "ns"),
  y.position = c(2.4,2.5,2.6)
)

stat_pvalue_manual(sig_results, label = "label", tip.length = 0.01)

ggplot(data = subset(df, ajcc_stage %in% c("Stage II","Stage III","Stage IV")), 
       mapping = aes(x = ajcc_stage, y = NC, color = BRAF)) +
  geom_boxplot() +
  theme_classic() +
  stat_pvalue_manual(sig_results, label = "label", tip.length = 0.01)





wilcox.test(x = df$NC[which(df$BRAF == "WT" & df$ajcc_stage == "Stage II")], 
            y = df$NC[which(df$BRAF == "WT" & df$ajcc_stage == "Stage III")])

sig_results_mut <- data.frame(
  group1 = c("Stage II","Stage II","Stage III"),
  group2 = c("Stage III","Stage IV","Stage IV"),
  p.adj = c(0.6412,0.7208,0.4804),
  label = c("ns","ns","ns"),
  y.position = c(2.3,2.35,2.4),
  BRAF = "Mut"
)

sig_results_wt <- data.frame(
  group1 = c("Stage II","Stage II","Stage III"),
  group2 = c("Stage III","Stage IV","Stage IV"),
  p.adj = c(0.8946,0.1278,0.09496),
  label = c("ns","ns","ns"),
  y.position = c(2.3,2.35,2.4),
  BRAF = "WT"
)
sig_results_wt$BRAF <- "WT"

sig_results_all <- rbind(sig_results_mut, sig_results_wt)

# 画图
ggplot(subset(df, ajcc_stage %in% c("Stage II","Stage III","Stage IV")), 
       aes(x = ajcc_stage, y = NC, color = BRAF)) +
  geom_boxplot() +
  facet_wrap(~BRAF) +
  theme_classic() +
  stat_pvalue_manual(sig_results_all, label = "label", tip.length = 0.01)


