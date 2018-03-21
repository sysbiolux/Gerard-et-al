getwd()

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

# load TEPIC predictions
ADay1 = read_delim("AHR_A_Day1_0.95_50000.dat", delim = "\t", col_names = FALSE) %>%
  rename(TF = X1, ensembl_gene_id = X2, input = X3, time = X4)
ADay1

ODay1 = read_delim("AHR_O_Day1_0.95_50000.dat", delim = "\t", col_names = FALSE) %>%
  rename(TF = X1, ensembl_gene_id = X2, input = X3, time = X4)
ODay1

Day0 = read_delim("AHR_ST2_0.95_50000.dat", delim = "\t", col_names = FALSE) %>%
  rename(TF = X1, ensembl_gene_id = X2, input = X3, time = X4)
Day0

# load DEGs from the KD per timepoint
KDADay1.DEG = read_delim("14042017-A-AHR-KD-D1-p0.1.txt", delim = "\t")
KDADay1.DEG

KDODay1.DEG = read_delim("14042017-O-AHR-KD-D1-p0.1.txt", delim = "\t")
KDODay1.DEG

KDDay0.DEG = read_delim("14042017-ST2-AHR-KD-D1-FDR0.1.txt", delim = "\t")
KDDay0.DEG

# load all GEnes per timepoint
KDADay1 = read_delim("08062017-A-AHR-KD-AllGenes.txt", delim = "\t")
KDADay1

KDODay1 = read_delim("08062017-O-AHR-KD-AllGenes.txt", delim = "\t")
KDODay1

KDDay0 = read_delim("08062017-ST2-AHR-KD-AllGenes.txt", delim = "\t")
KDDay0

# load affinitites per time point
AffDay0 = read_delim("ST2-D0-095_Decay_Affinity_Gene_View.txt.gz", delim = "\t")
AffDay0

# Select only affinities for AHR_ARNT and the top 1000
AffDay0.sort.AHR_ARNT.top1000 = AffDay0 %>%
  select(geneID, `AHR::ARNT`) %>%
  arrange(desc(`AHR::ARNT`)) %>%
  slice(1:1000)
  
AffDay0.sort.AHR_ARNT.top1000

# Select only affinities for AHR_ARNT and the top 200
AffDay0.sort.AHR_ARNT.top200 = AffDay0 %>%
  select(geneID, `AHR::ARNT`) %>%
  arrange(desc(`AHR::ARNT`)) %>%
  slice(1:200)

AffDay0.sort.AHR_ARNT.top200

# Do the same but this time for AHR alone (top 1000)

AffDay0.sort.AHR.top1000 = AffDay0 %>%
  select(geneID, AHR) %>%
  arrange(desc(AHR)) %>%
  slice(1:1000)

AffDay0.sort.AHR.top1000

# Do the same but this time for AHR alone (top 200)

AffDay0.sort.AHR.top200 = AffDay0 %>%
  select(geneID, AHR) %>%
  arrange(desc(AHR)) %>%
  slice(1:200)

AffDay0.sort.AHR.top200

## Affinities for adipo day1
AffA.D1 = read_delim("A-D1-095_Decay_Affinity_Gene_View.txt.gz", delim = "\t")
AffA.D1

## Take the top 200 for AHR::ARNT
AffA.D1.sort.AHR_ARNT.top200 = AffA.D1 %>%
  select(geneID, `AHR::ARNT`) %>%
  arrange(desc(`AHR::ARNT`)) %>%
  slice(1:200)

AffA.D1.sort.AHR_ARNT.top200

## Take the top 200 for AHR
AffA.D1.sort.AHR.top200 = AffA.D1 %>%
  select(geneID, AHR) %>%
  arrange(desc(AHR)) %>%
  slice(1:200)

AffA.D1.sort.AHR.top200

## Affinities for osteo day1
AffO.D1 = read_delim("O-D1-095_Decay_Affinity_Gene_View.txt.gz", delim = "\t")
AffO.D1

## Take the top 200 for AHR::ARNT
AffO.D1.sort.AHR_ARNT.top200 = AffO.D1 %>%
  select(geneID, `AHR::ARNT`) %>%
  arrange(desc(`AHR::ARNT`)) %>%
  slice(1:200)

AffO.D1.sort.AHR_ARNT.top200

## Take the top 200 for AHR
AffO.D1.sort.AHR.top200 = AffO.D1 %>%
  select(geneID, AHR) %>%
  arrange(desc(AHR)) %>%
  slice(1:200)

AffO.D1.sort.AHR.top200

#plot cumulative distribution dor ST2 day0 and calulcate p-value
pdf("14062017-CumPlot-ST2.pdf", width = 8.0, height = 8.0)
plot(ecdf(KDDay0$log2FoldChange),
     xlim = c(-0.3, 0.3), 
     verticals = TRUE, 
     do.p = FALSE,
     lwd = 2.7,
     main = "",
     ylab = "",xlab = "",
     cex.axis = 2.0)  ## Plot all genes

lines(as.list(environment(ecdf(KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0$ensembl_gene_id),]$log2FoldChange))),
      col = "blue",
      lwd = 2.7,
      lty = 1)  ## Plot genes predicted as AHR targets

# lines(as.list(environment(ecdf(KDDay0[match(unique(AffDay0.sort.AHR_ARNT.top200$geneID), KDDay0$ensembl_gene_id),]$log2FoldChange))),
#       col = "red",
#       lwd = 2.7,
#       lty = 1)  ## Plot genes with the top 200 affinities for AHR_ARNT

lines(as.list(environment(ecdf(KDDay0[match(unique(AffDay0.sort.AHR.top200$geneID), KDDay0$ensembl_gene_id),]$log2FoldChange))),
      col = "#8B1A1A",
      lwd = 2.7,
      lty = 1)  ## Plot genes with the top 200 affinities for AHR alone
dev.off()
# test = Day0 %>%
#   distinct(ensembl_gene_id) %>%
#   filter(!ensembl_gene_id %in% KDDay0.DEG$ensembl_gene_id)
  
#compute significance values for comparison of ecdfs
KS.Day0.AllGenesvsAHRpredicted = ks.test(KDDay0$log2FoldChange, KDDay0[match(unique(Day0$ensembl_gene_id), 
                                                               KDDay0$ensembl_gene_id),]$log2FoldChange)
KS.Day0.AllGenesvstop200AHR_ARNT = ks.test(KDDay0$log2FoldChange, 
                                           KDDay0[match(unique(AffDay0.sort.AHR_ARNT.top200$geneID), KDDay0$ensembl_gene_id),]$log2FoldChange)

KS.Day0.AllGenesvstop200AHR = ks.test(KDDay0$log2FoldChange, 
                                      KDDay0[match(unique(AffDay0.sort.AHR.top200$geneID), KDDay0$ensembl_gene_id),]$log2FoldChange)


#compute numbers for legend
enter = c(nrow(KDDay0),
          sum(!is.na(match(unique(Day0$ensembl_gene_id), KDDay0$ensembl_gene_id))),
          sum(!is.na(match(unique(AffDay0.sort.AHR_ARNT.top200$geneID), KDDay0$ensembl_gene_id))),
          sum(!is.na(match(unique(AffDay0.sort.AHR.top200$geneID), KDDay0$ensembl_gene_id))))

legend("topleft", legend = c(paste("all genes"), 
                             paste("AHR predicted targets genes"),
                             paste("Top 200 genes with highest affinities for AHR")),
                             # paste("p-value : ", KS.Day0.AllGenesvsAHRpredicted$p.value), sep = ""),
       lty = c(1, 1, 1), 
       col = c("black", "blue", "#8B1A1A"),
       lwd = 2.7, bty = "n")
# dev.off()
#plot cumulative distribution dor adipo day 1 and calulcate p-value
pdf("14062017-CumPlot-A.pdf", width = 8.0, height = 8.0)
plot(ecdf(KDADay1$log2FoldChange),
     xlim = c(-0.3, 0.3), 
     verticals = TRUE, 
     do.p = FALSE,
     lwd = 2.5,
     main = "",
     ylab = "",xlab = "",
     cex.axis = 2.0)

lines(as.list(environment(ecdf(KDADay1[match(unique(ADay1$ensembl_gene_id),KDADay1$ensembl_gene_id),]$log2FoldChange))),
      col = "blue",
      lwd = 2.5,
      lty = 1)

# lines(as.list(environment(ecdf(KDADay1[match(unique(AffA.D1.sort.AHR_ARNT.top200$geneID), KDADay1$ensembl_gene_id),]$log2FoldChange))),
#       col = "red",
#       lwd = 2.5,
#       lty = 1)  ## Plot genes with the top 200 affinities for AHR_ARNT

lines(as.list(environment(ecdf(KDADay1[match(unique(AffA.D1.sort.AHR.top200$geneID), KDADay1$ensembl_gene_id),]$log2FoldChange))),
      col = "#8B1A1A",
      lwd = 2.5,
      lty = 1)  ## Plot genes with the top 200 affinities for AHR alone
dev.off()
#compute significance values for comparison of ecdfs
KS.A.D1.AllGenesvsAHRpredicted = ks.test(KDADay1$log2FoldChange, KDADay1[match(unique(ADay1$ensembl_gene_id), 
                                                                               KDADay1$ensembl_gene_id),]$log2FoldChange)

KS.A.D1.AllGenesvsAHR_ARNTtop200 = ks.test(KDADay1$log2FoldChange, KDADay1[match(unique(AffA.D1.sort.AHR_ARNT.top200$geneID), 
                                                                                 KDADay1$ensembl_gene_id),]$log2FoldChange)
KS.A.D1.AllGenesvsAHRtop200 = ks.test(KDADay1$log2FoldChange, KDADay1[match(unique(AffA.D1.sort.AHR.top200$geneID), 
                                                                            KDADay1$ensembl_gene_id),]$log2FoldChange)

#compute numbers for legend
enter.A = c(nrow(KDADay1),
          sum(!is.na(match(unique(ADay1$ensembl_gene_id), KDADay1$ensembl_gene_id))),
          sum(!is.na(match(unique(AffA.D1.sort.AHR_ARNT.top200$geneID), KDADay1$ensembl_gene_id))),
          sum(!is.na(match(unique(AffA.D1.sort.AHR.top200$geneID), KDADay1$ensembl_gene_id))))

legend("topleft", legend = c(paste("all genes"), 
                             paste("AHR predicted targets"),
                             paste("Top 200 genes with highest affinities for AHR")),
                             # paste("p-value < 2.2e-16", sep = "")),
       lty = c(1, 1, 1), 
       col = c("black", "blue", "#8B1A1A"),
       lwd = 2.7, bty = "n")  ## Nedd to mark the pvalue manually because otherwise return 0

# dev.off()
#plot cumulative distribution dor osteo day 1 and calulcate p-value
pdf("14062017-CumPlot-O.pdf", width = 8.0, height = 8.0)
plot(ecdf(KDODay1$log2FoldChange),
     xlim = c(-0.3, 0.3), 
     verticals = TRUE, 
     do.p = FALSE,
     lwd = 2.5,
     main = "",
     ylab = "",xlab = "",
     cex.axis = 2.0)

lines(as.list(environment(ecdf(KDODay1[match(unique(ODay1$ensembl_gene_id),KDODay1$ensembl_gene_id),]$log2FoldChange))),
      col = "blue",
      lwd = 2.5,
      lty = 1)


# lines(as.list(environment(ecdf(KDODay1[match(unique(AffO.D1.sort.AHR_ARNT.top200$geneID), KDODay1$ensembl_gene_id),]$log2FoldChange))),
#       col = "red",
#       lwd = 2.5,
#       lty = 1)

lines(as.list(environment(ecdf(KDODay1[match(unique(AffO.D1.sort.AHR.top200$geneID), KDODay1$ensembl_gene_id),]$log2FoldChange))),
      col = "#8B1A1A",
      lwd = 2.5,
      lty = 1)
dev.off()
#compute significance values for comparison of ecdfs
KS.O.D1.AllGenesvsAHRpredicted = ks.test(KDODay1$log2FoldChange, 
                                         KDODay1[match(unique(ODay1$ensembl_gene_id), KDODay1$ensembl_gene_id),]$log2FoldChange)

KS.O.D1.AllGenesvsAHR_ARNTtop200 = ks.test(KDODay1$log2FoldChange, 
                                           KDODay1[match(unique(AffO.D1.sort.AHR_ARNT.top200$geneID), KDODay1$ensembl_gene_id),]$log2FoldChange)

KS.O.D1.AllGenesvsAHRtop200 = ks.test(KDODay1$log2FoldChange, 
                                      KDODay1[match(unique(AffO.D1.sort.AHR.top200$geneID), KDODay1$ensembl_gene_id),]$log2FoldChange)

#compute numbers for legend
enter.O = c(nrow(KDODay1),
            sum(!is.na(match(unique(ODay1$ensembl_gene_id), KDODay1$ensembl_gene_id))),
            sum(!is.na(match(unique(AffO.D1.sort.AHR_ARNT.top200$geneID), KDODay1$ensembl_gene_id))),
            sum(!is.na(match(unique(AffO.D1.sort.AHR.top200$geneID), KDODay1$ensembl_gene_id))))

legend("topleft", legend = c(paste("all genes"), 
                             paste("AHR predicted targets"),
                             paste("Top 200 genes with highest affinities for AHR")),
                             # paste("p-value : ", KS.O.D1.AllGenesvsAHRpredicted$p.value), sep = ""),
       lty = c(1, 1, 1), 
       col = c("black", "blue", "#8B1A1A"),
       lwd = 2.7, bty = "n")
# dev.off()

# lines(as.list(environment(ecdf(KDDay0.DEG[match(unique(Day0$ensembl_gene_id),KDDay0.DEG$ensembl_gene_id),]$log2FoldChange))),col="orange",lwd=2.7,lty=1) ## 444
# lines(as.list(environment(ecdf(KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0.DEG$ensembl_gene_id),]$log2FoldChange))),col="pink",lwd=2.7,lty=1) ## 444
# lines(as.list(environment(ecdf(genes170$log2FoldChange))),col="pink",lwd=2.7,lty=1) ## 170


## Plot DEG from all 3 KDs as boxplot => ugly on the cumulative plot ##
# ST2 Day0
sum(!is.na(match(unique(Day0$ensembl_gene_id),KDDay0.DEG$ensembl_gene_id))) ## 444 genes are correctly predicted by EpiC-DREM out of 614 DEG in in the KD for ST2 ##

Genes444 = KDDay0.DEG[match(unique(Day0$ensembl_gene_id),KDDay0.DEG$ensembl_gene_id),] %>% 
  filter(!is.na(ensembl_gene_id)) %>%
  arrange(log2FoldChange)## Check the 170 remaining ##
Genes444

# boxplot(Genes444$log2FoldChange, xlab = "AHR targets correctly predicted by EpiC-DREM (444) in day0", ylab = "log2FC")

Genes170 = KDDay0.DEG %>% filter(!ensembl_gene_id %in% Genes444$ensembl_gene_id)
Genes170
# boxplot(Genes170$log2FoldChange, xlab = "AHR targets not correctly predicted by EpiC-DREM (170) in day0", ylab = "log2FC")
# lmts = range(Genes444$log2FoldChange, Genes170$log2FoldChange)
# par(mfrow = c(1, 2))
# boxplot(Genes444$log2FoldChange, xlab = "AHR targets correctly predicted by EpiC-DREM (444) in day0", ylab = "log2FC", ylim = lmts)
# abline(h = 0, col = "red")
# boxplot(Genes170$log2FoldChange, xlab = "AHR targets not correctly predicted by EpiC-DREM (170) in day0", ylab = "log2FC", ylim = lmts)
# abline(h = 0, col = "red")

## Boxplot not informative, try scatterplot ##
#ST2 day0 - Combine absolute log2FC and affinities
# AffDay0 = AffDay0 %>% rename(ensembl_gene_id = geneID)
# 
# KDDay0.lin = KDDay0 %>%
#   inner_join(AffDay0, by = "ensembl_gene_id") %>%
#   select(ensembl_gene_id, log2FoldChange, `AHR::ARNT`, AHR) %>%
#   arrange(desc(`AHR::ARNT`)) 
# KDDay0.lin
# # ST2.AHR.pred.tar = KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0$ensembl_gene_id),]
# plot(KDDay0.lin$lin_FC, KDDay0.lin$`AHR::ARNT`)

# ggplot(KDDay0, aes(x = abs(log2FoldChange))) + 
#   geom_density(alpha = 0.5, fill = "black") +
#   geom_density(data = KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0$ensembl_gene_id),], aes(x = abs(log2FoldChange)), fill = "blue", alpha = 0.2)
#   

## Plot the mean of absolute log2 fold change as barplot ##
# ST2
mean.abs.KDDay0 = mean(abs(KDDay0$log2FoldChange))
mean.abs.KDDay0.AHR.predTar = mean(abs(!is.na(KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0$ensembl_gene_id),]$log2FoldChange)))
df.ST2 = rownames_to_column(as.data.frame(rbind(mean.abs.KDDay0, mean.abs.KDDay0.AHR.predTar))) %>%
  rename(mean = rowname, value = V1)

pdf("13062017-CumPlot-ST2-MEAN.pdf", width = 8.0, height = 4.0)
ST2.barP = ggplot(data = df.ST2, aes(x = mean, y = value, fill = mean)) +
  coord_flip() +
  geom_bar(stat = "identity", width = 1.0) +
  theme_classic() +
  scale_fill_manual(values = c("black", "red")) +
  # scale_x_discrete(breaks = c("mean.abs.KDDay0", "mean.abs.KDDay0.AHR.predTar"),
                   # labels = c("all genes", "AHR predicted target genes")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 35, colour = "black"),
        axis.text.y = element_blank(),
        legend.position = "none") +
  ylim(0, 1)


ST2.barP
dev.off()

# Adipo day 1
mean.abs.KDADay1 = mean(abs(KDADay1$log2FoldChange))
mean.abs.KDADay1.AHR.predTar = mean(abs(!is.na(KDADay1[match(unique(ADay1$ensembl_gene_id),KDADay1$ensembl_gene_id),]$log2FoldChange)))
df.AD1 = rownames_to_column(as.data.frame(rbind(mean.abs.KDADay1, mean.abs.KDADay1.AHR.predTar))) %>%
  rename(mean = rowname, value = V1)

pdf("13062017-CumPlot-A-MEAN.pdf", width = 8.0, height = 4.0)
AD1.barP = ggplot(data = df.AD1, aes(x = mean, y = value, fill = mean)) +
  coord_flip() +
  geom_bar(stat = "identity", width = 1.0) +
  theme_classic() +
  scale_fill_manual(values = c("black", "red")) +
  scale_x_discrete(breaks = c("mean.abs.KDADay1", "mean.abs.KDADay1.AHR.predTar"),
                   labels = c("all genes", "AHR predicted target genes")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 35, colour = "black"),
        axis.text.y = element_blank(),
        legend.position = "none") +
  ylim(0, 1)


AD1.barP
dev.off()

# osteo day 1
mean.abs.KDODay1 = mean(abs(KDODay1$log2FoldChange))
mean.abs.KDODay1.AHR.predTar = mean(abs(!is.na(KDODay1[match(unique(ODay1$ensembl_gene_id),KDODay1$ensembl_gene_id),]$log2FoldChange)))
df.OD1 = rownames_to_column(as.data.frame(rbind(mean.abs.KDODay1, mean.abs.KDODay1.AHR.predTar))) %>%
  rename(mean = rowname, value = V1)

pdf("13062017-CumPlot-O-MEAN.pdf", width = 8.0, height = 4.0)
OD1.barP = ggplot(data = df.OD1, aes(x = mean, y = value, fill = mean)) +
  coord_flip() +
  geom_bar(stat = "identity", width = 1.0) +
  theme_classic() +
  scale_fill_manual(values = c("black", "red")) +
  scale_x_discrete(breaks = c("mean.abs.KDODay1", "mean.abs.KDODay1.AHR.predTar"),
                   labels = c("all genes", "AHR predicted target genes")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 35, colour = "black"),
        axis.text.y = element_blank(),
        legend.position = "none") +
  ylim(0, 1)


OD1.barP
dev.off()

## Instead of barplots, try violin plots of log2FC ##
# ST2
KDDay0.mod = KDDay0 %>% mutate(cond = "All genes")   ## All genes from day0
Day0.AHR_predTar = KDDay0[match(unique(Day0$ensembl_gene_id),KDDay0$ensembl_gene_id),] %>%
  mutate(cond = "AHR predicted target genes")

Day0.AHR_ARNT.top200 = KDDay0[match(unique(AffDay0.sort.AHR_ARNT.top200$geneID), KDDay0$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR::ARNT")

Day0.AHR.top200 = KDDay0[match(unique(AffDay0.sort.AHR.top200$geneID), KDDay0$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR")

dat.ST2 = KDDay0.mod %>% bind_rows(Day0.AHR_predTar, Day0.AHR_ARNT.top200, Day0.AHR.top200)
dat.ST2$cond = as.factor(dat.ST2$cond)

ST2.violinP = ggplot(dat.ST2, aes(x = cond, y = log2FoldChange, fill = cond)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Violin plot - ST2", x = "Condition", y = "log2FC") +
  theme_classic() +
  ylim(-2.0, 2.0) +
  scale_fill_manual(values = c("blue", "black", "#8B1A1A", "red")) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = c("AHR predicted target genes" = "AHR predicted \n target genes",
                              "All genes" = "All \n genes",
                              "Top 200 genes with highest affinities for AHR::ARNT" = "Top 200 genes with highest affinities \n for AHR::ARNT",
                              "Top 200 genes with highest affinities for AHR" = "Top 200 genes with highest affinities \n for AHR"))
pdf("14062017-ST2-ViolinPlot.pdf", width = 8.0, height = 8.0)
ST2.violinP
dev.off()

# Adipo
KDADay1.mod = KDADay1 %>% mutate(cond = "All genes")
ADay1.AHR_predTar = KDADay1[match(unique(ADay1$ensembl_gene_id),KDADay1$ensembl_gene_id),] %>%
  mutate(cond = "AHR predicted target genes")

ADay1.AHR_ARNT.top200 = KDADay1[match(unique(AffA.D1.sort.AHR_ARNT.top200$geneID), KDADay1$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR::ARNT")

ADay1.AHR.top200 = KDADay1[match(unique(AffA.D1.sort.AHR.top200$geneID), KDADay1$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR")

dat.AD1 = KDADay1.mod %>% bind_rows(ADay1.AHR_predTar, ADay1.AHR_ARNT.top200, ADay1.AHR.top200)
dat.AD1$cond = as.factor(dat.AD1$cond)

AD1.violinP = ggplot(dat.AD1, aes(x = cond, y = log2FoldChange, fill = cond)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Violin plot - adipo day1", x = "Condition", y = "log2FC") +
  theme_classic() +
  ylim(-2.0, 2.0) +
  scale_fill_manual(values = c("blue", "black", "#8B1A1A", "red")) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = c("AHR predicted target genes" = "AHR predicted \n target genes",
                              "All genes" = "All \n genes",
                              "Top 200 genes with highest affinities for AHR::ARNT" = "Top 200 genes with highest affinities \n for AHR::ARNT",
                              "Top 200 genes with highest affinities for AHR" = "Top 200 genes with highest affinities \n for AHR"))

pdf("14062017-AD1-ViolinPlot.pdf", width = 8.0, height = 8.0)
AD1.violinP
dev.off()

# Osteo
KDODay1.mod = KDODay1 %>% mutate(cond = "All genes")
ODay1.AHR_predTar = KDODay1[match(unique(ODay1$ensembl_gene_id),KDODay1$ensembl_gene_id),] %>%
  mutate(cond = "AHR predicted target genes")

ODay1.AHR_ARNT.top200 = KDODay1[match(unique(AffO.D1.sort.AHR_ARNT.top200$geneID), KDODay1$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR::ARNT")

ODay1.AHR.top200 = KDODay1[match(unique(AffO.D1.sort.AHR.top200$geneID), KDODay1$ensembl_gene_id),] %>%
  mutate(cond = "Top 200 genes with highest affinities for AHR")

dat.OD1 = KDODay1.mod %>% bind_rows(ODay1.AHR_predTar, ODay1.AHR_ARNT.top200, ODay1.AHR.top200)
dat.OD1$cond = as.factor(dat.OD1$cond)

OD1.violinP = ggplot(dat.OD1, aes(x = cond, y = log2FoldChange, fill = cond)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Violin plot - osteo day1", x = "Condition", y = "log2FC") +
  theme_classic() +
  ylim(-2.0, 2.0) +
  scale_fill_manual(values = c("blue", "black", "#8B1A1A", "red")) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = c("AHR predicted target genes" = "AHR predicted \n target genes",
                              "All genes" = "All \n genes",
                              "Top 200 genes with highest affinities for AHR::ARNT" = "Top 200 genes with highest affinities \n for AHR::ARNT",
                              "Top 200 genes with highest affinities for AHR" = "Top 200 genes with highest affinities \n for AHR"))

pdf("14062017-OD1-ViolinPlot.pdf", width = 8.0, height = 8.0)
OD1.violinP
dev.off()