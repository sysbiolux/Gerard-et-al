## Plot qPCR results of marker genes in ST2 cells where Ahr was knocked down (ST2 siAHR) and Glis1 was knocked down (ST2 sGLIS1)##

library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)


# Ahr, Glis1 and marker gene expression in ST2 cells where GLIS1 is knocked-down

GLIS1.KD.ST2.qPCR = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/23022018-TF7_TF8_TF9_ST2_siGLIS1_qPCR.txt", 
                         delim = "\t") %>%
  filter(Samples == "ST2 siGLIS1 D1/ST2_NoDiff_siCTRL") %>%
  filter(Gene != "Ahr") %>%
  add_row(Samples = "ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", Gene = "empty1", Average = 0, SEM = 0, `Sig?` = "NS", .after = 1) %>%  #Add a dummy sample to make a space between bars
  add_row(Samples = "ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", Gene = "empty2", Average = 0, SEM = 0, `Sig?` = "NS", .after = 5)  #Add a dummy sample to make a space between bars

GLIS1.KD.ST2.qPCR$Gene = factor(GLIS1.KD.ST2.qPCR$Gene, levels = c("Glis1", "empty1","Cebpa", "Pparg", "Lpl", "empty2", "Runx2", "Sp7", "Bglap"))

GLIS1.KD.ST2.qPCR.all.pl = ggplot(data = GLIS1.KD.ST2.qPCR, aes(x = Gene, y = Average, fill = Samples)) +
  geom_bar(stat = "identity", size = 5, width = 0.9, position = position_dodge()) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.5) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 30, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "top", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6"), name = "Condition", breaks = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL"), 
                      labels = c("ST2 siGlis1")) +
  scale_x_discrete(breaks = c("Glis1", "empty1","Cebpa", "Pparg", "Lpl", "empty2", "Runx2", "Sp7", "Bglap"),
                   labels = c("Glis1", "","Cebpa", "Pparg", "Lpl", "", "Runx2", "Sp7", "Bglap")) +
  geom_text(x = 1, y = 0.75, label = "*", color = "black", size = 20)  +
  geom_text(x = 3.15, y = 1.4, label = "**", color = "black", size = 20, angle = 90) +
  geom_text(x = 5, y = 1.2, label = "*", color = "black", size = 20) +
  geom_text(x = 9, y = 1.48, label = "*", color = "black", size = 20)

print(GLIS1.KD.ST2.qPCR.all.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/15032018-ST2-Glis1.KD.Fig7A.MarkerGenesExp.pdf", width = 10, height = 10)
GLIS1.KD.ST2.qPCR.all.pl
dev.off()

# Ahr, Glis1 and marker gene expression in ST2 cells where AHR is knocked-down

AHR.KD.ST2.qPCR = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/23022018-TF7_TF8_TF9_ST2_siGLIS1_qPCR.txt", 
                               delim = "\t") %>%
  filter(Samples == "ST2 siAHR D1/ST2_NoDiff_siCTRL") %>%
  filter(Gene != "Glis1") %>%
  add_row(Samples = "ST2 siAHR D1/ST2_NoDiff_siCTRL", Gene = "empty1", Average = 0, SEM = 0, `Sig?` = "NS", .after = 1) %>%  #Add a dummy sample to make a space between bars
  add_row(Samples = "ST2 siAHR D1/ST2_NoDiff_siCTRL", Gene = "empty2", Average = 0, SEM = 0, `Sig?` = "NS", .after = 5)  #Add a dummy sample to make a space between bars

AHR.KD.ST2.qPCR$Gene = factor(AHR.KD.ST2.qPCR$Gene, levels = c("Ahr", "empty1","Cebpa", "Pparg", "Lpl", "empty2", "Runx2", "Sp7", "Bglap"))

AHR.KD.ST2.qPCR.all.pl = ggplot(data = AHR.KD.ST2.qPCR, aes(x = Gene, y = Average, fill = Samples)) +
  geom_bar(stat = "identity", size = 5, width = 0.9, position = position_dodge()) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 2.0) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 30, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "top", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#808080"), name = "Condition", breaks = c("ST2 siAHR D1/ST2_NoDiff_siCTRL"), 
                    labels = c("ST2 siAhr")) +
  scale_x_discrete(breaks = c("Ahr", "empty1","Cebpa", "Pparg", "Lpl", "empty2", "Runx2", "Sp7", "Bglap"),
                   labels = c("Ahr", "","Cebpa", "Pparg", "Lpl", "", "Runx2", "Sp7", "Bglap")) +
  geom_text(x = 1, y = 0.6, label = "*", color = "black", size = 20)

print(AHR.KD.ST2.qPCR.all.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/15032018-ST2-Ahr.KD.Fig7A.MarkerGenesExp.pdf", width = 10, height = 10)
AHR.KD.ST2.qPCR.all.pl
dev.off()

## Make the plot of Glis1 expression in Ahr knockdown (St2, adipocyte and osteoblast) ##

Glis1InAHR.KD = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/28022018-TF1_TF2_TF3_ST2&Ad&Ob_siAHR_Glis1Exp_qPCR.txt", 
                           delim = "\t", progress = TRUE)

Glis1InAHR.KD$Samples = factor(Glis1InAHR.KD$Samples, levels = c("ST2 siAHR D1/ST2_NoDiff_siCTRL", "Adipo siAHR D1/Adipo_NoDiff_siCTRL","Osteo siAHR D1/Osteo_NoDiff_siCTRL"))

Glis1InAHR.KD.pl = ggplot(data = Glis1InAHR.KD, aes(x = Samples, y = Average, fill = Samples)) +
  geom_bar(stat = "identity", size = 5, width = 0.6, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 3.5) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 40, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6", "#FF6666", "#0066CC"), breaks = c("ST2 siAHR D1/ST2_NoDiff_siCTRL", "Adipo siAHR D1/Adipo_NoDiff_siCTRL","Osteo siAHR D1/Osteo_NoDiff_siCTRL"), 
                    labels = c("ST2 + siAhr", "Ad + siAhr", "Ob + siAhr")) +
  scale_x_discrete(breaks = c("ST2 siAHR D1/ST2_NoDiff_siCTRL", "Adipo siAHR D1/Adipo_NoDiff_siCTRL","Osteo siAHR D1/Osteo_NoDiff_siCTRL"),
                   labels = c("ST2 + siAhr", "Ad + siAhr", "Ob + siAhr")) +
  geom_text(x = 2, y = 3.3, label = "*", color = "black", size = 20) +
  geom_text(x = 3, y = 1.8, label = "*", color = "black", size = 20)

print(Glis1InAHR.KD.pl)
pdf("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/02032018-TF1_TF2_TF3_ST2&Ad&Ob_siAHR_Glis1Exp_qPCR.pdf", width = 10, height = 10)
Glis1InAHR.KD.pl
dev.off()

## For Figure 7, plot the level of Notch1, Notch2 and Nocth3 from the knock-down of Ahr in ST2 cells, adipocyte and osteoblasts ##
Notch1 = "ENSMUSG00000026923"
Notch2 = "ENSMUSG00000027878"
Notch3 = "ENSMUSG00000038146"

Notches = c(Notch1, Notch2, Notch3)

# Load all genes from the knock-down of Ahr in ST2 cells, adipocyte and osteoblasts

ST2.AhrKD.AllGenes = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/08062017-ST2-AHR-KD-AllGenes.txt",
                                delim = "\t", progress = TRUE, col_names = TRUE)

ST2.AhrKD.Notches = ST2.AhrKD.AllGenes %>%
  filter(ensembl_gene_id %in% Notches) %>%
  mutate(linearFoldChange = 2^log2FoldChange,
         Condition = "ST2 siAhr",
         GeneName = c("Notch1", "Notch2", "Notch3"))
  

Ad.AhrKD.AllGenes = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/08062017-A-AHR-KD-AllGenes.txt",
                                delim = "\t", progress = TRUE, col_names = TRUE)

Ad.AhrKD.Notches = Ad.AhrKD.AllGenes %>%
  filter(ensembl_gene_id %in% Notches) %>%
  mutate(linearFoldChange = 2^log2FoldChange,
         Condition = "Ad siAhr",
         GeneName = c("Notch1", "Notch2", "Notch3"))

Ob.AhrKD.AllGenes = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/08062017-O-AHR-KD-AllGenes.txt",
                               delim = "\t", progress = TRUE, col_names = TRUE)

Ob.AhrKD.Notches = Ob.AhrKD.AllGenes %>%
  filter(ensembl_gene_id %in% Notches) %>%
  mutate(linearFoldChange = 2^log2FoldChange,
         Condition = "Ob siAhr",
         GeneName = c("Notch1", "Notch2", "Notch3"))

Notches.all = bind_rows(ST2.AhrKD.Notches, Ad.AhrKD.Notches, Ob.AhrKD.Notches)
Notches.all$Condition = factor(Notches.all$Condition, levels = c("ST2 siAhr", "Ad siAhr", "Ob siAhr"))

Notches.pl = ggplot(data = Notches.all, aes(x = GeneName, y = linearFoldChange, fill = Condition)) +
  geom_bar(stat = "identity", size = 5, width = 0.6, position = position_dodge()) +
  geom_errorbar(aes(ymin = 2^(log2FoldChange - lfcSE), ymax = 2^(log2FoldChange + lfcSE)), width = 0.2, size = 1.0, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 2.0) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 40, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "none", legend.text = element_text(size = 40, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6", "#FF6666", "#0066CC")) +
  geom_text(x = 2.8, y = 1.4, label = "*", color = "black", size = 20) +
  geom_text(x = 1, y = 1.4, label = "*", color = "black", size = 20) +
  geom_text(x = 2, y = 1.5, label = "0.07", color = "black", size = 10) +
  geom_text(x = 3.1, y = 1.8, label = "***", color = "black", size = 20, angle = 90) +
  geom_text(x = 3.3, y = 1.8, label = "***", color = "black", size = 20, angle = 90)

print(Notches.pl)
pdf("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/02032018-RNASeqST2&Ad&Ob_siAHR_Notches_Exp.pdf", width = 10, height = 10)
Notches.pl
dev.off()

## For Figure 7, plot the level of Notch1, Notch2 and Nocth3 from the knock-down of Glis1 in ST2 cells, adipocyte and osteoblasts ##



# Load all genes from the knock-down of Glis1 in ST2 cells, adipocyte and osteoblasts

Notches.siGlis1 = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/01032018-TF7_TF8_TF9_ST2&Ad&Ob_siGLIS1_Notches_qPCR.txt",
                                delim = "\t", progress = TRUE, col_names = TRUE)

Notches.siGlis1$Samples = factor(Notches.siGlis1$Samples, levels = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", "Ad siGLIS1 D1/Ad_siCTRL",
                                                                     "Ob siGLIS1 D1/Ob_siCTRL"))

Notches.siGLIS1.pl = ggplot(data = Notches.siGlis1, aes(x = Gene, y = Average, fill = Samples)) +
  geom_bar(stat = "identity", size = 5, width = 0.6, position = position_dodge()) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0, position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 2.0) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 40, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6", "#FF6666", "#0066CC"), breaks = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", "Ad siGLIS1 D1/Ad_siCTRL","Ob siGLIS1 D1/Ob_siCTRL"), 
                    labels = c("ST2 + siGlis1", "Ad + siGlis1", "Ob + siGlis1")) +
  geom_text(x = 0.8, y = 1.2, label = "*", color = "black", size = 20) +
  geom_text(x = 1, y = 1.6, label = "*", color = "black", size = 20) +
  geom_text(x = 1.9, y = 1.2, label = "**", color = "black", size = 20, angle = 90) +
  geom_text(x = 2.0, y = 1.8, label = "*", color = "black", size = 20) 

print(Notches.siGLIS1.pl)
pdf("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/02032018-RNASeqST2&Ad&Ob_siGLIS1_Notches_Exp.pdf", width = 10, height = 10)
Notches.siGLIS1.pl
dev.off()

## Plot Glis1 levels in ST2, adipocytes and osteoblasts samples where Glis1 is knocked down ##

Glis1.siGLIS1 = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/07032018-TF7_TF8_TF9_ST2&Ad&Ob_siGLIS1_Glis1_qPCR.txt",
                           delim = "\t", progress = TRUE)

Glis1.siGLIS1$Samples = factor(Glis1.siGLIS1$Samples, levels = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", 
                                                                 "Adipo siGLIS1/Adipo siCTRL", 
                                                                 "Osteo siGLIS1/Osteo siCTRL"))


Glis1.siGLIS1.pl = ggplot(data = Glis1.siGLIS1, aes(x = Samples, y = Average, fill = Samples)) +
  geom_bar(stat = "identity", size = 5, width = 0.3) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.0) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6", "#FF6666", "#0066CC"), breaks = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", "Adipo siGLIS1/Adipo siCTRL","Osteo siGLIS1/Osteo siCTRL")) +
  scale_x_discrete(breaks = c("ST2 siGLIS1 D1/ST2_NoDiff_siCTRL", "Adipo siGLIS1/Adipo siCTRL","Osteo siGLIS1/Osteo siCTRL"),
                   labels = c("ST2 + siGlis1", "Ad + siGlis1", "Ob + siGlis1")) +
  geom_text(x = 1.0, y = 0.75, label = "*", color = "black", size = 10) +
  geom_text(x = 3.0, y = 0.80, label = "*", color = "black", size = 10)

print(Glis1.siGLIS1.pl)
pdf("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/07032018-ST2&Ad&Ob_siGLIS1_Glis1_qPCRExp.pdf", width = 10, height = 10)
Glis1.siGLIS1.pl
dev.off()

## Plot Ahr levels in ST2, adipocytes and osteoblasts samples where Ahr is knocked down ##

Ahr.siAHR = read_delim("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/14042017-ST2.AhrKD-Ad.AhrKD-Ob.AhrKD_qPCR.txt",
                           delim = "\t", progress = TRUE)
View(Ahr.siAHR)

Ahr.siAHR$`Ahr NO DIFF` = factor(Ahr.siAHR$`Ahr NO DIFF`, levels = c("ST2siAhr/ST2siCTRL", 
                                                                 "AdsiAhr/AdsiCTRL", 
                                                                 "ObsiAhr/ObsiCTRL"))


Ahr.siAHR.pl = ggplot(data = Ahr.siAHR, aes(x = `Ahr NO DIFF`, y = Average_1, fill = `Ahr NO DIFF`)) +
  geom_bar(stat = "identity", size = 5, width = 0.3) +
  geom_errorbar(aes(ymin = Average_1 - SEM_1, ymax = Average_1 + SEM_1), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.0) +
  labs(y = "Relative expression") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold.italic", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_fill_manual(values = c("#A6A6A6", "#FF6666", "#0066CC"), breaks = c("ST2siAhr/ST2siCTRL", "AdsiAhr/AdsiCTRL","ObsiAhr/ObsiCTRL")) +
  scale_x_discrete(breaks = c("ST2siAhr/ST2siCTRL", "AdsiAhr/AdsiCTRL","ObsiAhr/ObsiCTRL"),
                   labels = c("ST2 + siAhr", "Ad + siAhr", "Ob + siAhr")) +
  geom_text(x = 1.0, y = 0.60, label = "*", color = "black", size = 10) +
  geom_text(x = 2.0, y = 0.62, label = "*", color = "black", size = 10) +
  geom_text(x = 3.0, y = 0.75, label = "0.055", color = "black", size = 10)

print(Ahr.siAHR.pl)
pdf("//atlas/FSTC_SYSBIO/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/SUPPLEMENT/Supplementary Figure S6/08032018-ST2&Ad&Ob_siAHR_Ahr_qPCRExp.pdf", width = 10, height = 10)
Ahr.siAHR.pl
dev.off()

