## Redo plots from the EnrichR analysis ##
library("readr")
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")

## Load the data for ChEA = TF ChIP-seq studies ##
chea.dat = read_delim("14042017-ST2.AhrKD-Ad.AhrKD-Ob.AhrKD_CommonGenes_ChEA_2016_table.txt", delim = "\t")
chea.dat.nee = select(chea.dat, Term, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  mutate(ChEA = 1:nrow(chea.dat))

## Make a horizontal barplot of the 1st 10 terms ranked by the combined score ##
chea.dat.nee$ChEA = factor(chea.dat.nee$ChEA, levels = chea.dat.nee$ChEA[order(chea.dat.nee$ChEA, decreasing = TRUE)])
# chea.dat.nee$`Adjusted P-value` = factor(chea.dat.nee$`Adjusted P-value`,
                                         # levels = chea.dat.nee$`Adjusted P-value`[order(chea.dat.nee$`Adjusted P-value`, 
                                                                                        # decreasing = TRUE)])

chea.pl = ggplot(data = chea.dat.nee %>% slice(1:5), aes(x = ChEA, y = `Combined Score`, fill = `Combined Score`)) +
  geom_bar(stat = "identity", width = 0.3) +
  theme_classic()  +
  coord_flip() +
  labs(title = "Transcription factor ChIP-Seq studies", y = "Combined enrichment score", x = "ChEA (Top 5)") +
  scale_x_discrete(breaks = c(1:5), labels = c("CEBPB_21427703_ChIP-Seq_3T3-L1_Mouse",
                                               "WT1_20215353_ChIP-ChIP_NEPHRON_PROGENITOR_Mouse",
                                               "NR1H3_23393188_ChIP-Seq_ATHEROSCLEROTIC-FOAM_Human",
                                               "PPARG_20176806_ChIP-Seq_3T3-L1_Mouse",
                                               "KDM2B_26808549_Chip-Seq_DND41_Human")) +
  scale_fill_gradient(low = "grey71", high = "chocolate") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "", legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(face = "bold", size = 40),
        plot.title = element_blank())
  #geom_text(x = 5.0, y = 2.2, label = "CEBPB - 3T3-L1", size = 10, color = "black") +
  #geom_text(x = 4.0, y = 2.9, label = "WT1 - Nephron prog.", size = 10, color = "black") +
  #geom_text(x = 3.0, y = 2.8, label = "NR1H3 - Foam cells", size = 10, color = "black") +
  #geom_text(x = 2.0, y = 2.3, label = "PPARG - 3T3-L1", size = 10, color = "black") +
  #geom_text(x = 1.0, y = 2.3, label = "KDM2B - DND41", size = 10, color = "black")
  

print(chea.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/15032018-Fig7F-ChEA.pdf", width = 10, height = 10)
chea.pl
dev.off()

## Load the data for the motifs ##
mot.dat = read_delim("14042017-ST2.AhrKD-Ad.AhrKD-Ob.AhrKD_Genome_Browser_PWMs_table.txt", delim = "\t")
mot.dat.nee = select(mot.dat, Term, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  mutate(mot = 1:nrow(mot.dat))

## Make a horizontal barplot of the 1st 10 terms ranked by the combined score ##
mot.dat.nee$mot = factor(mot.dat.nee$mot, levels = mot.dat.nee$mot[order(mot.dat.nee$mot, decreasing = TRUE)])


mot.pl = ggplot(data = mot.dat.nee %>% slice(1:5), aes(x = mot, y = `Combined Score`, fill = `Combined Score`)) +
  geom_bar(stat = "identity", width = 0.3) +
  theme_classic()  +
  coord_flip() +
  labs(title = "UCSC genome browser PWMs", y = "Combined enrichment score", x = "Transcription factor motifs (Top 5)") +
  scale_x_discrete(breaks = c(1:5), labels = c("UNKNOWN",
                                               "AHR-ARNT",
                                               "AP2GAMMA",
                                               "GABP_B",
                                               "UNKNOWN")) +
  labs(x = "Transcription factor motifs (Top 5)") +
  scale_fill_gradient(low = "grey71", high = "chocolate") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "", legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(face = "bold", size = 20),
        plot.title = element_text(hjust = 0.7))


print(mot.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/EnrichR/24042017-Fig5-motifs.pdf", width = 8, height = 8)
mot.pl
dev.off()

## Load the data for the human atlas ##
hg.at.dat = read_delim("14042017-ST2.Ahr.KD-Ad.AhrKD-Ob.AhrKD_Human_Gene_Atlas_table.txt", delim = "\t")
hg.at.dat.nee = select(hg.at.dat, Term, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  mutate(Tissue = 1:nrow(hg.at.dat))

## Make a horizontal barplot of the 1st 10 terms ranked by the combined score ##
hg.at.dat.nee$Tissue = factor(hg.at.dat.nee$Tissue, levels = hg.at.dat.nee$Tissue[order(hg.at.dat.nee$Tissue, decreasing = TRUE)])


hg.at.pl = ggplot(data = hg.at.dat.nee %>% slice(1:5), aes(x = Tissue, y = `Combined Score`, fill = `Combined Score`)) +
  geom_bar(stat = "identity", width = 0.3) +
  theme_classic()  +
  coord_flip() +
  labs(title = "Human gene atlas", y = "Combined enrichment score") +
  scale_x_discrete(breaks = c(1:5), labels = c("Smooth muscle",
                                               "Adipocyte",
                                               "CD33+ Myeloid",
                                               "Placenta",
                                               "Prostate")) +
  labs(x = "Human tissue (Top 5)") +
  scale_fill_gradient(low = "grey71", high = "chocolate") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "", legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(face = "bold", size = 40),
        plot.title = element_blank())


print(hg.at.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/15032018-Fig7F-HumanGeneAtlas.pdf", width = 10, height = 10)
hg.at.pl
dev.off()

## Load the data for the mouse atlas ##
mm.at.dat = read_delim("14042017-ST2.Ahr.KD-Ad.AhrKD-Ob.AhrKD_Mouse_Gene_Atlas_table.txt", delim = "\t")
mm.at.dat.nee = select(mm.at.dat, Term, `Combined Score`, `P-value`, `Adjusted P-value`) %>%
  mutate(Tissue = 1:nrow(mm.at.dat))

## Make a horizontal barplot of the 1st 10 terms ranked by the combined score ##
mm.at.dat.nee$Tissue = factor(mm.at.dat.nee$Tissue, levels = mm.at.dat.nee$Tissue[order(mm.at.dat.nee$Tissue, decreasing = TRUE)])

mm.at.pl = ggplot(data = mm.at.dat.nee %>% slice(1:5), aes(x = Tissue, y = `Combined Score`, fill = `Combined Score`)) +
  geom_bar(stat = "identity", width = 0.3) +
  theme_classic()  +
  coord_flip() +
  labs(title = "Mouse gene atlas", y = "Combined enrichment score") +
  scale_x_discrete(breaks = c(1:5), labels = c("Osteoblast Day 21",
                                               "Osteoblast Day 14",
                                               "Lung",
                                               "Adrenal gland",
                                               "Macrophage")) +
  labs(x = "Mouse tissue (Top 5)") +
  scale_fill_gradient(low = "grey71", high = "chocolate") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "", legend.text = element_text(size = 40, face = "bold"),
        legend.title = element_text(face = "bold", size = 40),
        plot.title = element_blank())


print(mm.at.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/15032018-Fig7F-MouseGeneAtlas.pdf", width = 10, height = 10)
mm.at.pl
dev.off()
