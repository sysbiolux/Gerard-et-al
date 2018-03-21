## Make the qPCR plot of adipogenic markers and osteoblastogenic markers ##

library("readr")
library("dplyr")
library("tibble")
library("tidyr")
library("ggplot2")
library("RColorBrewer")


dat = read_delim("09052017-TC-MarkerGenes-qPCR.txt", delim = "\t")
dat.fin = dat %>% 
  mutate(TP = rep(c(1, 3, 5, 9, 15), 12)) %>%
  mutate(Samples = rep(c("Adipo", "Adipo", "Adipo", "Adipo", "Adipo", "Osteo", "Osteo", "Osteo", "Osteo", "Osteo"), 6))
  # add_row(Samples = "Adipo", Gene = "Pparg", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
# dat.fin$TP = factor(dat.fin$TP)
dat.fin.Ad = filter(dat.fin, Samples == "Adipo")
dat.fin.Ob = filter(dat.fin, Samples == "Osteo")  

## Make the plots for the marher genes in the adipo lineage ##
## Pparg

dat.fin.Ad.Pparg = dat.fin.Ad %>%
  filter(Gene == "Pparg") %>%
  add_row(Samples = "ST2", Gene = "Pparg", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Pparg$TP = factor(dat.fin.Ad.Pparg$TP)
## Plot

Pparg.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Pparg, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Pparg, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 3.0) +
  geom_text(x = 3, y = 1.7, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 2.5, label = "**", color = "black", size = 20) +
  geom_text(x = 5, y = 2.2, label = "***", color = "black", size = 20) +
  geom_text(x = 6, y = 1.7, label = "*", color = "black", size = 20)
  


print(Pparg.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Pparg-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Pparg.qPCR.Ad.pl
dev.off()

## Cebpa

dat.fin.Ad.Cebpa = dat.fin.Ad %>%
  filter(Gene == "Cebpa") %>%
  add_row(Samples = "ST2", Gene = "Cebpa", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Cebpa$TP = factor(dat.fin.Ad.Cebpa$TP)

## Plot

Cebpa.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Cebpa, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Cebpa, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 4.5) +
  geom_text(x = 2, y = 1.6, label = "**", color = "black", size = 20) +
  geom_text(x = 4, y = 4.0, label = "**", color = "black", size = 20) +
  geom_text(x = 5, y = 4.2, label = "**", color = "black", size = 20) +
  geom_text(x = 6, y = 4.2, label = "*", color = "black", size = 20)



print(Cebpa.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Cebpa-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Cebpa.qPCR.Ad.pl
dev.off()

## Lpl

dat.fin.Ad.Lpl = dat.fin.Ad %>%
  filter(Gene == "Lpl") %>%
  add_row(Samples = "ST2", Gene = "Lpl", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Lpl$TP = factor(dat.fin.Ad.Lpl$TP)

## Plot

Lpl.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Lpl, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Lpl, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 7.5) +
  geom_text(x = 2, y = 1.2, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 7.4, label = "*", color = "black", size = 20) +
  geom_text(x = 5, y = 6.0, label = "**", color = "black", size = 20) +
  geom_text(x = 6, y = 4.0, label = "*", color = "black", size = 20)



print(Lpl.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Lpl-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Lpl.qPCR.Ad.pl
dev.off()

## Runx2

dat.fin.Ad.Runx2 = dat.fin.Ad %>%
  filter(Gene == "Runx2") %>%
  add_row(Samples = "ST2", Gene = "Runx2", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Runx2$TP = factor(dat.fin.Ad.Runx2$TP)

## Plot

Runx2.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Runx2, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Runx2, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 2, y = 1.3, label = "**", color = "black", size = 20) +
  geom_text(x = 3, y = 0.7, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 0.8, label = "**", color = "black", size = 20)



print(Runx2.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Runx2-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Runx2.qPCR.Ad.pl
dev.off()

## Sp7

dat.fin.Ad.Sp7 = dat.fin.Ad %>%
  filter(Gene == "Sp7") %>%
  add_row(Samples = "ST2", Gene = "Sp7", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Sp7$TP = factor(dat.fin.Ad.Sp7$TP)

## Plot

Sp7.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Sp7, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Sp7, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 2, y = 0.4, label = "**", color = "black", size = 20) +
  geom_text(x = 3, y = 0.5, label = "*", color = "black", size = 20)



print(Sp7.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Sp7-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Sp7.qPCR.Ad.pl
dev.off()

## Bglap

dat.fin.Ad.Bglap = dat.fin.Ad %>%
  filter(Gene == "Bglap") %>%
  add_row(Samples = "ST2", Gene = "Bglap", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ad.Bglap$TP = factor(dat.fin.Ad.Bglap$TP)

## Plot

Bglap.qPCR.Ad.pl = ggplot(data = dat.fin.Ad.Bglap, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = dat.fin.Ad.Bglap, aes(x = TP, y = Average, linetype = Gene), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 3.0) +
  geom_text(x = 2, y = 0.5, label = "***", color = "black", size = 20) +
  geom_text(x = 3, y = 0.3, label = "**", color = "black", size = 20) +
  geom_text(x = 4, y = 0.7, label = "**", color = "black", size = 20) +
  geom_text(x = 6, y = 2.9, label = "*", color = "black", size = 20) 



print(Bglap.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Bglap-Ad-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Bglap.qPCR.Ad.pl
dev.off()

## Make the plots for the marher genes in the osteo lineage ##
## Pparg

dat.fin.Ob.Pparg = dat.fin.Ob %>%
  filter(Gene == "Pparg") %>%
  add_row(Samples = "ST2", Gene = "Pparg", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Pparg$TP = factor(dat.fin.Ob.Pparg$TP)
## Plot

Pparg.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Pparg, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Pparg, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 5, y = 0.8, label = "*", color = "black", size = 20) +
  geom_text(x = 6, y = 0.5, label = "*", color = "black", size = 20)



print(Pparg.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Pparg-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Pparg.qPCR.Ob.pl
dev.off()

## Cebpa

dat.fin.Ob.Cebpa = dat.fin.Ob %>%
  filter(Gene == "Cebpa") %>%
  add_row(Samples = "ST2", Gene = "Cebpa", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Cebpa$TP = factor(dat.fin.Ob.Cebpa$TP)
## Plot

Cebpa.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Cebpa, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Cebpa, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 3, y = 0.8, label = "*", color = "black", size = 20) +
  geom_text(x = 6, y = 0.6, label = "*", color = "black", size = 20)



print(Cebpa.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Cebpa-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Cebpa.qPCR.Ob.pl
dev.off()

## Lpl

dat.fin.Ob.Lpl = dat.fin.Ob %>%
  filter(Gene == "Lpl") %>%
  add_row(Samples = "ST2", Gene = "Lpl", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Lpl$TP = factor(dat.fin.Ob.Lpl$TP)
## Plot

Lpl.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Lpl, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Lpl, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 2, y = 0.6, label = "**", color = "black", size = 20)



print(Lpl.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Lpl-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Lpl.qPCR.Ob.pl
dev.off()

## Runx2

dat.fin.Ob.Runx2 = dat.fin.Ob %>%
  filter(Gene == "Runx2") %>%
  add_row(Samples = "ST2", Gene = "Runx2", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Runx2$TP = factor(dat.fin.Ob.Runx2$TP)
## Plot

Runx2.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Runx2, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Runx2, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 1.5) +
  geom_text(x = 2, y = 1.3, label = "*", color = "black", size = 20) +
  geom_text(x = 5, y = 0.7, label = "**", color = "black", size = 20) +
  geom_text(x = 6, y = 0.4, label = "**", color = "black", size = 20)



print(Runx2.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Runx2-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Runx2.qPCR.Ob.pl
dev.off()

## Sp7

dat.fin.Ob.Sp7 = dat.fin.Ob %>%
  filter(Gene == "Sp7") %>%
  add_row(Samples = "ST2", Gene = "Sp7", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Sp7$TP = factor(dat.fin.Ob.Sp7$TP)
## Plot

Sp7.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Sp7, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Sp7, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 25) +
  geom_text(x = 2, y = 21, label = "*", color = "black", size = 20) +
  geom_text(x = 3, y = 22, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 21, label = "**", color = "black", size = 20) +
  geom_text(x = 5, y = 18, label = "*", color = "black", size = 20) +
  geom_text(x = 6, y = 9, label = "**", color = "black", size = 20)



print(Sp7.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Sp7-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Sp7.qPCR.Ob.pl
dev.off()

## Bglap

dat.fin.Ob.Bglap = dat.fin.Ob %>%
  filter(Gene == "Bglap") %>%
  add_row(Samples = "ST2", Gene = "Bglap", Average = 1, SEM = 0, `Sig?` = "NS", TP = 0, .before = 1)
dat.fin.Ob.Bglap$TP = factor(dat.fin.Ob.Bglap$TP)
## Plot

Bglap.qPCR.Ob.pl = ggplot(data = dat.fin.Ob.Bglap, aes(x = TP, y = Average, fill = Gene, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = dat.fin.Ob.Bglap, aes(x = TP, y = Average, linetype = Gene), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = "right", legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  ylim(0, 30) +
  geom_text(x = 2, y = 2.0, label = "*", color = "black", size = 20) +
  geom_text(x = 3, y = 5.5, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 11, label = "**", color = "black", size = 20) +
  geom_text(x = 5, y = 19, label = "***", color = "black", size = 20) +
  geom_text(x = 6, y = 30, label = "***", color = "black", size = 20)



print(Bglap.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig1/16032018-Bglap-Ob-qPCR-SuppFig1-biggerFont.pdf", width = 10, height = 10)
Bglap.qPCR.Ob.pl
dev.off()
