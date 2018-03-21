library("dplyr")
library("readr")
library("ggplot2")
library("tidyr")
library("tibble")


setwd("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/")
metaSE = read_delim("05042017-SEMapAdipoOsteoMERGED.bed", delim = "\t", col_names = FALSE)

## The last four lines are NA and it should be X #
metaSE[is.na(metaSE)] = "X"
View(metaSE)
j = c()
## Give names to the SE ##
for (i in 1:nrow(metaSE)) {
  j[i] = paste0("SE-",i)
}

## Assign names to the bed file ##
metaSE.fin = mutate(metaSE, SE_ID = j)
View(metaSE.fin)
write_delim(metaSE.fin, "05042017-SEMapAdipoOsteoMERGEDFin.bed", delim = "\t", col_names = FALSE)

## Quantification of the merged super-enhancers using annotatePeaks.pl from HOMER ##
# HOMER automatically normalizes each directory by the total number of mapped tags such that each directory contains 10 million tags.  
# This total can be changed by specifying "-norm <#>" or by specifying "-noadj" or "-raw", both of which will skip this normalization step and report integer read counts.
# The other important parameter when counting tags is to specify the size of the region you would like to count tags in with "-size <#>".  
# For example, "-size 1000" will count tags in the 1kb region centered on each peak, while "-size 50" will count tags in the 50 bp region centered on the peak (default is 200).  
# The number of tags is not normalized by the size of the region.
# -d = list of experiment directories
# -size given =  size given by actual coordinates in BED file
# -noann = skip genome annotation step, skip TSS annotation

# annotatePeaks.pl 05042017-SEMapAdipoOsteoMERGEDFin.bed mm10 -d ./../../TC1-H3K27-ST2-D0-TagDirectory/ ./../../TC1-H3K27-A-D1-TagDirectory/ ./../../TC1-H3K27-A-D3-TagDirectory/ ./../../TC1-H3K27-A-D5-TagDirectory/ ./../../TC1-H3K27-A-D15-TagDirectory/ ./../../TC1-H3K27-O-D1-TagDirectory/ ./../../TC1-H3K27-O-D3-TagDirectory/ ./../../TC1-H3K27-O-D5-TagDirectory/ ./../../TC1-H3K27-O-D9-TagDirectory/ ./../../TC1-H3K27-O-D15-TagDirectory/ -size given -noann)


## Make proper matrices for STEM ##
metaSE.Count = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoOsteoMERGEDFinCounts.txt", 
                          delim = "\t")
metaSE.Count = metaSE.Count[,c(1, 20:29)]
colnames(metaSE.Count) = c("SE_ID", "D0", "AD1", "AD3", "AD5", "AD15", "OD1", "OD3", "OD5", "OD9", "OD15")
metaSE.Count

## Matrix for SE adipo ##
metaSE.Count.Ad = select(metaSE.Count, SE_ID, D0, AD1, AD3, AD5, AD15)
metaSE.Count.Ad
write_delim(metaSE.Count.Ad, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt", 
            delim = "\t")

## MAtric for SE osteo ##
metaSE.Count.Ob = select(metaSE.Count, SE_ID, D0, OD1, OD3, OD5, OD9, OD15)
metaSE.Count.Ob
write_delim(metaSE.Count.Ob, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt", 
            delim = "\t")

## Replot the profiles in R (publication quality) for the adipo##
metaSE.STEM.Ad = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts_AllProf.txt", 
                            delim = "\t", skip = 1) %>% 
  arrange(Profile) %>%
  select(-Selected, -SPOT)
metaSE.STEM.Ad

## Extract SE by profile ##
# metaSE.STEM.Ad.Prof = gather(metaSE.STEM.Ad, D0, AD1, AD3, AD5, AD15, key = "TP", value = "FC")
# metaSE.STEM.Ad.Prof

id = factor(metaSE.STEM.Ad$Profile) %>% levels()

for (i in 1:length(id)) {
  assign(paste0("metaSE.STEM.Ad.Prof", id[i]), filter(metaSE.STEM.Ad, Profile == id[i]))
}

## Profile 1 from STEM for the adipo ##
metaSE.STEM.Ad.Prof1.fin = gather(metaSE.STEM.Ad.Prof1, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof1)) %>%
  spread(Profile, FC) %>%
  rename(FC = `1`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 24))

metaSE.STEM.Ad.Prof1.fin$TP = factor(metaSE.STEM.Ad.Prof1.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp1 = ggplot(metaSE.STEM.Ad.Prof1.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof1.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 5, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof1)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))

print(zp1)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ad-STEM-Prof1-Fig4C-biggerFont.pdf", width = 10, height = 10)
zp1
dev.off()

## Profile 3 from STEM for the adipo ##
metaSE.STEM.Ad.Prof3.fin = gather(metaSE.STEM.Ad.Prof3, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof3)) %>%
  spread(Profile, FC) %>%
  rename(FC = `3`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 16))

metaSE.STEM.Ad.Prof3.fin$TP = factor(metaSE.STEM.Ad.Prof3.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp3 = ggplot(metaSE.STEM.Ad.Prof3.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof3.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 5, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof3)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))

print(zp3)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ad-STEM-Prof3-Fig4C-biggerFont.pdf", width = 10, height = 10)
zp3
dev.off()

## Profile 4 from STEM for the adipo ##
metaSE.STEM.Ad.Prof4.fin = gather(metaSE.STEM.Ad.Prof4, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof4)) %>%
  spread(Profile, FC) %>%
  rename(FC = `4`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 16))

metaSE.STEM.Ad.Prof4.fin$TP = factor(metaSE.STEM.Ad.Prof4.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp4 = ggplot(metaSE.STEM.Ad.Prof4.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof4.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 5, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof4)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp4)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ad-STEM-Prof4-Fig4C-biggerFont.pdf", width = 10, height = 10)
zp4
dev.off()

## Profile 2 from STEM for the adipo ##
metaSE.STEM.Ad.Prof2.fin = gather(metaSE.STEM.Ad.Prof2, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof2)) %>%
  spread(Profile, FC) %>%
  rename(FC = `2`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 12))

metaSE.STEM.Ad.Prof2.fin$TP = factor(metaSE.STEM.Ad.Prof2.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp2 = ggplot(metaSE.STEM.Ad.Prof2.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof2.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 5, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof2)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp2)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ad-STEM-Prof2-Fig4C-biggerFont.pdf", width = 10, height = 10)
zp2
dev.off()

## Profile 0 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof0.fin = gather(metaSE.STEM.Ad.Prof0, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof0)) %>%
  spread(Profile, FC) %>%
  rename(FC = `0`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof0.fin$TP = factor(metaSE.STEM.Ad.Prof0.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp0 = ggplot(metaSE.STEM.Ad.Prof0.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof0.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof0)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp0)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof0-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp0
dev.off()

## Profile 6 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof6.fin = gather(metaSE.STEM.Ad.Prof6, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof6)) %>%
  spread(Profile, FC) %>%
  rename(FC = `6`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 3))

metaSE.STEM.Ad.Prof6.fin$TP = factor(metaSE.STEM.Ad.Prof6.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp6 = ggplot(metaSE.STEM.Ad.Prof6.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof6.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof6)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp6)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof6-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp6
dev.off()

## Profile 7 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof7.fin = gather(metaSE.STEM.Ad.Prof7, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof7)) %>%
  spread(Profile, FC) %>%
  rename(FC = `7`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 8))

metaSE.STEM.Ad.Prof7.fin$TP = factor(metaSE.STEM.Ad.Prof7.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp7 = ggplot(metaSE.STEM.Ad.Prof7.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof7.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof7)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp7)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof7-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp7
dev.off()

## Profile 9 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof9.fin = gather(metaSE.STEM.Ad.Prof9, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof9)) %>%
  spread(Profile, FC) %>%
  rename(FC = `9`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof9.fin$TP = factor(metaSE.STEM.Ad.Prof9.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp9 = ggplot(metaSE.STEM.Ad.Prof9.fin,
             aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof9.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof9)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp9)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof9-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp9
dev.off()

## Profile 11 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof11.fin = gather(metaSE.STEM.Ad.Prof11, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof11)) %>%
  spread(Profile, FC) %>%
  rename(FC = `11`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof11.fin$TP = factor(metaSE.STEM.Ad.Prof11.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp11 = ggplot(metaSE.STEM.Ad.Prof11.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof11.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof11)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp11)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof11-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp11
dev.off()

## Profile 12 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof12.fin = gather(metaSE.STEM.Ad.Prof12, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof12)) %>%
  spread(Profile, FC) %>%
  rename(FC = `12`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 2))

metaSE.STEM.Ad.Prof12.fin$TP = factor(metaSE.STEM.Ad.Prof12.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp12 = ggplot(metaSE.STEM.Ad.Prof12.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof12.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof12)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp12)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof12-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp12
dev.off()

## Profile 13 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof13.fin = gather(metaSE.STEM.Ad.Prof13, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof13)) %>%
  spread(Profile, FC) %>%
  rename(FC = `13`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof13.fin$TP = factor(metaSE.STEM.Ad.Prof13.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp13 = ggplot(metaSE.STEM.Ad.Prof13.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof13.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof13)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp13)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof13-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp13
dev.off()

## Profile 15 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof15.fin = gather(metaSE.STEM.Ad.Prof15, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof15)) %>%
  spread(Profile, FC) %>%
  rename(FC = `15`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 2))

metaSE.STEM.Ad.Prof15.fin$TP = factor(metaSE.STEM.Ad.Prof15.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp15 = ggplot(metaSE.STEM.Ad.Prof15.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof15.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 9, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof15)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2, 2) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp15)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof15-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp15
dev.off()

## Profile 18 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof18.fin = gather(metaSE.STEM.Ad.Prof18, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof18)) %>%
  spread(Profile, FC) %>%
  rename(FC = `18`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 7))

metaSE.STEM.Ad.Prof18.fin$TP = factor(metaSE.STEM.Ad.Prof18.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp18 = ggplot(metaSE.STEM.Ad.Prof18.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof18.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 6, y = 2.5, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof18)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.0, 3.0) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp18)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof18-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp18
dev.off()

## Profile 19 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof19.fin = gather(metaSE.STEM.Ad.Prof19, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof19)) %>%
  spread(Profile, FC) %>%
  rename(FC = `19`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof19.fin$TP = factor(metaSE.STEM.Ad.Prof19.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp19 = ggplot(metaSE.STEM.Ad.Prof19.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof19.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof19)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2.0, 2.0) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp19)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof19-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp19
dev.off()

## Profile 20 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof20.fin = gather(metaSE.STEM.Ad.Prof20, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof20)) %>%
  spread(Profile, FC) %>%
  rename(FC = `20`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof20.fin$TP = factor(metaSE.STEM.Ad.Prof20.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp20 = ggplot(metaSE.STEM.Ad.Prof20.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof20.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof20)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-2.0, 2.0) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp20)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof20-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp20
dev.off()

## Profile 23 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof23.fin = gather(metaSE.STEM.Ad.Prof23, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof23)) %>%
  spread(Profile, FC) %>%
  rename(FC = `23`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 3))

metaSE.STEM.Ad.Prof23.fin$TP = factor(metaSE.STEM.Ad.Prof23.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp23 = ggplot(metaSE.STEM.Ad.Prof23.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof23.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof23)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.0, 3.0) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp23)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof23-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp23
dev.off()

## Profile 24 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof24.fin = gather(metaSE.STEM.Ad.Prof24, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof24)) %>%
  spread(Profile, FC) %>%
  rename(FC = `24`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 3))

metaSE.STEM.Ad.Prof24.fin$TP = factor(metaSE.STEM.Ad.Prof24.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp24 = ggplot(metaSE.STEM.Ad.Prof24.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof24.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof24)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.0, 3.0) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp24)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof24-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp24
dev.off()

## Profile 25 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof25.fin = gather(metaSE.STEM.Ad.Prof25, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof25)) %>%
  spread(Profile, FC) %>%
  rename(FC = `25`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 8))

metaSE.STEM.Ad.Prof25.fin$TP = factor(metaSE.STEM.Ad.Prof25.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp25 = ggplot(metaSE.STEM.Ad.Prof25.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof25.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 3.5, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof25)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.5, 3.5) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp25)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof25-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp25
dev.off()

## Profile 26 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof26.fin = gather(metaSE.STEM.Ad.Prof26, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof26)) %>%
  spread(Profile, FC) %>%
  rename(FC = `26`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 1))

metaSE.STEM.Ad.Prof26.fin$TP = factor(metaSE.STEM.Ad.Prof26.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp26 = ggplot(metaSE.STEM.Ad.Prof26.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof26.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof26)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.5, 3.5) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp26)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof26-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp26
dev.off()

## Profile 28 from TSEM for the adipo ##

metaSE.STEM.Ad.Prof28.fin = gather(metaSE.STEM.Ad.Prof28, TP, FC, 3:ncol(metaSE.STEM.Ad.Prof28)) %>%
  spread(Profile, FC) %>%
  rename(FC = `28`) %>%
  mutate(Time = rep(c(1, 15, 3, 5, 0), 9))

metaSE.STEM.Ad.Prof28.fin$TP = factor(metaSE.STEM.Ad.Prof28.fin$TP, levels = c("D0", "AD1", "AD3", "AD5", "AD15"))

zp28 = ggplot(metaSE.STEM.Ad.Prof28.fin,
              aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ad.Prof28.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#FF6666",
            size = 1.0) +
  geom_text(x = 13, y = 2, label = paste0("n = ", nrow(metaSE.STEM.Ad.Prof28)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3.5, 3.5) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 15),
                     labels = c("0", "1", "3", "5", "15"))
print(zp28)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ad-STEM-Prof28-SuppFig5-biggerFont.pdf", width = 10, height = 10)
zp28
dev.off()

## Replot the profiles in R (publication quality) for the osteo##
metaSE.STEM.Ob = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts_AllProf.txt", 
                            delim = "\t", skip = 1) %>% 
  arrange(Profile) %>%
  select(-Selected, -SPOT)
metaSE.STEM.Ob

## Extract SE by profile ##
# metaSE.STEM.Ad.Prof = gather(metaSE.STEM.Ad, D0, AD1, AD3, AD5, AD15, key = "TP", value = "FC")
# metaSE.STEM.Ad.Prof

id = factor(metaSE.STEM.Ob$Profile) %>% levels()

for (i in 1:length(id)) {
  assign(paste0("metaSE.STEM.Ob.Prof", id[i]), filter(metaSE.STEM.Ob, Profile == id[i]))
}

## Profile 21 from STEM for the Osteo ##
metaSE.STEM.Ob.Prof21.fin = gather(metaSE.STEM.Ob.Prof21, TP, FC, 3:ncol(metaSE.STEM.Ob.Prof21)) %>%
  spread(Profile, FC) %>%
  rename(FC = `21`) %>%
  mutate(Time = rep(c(0, 1, 15, 3, 5, 9), 38))
metaSE.STEM.Ob.Prof21.fin$TP = factor(metaSE.STEM.Ob.Prof21.fin$TP, levels = c("D0", "OD1", "OD3", "OD5", "OD9", "OD15"))

Ob.zp1 = ggplot(metaSE.STEM.Ob.Prof21.fin,
                aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ob.Prof21.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#0066CC",
            size = 1.0) +
  geom_text(x = 5, y = 3, label = paste0("n = ", nrow(metaSE.STEM.Ob.Prof21)), color = "black", size = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3, 3) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 9, 15),
                     labels = c("0", "1", "3", "5", "9","15"))

print(Ob.zp1)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ob-STEM-Prof21-Fig4C-biggerFont.pdf", width = 10, height = 10)
Ob.zp1
dev.off()

## Profile 4 from STEM for the Osteo ##
metaSE.STEM.Ob.Prof4.fin = gather(metaSE.STEM.Ob.Prof4, TP, FC, 3:ncol(metaSE.STEM.Ob.Prof4)) %>%
  spread(Profile, FC) %>%
  rename(FC = `4`) %>%
  mutate(Time = rep(c(0, 1, 15, 3, 5, 9), 24))

metaSE.STEM.Ob.Prof4.fin$TP = factor(metaSE.STEM.Ob.Prof4.fin$TP, levels = c("D0", "OD1", "OD3", "OD5", "OD9", "OD15"))

Ob.zp4 = ggplot(metaSE.STEM.Ob.Prof4.fin,
                aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ob.Prof4.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#0066CC",
            size = 1.0) +
  geom_text(x = 5, y = 3, label = paste0("n = ", nrow(metaSE.STEM.Ob.Prof4)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3, 3) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 9, 15),
                     labels = c("0", "1", "3", "5", "9","15"))

print(Ob.zp4)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure4/15032018-Ob-STEM-Prof4-Fig4C-biggerFont.pdf", width = 10, height = 10)
Ob.zp4
dev.off()

## Profile 5 from STEM for the Osteo ##
metaSE.STEM.Ob.Prof5.fin = gather(metaSE.STEM.Ob.Prof5, TP, FC, 3:ncol(metaSE.STEM.Ob.Prof5)) %>%
  spread(Profile, FC) %>%
  rename(FC = `5`) %>%
  mutate(Time = rep(c(0, 1, 15, 3, 5, 9), 9))

metaSE.STEM.Ob.Prof5.fin$TP = factor(metaSE.STEM.Ob.Prof5.fin$TP, levels = c("D0", "OD1", "OD3", "OD5", "OD9", "OD15"))

Ob.zp5 = ggplot(metaSE.STEM.Ob.Prof5.fin,
                aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ob.Prof5.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#0066CC",
            size = 1.0) +
  geom_text(x = 13, y = 3, label = paste0("n = ", nrow(metaSE.STEM.Ob.Prof5)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3, 3) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 9, 15),
                     labels = c("0", "1", "3", "5", "9","15"))

print(Ob.zp5)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ob-STEM-Prof5-SuppFig5-biggerFont.pdf", width = 10, height = 10)
Ob.zp5
dev.off()

## Profile 22 from STEM for the Osteo ##
metaSE.STEM.Ob.Prof22.fin = gather(metaSE.STEM.Ob.Prof22, TP, FC, 3:ncol(metaSE.STEM.Ob.Prof22)) %>%
  spread(Profile, FC) %>%
  rename(FC = `22`) %>%
  mutate(Time = rep(c(0, 1, 15, 3, 5, 9), 2))

metaSE.STEM.Ob.Prof22.fin$TP = factor(metaSE.STEM.Ob.Prof22.fin$TP, levels = c("D0", "OD1", "OD3", "OD5", "OD9", "OD15"))

Ob.zp22 = ggplot(metaSE.STEM.Ob.Prof22.fin,
                 aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ob.Prof22.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#0066CC",
            size = 1.0) +
  geom_text(x = 13, y = 3, label = paste0("n = ", nrow(metaSE.STEM.Ob.Prof22)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3, 3) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 9, 15),
                     labels = c("0", "1", "3", "5", "9","15"))

print(Ob.zp22)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ob-STEM-Prof22-SuppFig5-biggestFont.pdf", width = 10, height = 10)
Ob.zp22
dev.off()

## Profile 27 from STEM for the Osteo ##
metaSE.STEM.Ob.Prof27.fin = gather(metaSE.STEM.Ob.Prof27, TP, FC, 3:ncol(metaSE.STEM.Ob.Prof27)) %>%
  spread(Profile, FC) %>%
  rename(FC = `27`) %>%
  mutate(Time = rep(c(0, 1, 15, 3, 5, 9), 6))

metaSE.STEM.Ob.Prof27.fin$TP = factor(metaSE.STEM.Ob.Prof27.fin$TP, levels = c("D0", "OD1", "OD3", "OD5", "OD9", "OD15"))

Ob.zp27 = ggplot(metaSE.STEM.Ob.Prof27.fin,
                 aes(x = Time, y = FC)) + 
  geom_line(data = metaSE.STEM.Ob.Prof27.fin, aes(x = Time, y = FC, group = SE_ID, linetype = "SE_ID"), color = "#0066CC",
            size = 1.0) +
  geom_text(x = 13, y = 3, label = paste0("n = ", nrow(metaSE.STEM.Ob.Prof27)), color = "black", size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(-3, 3) +
  ylab(quote(bold("Signal change ") ~(bold(log[bold("2")]~"H3K27Ac"))*"")) +
  xlab("Time [Days]") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 40, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 40, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 40),
        axis.title.y = element_text(face = "bold", size = 40),
        legend.position = "right", legend.text = element_text(size = 40, face = "bold")) +
  scale_linetype_manual(name = "", values = "solid") +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 9, 15),
                     labels = c("0", "1", "3", "5", "9","15"))

print(Ob.zp27)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Supplementary figures/SuppFig5/16032018-Ob-STEM-Prof27-SuppFig5-biggerFont.pdf", width = 10, height = 10)
Ob.zp27
dev.off()