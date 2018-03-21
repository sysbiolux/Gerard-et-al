## RNA-seq analysis uising DEseq2 of Ahr knockdown (KD) samples ##

library(DESeq2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
setwd("G:/RNASeq-Ahr-Gzf1-KD/CountMat/")

## COunt data  - Program:featureCounts v1.4.6-p3; Command:"featureCounts" "-a" 
##"./../../../RNA-seq-rep1/mm10fromENSEMBL/Mus_musculus.GRCm38.79.gtf" "-o" "03022017-Ahr.Gzf1.ReadsCount.try2txt" "-F" "GTF" "-t" "exon" "-g" 
##"gene_id" "-s" "2" "-Q" "1" "-T" "12" "--minReadOverlap" "1" "-M" "A-siAHR-D1-rep1.bam" "A-siAHR-D1-rep2.bam" "A-siAHR-D1-rep3.bam" 
##"A-siCTRL-D1-rep1.bam" "A-siCTRL-D1-rep2.bam" "A-siCTRL-D1-rep3.bam" "A-siCTRL-D9-rep1.bam" "A-siCTRL-D9-rep2.bam" "A-siCTRL-D9-rep3.bam" 
##"A-siGZF1-D9-rep1.bam" "A-siGZF1-D9-rep2.bam" "A-siGZF1-D9-rep3.bam" "O-siAHR-D1-rep1.bam" "O-siAHR-D1-rep2.bam" "O-siAHR-D1-rep3.bam" 
##"O-siCTRL-D1-rep1.bam" "O-siCTRL-D1-rep2.bam" "O-siCTRL-D1-rep3.bam" "ST2-siAHR-D1-rep1.bam" "ST2-siAHR-D1-rep2.bam" "ST2-siAHR-D1-rep3.bam" 
##"ST2-siCTRL-D1-rep1.bam" "ST2-siCTRL-D1-rep2.bam" "ST2-siCTRL-D1-rep3.bam" ##

library(readr)
data = read_delim("03022017-Ahr.Gzf1.ReadsCount.try2.DESEQ.txt", delim = "\t", progress = TRUE) ## each row is a gene and each column is a 

#RNA library, values give the raw numbers of sequencing reads mapped to the genes in each library ##
data
dataCount = data

## Only keep the counts for the siAHR (and their matched siCTRL) samples ##
dataCount = as.data.frame(select(dataCount, starts_with("ST2-siAHR"), starts_with("ST2-siCTRL"), starts_with("A-siAHR"),
                                 starts_with("A-siCTRL-D1"), starts_with("O-siAHR"), starts_with("O-siCTRL")))
rownames(dataCount) = data$Geneid
head(dataCount)

## Provide metadata ##

mData = data.frame(row.names = colnames(dataCount), Cell = factor(c(rep("Msc",6), rep("Adipo",6), rep("Osteo", 6))), 
                   type = factor(c(rep("siAHR",3), rep("siCTRL",3), rep("siAHR",3), rep("siCTRL",3), rep("siAHR",3), 
                                   rep("siCTRL", 3))), 
                   Time = factor(c(rep("Day1", 6), rep("Day1",6), rep("Day1", 6))))
mData


## DESeq matrix ##
dds = DESeq2::DESeqDataSetFromMatrix(dataCount, mData, design = ~type)
dds
# Make sure that CTRL is the 1st level in the condition factor (R chose them alphabetically otherwise)
dds$type = relevel(dds$type, "siCTRL")
dds$group = factor(paste0(dds$Cell, dds$type, dds$Time))
design(dds) = ~group
dds$group
dds

## Remove rows that have only 0 or 1 read ##
dds = dds[rowSums(counts(dds)) > 1,]

## Visually explore the dataset ##
rld = rlogTransformation(dds, blind = FALSE) ##genes with high counts, rlog =~ log2 transformation / genes with low counts, values are shrunkrn towards the 
                                             # gene's average across all samples. NOT FOR DIFFERENTIAL TESTING!!!
head(assay(rld))

# => Genes with low counts are high variable using log2 transformation (variance increases with the mean because of poisson distribution) 
# while rlog transformation compresses differences for genes

## Run Deseq for ST2 ##
dds2 = DESeq(dds)
resultsNames(dds2)  #Check the result name to apply fold change

## DEGs in ST2 siAhr vs ST2 siCTRL ##
res.ST2.D1 = results(dds2, contrast = c("group", "MscsiAHRDay1", "MscsiCTRLDay1"))
write_delim(as.data.frame(res.ST2.D1), "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/28022017-ST2-AhrKD-AllGenes.txt", delim = "\t", col_names = TRUE)
resSig.ST2.D1 = subset(res.ST2.D1, padj < 0.1)
resSig.ST2.D1.df = as.data.frame(resSig.ST2.D1[order(resSig.ST2.D1$log2FoldChange), ])

## Add the gene symbol to the data => easier to read
library("biomaRt")
listMarts(host = 'mar2015.archive.ensembl.org') ## TO HAVE THE SAME ANNOTATION THAN THE ONE USED FOR THE TIME COURSE OF DIFFERENTIATION ##
ensembl79 = useMart(host = 'mar2015.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'mmusculus_gene_ensembl')
filters = listFilters(ensembl79)
head(filters)
attributes = listAttributes(ensembl79)
head(attributes)

ensID.ST2 = rownames(resSig.ST2.D1.df)
qry.ST2 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                values = ensID.ST2, mart = ensembl79)
head(qry.ST2)
resSig.ST2.D1.df = rownames_to_column(resSig.ST2.D1.df, var = "ensembl_gene_id")

res.ST2.D1.GS = full_join(resSig.ST2.D1.df, qry.ST2, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  dplyr::rename(geneSymbol = external_gene_name)
head(res.ST2.D1.GS)
write.table(res.ST2.D1.GS, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-ST2-AHR-KD-D1-FDR0.1.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## DEGs in A siAhr vs A siCTRL ##
res.A.D1 = results(dds2, contrast = c("group", "AdiposiAHRDay1", "AdiposiCTRLDay1"))
write_delim(as.data.frame(res.A.D1), "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/28022017-Ad-AhrKD-AllGenes.txt", delim = "\t", col_names = TRUE)
resSig.A.D1 = subset(res.A.D1, padj < 0.1)
resSig.A.D1.df = as.data.frame(resSig.A.D1[order(resSig.A.D1$log2FoldChange), ])
head(resSig.A.D1.df)

## Add the gene symbol to the data => easier to read
ensID.A = rownames(resSig.A.D1.df)
qry.A = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                values = ensID.A, mart = ensembl79)
head(qry.A)
resSig.A.D1.df = rownames_to_column(resSig.A.D1.df, var = "ensembl_gene_id")

res.A.D1.GS = full_join(resSig.A.D1.df, qry.A, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  dplyr::rename(geneSymbol = external_gene_name)
head(res.A.D1.GS)
write.table(res.A.D1.GS, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-A-AHR-KD-D1-p0.1.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## DEGs in O siAhr vs O siCTRL ##
res.O.D1 = results(dds2, contrast = c("group", "OsteosiAHRDay1", "OsteosiCTRLDay1"))
write_delim(as.data.frame(res.O.D1), "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/28022017-Ob-AhrKD-AllGenes.txt", delim = "\t", col_names = TRUE)

resSig.O.D1 = subset(res.O.D1, padj < 0.1)
resSig.O.D1.df = as.data.frame(resSig.O.D1[order(resSig.O.D1$log2FoldChange), ])
head(resSig.O.D1.df)

ensID.O = rownames(resSig.O.D1.df)
qry.O = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
              values = ensID.O, mart = ensembl79)
head(qry.O)
resSig.O.D1.df = rownames_to_column(resSig.O.D1.df, var = "ensembl_gene_id")

res.O.D1.GS = full_join(resSig.O.D1.df, qry.O, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  dplyr::rename(geneSymbol = external_gene_name)
head(res.O.D1.GS)
write.table(res.O.D1.GS, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-O-AHR-KD-D1-p0.1.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Make a scatter plot on DEG between adipo diff for 15 days and Ahr KD in ST2##
head(res.ST2.D1)  ## DEG in ST2 KD for Ahr
res.ST2.D1$ensembl_gene_id = rownames(res.ST2.D1)
res.ST2.D1.sc = as.data.frame(res.ST2.D1)
head(res.ST2.D1.sc)

load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-DEG-AD1vsD0.rda")  ## DEG in diff AD15 vs D0
ComGenes266 = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/EnrichR/18042017-CommonGenesBetweenAhrKG_266.txt", delim = "\t", col_names = FALSE) %>%
  transmute(ensembl_gene_id = X1, external_gene_name = X2)
resAd15vsD0$ensembl_gene_id = rownames(resAd15vsD0)
resAd15vsD0.sc = as.data.frame(resAd15vsD0)
head(resAd15vsD0.sc)

# Table combining fold change and FDR for DEG in both sets ##
matScST2.AhrKD = full_join(res.ST2.D1.sc, resAd15vsD0.sc, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange.x, padj.x, log2FoldChange.y, padj.y) %>%
  dplyr::rename(log2FCAhrKD.ST2 = log2FoldChange.x, padjAhrKD.ST2 = padj.x, log2FCAd15vsD0 = log2FoldChange.y,
                padjAd15vsD0 = padj.y) %>%
  mutate(log2FCAhrKD.ST2 = ifelse(is.na(log2FCAhrKD.ST2), "0", log2FCAhrKD.ST2),
         log2FCAd15vsD0 = ifelse(is.na(log2FCAd15vsD0), "0", log2FCAd15vsD0),
         padjAhrKD.ST2 = ifelse(is.na(padjAhrKD.ST2), "1", padjAhrKD.ST2),
         padjAd15vsD0 = ifelse(is.na(padjAd15vsD0), "1", padjAd15vsD0)) %>%
  mutate(Sig = ifelse(as.numeric(padjAhrKD.ST2) < 0.1 & as.numeric(padjAd15vsD0) < 0.1, "Sig in both (FDR < 0.1)",
                      ifelse(as.numeric(padjAhrKD.ST2) < 0.1 & as.numeric(padjAd15vsD0) > 0.1, "Sig in KD",
                             ifelse(as.numeric(padjAd15vsD0) < 0.1 & as.numeric(padjAhrKD.ST2) > 0.1, "Sig in diff", "No Sig"))))

head(matScST2.AhrKD)
# Take the gene symbol instead of the ensembl id
ensG = dplyr::select(matScST2.AhrKD, ensembl_gene_id)
qry.ensG = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id',
                 values = ensG, mart = ensembl79)
matScST2.AhrKD = full_join(matScST2.AhrKD, qry.ensG, by = "ensembl_gene_id")
head(matScST2.AhrKD)
## Scatterplot
library(ggplot2)
library(viridis)
library(ggrepel)

# Highlight the genes name that are significantly UP in the KD and DOWN in diff

scST2vsAD15 = ggplot(matScST2.AhrKD %>% filter(Sig != "No Sig") %>% filter(Sig != "Sig in diff") , aes(x = as.numeric(log2FCAhrKD.ST2), y = as.numeric(log2FCAd15vsD0),
                                                                                                       color = Sig)) +
  geom_point() +
  # geom_hex(bins = 100) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 1.0, color = "#A6A6A6", linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1.0,color = "#A6A6A6", linetype = "dashed") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic")) +
  scale_color_manual(values = c("red", "#E69F00")) +
  # scale_color_manual(values = c("black", "red", "#CCCCCC", "#E69F00")) +
  # viridis::scale_fill_viridis(option = "B") +
  # scale_x_continuous(limits = c(-2, 2)) +
  # scale_y_continuous(limits = c(-6, 6)) +
  labs(x = "Fold change ST2 siAhr/ST2 siCtrl (log2)",
       y = "Fold change adipocytes day 15/ST2 day 0 (log2)") +
  geom_point(data = matScST2.AhrKD %>% filter(Sig == "Sig in both (FDR < 0.1)"), colour = "red") +
  geom_point(data = matScST2.AhrKD %>% filter(Sig == "Sig in KD"), colour = "#E69F00") +
  geom_text_repel(data = matScST2.AhrKD %>% filter(Sig == "Sig in both (FDR < 0.1)") %>%
                    filter(external_gene_name == "Notch3"), size = 10, fontface = "bold.italic", color = "black", aes(label = external_gene_name), nudge_x = 0.5, nudge_y = 5, 
                  box.padding = unit(0.75, "lines")) +
  geom_text_repel(data = matScST2.AhrKD %>% filter(Sig == "Sig in both (FDR < 0.1)") %>%
                    filter(external_gene_name == "Ahr"), size = 10, fontface = "bold.italic", color = "black", aes(label = external_gene_name), nudge_x = -0.5, 
                  nudge_y = -0.5, box.padding = unit(0.75, "lines"))
# geom_point(data = matScST2.AhrKD %>% inner_join(ComGenes266, by = "external_gene_name") , color = "blue")
scST2vsAD15
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/02032018-FlyPlot-AD15vsST2KD.pdf", width = 8, height = 8)
scST2vsAD15
dev.off()

# Make a scatter plot on DEG between osteo diff for 15 days  and Ahr KD in ST2##
# head(res.ST2.D1)  ## DEG in ST2 KD for Ahr
# res.ST2.D1$ensembl_gene_id = rownames(res.ST2.D1)
# res.ST2.D1.sc = as.data.frame(res.ST2.D1)
# head(res.ST2.D1.sc)

load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-DEG-OD15vsD0.rda")  ## DEG in diff OB15 vs D0
resOd15vsD0$ensembl_gene_id = rownames(resOd15vsD0)
resOd15vsD0.sc = as.data.frame(resOd15vsD0)
head(resOd15vsD0.sc)

# Table combining fold change and FDR for DEG in both sets ##
matScST2.AhrKD.Ob = full_join(res.ST2.D1.sc, resOd15vsD0.sc, by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, log2FoldChange.x, padj.x, log2FoldChange.y, padj.y) %>%
  dplyr::rename(log2FCAhrKD.ST2 = log2FoldChange.x, padjAhrKD.ST2 = padj.x, log2FCOd15vsD0 = log2FoldChange.y,
                padjOd15vsD0 = padj.y) %>%
  mutate(log2FCAhrKD.ST2 = ifelse(is.na(log2FCAhrKD.ST2), "0", log2FCAhrKD.ST2),
         log2FCOd15vsD0 = ifelse(is.na(log2FCOd15vsD0), "0", log2FCOd15vsD0),
         padjAhrKD.ST2 = ifelse(is.na(padjAhrKD.ST2), "1", padjAhrKD.ST2),
         padjOd15vsD0 = ifelse(is.na(padjOd15vsD0), "1", padjOd15vsD0)) %>%
  mutate(Sig = ifelse(as.numeric(padjAhrKD.ST2) < 0.1 & as.numeric(padjOd15vsD0) < 0.1, "Sig in both (FDR < 0.1)",
                      ifelse(as.numeric(padjAhrKD.ST2) < 0.1 & as.numeric(padjOd15vsD0) > 0.1, "Sig in KD",
                             ifelse(as.numeric(padjOd15vsD0) < 0.1 & as.numeric(padjAhrKD.ST2) > 0.1, "Sig in diff", "No Sig"))))

head(matScST2.AhrKD.Ob)
# Take the gene symbol instead of the ensembl id
ensG = dplyr::select(matScST2.AhrKD.Ob, ensembl_gene_id)
qry.ensG = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id',
                 values = ensG, mart = ensembl79)
matScST2.AhrKD.Ob = full_join(matScST2.AhrKD.Ob, qry.ensG, by = "ensembl_gene_id")
head(matScST2.AhrKD.Ob)
## Scatterplot
library(ggplot2)
library(viridis)
library(ggrepel)

# Highlight the genes name that are significantly UP in the KD and DOWN in diff

scST2vsOb15 = ggplot(matScST2.AhrKD.Ob %>% filter(Sig != "No Sig") %>% filter(Sig != "Sig in diff"), aes(x = as.numeric(log2FCAhrKD.ST2), y = as.numeric(log2FCOd15vsD0),
                                                                                                         color = Sig)) +
  geom_point() +
  # geom_hex(bins = 100) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 1.0, color = "#A6A6A6", linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1.0, color = "#A6A6A6", linetype = "dashed") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "none", legend.text = element_text(size = 20, face = "bold.italic")) +
  scale_color_manual(values = c("red", "#E69F00")) +
  # viridis::scale_fill_viridis(option = "B") +
  # scale_x_continuous(limits = c(-2, 2)) +
  # scale_y_continuous(limits = c(-6, 6)) +
  labs(x = "Fold change ST2 siAhr/ST2 siCtrl (log2)",
       y = "Fold change osteoblasts day 15/ST2 day 0 (log2)") +
  geom_point(data = matScST2.AhrKD.Ob %>% filter(Sig == "Sig in both (FDR < 0.1)"), colour = "red") +
  geom_point(data = matScST2.AhrKD.Ob %>% filter(Sig == "Sig in KD"), colour = "#E69F00") +
  geom_text_repel(data = matScST2.AhrKD.Ob %>% filter(Sig == "Sig in both (FDR < 0.1)") %>%
                    filter(external_gene_name == "Notch3"), size = 10, fontface = "bold.italic", color = "black", 
                  aes(label = external_gene_name), nudge_x = 0.5, nudge_y = 3.0, box.padding = unit(0.75, "lines")) +
  geom_text_repel(data = matScST2.AhrKD.Ob %>% filter(Sig == "Sig in both (FDR < 0.1)") %>%
                    filter(external_gene_name == "Ahr"), size = 10, fontface = "bold.italic", color = "black", 
                  aes(label = external_gene_name), nudge_x = -0.5, nudge_y = -0.5, box.padding = unit(0.75, "lines"))
# geom_point(data = matScST2.AhrKD.Ob %>% inner_join(ComGenes266, by = "external_gene_name") , color = "blue")
# geom_label_repel(data = matScST2.AhrKD %>% filter(Sig == "Sig in both (FDR < 0.1)") %>%
# filter(as.numeric(log2FCAhrKD.ST2) >= 0.27) %>%
# filter(as.numeric(log2FCAd1vsD0) <= -1),
# aes(label = external_gene_name))
scST2vsOb15
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure7/02032018-FlyPlot-OB15vsST2KD.pdf", width = 8, height = 8)
scST2vsOb15
dev.off()
















