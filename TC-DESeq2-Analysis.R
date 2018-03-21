## DESeq2 - Time course and knockdown analysis ##
library(readr)
library(dplyr)
library(DESeq2)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(genefilter)
library(pheatmap)
library(biomaRt)

## Count reads with featureCounts version 1.4.6-p3 ##
#featureCounts -a ./../RNA-seq-rep1/mm10fromENSEMBL/Mus_musculus.GRCm38.79.gtf -o 24022017-AllRNASeqData-CountMat.txt -F GTF -t exon -g gene_id -s 0 -Q 1 -T 12 --minReadOverlap 1 -M Mapped.Trim.H0Y3NBGXX_TC1-O-D3-2_14s006637-1-1_Sinkkonen_lane114s006637_1_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-A-D1-1_15s009483-1-1_Sinkkonen_lane115s009483_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-A-D15-3_15s009491-1-1_Sinkkonen_lane115s009491_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-A-D3-2_15s009485-1-1_Sinkkonen_lane115s009485_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-A-D5-3_15s009487-1-1_Sinkkonen_lane115s009487_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-A-D9-1_15s009489-1-1_Sinkkonen_lane115s009489_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-O-D1-3_15s009484-1-1_Sinkkonen_lane115s009484_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-O-D15-3_15s009492-1-1_Sinkkonen_lane115s009492_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-O-D3-3_15s009486-1-1_Sinkkonen_lane115s009486_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-O-D5-2_15s009488-1-1_Sinkkonen_lane115s009488_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-O-D9-3_15s009490-1-1_Sinkkonen_lane115s009490_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H5TCCBGXX_TC2-ST2-D0-3_15s009482-1-1_Sinkkonen_lane115s009482_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-A-D1-3_15s009494-1-1_Sinkkonen_lane115s009494_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-A-D15-3_15s009502-1-1_Sinkkonen_lane115s009502_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-A-D3-3_15s009496-1-1_Sinkkonen_lane115s009496_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-A-D5-2_15s009498-1-1_Sinkkonen_lane115s009498_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-A-D9-2_15s009500-1-1_Sinkkonen_lane115s009500_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-O-D1-1_15s009495-1-1_Sinkkonen_lane115s009495_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-O-D15-2_15s009503-1-1_Sinkkonen_lane115s009503_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-O-D3-2_15s009497-1-1_Sinkkonen_lane115s009497_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-O-D5-2_15s009499-1-1_Sinkkonen_lane115s009499_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-O-D9-1_15s009501-1-1_Sinkkonen_lane115s009501_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.H7LKGBGXX_TC3-ST2-D0-1_15s009493-1-1_Sinkkonen_lane115s009493_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_A_D15_2_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_A_D1_3_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_A_D3_2_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_A_D5_1_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_A_D9_2_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_O_D15_3_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_O_D1_2_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_O_D5_1_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_O_D9_1_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam Mapped.Trim.TC1_ST2_D0_1_sequence.txt.worRNAreads.txt.gzAligned.sortedByCoord.out.bam 

## Load data ##
TCData = read_delim("Counts/24022017-AllRNASeqData-CountMat.txt", delim = "\t", progress = TRUE, skip = 1)
View(TCData)

## Rearrange data to have first the time course ##
spleNamesTC = read_delim("G:/RNA-Seq-NewAnalysis/Counts/SamplesNamesTC.txt", progress = TRUE, delim = "\t", col_names = FALSE)            
TCData.Mat = column_to_rownames(TCData, var = "Geneid") %>% 
  select(-Chr, -Start, -End, -Strand, -Length)
colnames(TCData.Mat) = t(spleNamesTC)

# TCData.Mat = mutate(TCData.Mat, `ST2-D0-rep1_c` = `ST2-D0-rep1`, `ST2-D0-rep2_c` = `ST2-D0-rep2`, `ST2-D0-rep3_c` = `ST2-D0-rep3`)
TCData.Mat
TCData.Mat.Fin = select(TCData.Mat, ends_with("ST2-D0-rep1"), ends_with("ST2-D0-rep2"), ends_with("ST2-D0-rep3"), contains("A-D1-rep1"), 
                        contains("A-D1-rep2"), contains("A-D1-rep3"), contains("A-D3-rep1"), contains("A-D3-rep2"), contains("A-D3-rep3"),
                        contains("A-D5-rep1"), contains("A-D5-rep2"), contains("A-D5-rep3"), contains("A-D9-rep1"), contains("A-D9-rep2"), 
                        contains("A-D9-rep3"), contains("A-D15-rep1"), contains("A-D15-rep2"), contains("A-D15-rep3"),contains("O-D1-rep1"),
                        # matches("ST2-D0-rep1_c"), matches("ST2-D0-rep2_c"), matches("ST2-D0-rep3_c"),
                        contains("O-D1-rep2"), contains("O-D1-rep3"), contains("O-D3-rep1"), contains("O-D3-rep2"), contains("O-D3-rep3"),
                        contains("O-D5-rep1"), contains("O-D5-rep2"), contains("O-D5-rep3"), contains("O-D9-rep1"), contains("O-D9-rep2"), 
                        contains("O-D9-rep3"), contains("O-D15-rep1"), contains("O-D15-rep2"), contains("O-D15-rep3"))
TCData.Mat.Fin
colnames(TCData.Mat.Fin)

## Cretae the metadata ##
mData = data.frame(row.names = colnames(TCData.Mat.Fin), Cell = factor(c(rep("ST2", 3), rep("Ad", 15), rep("Ob", 15))),
                   Days = factor(c(rep("0", 3), rep("1", 3), rep("3", 3), rep("5", 3), rep("9", 3), rep("15", 3),
                                   rep("1", 3), rep("3", 3), rep("5", 3), rep("9", 3), rep("15", 3))),
                   replicate = factor(rep(1:3, 11)), 
                   SampleID = factor(c(rep("ST2-D0", 3), rep("A-D1", 3), rep("A-D3", 3), rep("A-D5", 3), rep("A-D9", 3), rep("A-D15", 3), 
                                       rep("O-D1", 3), rep("O-D3", 3), rep("O-D5", 3), rep("O-D9", 3), rep("O-D15", 3))))

mData


## DESeq matrix ##

ddsTC = DESeq2::DESeqDataSetFromMatrix(TCData.Mat.Fin, mData, design = ~SampleID)
colData(ddsTC)
ddsTC$SampleID = relevel(ddsTC$SampleID, "ST2-D0")
colData(ddsTC)

## Remove rows that have only 0 or 1 read ##
ddsTC = ddsTC[rowSums(counts(ddsTC)) > 1,]
head(assay(ddsTC))

## R-Log transformation ##
rldTC = rlogTransformation(ddsTC, blind = FALSE)
colData(rldTC)
head(assay(rldTC))

## PCA plot for Figure 1B ##
# PCAPlot = plotPCA(rldTC, intgroup = c("Cell", "type"))

PCAPlot.HandMade = plotPCA(rldTC, intgroup = c("Cell", "Days", "replicate"), returnData = TRUE)
percentVar = round(100 * attr(PCAPlot.HandMade, "percentVar"))

PCAPlot.HandMade$Days =  factor(PCAPlot.HandMade$Days, levels = c(0, 1, 3, 5, 9, 15))
PCAPlot.HandMade$Cell =  factor(PCAPlot.HandMade$Cell, levels = c("ST2", "Ad", "Ob"))
PCA.Diff = ggplot(PCAPlot.HandMade, aes(x = PC1, y = PC2, shape = Cell, color = Days)) + 
  geom_point(size = 5.0) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance"),
       fill = "Time") +
  theme_classic() +
  scale_shape_manual(values = c(17, 16, 15)) +
  coord_cartesian(ylim = c(-40, 30)) +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(face = "bold", size = 16, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 16, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 16),
        axis.title.y = element_text(face = "bold", size = 16),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16, face = "bold")) 

pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/24032017-PCA.Diff.pdf")
PCA.Diff
dev.off()

## Pheatmap for Figure 1D ##
# Get the gene name instead of Ensembl ID 
# Use the same reference as the one that has been used during the mapping (Mus_musculus.GRCm38.79.gtf)

listMarts(host = 'mar2015.archive.ensembl.org')
ensembl79 = useMart(host = 'mar2015.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
filters = listFilters(ensembl79)
head(filters)
attributes = listAttributes(ensembl79)
head(attributes)

topVarGenes = head(order(rowVars(assay(rldTC)),decreasing = TRUE), 100) ## Take top 100 variable genes 
mat = as.data.frame(assay(rldTC)[topVarGenes, ])
mat = mat - rowMeans(mat) ## how much a gene deviates in a specific sample from the gene average across all samples
head(mat)
ensGId.mat = rownames(mat)
qry.mat = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = ensGId.mat, mart = ensembl79)
head(qry.mat)
mat$ensembl_gene_id = ensGId.mat

NameMat.mat = dplyr::full_join(mat, qry.mat, by = "ensembl_gene_id")
head(NameMat.mat)
rownames(NameMat.mat) = NameMat.mat$external_gene_name

NameMat.mat = NameMat.mat[,-34:-35]
head(NameMat.mat, 50)
df = as.data.frame(colData(rldTC)["Days"])
df$Days = factor(df$Days, levels = c(0, 1, 3, 5, 9, 15))
Days = c("#A6CEE3" ,"#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99" ,"#E31A1C")
names(Days) = c("0", "1", "3", "5", "9", "15")
ann_colors = list(Days = Days)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/24032017-PheatmapTop100VarGenes.pdf", width = 8.27, height = 11.69,
    onefile = FALSE)
ph = pheatmap(NameMat.mat, annotation_col = df, fontsize_row = 9, cluster_cols = FALSE, scale = "row",
              color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
              annotation_colors = ann_colors, annotation_legend = TRUE, show_colnames = FALSE, cellwidth = 8,
              cellheight = 8, gaps_col = c(3, 18, 33))
dev.off()

## DEG ##
ddsTC2 = DESeq(ddsTC)
resultsNames(ddsTC2)
TC.NormCountMat = counts(ddsTC2, normalized = TRUE)
View(TC.NormCountMat)
write.table(TC.NormCountMat, "31032017-TC-NormCountMat.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Calculate FPKM for TEPIC-DREM pipeline ##
# Take the gene length from the featureCounts output
gene_Length = dplyr::select(TCData, Geneid, Length)
gene_Length
GeneEnsId = tbl_df(rownames(ddsTC)) %>%
  dplyr::rename(Geneid = value)
GeneEnsId

# Combine the gene name from above with the gene length
geneID.Length = left_join(GeneEnsId, gene_Length, by = "Geneid")

# Calculate the FPKM 
mcols(ddsTC)$basepairs = as.vector(geneID.Length$Length)
all(rownames(ddsTC) == geneID.Length$Geneid)  ## make sure that these new columns line up row for row with the DESeqDataSet
TC.FPKM.mat = fpkm(ddsTC, robust = TRUE)
View(TC.FPKM.mat)
TC.FPKM.mat = rownames_to_column(as.data.frame(TC.FPKM.mat)) %>%
  dplyr::rename(ensembl_gene_id = rowname)
head(TC.FPKM.mat)
write.table(TC.FPKM.mat, "18052017_TC_FPKM_Mat.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in adipo day 1 vs ST2 day 0
FCAD1vsD0 = list("SampleIDA.D1", "SampleIDST2.D0")
resAd1vsD0 = results(ddsTC2, contrast = FCAD1vsD0)
# save(resAd1vsD0, file = "G:/RNASeq-Ahr-Gzf1-KD/CountMat/13032017-DEG-AD1vsD0.rda")
resAd1vsD0.Sig = as.data.frame(subset(resAd1vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Ad1vsD0 = select(resAd1vsD0.Sig, ensembl_gene_id)
qry.Ad1vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Ad1vsD0, mart = ensembl79)
head(qry.Ad1vsD0)

resAd1vsD0.Sig.GS = full_join(resAd1vsD0.Sig, qry.Ad1vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resAd1vsD0.Sig.GS)

write.table(resAd1vsD0.Sig.GS, "./DEG/24032017-AD1vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in adipo day 3 vs ST2 day 0
FCAD3vsD0 = list("SampleIDA.D3", "SampleIDST2.D0")
resAd3vsD0 = results(ddsTC2, contrast = FCAD3vsD0)
resAd3vsD0.Sig = as.data.frame(subset(resAd3vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Ad3vsD0 = select(resAd3vsD0.Sig, ensembl_gene_id)
qry.Ad3vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Ad3vsD0, mart = ensembl79)
head(qry.Ad3vsD0)

resAd3vsD0.Sig.GS = full_join(resAd3vsD0.Sig, qry.Ad3vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resAd3vsD0.Sig.GS)

write.table(resAd3vsD0.Sig.GS, "./DEG/24032017-AD3vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in adipo day 5 vs ST2 day 0
FCAD5vsD0 = list("SampleIDA.D5", "SampleIDST2.D0")
resAd5vsD0 = results(ddsTC2, contrast = FCAD5vsD0)
resAd5vsD0.Sig = as.data.frame(subset(resAd5vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Ad5vsD0 = select(resAd5vsD0.Sig, ensembl_gene_id)
qry.Ad5vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Ad5vsD0, mart = ensembl79)
head(qry.Ad5vsD0)

resAd5vsD0.Sig.GS = full_join(resAd5vsD0.Sig, qry.Ad5vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resAd5vsD0.Sig.GS)

write.table(resAd5vsD0.Sig.GS, "./DEG/24032017-AD5vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in adipo day 9 vs ST2 day 0
FCAD9vsD0 = list("SampleIDA.D9", "SampleIDST2.D0")
resAd9vsD0 = results(ddsTC2, contrast = FCAD9vsD0)
resAd9vsD0.Sig = as.data.frame(subset(resAd9vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Ad9vsD0 = select(resAd9vsD0.Sig, ensembl_gene_id)
qry.Ad9vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Ad9vsD0, mart = ensembl79)
head(qry.Ad9vsD0)

resAd9vsD0.Sig.GS = full_join(resAd9vsD0.Sig, qry.Ad9vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resAd9vsD0.Sig.GS)

write.table(resAd9vsD0.Sig.GS, "./DEG/24032017-AD9vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in adipo day 15 vs ST2 day 0
FCAD15vsD0 = list("SampleIDA.D15", "SampleIDST2.D0")
resAd15vsD0 = results(ddsTC2, contrast = FCAD15vsD0)
save(resAd15vsD0, file = "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-DEG-AD15vsD0.rda")
resAd15vsD0.Sig = as.data.frame(subset(resAd15vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Ad15vsD0 = select(resAd15vsD0.Sig, ensembl_gene_id)
qry.Ad15vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                     values = ensID.Ad15vsD0, mart = ensembl79)
head(qry.Ad15vsD0)

resAd15vsD0.Sig.GS = full_join(resAd15vsD0.Sig, qry.Ad15vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resAd15vsD0.Sig.GS)

write.table(resAd15vsD0.Sig.GS, "./DEG/24032017-AD15vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in osteoday 1 vs ST2 day 0

FCOD1vsD0 = list("SampleIDO.D1", "SampleIDST2.D0")
resOd1vsD0 = results(ddsTC2, contrast = FCOD1vsD0)
# save(resOd1vsD0, file = "G:/RNASeq-Ahr-Gzf1-KD/CountMat/15032017-DEG-OD1vsD0.rda")
resOd1vsD0.Sig = as.data.frame(subset(resOd1vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Od1vsD0 = select(resOd1vsD0.Sig, ensembl_gene_id)
qry.Od1vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Od1vsD0, mart = ensembl79)
head(qry.Od1vsD0)

resOd1vsD0.Sig.GS = full_join(resOd1vsD0.Sig, qry.Od1vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resOd1vsD0.Sig.GS)

write.table(resOd1vsD0.Sig.GS, "./DEG/24032017-OD1vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in osteoday 3 vs ST2 day 0
FCOD3vsD0 = list("SampleIDO.D3", "SampleIDST2.D0")
resOd3vsD0 = results(ddsTC2, contrast = FCOD3vsD0)
resOd3vsD0.Sig = as.data.frame(subset(resOd3vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Od3vsD0 = select(resOd3vsD0.Sig, ensembl_gene_id)
qry.Od3vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Od3vsD0, mart = ensembl79)
head(qry.Od3vsD0)

resOd3vsD0.Sig.GS = full_join(resOd3vsD0.Sig, qry.Od3vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resOd3vsD0.Sig.GS)

write.table(resOd3vsD0.Sig.GS, "./DEG/24032017-OD3vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in osteoday 5 vs ST2 day 0
FCOD5vsD0 = list("SampleIDO.D5", "SampleIDST2.D0")
resOd5vsD0 = results(ddsTC2, contrast = FCOD5vsD0)
resOd5vsD0.Sig = as.data.frame(subset(resOd5vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Od5vsD0 = select(resOd5vsD0.Sig, ensembl_gene_id)
qry.Od5vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Od5vsD0, mart = ensembl79)
head(qry.Od5vsD0)

resOd5vsD0.Sig.GS = full_join(resOd5vsD0.Sig, qry.Od5vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resOd5vsD0.Sig.GS)

write.table(resOd5vsD0.Sig.GS, "./DEG/24032017-OD5vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in osteoday 9 vs ST2 day 0
FCOD9vsD0 = list("SampleIDO.D9", "SampleIDST2.D0")
resOd9vsD0 = results(ddsTC2, contrast = FCOD9vsD0)
resOd9vsD0.Sig = as.data.frame(subset(resOd9vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Od9vsD0 = select(resOd9vsD0.Sig, ensembl_gene_id)
qry.Od9vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                    values = ensID.Od9vsD0, mart = ensembl79)
head(qry.Od9vsD0)

resOd9vsD0.Sig.GS = full_join(resOd9vsD0.Sig, qry.Od9vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resOd9vsD0.Sig.GS)

write.table(resOd9vsD0.Sig.GS, "./DEG/24032017-OD9vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# DEG in osteoday 15 vs ST2 day 0
FCOD15vsD0 = list("SampleIDO.D15", "SampleIDST2.D0")
resOd15vsD0 = results(ddsTC2, contrast = FCOD15vsD0)
save(resOd15vsD0, file = "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/DEG/14042017-DEG-OD15vsD0.rda")
resOd15vsD0.Sig = as.data.frame(subset(resOd15vsD0, padj < 0.05)) %>%
  rownames_to_column(var = "ensembl_gene_id")

ensID.Od15vsD0 = select(resOd15vsD0.Sig, ensembl_gene_id)
qry.Od15vsD0 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', 
                     values = ensID.Od15vsD0, mart = ensembl79)
head(qry.Od15vsD0)

resOd15vsD0.Sig.GS = full_join(resOd15vsD0.Sig, qry.Od15vsD0, by = "ensembl_gene_id") %>%
  select(ensembl_gene_id, external_gene_name, baseMean, log2FoldChange,lfcSE, stat, pvalue, padj)
head(resOd15vsD0.Sig.GS)

write.table(resOd15vsD0.Sig.GS, "./DEG/24032017-OD15vsD0.Sig.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## Keep genes that are significicant (FDR < 5%) at any timepoints for adipo for Figure 1C ##
resAd1.any = as.data.frame(resAd1vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeAD1 = log2FoldChange, padjAD1 = padj)

resAd3.any = as.data.frame(resAd3vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeAD3 = log2FoldChange, padjAD3 = padj)

resAd5.any = as.data.frame(resAd5vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeAD5 = log2FoldChange, padjAD5 = padj)

resAd9.any = as.data.frame(resAd9vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeAD9 = log2FoldChange, padjAD9 = padj)

resAd15.any = as.data.frame(resAd15vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeAD15 = log2FoldChange, padjAD15 = padj)

resOb1.any = as.data.frame(resOd1vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeOB1 = log2FoldChange, padjOB1 = padj)

resOb3.any = as.data.frame(resOd3vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeOB3 = log2FoldChange, padjOB3 = padj)

resOb5.any = as.data.frame(resOd5vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeOB5 = log2FoldChange, padjOB5 = padj)

resOb9.any = as.data.frame(resOd9vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeOB9 = log2FoldChange, padjOB9 = padj)

resOb15.any = as.data.frame(resOd15vsD0) %>% 
  rownames_to_column(var = "ensembl_gene_id") %>% 
  select(ensembl_gene_id, log2FoldChangeOB15 = log2FoldChange, padjOB15 = padj)

DEGAll.Adipo = left_join(resAd1.any, resAd3.any, by = "ensembl_gene_id") %>% 
  left_join(resAd5.any, by = "ensembl_gene_id") %>%
  left_join(resAd9.any, by = "ensembl_gene_id") %>%
  left_join(resAd15.any, by = "ensembl_gene_id")
head(DEGAll.Adipo)

DEGAll.Adipo.padj = select(DEGAll.Adipo, ensembl_gene_id, starts_with("padj"))
head(DEGAll.Adipo.padj)  


DEGAll.Adipo.any = DEGAll.Adipo.padj[apply(DEGAll.Adipo.padj[, -1], MARGIN = 1, function(x) any(x < 0.05)), ] %>%
  na.omit()
head(DEGAll.Adipo.any)
write_delim(DEGAll.Adipo.any, "./DEG/24032017-Ad-DEGatAnyTP-FDR0.05.txt", delim = "\t", col_names = TRUE)

## Apply a cut-off on the threshold ##

DEGAll.Adipo.FCThr = filter(DEGAll.Adipo, DEGAll.Adipo$ensembl_gene_id %in% DEGAll.Adipo.any$ensembl_gene_id) %>%
  select(ensembl_gene_id,starts_with("log2"))
head(DEGAll.Adipo.FCThr) 
DEGAll.Adipo.FCThr.abs1UP = DEGAll.Adipo.FCThr[apply(DEGAll.Adipo.FCThr[, -1], MARGIN = 1, 
                                                     function(x) any(x >= 1)), ]
head(DEGAll.Adipo.FCThr.abs1UP)
write_delim(DEGAll.Adipo.FCThr.abs1UP, "./DEG/27032017-Ad-DEGatAnyTP-FDR0.05-UPFClin2.txt", delim = "\t", col_names = TRUE)


DEGAll.Adipo.FCThr.abs1DOWN = DEGAll.Adipo.FCThr[apply(DEGAll.Adipo.FCThr[, -1], MARGIN = 1, 
                                                       function(x) any(x <= -1)), ]
head(DEGAll.Adipo.FCThr.abs1DOWN)
write_delim(DEGAll.Adipo.FCThr.abs1DOWN, "./DEG/27032017-Ad-DEGatAnyTP-FDR0.05-DOWNFClin2.txt", delim = "\t", col_names = TRUE)


## Keep genes that are significicant (FDR < 5%) at any timepoints for osteo for Figure 1C ##
DEGAll.Osteo = left_join(resOb1.any, resOb3.any, by = "ensembl_gene_id") %>%
  left_join(resOb5.any, by = "ensembl_gene_id") %>%
  left_join(resOb9.any, by = "ensembl_gene_id") %>%
  left_join(resOb15.any, by = "ensembl_gene_id")
head(DEGAll.Osteo)

DEGAll.Osteo.padj = select(DEGAll.Osteo, ensembl_gene_id, starts_with("padj"))
head(DEGAll.Osteo.padj)  

DEGAll.Osteo.any = DEGAll.Osteo.padj[apply(DEGAll.Osteo.padj[, -1], MARGIN = 1, function(x) any(x < 0.05)), ] %>%
  na.omit()
head(DEGAll.Osteo.any)
write_delim(DEGAll.Osteo.any, "./DEG/24032017-Ob-DEGatAnyTP-FDR0.05.txt", delim = "\t", col_names = TRUE)

## Apply a cut-off on the threshold osteo
DEGAll.Osteo.FCThr = filter(DEGAll.Osteo, DEGAll.Osteo$ensembl_gene_id %in% DEGAll.Osteo.any$ensembl_gene_id) %>%
  select(ensembl_gene_id,starts_with("log2"))
head(DEGAll.Osteo.FCThr) 
DEGAll.Osteo.FCThr.abs1UP = DEGAll.Osteo.FCThr[apply(DEGAll.Osteo.FCThr[, -1], MARGIN = 1, 
                                                     function(x) any(x >= 1)), ]
head(DEGAll.Osteo.FCThr.abs1UP)
write_delim(DEGAll.Osteo.FCThr.abs1UP, "./DEG/27032017-Ob-DEGatAnyTP-FDR0.05-UPFClin2.txt", delim = "\t", col_names = TRUE)

DEGAll.Osteo.FCThr.abs1DOWN = DEGAll.Osteo.FCThr[apply(DEGAll.Osteo.FCThr[, -1], MARGIN = 1, 
                                                       function(x) any(x <= -1)), ]
write_delim(DEGAll.Osteo.FCThr.abs1DOWN, "./DEG/27032017-Ob-DEGatAnyTP-FDR0.05-DOWNFClin2.txt", delim = "\t", col_names = TRUE)

head(DEGAll.Osteo.FCThr.abs1DOWN)
