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

## Plot time-series profile of Ahr, Glis1, Cxcl12, Notch1, Notch2, Notch3, Notch4

Ahr = "ENSMUSG00000019256"
Notch1 = "ENSMUSG00000026923"
Notch2 = "ENSMUSG00000027878"
Notch3 = "ENSMUSG00000038146"
Notch4 = "ENSMUSG00000015468"
Cxcl12 = "ENSMUSG00000061353"
Glis1 = "ENSMUSG00000034762"

# Retrieve the normalized counts for Ahr #
AhrCounts = plotCounts(ddsTC2, gene = Ahr, intgroup = c("SampleID", "Cell"), returnData = TRUE)
AhrCounts$SampleID = factor(AhrCounts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15",
                                                           "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Glis1 #
Glis1Counts = plotCounts(ddsTC2, gene = Glis1, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Glis1Counts$SampleID = factor(Glis1Counts$SampleID, levels = c("ST2-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15",
                                                               "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Notch1 #
Notch1Counts = plotCounts(ddsTC2, gene = Notch1, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Notch1Counts$SampleID = factor(Notch1Counts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15", 
                                                                 "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Notch2 #
Notch2Counts = plotCounts(ddsTC2, gene = Notch2, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Notch2Counts$SampleID = factor(Notch2Counts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15", 
                                                                 "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Notch3 #
Notch3Counts = plotCounts(ddsTC2, gene = Notch3, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Notch3Counts$SampleID = factor(Notch3Counts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15", 
                                                                 "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Notch4 #
Notch4Counts = plotCounts(ddsTC2, gene = Notch4, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Notch4Counts$SampleID = factor(Notch4Counts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15", 
                                                                 "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

# Retrieve the normalized counts for Cxcl12 #
Cxcl12Counts = plotCounts(ddsTC2, gene = Cxcl12, intgroup = c("SampleID", "Cell"), returnData = TRUE)
Cxcl12Counts$SampleID = factor(Cxcl12Counts$SampleID, levels = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15", 
                                                                 "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))

## Cxcl12 plot (mRNA and SE) only for Ad for Fig4 in the paper##
# library(scales) 
# show_col(hue_pal()(8))  ## Apply the same color at the line and the data points
Cxcl12Counts.Ad = Cxcl12Counts[c(1:18),] ## Select only the count for adipo
rel = "Msc-D0"
Cxcl12Counts.Ad.rel = Cxcl12Counts.Ad %>% 
  mutate(countRel = count/count[SampleID == rel]) %>% 
  mutate(type = "Cxcl12 mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

Cxcl12SE = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt", 
                      delim = "\t") ## Load the count for the meta SE

AdipoSENm = c("SE_ID", "Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")
colnames(Cxcl12SE) = AdipoSENm
Cxcl12SE.SE921 = filter(Cxcl12SE, SE_ID == "SE-921") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame()## Take only the counts for SE-921 which targets Cxcl12

head(Cxcl12SE.SE921)

Cxcl12SE.rel = Cxcl12SE.SE921 %>% 
  mutate(countRel.SE = `SE-921`/`SE-921`[1]) %>%
  mutate(SampleID = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = "merged SE") %>%
  as.data.frame()

# Load the Pearson correaltion for SE-921 and Cxcl12
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/10042017-PearsonCor-Cxcl12-SE921.RData")
Cxcl12_SE.921
Cxcl12_Plot_Ad = ggplot(Cxcl12Counts.Ad.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#FF6666") + 
  geom_smooth(data = Cxcl12Counts.Ad.rel, aes(x = SampleID, y = countRel,
                                              linetype = "Cxcl12 mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") + 
  geom_point(data = Cxcl12SE.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, shape = 17) +
  geom_smooth(data = Cxcl12SE.rel, aes(x = SampleID, y = countRel.SE, linetype = "merged SE"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") +
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "right", legend.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.8)) +
  scale_linetype_manual(name = "", values = c("solid", "dotdash")) +
  
  scale_x_discrete(breaks = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

Cxcl12_Plot_Ad = Cxcl12_Plot_Ad + 
  geom_text(x = 5.5, y = 1.5, label = paste0("r = ", round(Cxcl12_SE.921, digits = 2)), color = "black", size = 10)

print(Cxcl12_Plot_Ad)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Cxcl12-Ad-Fig3B.pdf", width = 8, height = 8)
Cxcl12_Plot_Ad
dev.off()

## Cxcl12 plot (mRNA and SE) only for Ob for Fig3 in the paper##
## Apply the same color at the line and the data points
Cxcl12Counts.Ob = Cxcl12Counts[c(1:3, 19:33),] ## Select only the count for adipo
rel.Ob = "Msc-D0"
Cxcl12Counts.Ob.rel = Cxcl12Counts.Ob %>% 
  mutate(countRel = count/count[SampleID == rel.Ob]) %>% 
  mutate(type = "Cxcl12 mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

Cxcl12SE.Ob = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt", 
                         delim = "\t") ## Load the count for the meta SE

OsteoSENm = c("SE_ID", "Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "O-D15")
colnames(Cxcl12SE.Ob) = OsteoSENm
Cxcl12SE.Ob.SE921 = filter(Cxcl12SE.Ob, SE_ID == "SE-921") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame()## Take only the counts for SE-889 which targets Cxcl12

head(Cxcl12SE.Ob.SE921)

Cxcl12SE.Ob.rel = Cxcl12SE.Ob.SE921 %>% 
  mutate(countRel.SE = `SE-921`/`SE-921`[1]) %>%
  mutate(SampleID = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "O-D15")) %>%
  mutate(type.SE = "merged SE") %>%
  as.data.frame()

load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/10042017-PearsonCor-Cxcl12-SE921.Ob.RData")
Cxcl12_Ob_SE.921
Cxcl12_Plot_Ob = ggplot(Cxcl12Counts.Ob.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#0066CC") + 
  geom_smooth(data = Cxcl12Counts.Ob.rel, aes(x = SampleID, y = countRel,
                                              linetype = "Cxcl12 mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") + 
  geom_point(data = Cxcl12SE.Ob.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, shape = 17) +
  geom_smooth(data = Cxcl12SE.Ob.rel, aes(x = SampleID, y = countRel.SE, linetype = "merged SE"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") +
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells",
       title = "Osteoblast differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = "right", legend.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.8)) +
  scale_linetype_manual(name = "", values = c("solid", "dotdash")) +
  
  scale_x_discrete(breaks = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

Cxcl12_Plot_Ob = Cxcl12_Plot_Ob +
  geom_text(x = 5.5, y = 1.5, label = paste0("r = ", round(Cxcl12_Ob_SE.921, digits = 2)), color = "black", size = 10)

print(Cxcl12_Plot_Ob)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Cxcl12-Ob-Fig3B.pdf", width = 8, height = 8)
Cxcl12_Plot_Ob
dev.off()

#################################### Ahr figure 4 - se and mRNA plot - Adipo ############################################
## Ahr plot (mRNA and SE) only for Ad for Fig4 in the paper##

AhrCounts.Ad = AhrCounts[c(1:18),] ## Select only the count for adipo
rel = "Msc-D0"
AhrCounts.Ad.rel = AhrCounts.Ad %>% 
  mutate(countRel = count/count[SampleID == rel]) %>% 
  mutate(type = "Ahr mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

AhrSE = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt", 
                   delim = "\t") ## Load the count for the meta SE

AdipoSENm = c("SE_ID", "Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")
colnames(AhrSE) = AdipoSENm


AhrSE.SE283 = filter(AhrSE, SE_ID == "SE-283") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-283 which targets Ahr

AhrSE.SE284 = filter(AhrSE, SE_ID == "SE-284") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-284 which targets Ahr

AhrSE.SE285 = filter(AhrSE, SE_ID == "SE-285") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-285 which targets Ahr

AhrSE.SE286 = filter(AhrSE, SE_ID == "SE-286") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-286 which targets Ahr

## Load the Pearson correlation between Ahr mRNA and their SEs ##
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-SE283.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-SE284.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-SE285.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-SE286.RData")
Ahr_SE.283
Ahr_SE.284
Ahr_SE.285
Ahr_SE.286

AhrSE.283.rel = AhrSE.SE283 %>% 
  mutate(countRel.SE = `SE-283`/`SE-283`[1]) %>%
  mutate(SampleID = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = paste("SE-283 | r =", round(Ahr_SE.283, digits = 2))) %>%
  as.data.frame()

AhrSE.284.rel = AhrSE.SE284 %>% 
  mutate(countRel.SE = `SE-284`/`SE-284`[1]) %>%
  mutate(SampleID = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = paste("SE-284 | r =", round(Ahr_SE.284, digits = 2))) %>%
  as.data.frame()

AhrSE.285.rel = AhrSE.SE285 %>% 
  mutate(countRel.SE = `SE-285`/`SE-285`[1]) %>%
  mutate(SampleID = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = paste("SE-285 | r =", round(Ahr_SE.285, digits = 2))) %>%
  as.data.frame()

AhrSE.286.rel = AhrSE.SE286 %>% 
  mutate(countRel.SE = `SE-286`/`SE-286`[1]) %>%
  mutate(SampleID = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = paste("SE-286 | r =", round(Ahr_SE.286, digits = 2))) %>%
  as.data.frame()

## Plotting

Ahr_SE_Plot_Ad = ggplot(AhrCounts.Ad.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#FF6666") + 
  geom_smooth(data = AhrCounts.Ad.rel, aes(x = SampleID, y = countRel,
                                           linetype = "Ahr mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") + 
  geom_point(data = AhrSE.283.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, 
             shape = 15) +
  geom_point(data = AhrSE.284.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, 
             shape = 16) +
  geom_point(data = AhrSE.285.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, 
             shape = 17) +
  geom_point(data = AhrSE.286.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, 
             shape = 19) +
  geom_smooth(data = AhrSE.283.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-283 | r = 0.99"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") + 
  geom_smooth(data = AhrSE.284.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-284 | r = 0.98"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") +
  geom_smooth(data = AhrSE.285.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-285 | r = 0.98"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") +
  geom_smooth(data = AhrSE.286.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-286 | r = 0.99"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") +
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5)) +
  scale_linetype_manual(name = "", values = c("solid", "dashed", "dotted", "dotdash", "longdash")) +
  
  scale_x_discrete(breaks = c("Msc-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

print(Ahr_SE_Plot_Ad)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Ahr-Ad-SE-Fig4B.pdf", width = 8, height = 8)
Ahr_SE_Plot_Ad
dev.off()

#################################### Ahr figure 4 - se and mRNA plot - Osteo ############################################
## Ahr plot (mRNA and SE) only for Ad for Fig3 in the paper##

AhrCounts.Ob = AhrCounts[c(1:3, 19:33),] ## Select only the count for adipo
rel = "Msc-D0"
AhrCounts.Ob.rel = AhrCounts.Ob %>% 
  mutate(countRel = count/count[SampleID == rel]) %>% 
  mutate(type = "Ahr mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

AhrSE = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt", 
                   delim = "\t") ## Load the count for the meta SE

OsteoSENm = c("SE_ID", "Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "A-D15")
colnames(AhrSE) = OsteoSENm


AhrSE.Ob.SE283 = filter(AhrSE, SE_ID == "SE-283") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-283 which targets Ahr

AhrSE.Ob.SE284 = filter(AhrSE, SE_ID == "SE-284") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-284 which targets Ahr

AhrSE.Ob.SE285 = filter(AhrSE, SE_ID == "SE-285") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-285 which targets Ahr

AhrSE.Ob.SE286 = filter(AhrSE, SE_ID == "SE-286") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-286 which targets Ahr

## Load the Pearson correlation between Ahr mRNA and their SEs ##
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-Ob-SE283.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-Ob-SE284.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-Ob-SE285.RData")
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/11042017-PearsonCor-Ahr-Ob-SE286.RData")
Ahr_Ob_SE.283
Ahr_Ob_SE.284
Ahr_Ob_SE.285
Ahr_Ob_SE.286

AhrSE.Ob.283.rel = AhrSE.Ob.SE283 %>% 
  mutate(countRel.SE = `SE-283`/`SE-283`[1]) %>%
  mutate(SampleID = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9","O-D15")) %>%
  mutate(type.SE = paste("SE-283 | r =", round(Ahr_Ob_SE.283, digits = 2))) %>%
  as.data.frame()

AhrSE.Ob.284.rel = AhrSE.Ob.SE284 %>% 
  mutate(countRel.SE = `SE-284`/`SE-284`[1]) %>%
  mutate(SampleID = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9","O-D15")) %>%
  mutate(type.SE = paste("SE-284 | r =", round(Ahr_Ob_SE.284, digits = 2))) %>%
  as.data.frame()

AhrSE.Ob.285.rel = AhrSE.Ob.SE285 %>% 
  mutate(countRel.SE = `SE-285`/`SE-285`[1]) %>%
  mutate(SampleID = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9","O-D15")) %>%
  mutate(type.SE = paste("SE-285 | r =", round(Ahr_Ob_SE.285, digits = 2))) %>%
  as.data.frame()

AhrSE.Ob.286.rel = AhrSE.Ob.SE286 %>% 
  mutate(countRel.SE = `SE-286`/`SE-286`[1]) %>%
  mutate(SampleID = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9","O-D15")) %>%
  mutate(type.SE = paste("SE-286 | r =", round(Ahr_Ob_SE.286, digits = 2))) %>%
  as.data.frame()

## Plotting

Ahr_SE_Plot_Ob = ggplot(AhrCounts.Ob.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#0066CC") + 
  geom_smooth(data = AhrCounts.Ob.rel, aes(x = SampleID, y = countRel,
                                           linetype = "Ahr mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") + 
  geom_point(data = AhrSE.Ob.283.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, 
             shape = 15) +
  geom_point(data = AhrSE.Ob.284.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, 
             shape = 16) +
  geom_point(data = AhrSE.Ob.285.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, 
             shape = 17) +
  geom_point(data = AhrSE.Ob.286.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, 
             shape = 19) +
  geom_smooth(data = AhrSE.Ob.283.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-283 | r = 0.98"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") + 
  geom_smooth(data = AhrSE.Ob.284.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-284 | r = 0.96"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") +
  geom_smooth(data = AhrSE.Ob.285.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-285 | r = 0.98"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") +
  geom_smooth(data = AhrSE.Ob.286.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-286 | r = 0.95"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") +
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5)) +
  scale_linetype_manual(name = "", values = c("solid", "dashed", "dotted", "dotdash", "longdash")) +
  
  scale_x_discrete(breaks = c("Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

print(Ahr_SE_Plot_Ob)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Ahr-Ob-SE-Fig4B.pdf", width = 8, height = 8)
Ahr_SE_Plot_Ob
dev.off()

#################################### Glis1 figure 4 - se and mRNA plot - Adipo ############################################
## Glis1 plot (mRNA and SE) only for Ad for Fig4 in the paper##

Glis1Counts.Ad = Glis1Counts[c(1:18),] ## Select only the count for adipo
rel = "ST2-D0"
Glis1Counts.Ad.rel = Glis1Counts.Ad %>% 
  mutate(countRel = count/count[SampleID == rel]) %>% 
  mutate(type = "Glis1 mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

Glis1SE = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt", 
                     delim = "\t") ## Load the count for the meta SE

AdipoSENm = c("SE_ID", "Msc-D0", "A-D1", "A-D3", "A-D5", "A-D15")
colnames(Glis1SE) = AdipoSENm


Glis1SE.SE831 = filter(Glis1SE, SE_ID == "SE-831") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-831 which targets Glis1

## Load the Pearson correlation between Glis1 mRNA and its SE ##
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/02022018-PearsonCor-Glis1-SE831.RData")

Glis1_SE.831

Glis1SE.831.rel = Glis1SE.SE831 %>% 
  mutate(countRel.SE = `SE-831`/`SE-831`[1]) %>%
  mutate(SampleID = c("ST2-D0", "A-D1", "A-D3", "A-D5", "A-D15")) %>%
  mutate(type.SE = paste("SE-831 | r =", round(Glis1_SE.831, digits = 2))) %>%
  as.data.frame()

## Plot together Glis1 mRNA (normalized count) and its SE

Glis1_SE_Plot_Ad = ggplot(Glis1Counts.Ad.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#FF6666") + 
  geom_smooth(data = Glis1Counts.Ad.rel, aes(x = SampleID, y = countRel,
                                             linetype = "Glis1 mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") + 
  geom_point(data = Glis1SE.831.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#FF6666", size = 5, 
             shape = 15) +
  geom_smooth(data = Glis1SE.831.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-831 | r = 0.99"), 
              method = "loess", se = FALSE, size = 1.5, color = "#FF6666") + 
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 30, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_linetype_manual(name = "", values = c("solid", "dotted")) +
  
  scale_x_discrete(breaks = c("ST2-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

print(Glis1_SE_Plot_Ad)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure6/15032018-Glis1-Ad-SE-Fig6C-biggerFont.pdf", width = 10, height = 10)
Glis1_SE_Plot_Ad
dev.off()

#################################### Glis1 figure 4 - se and mRNA plot - Osteo ############################################
## Glis1 plot (mRNA and SE) only for Ob for Fig4 in the paper##

Glis1Counts.Ob = Glis1Counts[c(1:3, 19:33),] ## Select only the count for osteo
rel = "ST2-D0"
Glis1Counts.Ob.rel = Glis1Counts.Ob %>% 
  mutate(countRel = count/count[SampleID == rel]) %>% 
  mutate(type = "Glis1 mRNA")  ## Normalize the count to Day0 and adda column specifying that it is mRNA samples

Glis1SE = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt", 
                     delim = "\t") ## Load the count for the meta SE

OsteoSENm = c("SE_ID", "Msc-D0", "O-D1", "O-D3", "O-D5", "O-D9", "A-D15")
colnames(Glis1SE) = OsteoSENm


Glis1SE.Ob.SE831 = filter(Glis1SE, SE_ID == "SE-831") %>%
  column_to_rownames("SE_ID") %>%
  t() %>%
  as_data_frame() %>% 
  head()  ## Take only the counts for SE-831 which targets Glis1


## Load the Pearson correlation between Ahr mRNA and their SEs ##
load("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/02022018-PearsonCor-Glis1-Ob-SE831.RData")

Glis1_SE.831_Ob

Glis1SE.Ob.831.rel = Glis1SE.Ob.SE831 %>% 
  mutate(countRel.SE = `SE-831`/`SE-831`[1]) %>%
  mutate(SampleID = c("ST2-D0", "O-D1", "O-D3", "O-D5", "O-D9","O-D15")) %>%
  mutate(type.SE = paste("SE-831 | r =", round(Glis1_SE.831_Ob, digits = 2))) %>%
  as.data.frame()


## Plotting


Glis1_SE_Plot_Ob = ggplot(Glis1Counts.Ob.rel, aes(x = SampleID, y = countRel, group = 1)) +
  geom_point(size = 5, color = "#0066CC") + 
  geom_smooth(data = Glis1Counts.Ob.rel, aes(x = SampleID, y = countRel,
                                             linetype = "Glis1 mRNA"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") + 
  geom_point(data = Glis1SE.Ob.831.rel, aes(x = SampleID, y = countRel.SE, group = 1), color = "#0066CC", size = 5, 
             shape = 15) +
  geom_smooth(data = Glis1SE.Ob.831.rel, aes(x = SampleID, y = countRel.SE, linetype = "SE-831 | r = 0.95"), 
              method = "loess", se = FALSE, size = 1.5, color = "#0066CC") + 
  labs(x = "Time [Days]", y = "Signal relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme_classic() +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 30, face = "bold"),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_linetype_manual(name = "", values = c("solid", "dotted")) +
  
  scale_x_discrete(breaks = c("ST2-D0", "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) + 
  ylim(0, 1.5)

print(Glis1_SE_Plot_Ob)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure6/15032018-Glis1-Ob-SE-Fig6C-biggerFont.pdf", width = 10, height = 10)
Glis1_SE_Plot_Ob
dev.off()

## qPCR plot figure 4 ##
Ahr.qPCR = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/03042017-Ahr-qPCR.txt", delim = "\t")


Ahr.qPCR.Ad = Ahr.qPCR[1:5,] %>% 
  mutate(Time = c(1, 3, 5, 9, 15), Ahr = "Ahr") %>%
  add_row(Ahr = "Ahr", Average = 1, SD = 0, SEM = 0, Time = 0, .before = 1)

Ahr.qPCR.Ad$Time = factor(Ahr.qPCR.Ad$Time)
Ahr.qPCR.Ob = Ahr.qPCR[6:10,] %>% 
  mutate(Time = c(1, 3, 5, 9, 15), Ahr = "Ahr") %>%
  add_row(Ahr = "Ahr", Average = 1, SD = 0, SEM = 0, Time = 0, .before = 1)

Ahr.qPCR.Ob$Time = factor(Ahr.qPCR.Ob$Time)

Ahr.qPCR.Ad.pl = ggplot(data = Ahr.qPCR.Ad, aes(x = Time, y = Average, fill = Ahr, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = Ahr.qPCR.Ad, aes(x = Time, y = Average, linetype = Ahr), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  geom_text(x = 2, y = 0.35, label = "**", color = "black", size = 10) +
  geom_text(x = 3, y = 0.35, label = "***", color = "black", size = 10) +
  geom_text(x = 4, y = 0.20, label = "**", color = "black", size = 10) +
  geom_text(x = 5, y = 0.25, label = "**", color = "black", size = 10) +
  geom_text(x = 6, y = 0.25, label = "**", color = "black", size = 10)


print(Ahr.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Ahr-Ab-qPCR-Fig4B.pdf", width = 8, height = 8)
Ahr.qPCR.Ad.pl
dev.off()

Ahr.qPCR.Ob.pl = ggplot(data = Ahr.qPCR.Ob, aes(x = Time, y = Average, fill = Ahr, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = Ahr.qPCR.Ob, aes(x = Time, y = Average, linetype = Ahr), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.5) +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 20, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 20, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 20, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  geom_text(x = 4, y = 0.7, label = "*", color = "black", size = 10) +
  geom_text(x = 5, y = 0.5, label = "*", color = "black", size = 10) +
  geom_text(x = 6, y = 0.3, label = "*", color = "black", size = 10) +
  scale_fill_manual(values = c("#0066CC"))


print(Ahr.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/12042017-Ahr-Ob-qPCR-Fig4B.pdf", width = 8, height = 8)
Ahr.qPCR.Ob.pl
dev.off()

## qPCR plot of Glis1 in adipo for figure 4 ##
Glis1.qPCR = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/qPCR/02022018-Glis1-qPCR.txt", delim = "\t")


Glis1.qPCR.Ad = Glis1.qPCR[1:5,] %>% 
  mutate(Time = c(1, 3, 5, 9, 15), Glis1 = "Glis1") %>%
  add_row(Glis1 = "Glis1", Average = 1, SD = 0, SEM = 0, Time = 0, .before = 1)

Glis1.qPCR.Ad$Time = factor(Glis1.qPCR.Ad$Time)
Glis1.qPCR.Ob = Glis1.qPCR[6:10,] %>% 
  mutate(Time = c(1, 3, 5, 9, 15), Glis1 = "Glis1") %>%
  add_row(Glis1 = "Glis1", Average = 1, SD = 0, SEM = 0, Time = 0, .before = 1)

Glis1.qPCR.Ob$Time = factor(Glis1.qPCR.Ob$Time)

Glis1.qPCR.Ad.pl = ggplot(data = Glis1.qPCR.Ad, aes(x = Time, y = Average, fill = Glis1, group = 1)) +
  geom_point(color = "#FF6666", size = 5) +
  geom_line(data = Glis1.qPCR.Ad, aes(x = Time, y = Average, linetype = Glis1), color = "#FF6666", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.0) +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Adipocyte differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  geom_text(x = 2, y = 0.50, label = "*", color = "black", size = 20) +
  geom_text(x = 3, y = 0.35, label = "**", color = "black", size = 20) +
  geom_text(x = 4, y = 0.60, label = "**", color = "black", size = 20) +
  geom_text(x = 5, y = 0.30, label = "**", color = "black", size = 20) +
  geom_text(x = 6, y = 0.25, label = "**", color = "black", size = 20)


print(Glis1.qPCR.Ad.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure6/15032018-Glis1-Ad-qPCR-Fig6C-biggerFont.pdf", width = 10, height = 10)
Glis1.qPCR.Ad.pl
dev.off()

Glis1.qPCR.Ob.pl = ggplot(data = Glis1.qPCR.Ob, aes(x = Time, y = Average, fill = Glis1, group = 1)) +
  geom_point(color = "#0066CC", size = 5) +
  geom_line(data = Glis1.qPCR.Ob, aes(x = Time, y = Average, linetype = Glis1), color = "#0066CC", size = 1.5) +
  geom_errorbar(aes(ymin = Average - SEM, ymax = Average + SEM), width = 0.2, size = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 2.0) +
  theme_classic() +
  ylim(0, 1.0) +
  labs(x = "Time [Days]", y = "Expression relative to undifferentiated ST2 cells", 
       title = "Osteoblast differentiation") +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.text.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.text.y = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.x = element_text(face = "bold", size = 30),
        axis.title.y = element_text(face = "bold", size = 30),
        legend.position = c(0.8, 0.8), legend.text = element_text(size = 30, face = "bold.italic"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 30, hjust = 0.5)) +
  scale_x_discrete(breaks = c("0", "1", "3", "5", "9", "15"),
                   labels = c("D0", "D1", "D3", "D5", "D9", "D15")) +
  geom_text(x = 2, y = 0.85, label = "*", color = "black", size = 20) +
  geom_text(x = 3, y = 0.65, label = "*", color = "black", size = 20) +
  geom_text(x = 4, y = 0.65, label = "*", color = "black", size = 20) +
  geom_text(x = 5, y = 0.4, label = "***", color = "black", size = 20) +
  geom_text(x = 6, y = 0.3, label = "***", color = "black", size = 20) +
  scale_fill_manual(values = c("#0066CC"))


print(Glis1.qPCR.Ob.pl)
pdf("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Figures/Figure6/15032018-Glis1-Ob-qPCR-Fig6D-biggerFont.pdf", width = 10, height = 10)
Glis1.qPCR.Ob.pl
dev.off()
