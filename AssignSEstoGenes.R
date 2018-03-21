## Assign super-enhancers (SEs) to their target genes ##
getwd()

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(stringr)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)

## Load the Rdata from the time-course analysed by DEseq2
TC = read.table("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/31032017-TC-NormCountMat.txt", sep = "\t",
                header = TRUE) %>%
  rownames_to_column(var = "ensembl_gene_id")

## Retrieve Ensembl ID of the time-course analysed by DEseq2 and assign the gene coordinates based on the mm10-GRCm38.79 annotation from Ensembl
head(TC)
ddsTC2.ensID = TC %>%
  dplyr::select(ensembl_gene_id)
head(ddsTC2.ensID)

## Load the annotation in a txdb variable ##
txdb.new = loadDb("Y:/Deborah.GERARD/CompBetweenOlfCUffAndNewCUff/Mus_musculus.GRCm38.79.sqlite")

## Extract information about the genes ##
txdb.genes = genes(txdb.new)
head(txdb.genes)
txdb.genes.df = as.data.frame(txdb.genes) %>%
  dplyr::rename(ensembl_gene_id = gene_id)
head(txdb.genes.df)

# Extract coordinates from txdb.genes.df and append them to ddsTC2.ensID
ddsTC2.ensID.coord = ddsTC2.ensID %>%
  left_join(txdb.genes.df, by = "ensembl_gene_id") %>%
  dplyr::select(seqnames, start, end, ensembl_gene_id) %>%
  dplyr::rename(chr = seqnames, Start = start, End = end)

write.table(ddsTC2.ensID.coord, "19062017-Mus_musculus.GRCm38.79.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Associate the merged SEs to genes within 1mb window ##
#bedmap --range 500000 --echo --echo-map-id --delim '|' 19062017-Mus_musculus.GRCm38.79.sorted.bed 05042017-SEMapAdipoOsteoMERGEDFin_sorted.bed | awk -F'|' '{ reference = $1; mappings = $2; n = split(mappings, m, "|"); for(i = 1; i <= n; i++) { print reference"|"m[i]; } }' > 19062017-SEsAssignedToMm10Genes_TSS_Genes_1000000kb.bed

## Load the data - here the merged SEs have been assigned to genes from the TSS of the genes +- 50kb ##
dat.1mb = read_delim("19062017-SEsAssignedToMm10Genes_TSS_Genes_1000000kb.bed", delim = "\t", progress = TRUE, col_names = FALSE)
dat.1mb.clean = dat.1mb %>%
  dplyr::rename(chr = X1, Start = X2, End = X3, ensembl_gene_id = X4)

# Extract the SEs information
Col.SEs = dat.1mb.clean %>%
  dplyr::select(chr) %>%
  filter(str_detect(chr, "SE")) %>%
  dplyr::rename(SE_ID = chr) #11895 rows

# Extract the gene information
dat.1mb.clean.fin = dat.1mb.clean %>%
  filter(!str_detect(chr, "SE"))

# Combine gene and SE information in a proper way
dat.SE.Genes = bind_cols(dat.1mb.clean.fin, Col.SEs) %>%
  dplyr::select(-chr, -Start, -End) 
dat.SE.Genes$SE_ID = gsub("\\|", "", dat.SE.Genes$SE_ID)
dat.SE.Genes$SE_ID = gsub("\\;", ":", dat.SE.Genes$SE_ID)
View(dat.SE.Genes)

# Multiples SEs for one gene, split SES on multiple rows 
dat.SE.Genes = dat.SE.Genes %>%
  separate_rows(SE_ID, sep = ":")
View(dat.SE.Genes)

###################################################################################### Assign genes to the 25 SEs that are shared between adipocyte and osteoblast ############################################################
# Extract the 25 genes from the venny Fig4 in Gerard et al. 2017
Dyn.SEs = read_delim("15062017-GerardCell-Fig4-Venny-25SEs.txt", delim = "\t", progress = TRUE, col_names = TRUE)
Dyn.SEs.25 = dat.SE.Genes %>%
  filter(SE_ID %in% Dyn.SEs$SE_ID) %>%
  unite(ensembl_gene_id, SE_ID, col = "gene_SEID", sep = "_", remove = FALSE)
head(Dyn.SEs.25)

## Retrieve the normalized counts for all genes which might be regulated by SEs ##

spl_Dyn.SEs.25.counts = TC %>%
  filter(ensembl_gene_id %in% Dyn.SEs.25$ensembl_gene_id) %>%
  tbl_df()

# spl_Dyn.SEs.25.counts = spl_Dyn.SEs.25 %>%
#   map(function(df) plotCounts(ddsTC2, gene = df, intgroup = c("SampleID", "Cell"), returnData = TRUE))
  
## Calculate the mean of genes per sample (remove the adipo day9 (because the IP for Adipo day9 was bad, so not used) ##
spl_Dyn.SEs.25.mean = spl_Dyn.SEs.25.counts %>% 
  mutate(ST2.D0.Av = ((ST2.D0.rep1 + ST2.D0.rep2 + ST2.D0.rep3)/3), 
                                A.D1.Av = ((A.D1.rep1 + A.D1.rep2 + A.D1.rep3)/3),
                                A.D3.Av = ((A.D3.rep1 + A.D3.rep2 + A.D3.rep3)/3),
                                A.D5.Av = ((A.D5.rep1 + A.D5.rep2 + A.D5.rep3)/3),
                                A.D15.Av = ((A.D15.rep1 + A.D15.rep2 + A.D15.rep3)/3),
                                O.D1.Av = ((O.D1.rep1 + O.D1.rep2 + O.D1.rep3)/3),
                                O.D3.Av = ((O.D3.rep1 + O.D3.rep2 + O.D3.rep3)/3),
                                O.D5.Av = ((O.D5.rep1 + O.D5.rep2 + O.D5.rep3)/3),
                                O.D9.Av = ((O.D9.rep1 + O.D9.rep2 + O.D9.rep3)/3),
                                O.D15.Av = ((O.D15.rep1 + O.D15.rep2 + O.D15.rep3)/3)) %>%
  dplyr::select(-contains("rep"), ensembl_gene_id)

spl_Dyn.SEs.25.mean.Ad = spl_Dyn.SEs.25.mean %>%
  dplyr::select(ensembl_gene_id, ST2.D0.Av:A.D15.Av)

spl_Dyn.SEs.25.mean.Ob = spl_Dyn.SEs.25.mean %>%
  dplyr::select(ensembl_gene_id, ST2.D0.Av, O.D1.Av:O.D15.Av)  
  
## Load the count for the merged SE in adipo and osteo
mergeSE.Ad.Counts = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt",
                               delim = "\t")

mergeSE.Ob.COunts = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt",
                               delim = "\t")

## Retrieve the normalized counts for the dynamic SEs ##
Dyn.SEs.25.Ad = mergeSE.Ad.Counts %>%
  filter(SE_ID %in% Dyn.SEs$SE_ID) %>%
  right_join(Dyn.SEs.25, by = "SE_ID") %>%
  right_join(spl_Dyn.SEs.25.mean.Ad, by = "ensembl_gene_id")
View(Dyn.SEs.25.Ad)

Dyn.SEs.25.Ob = mergeSE.Ob.COunts %>%
  filter(SE_ID %in% Dyn.SEs$SE_ID) %>%
  right_join(Dyn.SEs.25, by = "SE_ID") %>%
  right_join(spl_Dyn.SEs.25.mean.Ob, by = "ensembl_gene_id")
View(Dyn.SEs.25.Ob)

## calculate the correlation between the dyn SE and their potentially target genes within a 1000 kb window (gene TSS +- 500 kb)
# adipo
cor.mat.Ad = matrix(ncol = 1, nrow = 327)
for (i in 1:nrow(Dyn.SEs.25.Ad)) {
  paste("i = ", i)
  cor.mat.Ad[i,] = cor(as.numeric(Dyn.SEs.25.Ad[i,2:6]), as.numeric(Dyn.SEs.25.Ad[i,9:13]), method = "pearson")
}
  
cor.mat.Ad.fin = cor.mat.Ad %>%
  tbl_df() %>%
  dplyr::rename(PearsonCor = V1) %>%
  cbind(Dyn.SEs.25.Ad) %>%
  dplyr::select(gene_SEID, SE_ID, ensembl_gene_id, PearsonCor)

## calculate the correlation between the dyn SE and their potentially target genes within a 1000 kb window (gene TSS +- 500 kb)
# osteo
cor.mat.Ob = matrix(ncol = 1, nrow = 327)
for (i in 1:nrow(Dyn.SEs.25.Ob)) {
  paste("i = ", i)
  cor.mat.Ob[i,] = cor(as.numeric(Dyn.SEs.25.Ob[i,2:7]), as.numeric(Dyn.SEs.25.Ob[i,10:15]), method = "pearson")
}

cor.mat.Ob.fin = cor.mat.Ob %>%
  tbl_df() %>%
  dplyr::rename(PearsonCor = V1) %>%
  cbind(Dyn.SEs.25.Ob) %>%
  dplyr::select(gene_SEID, SE_ID, ensembl_gene_id, PearsonCor)

## Get the gene name using the annotation from Ensembl (mar2015.archive.ensembl.org)
ensembl79 = useMart(host = 'mar2015.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
filters = listFilters(ensembl79)
head(filters)
attributes = listAttributes(ensembl79)
head(attributes)

ensGId.Ad = cor.mat.Ad.fin %>%
  dplyr::select(ensembl_gene_id)

# ensGId.Ob = cor.mat.Ob.fin %>%
#   dplyr::select(ensembl_gene_id)

qry.SEs = getBM(attributes = c('external_gene_name', 'ensembl_gene_id', 'gene_biotype'), filters = 'ensembl_gene_id', values = ensGId.Ad, mart = ensembl79)

# Add the gene name to the matrix => easier than Ensembl IDs
cor.mat.Ad.fin.GN = cor.mat.Ad.fin %>%
  dplyr::full_join(qry.SEs, by = "ensembl_gene_id") %>%
  tbl_df()

cor.mat.Ob.fin.GN = cor.mat.Ob.fin %>%
  dplyr::full_join(qry.SEs, by = "ensembl_gene_id") %>%
  tbl_df()

## Create correlation matrix for every dyn SEs in adipo and osteo
id.Ad = factor(cor.mat.Ad.fin.GN$SE_ID) %>% levels()
id.Ob = factor(cor.mat.Ob.fin.GN$SE_ID) %>% levels()

for (i in 1:length(id.Ad)) {
  assign(paste0("PearCor.mat.Ad_", id.Ad[i]), filter(cor.mat.Ad.fin.GN, SE_ID == id.Ad[i]) %>%
           arrange(desc(PearsonCor)))
}

`PearCor.mat.Ad_SE-286` %>% View

# Extract top gene associated to SE 
## Adipo
Shared.Ad.DynSE.BestIt = matrix(nrow = 25, ncol = 5)
Shared.Ad.datalist = list()

for (i in 1:length(id.Ad)) {
  Shared.Ad.DynSE.BestIt$i = filter(cor.mat.Ad.fin.GN, SE_ID == id.Ad[i]) %>%
           arrange(desc(PearsonCor)) %>%
    filter(gene_biotype == "protein_coding") %>%
    top_n(1, PearsonCor)
  Shared.Ad.datalist[[i]] = Shared.Ad.DynSE.BestIt$i
}

big_Shared.Ad.DynSE.BestIt = bind_rows(Shared.Ad.datalist) 
big_Shared.Ad.DynSE.BestIt = big_Shared.Ad.DynSE.BestIt %>%
  rename(gene_SEID_Ad = gene_SEID, 
         PearsonCor_Ad = PearsonCor,
         ensembl_gene_id_Ad = ensembl_gene_id,
         external_gene_name_Ad = external_gene_name,
         gene_biotype_Ad = gene_biotype)

## Osteo
for (i in 1:length(id.Ob)) {
  assign(paste0("PearCor.mat.Ob_", id.Ob[i]), filter(cor.mat.Ob.fin.GN, SE_ID == id.Ob[i]) %>%
           arrange(desc(PearsonCor)))
}

`PearCor.mat.Ob_SE-921` %>% View

Shared.Ob.DynSE.BestIt = matrix(nrow = 25, ncol = 5)
Shared.Ob.datalist = list()

for (i in 1:length(id.Ob)) {
  Shared.Ob.DynSE.BestIt$i = filter(cor.mat.Ob.fin.GN, SE_ID == id.Ob[i]) %>%
           arrange(desc(PearsonCor)) %>%
    filter(gene_biotype == "protein_coding") %>%
    top_n(1, PearsonCor)
  Shared.Ob.datalist[[i]] = Shared.Ob.DynSE.BestIt$i
}

big_Shared.Ob.DynSE.BestIt = bind_rows(Shared.Ob.datalist)
big_Shared.Ob.DynSE.BestIt = big_Shared.Ob.DynSE.BestIt %>%
  rename(gene_SEID_Ob = gene_SEID, 
         PearsonCor_Ob = PearsonCor,
         ensembl_gene_id_Ob = ensembl_gene_id,
         external_gene_name_Ob = external_gene_name,
         gene_biotype_Ob = gene_biotype)

# Combine the Paerson correlation from the Shared SEs in adipo and osteo together
big_Shared.All.DynSE.BestIt = big_Shared.Ad.DynSE.BestIt %>%
  full_join(big_Shared.Ob.DynSE.BestIt, by = "SE_ID") %>%
  select(-gene_biotype_Ad, -gene_biotype_Ob)
big_Shared.All.DynSE.BestIt %>% View

write_delim(big_Shared.All.DynSE.BestIt, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Supplement/08032018-Shared_DynSEs_PearsonCor_Ad&Ob.txt",
            delim = "\t", col_names = TRUE)




###################################################################################### Assign genes to the 95 SEs that are unique to adipocyte ############################################################
# Extract the 95 genes from the venny Fig4 in Gerard et al. 2017
Dyn.SEs.Ad = read_delim("21022018-GerardCell-Fig4-Venny-95SEs-onlyAdipo.txt", delim = "\t", progress = TRUE, col_names = TRUE)
Dyn.SEs.95 = dat.SE.Genes %>%
  filter(SE_ID %in% Dyn.SEs.Ad$SE_ID) %>%
  unite(ensembl_gene_id, SE_ID, col = "gene_SEID", sep = "_", remove = FALSE)
Dyn.SEs.95

## Retrieve the normalized counts for all genes which might be regulated by SEs ##

spl_Dyn.SEs.95.Ad.counts = TC %>%
  filter(ensembl_gene_id %in% Dyn.SEs.95$ensembl_gene_id) %>%
  tbl_df()

## Calculate the mean of genes per sample (remove the adipo day9 (because the IP for Adipo day9 was bad, so not used) ##
spl_Dyn.SEs.95.Ad.mean = spl_Dyn.SEs.95.Ad.counts %>% 
  mutate(ST2.D0.Av = ((ST2.D0.rep1 + ST2.D0.rep2 + ST2.D0.rep3)/3), 
         A.D1.Av = ((A.D1.rep1 + A.D1.rep2 + A.D1.rep3)/3),
         A.D3.Av = ((A.D3.rep1 + A.D3.rep2 + A.D3.rep3)/3),
         A.D5.Av = ((A.D5.rep1 + A.D5.rep2 + A.D5.rep3)/3),
         A.D15.Av = ((A.D15.rep1 + A.D15.rep2 + A.D15.rep3)/3)) %>%
  dplyr::select(-contains("rep"), ensembl_gene_id)


## Load the count for the merged SE in adipo
mergeSE.Ad.Counts = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapAdipoMERGEDFinCounts.txt",
                               delim = "\t")

## Retrieve the normalized counts for the dynamic SEs ##
Dyn.SEs.95.Ad = mergeSE.Ad.Counts %>%
  filter(SE_ID %in% Dyn.SEs.Ad$SE_ID) %>%
  right_join(Dyn.SEs.95, by = "SE_ID") %>%
  right_join(spl_Dyn.SEs.95.Ad.mean, by = "ensembl_gene_id")
View(Dyn.SEs.95.Ad)


## calculate the correlation between the dyn SE and their potentially target genes within a 1000 kb window (gene TSS +- 500 kb)
# adipo
cor.mat.Ad.95 = matrix(ncol = 1, nrow = 1398)
for (i in 1:nrow(Dyn.SEs.95.Ad)) {
  paste("i = ", i)
  cor.mat.Ad.95[i,] = cor(as.numeric(Dyn.SEs.95.Ad[i,2:6]), as.numeric(Dyn.SEs.95.Ad[i,9:13]), method = "pearson")
}

cor.mat.Ad.95.fin = cor.mat.Ad.95 %>%
  tbl_df() %>%
  dplyr::rename(PearsonCor = V1) %>%
  cbind(Dyn.SEs.95.Ad) %>%
  dplyr::select(gene_SEID, SE_ID, ensembl_gene_id, PearsonCor)


## Get the gene name using the annotation from Ensembl (mar2015.archive.ensembl.org)
ensembl79 = useMart(host = 'mar2015.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
filters = listFilters(ensembl79)
head(filters)
attributes = listAttributes(ensembl79)
head(attributes)

ensGId.Ad.95 = cor.mat.Ad.95.fin %>%
  dplyr::select(ensembl_gene_id)

qry.SEs.Ad.95 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = ensGId.Ad.95, mart = ensembl79)

# Add the gene name to the matrix => easier than Ensembl IDs
cor.mat.Ad.95.fin.GN = cor.mat.Ad.95.fin %>%
  dplyr::full_join(qry.SEs.Ad.95, by = "ensembl_gene_id") %>%
  tbl_df()
write_delim(cor.mat.Ad.95.fin.GN, "21022018-UniqDynSEAdipo.txt", delim = "\t", col_names = TRUE)

## Create correlation matrix for every dyn SEs in adipo 
id.Ad.95 = factor(cor.mat.Ad.95.fin.GN$SE_ID) %>% levels()

for (i in 1:length(id.Ad.95)) {
  assign(paste0("PearCor.mat.Ad_Unique_SEs-", id.Ad.95[i]), filter(cor.mat.Ad.95.fin.GN, SE_ID == id.Ad.95[i]) %>%
           arrange(desc(PearsonCor)))
}

View(`PearCor.mat.Ad_Unique_SEs-SE-1032`) 

# Extract top gene associated to SE
Dyn.Ad.BestIt = matrix(nrow = 93, ncol = 5)
datalist = list()
for (i in 1:length(id.Ad.95)) {
  Dyn.Ad.BestIt$i = filter(cor.mat.Ad.95.fin.GN, SE_ID == id.Ad.95[i]) %>%
           arrange(desc(PearsonCor)) %>%
             top_n(1, PearsonCor)
  datalist[[i]] = Dyn.Ad.BestIt$i

  }
big_Dyn.Ad.BestIt = bind_rows(datalist)
View(big_Dyn.Ad.BestIt)
write_delim(big_Dyn.Ad.BestIt, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Supplement/22022018-SupplTab3-DynSEOnlyAdipo-BestCor.txt", col_names = TRUE, delim = "\t")


## Check which genes that are associated to dynamic SEs in adipo are TFs ##
mm_Ortho_TF = read_delim("MOUSE_TFs_CURATED.txt", delim = "\t", progress = TRUE) %>%
  rename(external_gene_name = `Ortholog Symbol`) %>%
  select(external_gene_name)

# Overlap mouse TFs with genes associated to dynamic SEs in adipo
Adipo.Dyn.SEs.BI.TF = mm_Ortho_TF %>%
  filter(external_gene_name %in% big_Dyn.Ad.BestIt$external_gene_name) %>%
  left_join(big_Dyn.Ad.BestIt, by = "external_gene_name")
Adipo.Dyn.SEs.BI.TF
write_delim(Adipo.Dyn.SEs.BI.TF, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Supplement/22022018-SupplTab3-DynSEOnlyAdipo-BestCor-OnlyTFs.txt", col_names = TRUE, delim  ="\t")

###################################################################################### Assign genes to the 54 SEs that are unique to osteoblast ############################################################
# Extract the 54 genes from the venny Fig4 in Gerard et al. 2017
Dyn.SEs.Ob = read_delim("22022018-GerardCell-Fig4-Venny-54SEs-onlyAdipo.txt", delim = "\t", progress = TRUE, col_names = TRUE)
Dyn.SEs.54 = dat.SE.Genes %>%
  filter(SE_ID %in% Dyn.SEs.Ob$SE_ID) %>%
  unite(ensembl_gene_id, SE_ID, col = "gene_SEID", sep = "_", remove = FALSE)
Dyn.SEs.54

## Retrieve the normalized counts for all genes which might be regulated by SEs ##

spl_Dyn.SEs.54.Ob.counts = TC %>%
  filter(ensembl_gene_id %in% Dyn.SEs.54$ensembl_gene_id) %>%
  tbl_df()

# spl_Dyn.SEs.25.counts = spl_Dyn.SEs.25 %>%
#   map(function(df) plotCounts(ddsTC2, gene = df, intgroup = c("SampleID", "Cell"), returnData = TRUE))

# head(spl_Dyn.SEs.25.counts)
# spl_Dyn.SEs.25.counts.fin = matrix(unlist(spl_Dyn.SEs.25.counts), ncol = 327, byrow = FALSE)
# spl_Dyn.SEs.25.counts.fin = as_data_frame(spl_Dyn.SEs.25.counts.fin[1:33,]) %>%
#   mutate(Sample = spl_Dyn.SEs.25.counts[[1]]$SampleID) %>%
#   dplyr::select(Sample, everything())
# 
# colnames(spl_Dyn.SEs.25.counts.fin)[2:length(spl_Dyn.SEs.25.counts.fin)] = Dyn.SEs.25$gene_SEID
# spl_Dyn.SEs.25.counts.fin$Sample = factor(spl_Dyn.SEs.25.counts.fin$Sample, levels = c("ST2-D0", "A-D1", "A-D3", "A-D5", "A-D9", "A-D15",
#                                                                                        "O-D1", "O-D3", "O-D5", "O-D9", "O-D15"))
# View(spl_Dyn.SEs.25.counts.fin)

## Calculate the mean of genes per sample  ##
spl_Dyn.SEs.54.Ob.mean = spl_Dyn.SEs.54.Ob.counts %>% 
  mutate(ST2.D0.Av = ((ST2.D0.rep1 + ST2.D0.rep2 + ST2.D0.rep3)/3), 
         O.D1.Av = ((O.D1.rep1 + O.D1.rep2 + O.D1.rep3)/3),
         O.D3.Av = ((O.D3.rep1 + O.D3.rep2 + O.D3.rep3)/3),
         O.D5.Av = ((O.D5.rep1 + O.D5.rep2 + O.D5.rep3)/3),
         O.D9.Av = ((O.D9.rep1 + O.D9.rep2 + O.D9.rep3)/3),
         O.D15.Av = ((O.D15.rep1 + O.D15.rep2 + O.D15.rep3)/3)) %>%
  dplyr::select(-contains("rep"), ensembl_gene_id)


## Load the count for the merged SE in osteo
mergeSE.Ob.COunts = read_delim("Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/STEM/05042017-SEMapOsteoMERGEDFinCounts.txt",
                               delim = "\t")

## Retrieve the normalized counts for the dynamic SEs ##
Dyn.SEs.54.Ob = mergeSE.Ob.COunts %>%
  filter(SE_ID %in% Dyn.SEs.Ob$SE_ID) %>%
  right_join(Dyn.SEs.54, by = "SE_ID") %>%
  right_join(spl_Dyn.SEs.54.Ob.mean, by = "ensembl_gene_id")
View(Dyn.SEs.54.Ob)


## calculate the correlation between the dyn SE and their potentially target genes within a 1000 kb window (gene TSS +- 500 kb)
# osteo
cor.mat.Ob.54 = matrix(ncol = 1, nrow = 867)
for (i in 1:nrow(Dyn.SEs.54.Ob)) {
  paste("i = ", i)
  cor.mat.Ob.54[i,] = cor(as.numeric(Dyn.SEs.54.Ob[i,2:7]), as.numeric(Dyn.SEs.54.Ob[i,10:15]), method = "pearson")
}

cor.mat.Ob.54.fin = cor.mat.Ob.54 %>%
  tbl_df() %>%
  dplyr::rename(PearsonCor = V1) %>%
  cbind(Dyn.SEs.54.Ob) %>%
  dplyr::select(gene_SEID, SE_ID, ensembl_gene_id, PearsonCor)


## Get the gene name using the annotation from Ensembl (mar2015.archive.ensembl.org)
ensembl79 = useMart(host = 'mar2015.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
filters = listFilters(ensembl79)
head(filters)
attributes = listAttributes(ensembl79)
head(attributes)

ensGId.Ob.54 = cor.mat.Ob.54.fin %>%
  dplyr::select(ensembl_gene_id)

qry.SEs.Ob.54 = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), filters = 'ensembl_gene_id', values = ensGId.Ob.54, mart = ensembl79)

# Add the gene name to the matrix => easier than Ensembl IDs
cor.mat.Ob.54.fin.GN = cor.mat.Ob.54.fin %>%
  dplyr::full_join(qry.SEs.Ob.54, by = "ensembl_gene_id") %>%
  tbl_df()
write_delim(cor.mat.Ob.54.fin.GN, "22022018-UniqDynSEOsteo.txt", delim = "\t", col_names = TRUE)

## Create correlation matrix for every dyn SEs in adipo 
id.Ob.54 = factor(cor.mat.Ob.54.fin.GN$SE_ID) %>% levels()

for (i in 1:length(id.Ob.54)) {
  assign(paste0("PearCor.mat.Ob_Unique_SEs-", id.Ob.54[i]), filter(cor.mat.Ob.54.fin.GN, SE_ID == id.Ob.54[i]) %>%
           arrange(desc(PearsonCor)))
}

View(`PearCor.mat.Ob_Unique_SEs-SE-913`) 

# Extract top gene associated to SE
Dyn.Ob.BestIt = matrix(nrow = 54, ncol = 5)
datalist.Ob = list()
for (i in 1:length(id.Ob.54)) {
  Dyn.Ob.BestIt$i = filter(cor.mat.Ob.54.fin.GN, SE_ID == id.Ob.54[i]) %>%
    arrange(desc(PearsonCor)) %>%
    top_n(1, PearsonCor)
  datalist.Ob[[i]] = Dyn.Ob.BestIt$i
  
}
big_Dyn.Ob.BestIt = bind_rows(datalist.Ob)
View(big_Dyn.Ob.BestIt)
write_delim(big_Dyn.Ob.BestIt, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Supplement/22022018-SupplTab3-DynSEOnlyOsteo-BestCor.txt", col_names = TRUE, delim = "\t")


## Check which genes that are associated to dynamic SEs in adipo are TFs ##
mm_Ortho_TF = read_delim("MOUSE_TFs_CURATED.txt", delim = "\t", progress = TRUE) %>%
  rename(external_gene_name = `Ortholog Symbol`) %>%
  select(external_gene_name)

# Overlap mouse TFs with genes associated to dynamic SEs in adipo
Osteo.Dyn.SEs.BI.TF = mm_Ortho_TF %>%
  filter(external_gene_name %in% big_Dyn.Ob.BestIt$external_gene_name) %>%
  left_join(big_Dyn.Ob.BestIt, by = "external_gene_name")
Osteo.Dyn.SEs.BI.TF
write_delim(Osteo.Dyn.SEs.BI.TF, "Y:/Deborah.GERARD/Gerard et al. - Manuscript 1/Supplement/22022018-SupplTab3-DynSEOnlyOsteo-BestCor-OnlyTFs.txt", col_names = TRUE, delim = "\t")
