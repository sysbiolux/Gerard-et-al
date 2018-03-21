#!/bin/bash

## CONSTRUCTION OF THE SE METAMAP TO BE ABLE TO TRACK THE SEs OVER TIME ##

export PATH=/work/users/dgerard/programs/bedtools2/bin:$PATH

cd /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/SuperEnhancers/

# Use sort (Unix) to sort the SE BED files by chromosome
## ADIPOGENESIS ##
sort -k1,1 TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.bed > TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-A-D1-minDist10kb-L0.bed > TC1-SuperEnhancers-A-D1-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-A-D3-minDist10kb-L0.bed > TC1-SuperEnhancers-A-D3-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-A-D5-minDist10kb-L0.bed > TC1-SuperEnhancers-A-D5-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-A-D15-minDist10kb-L0.bed > TC1-SuperEnhancers-A-D15-minDist10kb-L0.sorted.bed

## OSTEOBLASTOGENESIS ##
sort -k1,1 TC1-SuperEnhancers-O-D1-minDist10kb-L0.bed > TC1-SuperEnhancers-O-D1-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-O-D3-minDist10kb-L0.bed > TC1-SuperEnhancers-O-D3-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-O-D5-minDist10kb-L0.bed > TC1-SuperEnhancers-O-D5-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-O-D9-minDist10kb-L0.bed > TC1-SuperEnhancers-O-D9-minDist10kb-L0.sorted.bed
sort -k1,1 TC1-SuperEnhancers-O-D15-minDist10kb-L0.bed > TC1-SuperEnhancers-O-D15-minDist10kb-L0.sorted.bed


# Then, use genomeCoverageBed (version V2.24.0) from bedtools2 to compute BEDGRAPH (-bg) summaries of SEs coverage (e.g., aligned sequences) for a given genome (mm10).
## ADIPOGENESIS ##
genomeCoverageBed -i TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-A-D1-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-A-D1-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-A-D3-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-A-D3-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-A-D5-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-A-D5-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-A-D15-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-A-D15-minDist10kb-L0.sorted.bg

## OSTEOBLASTOGENESIS ##
genomeCoverageBed -i TC1-SuperEnhancers-O-D1-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-O-D1-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-O-D3-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-O-D3-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-O-D5-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-O-D5-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-O-D9-minDist10kb-L0.sorted.bed -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-O-D9-minDist10kb-L0.sorted.bg
genomeCoverageBed -i TC1-SuperEnhancers-O-D15-minDist10kb-L0.sorted.bed. -g ./../../mm10.chr.sizesNOPREFIX -bg > TC1-SuperEnhancers-O-D15-minDist10kb-L0.sorted.bg

# Then, use unionBedGraphs (version V2.24.0) from bedtools2 to combine multiple BEDGRAPH files into a single file such that one can directly compare coverage (and other text-values such as genotypes) across multiple sample
unionBedGraphs -i TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-A-D1-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-A-D3-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-A-D5-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-A-D15-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-O-D1-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-O-D3-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-O-D5-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-O-D9-minDist10kb-L0.sorted.bg TC1-SuperEnhancers-O-D15-minDist10kb-L0.sorted.bg > 05042017-SEMapAdipoOsteo.bg

# Then, use mergeBed (version V2.24.0) from bedtools2 to combine overlapping (>= 1 bp) or “book-ended” features in an interval file into a single feature which spans all of the combined features.
mergeBed -i 05042017-SEMapAdipoOsteo.bg > 05042017-SEMapAdipoOsteoMERGED.bed