#!/bin/bash

## Peak calling for super-enhancers and enhancers (H3K27Ac) usinh HOMER
## -minDist 10000 = maximun distance used to stitch peaks together, here is 10kb (12500 bp is the default)
## -L 0 = the local fold change is used to identify initial peak by default for super enhancers, you may want to disable this set by setting it to 0 or at least 1.

export PATH=/work/users/dgerard/programs/HOMER/bin:$PATH

## ADIPOGENESIS ##
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-ST2-D0-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-ST2-D0-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-ST2-D0-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-ST2-D0-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-A-D1-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-A-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-A-D1-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-A-D1-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-A-D3-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-A-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-A-D3-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-A-D3-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-A-D5-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-A-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-A-D5-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-A-D5-minDist10kb-L0 -minDist 10000 -L 0
# IP for adipocytes day 9 did not work => Not useable
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-A-D15-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-A-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-A-D15-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-A-D15-minDist10kb-L0 -minDist 10000 -L 0

## OSTEOBLASTOGENESIS ##
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-O-D1-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-O-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-O-D1-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-O-D1-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-O-D3-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-O-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-O-D3-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-O-D3-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-O-D5-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-O-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-O-D5-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-O-D5-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-O-D9-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-O-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-O-D9-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-O-D9-minDist10kb-L0 -minDist 10000 -L 0
findPeaks /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-H3K27-O-D15-TagDirectory/ -i /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-I-O-D5-TagDirectory/ -style super -o /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-SuperEnhancers-O-D15-minDist10kb-L0 -typical /scratch/users/dgerard/TC1-ChIP-seq/H3K27Ac/TC1-TypicalEnhancers-O-D15-minDist10kb-L0 -minDist 10000 -L 0

# Then, convert the HOMER SE files to BED format using pos2bed.pl from the HOMER suite
## ADIPOGENESIS ##
pos2bed.pl TC1-SuperEnhancers-ST2-D0-minDist10kb-L0 > TC1-SuperEnhancers-ST2-D0-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-A-D1-minDist10kb-L0 > TC1-SuperEnhancers-A-D1-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-A-D3-minDist10kb-L0 > TC1-SuperEnhancers-A-D3-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-A-D5-minDist10kb-L0 > TC1-SuperEnhancers-A-D5-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-A-D15-minDist10kb-L0 > TC1-SuperEnhancers-A-D15-minDist10kb-L0.bed

## OSTEOBLASTOGENESIS ##
pos2bed.pl TC1-SuperEnhancers-O-D1-minDist10kb-L0 > TC1-SuperEnhancers-O-D1-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-O-D3-minDist10kb-L0 > TC1-SuperEnhancers-O-D3-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-O-D5-minDist10kb-L0 > TC1-SuperEnhancers-O-D5-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-O-D9-minDist10kb-L0 > TC1-SuperEnhancers-O-D9-minDist10kb-L0.bed
pos2bed.pl TC1-SuperEnhancers-O-D15-minDist10kb-L0 > TC1-SuperEnhancers-O-D15-minDist10kb-L0.bed