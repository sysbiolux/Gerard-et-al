#!/bin/bash

# DG - H3K36me3 peak detection using SICER
# window size = 200 bp
# fragment size = 300 bp
# gap size = 600 bp
# FDR = 0.05
# redundancy threshold = 1
# effective genome fraction = 0.8
# species = Mus Musculus mm10
 


cd /work/users/dgerard/programs/SICER/SICER_V1.1/SICER
export PYTHONPATH=/work/users/dgerard/programs/SICER/SICER_V1.1/SICER/lib/python2.7/site-packages/:$PYTHONPATH
echo $PYTHONPATH
source ~/.bashrc

sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ TC1-H3K36-ST2-D0.GRCm38.p3.q30.bed TC1-I-ST2-D0.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-A-D1.GRCm38.p3.q30.bed TC1-I-A-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05 
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-A-D3.GRCm38.p3.q30.bed TC1-I-A-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-A-D5.GRCm38.p3.q30.bed TC1-I-A-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-A-D9.GRCm38.p3.q30.bed TC1-I-A-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-A-D15.GRCm38.p3.q30.bed TC1-I-A-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-O-D1.GRCm38.p3.q30.bed TC1-I-O-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-O-D3.GRCm38.p3.q30.bed TC1-I-O-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-O-D5.GRCm38.p3.q30.bed TC1-I-O-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-O-D9.GRCm38.p3.q30.bed TC1-I-O-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
sh SICER.sh /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/BED_files_SICER TC1-H3K36-O-D15.GRCm38.p3.q30.bed TC1-I-O-D5.GRCm38.p3.q30.bed /scratch/users/dgerard/TC1-ChIP-seq/H3K36me3/ mm10 1 200 300 0.8 600 0.05
