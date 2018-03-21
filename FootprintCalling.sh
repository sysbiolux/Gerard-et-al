#!/bin/bash

# Concatenate typical enhancers and super-enhancers per time-point in order to call footprints
cat TC1-TypicalEnhancers-ST2-D0-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-ST2-D0-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-ST2-D0.sorted.bed
cat TC1-TypicalEnhancers-A-D1-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-A-D1-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-A-D1.sorted.bed
cat TC1-TypicalEnhancers-A-D3-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-A-D3-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-A-D3.sorted.bed
cat TC1-TypicalEnhancers-A-D5-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-A-D5-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-A-D5.sorted.bed
cat TC1-TypicalEnhancers-A-D15-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-A-D15-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-A-D15.sorted.bed
cat TC1-TypicalEnhancers-O-D1-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-O-D1-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-O-D1.sorted.bed
cat TC1-TypicalEnhancers-O-D3-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-O-D3-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-O-D3.sorted.bed
cat TC1-TypicalEnhancers-O-D5-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-O-D5-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-O-D5.sorted.bed
cat TC1-TypicalEnhancers-O-D9-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-O-D9-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-O-D9.sorted.bed
cat TC1-TypicalEnhancers-O-D15-minDist10kb-L0-HINT.bed TC1-SuperEnhancers-O-D15-minDist10kb-L0-HINT.bed | sort -k1,1 > 06042017-SETE-O-D15.sorted.bed

# Call the footprints using HINT 0.9.9

parallel -j5 'rgt-hint --organism mm10 --histone-norm-per 98 --histone-slope-per 98 --default-bias-correction --output-fname {}.bed {} --output-location ../Output/FP-08052017/
