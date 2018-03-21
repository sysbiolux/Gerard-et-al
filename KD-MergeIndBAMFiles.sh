#!/bin/bash

## Use Samtools to merge individual BAM files into one BAM file 
export PATH=/work/users/dgerard/programs/samtools-1.2:$PATH

cd /scratch/users/dgerard/Ahr-Gzf1-KD/MappedReads/BAMFiles

echo "Merging BAM files for A-siAHR-D1"

samtools merge -r -@ 12 A-siAHR-D1-rep1.bam Mapped.Trim.AsiAHRD11*
samtools merge -r -@ 12 A-siAHR-D1-rep2.bam Mapped.Trim.AsiAHRD12*
samtools merge -r -@ 12 A-siAHR-D1-rep3.bam Mapped.Trim.AsiAHRD13*

echo "Merging BAM files for A-siCTRL-D1"

samtools merge -r -@ 12 A-siCTRL-D1-rep1.bam Mapped.Trim.AsiCTRLD11*
samtools merge -r -@ 12 A-siCTRL-D1-rep2.bam Mapped.Trim.AsiCTRLD12*
samtools merge -r -@ 12 A-siCTRL-D1-rep3.bam Mapped.Trim.AsiCTRLD13*

echo "Merging BAM files for A-siCTRL-D9"

samtools merge -r -@ 12 A-siCTRL-D9-rep1.bam Mapped.Trim.AsiCTRLD91*
samtools merge -r -@ 12 A-siCTRL-D9-rep2.bam Mapped.Trim.AsiCTRLD92*
samtools merge -r -@ 12 A-siCTRL-D9-rep3.bam Mapped.Trim.AsiCTRLD93*

echo "Merging BAM files for A-siGZF1-D9"

samtools merge -r -@ 12 A-siGZF1-D9-rep1.bam Mapped.Trim.AsiGZF1D91*
samtools merge -r -@ 12 A-siGZF1-D9-rep2.bam Mapped.Trim.AsiGZF1D92*
samtools merge -r -@ 12 A-siGZF1-D9-rep3.bam Mapped.Trim.AsiGZF1D93*

echo "Merging BAM files for O-siAHR-D1"

samtools merge -r -@ 12 O-siAHR-D1-rep1.bam Mapped.Trim.OsiAHRD11*
samtools merge -r -@ 12 O-siAHR-D1-rep2.bam Mapped.Trim.OsiAHRD12*
samtools merge -r -@ 12 O-siAHR-D1-rep3.bam Mapped.Trim.OsiAHRD13*

echo "Merging BAM files for O-siCTRL-D1"

samtools merge -r -@ 12 O-siCTRL-D1-rep1.bam Mapped.Trim.OsiCTRLD11*
samtools merge -r -@ 12 O-siCTRL-D1-rep2.bam Mapped.Trim.OsiCTRLD12*
samtools merge -r -@ 12 O-siCTRL-D1-rep3.bam Mapped.Trim.OsiCTRLD13*

echo "Merging BAM files for ST2-siAHR-D1"

samtools merge -r -@ 12 ST2-siAHR-D1-rep1.bam Mapped.Trim.ST2siAHRD11*
samtools merge -r -@ 12 ST2-siAHR-D1-rep2.bam Mapped.Trim.ST2siAHRD12*
samtools merge -r -@ 12 ST2-siAHR-D1-rep3.bam Mapped.Trim.ST2siAHRD13*

echo "Merging BAM files for ST2-siCTRL-D1"

samtools merge -r -@ 12 ST2-siCTRL-D1-rep1.bam Mapped.Trim.ST2siCTRLD11*
samtools merge -r -@ 12 ST2-siCTRL-D1-rep2.bam Mapped.Trim.ST2siCTRLD12*
samtools merge -r -@ 12 ST2-siCTRL-D1-rep3.bam Mapped.Trim.ST2siCTRLD13*









