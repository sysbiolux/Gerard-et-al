#!/bin/bash

## Remove Illumina adapters from the data and do the trimming based on quality (phred score of 30) using Cutadapts version 1.8.1
# -q 30 = trim low quality ends from reads before adapter removal
# -m 20 = Throw away processed reads shorter than 20 bases
# -e = error tolerance for adapter. 
# -a = Removing Illumina adapters. Prefix of the adapter sequence that is common to all indexed adaptors
#####
export PATH=/work/users/dgerard/programs/cutadapt-1.8.1/bin:$PATH

## Adapters removal and quality trimming for the first biological replicate

cd /scratch/users/dgerard/RNA-seq-rep1/RawData/

for i in *.gz
do
	echo "Cutadapt on file $i"
	# Create a name for the output
	outTrim=Trim.$i
	
	# Create a name for the output stats
	outStats=cutadaptLog.$i

cutadapt -q 30 -m 20 -e 0.1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $i 2> $outStats > $outTrim
done 

## Adapters removal and quality trimming for the second biological replicate

cd /scratch/users/dgerard/RNA-seq-rep2/RawData/

for i in *.gz
do
        echo "Cutadapt on file $i"
        # Create a name for the output
        outTrim=Trim.$i

        # Create a name for the output stats
        outStats=cutadaptLog.$i

cutadapt -q 30 -m 20 -e 0.1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $i 2> $outStats > $outTrim
done

## Adapters removal and quality trimming for the third biological replicate
## Sample TC3-A-D5 contains contaminants that have been detected as mitochondrial sequences (via a blast)
## DO NOT REMOVE THEM, they will be mapped to the mitochondrial chromosome anyway.

cd /scratch/users/dgerard/RNA-seq-rep3/RawData/

for i in *.gz
do
        echo "Cutadapt on file $i"
        # Create a name for the output
        outTrim=Trim.$i

        # Create a name for the output stats
        outStats=cutadaptLog.$i

cutadapt -q 30 -m 20 -e 0.1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $i 2> $outStats > $outTrim
done
