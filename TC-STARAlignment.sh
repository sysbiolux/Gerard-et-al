#!/bin/bash

## Use STAR 2.5.2b to align RNA-Seq reads to the genome 
export PATH=/work/users/dgerard/programs/STAR-2.5.2b/bin/Linux_x86_64_static:$PATH

## Align 1st biological replicate ##

cd /scratch/users/dgerard/RNA-seq-rep1/ReadsworRNA
echo "Alignment to the genome with STAR 2.5.2b"

for i in  *worRNAreads*
do
	echo "STAR alignment for file $i"
	# Create names for the output
	outSTAR=Mapped.$i

	# Run the STAR command
	
	STAR --genomeDir /scratch/users/dgerard/Ahr-Gzf1-KD/STAR.GRCm38mm10.anno79 --runThreadN 12 --readFilesIn \
	$i --readFilesCommand zcat --outFileNamePrefix /scratch/users/dgerard/RNA-seq-rep1/MappedReads/$outSTAR --twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed \
	1000000 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --outFilterMismatchNoverLmax 0.3 --seedSearchStartLmax 12 --alignSJoverhangMin \
	15 --alignEndsType Local --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --winAnchorMultimapNmax 50 --alignSJDBoverhangMin 3 --outFilterType \
	BySJout --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
done

## Align 2nd  biological replicate ##

cd /scratch/users/dgerard/RNA-seq-rep2/ReadsworRNA

for i in  *worRNAreads*
do
       echo "STAR alignment for file $i"
        # Create names for the output
        outSTAR=Mapped.$i

        # Run the STAR command

	STAR --genomeDir /scratch/users/dgerard/Ahr-Gzf1-KD/STAR.GRCm38mm10.anno79 --runThreadN 12 --readFilesIn \
        $i --readFilesCommand zcat --outFileNamePrefix /scratch/users/dgerard/RNA-seq-rep2/MappedReads/$outSTAR --twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed \
        1000000 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --outFilterMismatchNoverLmax 0.3 --seedSearchStartLmax 12 --alignSJoverhangMin \
        15 --alignEndsType Local --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --winAnchorMultimapNmax 50 --alignSJDBoverhangMin 3 --outFilterType \
        BySJout --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif

done

## Align 3rd  biological replicate ##

cd /scratch/users/dgerard/RNA-seq-rep3/ReadsworRNA

for i in  *worRNAreads*
do
        echo "STAR alignment for file $i"
        # Create names for the output
        outSTAR=Mapped.$i

        # Run the STAR command

        STAR --genomeDir /scratch/users/dgerard/Ahr-Gzf1-KD/STAR.GRCm38mm10.anno79 --runThreadN 12 --readFilesIn \
        $i --readFilesCommand zcat --outFileNamePrefix /scratch/users/dgerard/RNA-seq-rep3/MappedReads/$outSTAR --twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed \
        1000000 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 100 --outFilterMismatchNmax 33 --outFilterMismatchNoverLmax 0.3 --seedSearchStartLmax 12 --alignSJoverhangMin \
        15 --alignEndsType Local --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0.3 --winAnchorMultimapNmax 50 --alignSJDBoverhangMin 3 --outFilterType \
        BySJout --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
done
