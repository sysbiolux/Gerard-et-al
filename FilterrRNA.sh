#!/bin/bash

## Use SortMeRNA to filter rRNA reads from RNA-Seq samples 
export PATH=/work/users/dgerard/programs/sortmerna-2.0:$PATH

## Remove any ribosomal reads from the 1st biological replicate

cd /scratch/users/dgerard/RNA-seq-rep1/TrimmedData

echo "Running Sortmerna-2.0"
for i in Trim*
do
	echo "rRNA reads removal from file $i"
	# Create names for the file containing ribosomal RNA reads
	outwrRNAreads=$i.wrRNAreads

	# Create names for the file containing reads WITHOUT rRNA reads
	outwoRNAreads=$i.worRNAreads

	# Run the sortmerna command	
	sortmerna --ref /work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-18s-id95.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk18S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-28s-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk28S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5s-db:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5.8s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5.8s-db: --reads \
	$i --aligned ./../ReadswrRNA/$outwrRNAreads --other ./../ReadsworRNA/$outwoRNAreads --log -a 4 -v --fastx

	# gzip files to save space
	gzip $i

done

## Remove any ribosomal reads from the 2nd biological replicate

cd /scratch/users/dgerard/RNA-seq-rep2/TrimmedData
gunzip *.gz
echo "Running Sortmerna-2.0"
for i in Trim*
do
	# Unzip files
	#gunzip $i

	echo "rRNA reads removal from file $i"
	# Create names for the file containing ribosomal RNA reads
	outwrRNAreads=$i.wrRNAreads

	# Create names for the file containing reads WITHOUT rRNA reads
	outwoRNAreads=$i.worRNAreads

	# Run the sortmerna command
	sortmerna --ref /work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-18s-id95.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk18S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-28s-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk28S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5s-db:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5.8s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5.8s-db: --reads \
	$i --aligned ./../ReadswrRNA/$outwrRNAreads --other ./../ReadsworRNA/$outwoRNAreads --log -a 12 -v --fastx

	# gzip files to save space
	gzip $i

done

## Remove any ribosomal reads from the 3rd biological replicate

cd /scratch/users/dgerard/RNA-seq-rep3/TrimmedData
gunzip *.gz
echo "Running Sortmerna-2.0"
for i in Trim*
do
	# Unzip files
	#gunzip $i

	echo "rRNA reads removal from file $i"
	# Create names for the file containing ribosomal RNA reads
	outwrRNAreads=$i.wrRNAreads

	# Create names for the file containing reads WITHOUT rRNA reads
	outwoRNAreads=$i.worRNAreads

	# Run the sortmerna command
	sortmerna --ref /work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-18s-id95.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk18S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/silva-euk-28s-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/silva-euk28S_DB:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5s-db:/work/users/dgerard/programs/sortmerna-2.0/rRNA_databases/rfam-5.8s-database-id98.fasta,/work/users/dgerard/programs/sortmerna-2.0/index/rfam-5.8s-db: --reads \
	$i --aligned ./../ReadswrRNA/$outwrRNAreads --other ./../ReadsworRNA/$outwoRNAreads --log -a 12 -v --fastx

	# gzip files to save space
	gzip $i

done
