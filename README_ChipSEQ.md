## 150601 Chip-seq mouse time series

Enquirer: Déborah Gérard / Lasse Sinkkonen  
Date: 20150112

### Input files
TC1 = Time course 1  
I = Input  
ST2 = Cell line  
D = day

3 histone modifications:  
1. H3K27ac
2. H3K4m3
3. H3K36m3

TruSeq Illumina sequencing


```
# copy files
rsync -av --progress */*gz /home/ginolhac/Projects/150601/TC1_ST2_ChIP-seq/
# check copying process
rsync -av --progress */*.md5 /home/ginolhac/Projects/150601/TC1_ST2_ChIP-seq/
parallel "md5sum {} > {.}.md5" ::: */*gz
cd /home/ginolhac/Projects/150601/TC1_ST2_ChIP-seq/
parallel "md5sum -c {}" ::: *.md5
```

### Reference, mouse genome
GRCm38.p3 (mm10, patch 3)


```
wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/Primary_Assembly/assembled_chromosomes/FASTA/*fa.gz
# concatenate chrs in correct order
for i in `seq 19; echo X; echo Y`
    do echo $i 
    zcat chr${i}.fa.gz >> GRCm38.p3.fasta
done 
sed -i 's/^>.*chromosome.\([0-9YX]*\),/>\1/' GRCm38.p3.fasta
```

Mitochondrial genome
```
http://www.ncbi.nlm.nih.gov/nuccore/NC_005089.1
```

### Run paleomix


#### Makefile

```
for f in *gz
  do IFS='_' read -a array <<< "$f"
  echo -en "${array[1]}:\n  ${array[1]}:\n    ${array[1]}:\n" >> tmp
  echo -en "      \"${array[2]}\": /home/ginolhac/Projects/150601/TC1_ST2_ChIP-seq/${f}\n" >> tmp
done
cat tmp >> ../mouse.makefile
```


```
nice -n +19 paleomix bam_pipeline run --bwa-max-threads=4 --max-threads=6 mouse.makefile
```

last for 28 hours


#### check for correct trimming

```
AdapRem AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
Fastqc  AGATCGGAAGAGCACACGTCTGAACTCCAGTCACG
TruSeq  GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCC
Index_input                              GATCAG

```

### Remove MT reads in the nuclear genome

```
for f in *q30.bam; 
  do ./filterConta2.py \
  $f ${f/GRCm38.p3.q30.bam/MT.bam} strict \
  > ${f/GRCm38.p3.q30.bam/GRCm38.p3.q30.noMT.bam} ;
done
```


### Peak calling

#### MACS2

mm for mus musculus!

Complete version

```
# process reads
for FILE in ../*H3K4*q30.bam
  do IFS='.' read -a array <<< "$FILE"
  TREAT=${array[2]}; TRE=${TREAT/\//}
  macs2 callpeak -t $FILE -c ${FILE/H3K4/I} -f BAM -g mm -n $TRE -B -q 0.01 --outdir $TRE
done 

```

