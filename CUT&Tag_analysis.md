# Raw data processing 

## Merging technical replicates 

```
for name in $sample
do
cat $name*.R1.fq.gz > merged/$name.R1.fq.gz
cat $name*.R2.fq.gz > merged/$name.R2.fq.gz
done
```

## Fastqc 
```
module load fastqc
for name in $sample
do
fastqc $name
done
```

## Trimming adapter sequence 
```
for name in $sample
do
fastp -i $name.R1.fq.gz -I $name.R2.fq.gz -o trimmed_$name.R1.fq.gz -O trimmed_$name.R2.fq.gz --detect_adapter_for_pe -l 20 -j $name.fastp.json -h $name.fastp.html
done
```
## Aligning reads 

```
#Download reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

#Build genome index
module load anaconda3/personal
module load bowtie2
bowtie2-build hg38.fa hg38

#Align reads
for name in $sample
do 
bowtie2 -x /genome_directory/hg38 -1 trimmed_$name.R1.fq.gz -2 trimmed_$name.R2.fq.gz  -S $name.sam  2> summary_$name.txt
done

#Filter unmapped, secondary, chimeric, low quality (q<10), blacklisted and duplicated reads
module load samtools
module load bedtools
module load picard
for file in *.sam
do
name=`basename $file .sam` 
samtools view -F 2308 -b -q 10 $name.sam > $name.bam
samtools sort -o $name.sort.bam -O 'bam' -T temp_$name $name.bam
picard MarkDuplicates quiet=TRUE REMOVE_DUPLICATES=TRUE I=$name.sort.bam O=$name.dedup.sort.bam M=$name.metrics.txt
bedtools intersect -v -abam $name.dedup.sort.bam -b 'blacklisted_regions/hg38-blacklist.v2.bed' > $name.bl.dedup.sort.bam
samtools index $name.bl.dedup.sort.bam
rm $name.bam $name.sort.bam $name.dedup.sort.bam
done

```

# Down-stream processing 

## Sample clustering 

```
library(DiffBind)

g4_DE <- dba(sampleSheet="g4_sample_sheet.csv")
plot(G4_DE)
```

## Peak calling 

```
#MACS-default
for file in *.sort.bam 
do
name=`basename $file .bl.dedup.sort.bam` 
macs2 callpeak -f BAMPE -t $file -n $name
done

#MACS-global
for file in *.sort.bam 
do
name=`basename $file .bl.dedup.sort.bam` 
macs2 callpeak -f BAMPE --nolambda -t $file -n $name
done

#SEARC
#Make bedgraph files from peak file
for name in $sample
do
bamCoverage -of bedgraph -b $name.dedup.bam -o $name.bedgraph
done

#peak calling with SEARC
module load bedtools 
for name in $sample
do
searc='SEACR_1.3.sh'
$searc $name.bedgraph 0.01 non stringent $name.searc
done

```

## Peak intersection 

```
#Find common peaks
module load bedtools
bedtools intersect -u -a file_1.bed -b file_2.bed > file_intersect.bed

#Find distinct peaks
module load bedtools
bedtools intersect -v -a file_1.bed -b file_2.bed > file_intersect.bed

## Genomic annotation
module load homer
annotatePeaks.pl hg38 file.bed > annotated_file.bed 
```

## Super-enhancer calling 
module load anaconda3/personal
module load bedtools
module load samtools
mkdir annotation 
cp hg38_refseq.ucsc annotation 
ROSE_main.py -g hg38 -i file.bed -r file.bam -o sample_SE 

