## Quality control 


## Aligning 

```
#Build genome index 
module load star
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir path_to_genome --genomeFastaFiles hg38_ucsc/hg38.fa --sjdbGTFfile hg38_ucsc/hg38.refGene.gtf --sjdbOverhang 149

#Align reads
module load star
module load samtools 

for name in $sample
do 
STAR --runThreadN 12 \
--genomeDir path_to_genome \
--sjdbGTFfile hg38_ucsc/hg38.refGene.gtf \
--sjdbOverhang 149 \
--outFileNamePrefix $name \
--readFilesCommand gunzip -c \
--readFilesIn $name.R1.fq.gz $name.R2.fq.gz \
--outSAMtype BAM SortedByCoordinate
done

```



## Differential gene expression 
