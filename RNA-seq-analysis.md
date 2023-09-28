## Quality control 

```
#Fastqc
module load fastqc
for name in $sample
do
fastqc $name
done
```

```
#Adapter trimming
for name in $sample
do
fastp -i $name.R1.fq.gz -I $name.R2.fq.gz -o trimmed_$name.R1.fq.gz -O trimmed_$name.R2.fq.gz --detect_adapter_for_pe -l 20 -j $name.fastp.json -h $name.fastp.html
done
```

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

```
#Feauture count 
nice R 
library(Rsubread)
fc_PEO_RNA <- featureCounts(files=list.files(pattern='.bam$'),annot.ext='hg38_ucsc/hg38.refGene.gtf',isGTFAnnotationFile=TRUE, nthread=32, GTF.featureType="exon", GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=T)
head(fc_PEO_RNA$counts)
write.table(x=data.frame(fc_PEO_RNA$counts,stringsAsFactors=FALSE), file="fc_PEO_RNA.txt",quote=FALSE,sep="\t",row.names=TRUE)

#DESeq2
library(DESeq2)
fc_PEO_RNA <- read.delim("fc_PEO_RNA.txt", header = TRUE, sep = "\t")
condition <- c("PEO1", "PEO1", "PEO1", "PEO4", "PEO4", "PEO4")
dds <- DESeqDataSetFromMatrix(countData=fc_PEO_RNA, colData=DataFrame(condition), design=~condition)

#remove rows that have <10 counts 
dds <- dds[rowSums(counts(dds))>=10,]

#conduct differential analysis
dds <- DESeq(dds)
res <- results(dds)

#get the name of the comparisons made aka coefficients
resultsNames(dds) 
 
#Shrink fold change - coef = name from above step 
res.shr <- lfcShrink(dds,coef="condition_PEO4_vs_PEO1")

#Remove rows with missing adj p value
resFix <- res.shr[!is.na(res.shr$padj),]

#Export data 
write.csv(as.data.frame(resFix), file="geneID_RNA_PEO4_v_PEO1.csv")

```

