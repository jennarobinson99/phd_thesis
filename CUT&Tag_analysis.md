# Raw data processing 

## Merging technical replicates 

for name in $sample
do
cat $name*.R1.fq.gz > merged/$name.R1.fq.gz
cat $name*.R2.fq.gz > merged/$name.R2.fq.gz
done

## Fastqc 
for name in $sample
do
fastqc $name
done

## Trimming adapter sequence 
for name in $sample
do
fastp -i $name.R1.fq.gz -I $name.R2.fq.gz -o trimmed_$name.R1.fq.gz -O trimmed_$name.R2.fq.gz --detect_adapter_for_pe -l 20 -j $name.fastp.json -h $name.fastp.html
done



