# nuclearEDNA
Make a conda environment for read trimming
```
conda create -n fastp fastp -c bioconda
conda activate fastp 
 ```
list and make unique reads for parallelizing 

```
ls reads_salmon/*gz | sed 's/\_R..fastq.gz//' | sed 's/reads_salmon\///' | sort | uniq > SalmonSeqIds
ls reads_charr/*gz | sed 's/\_R..fastq.gz//' | sed 's/reads_salmon\///' | sort | uniq > CharrSeqIds
```

Trim reads
```
mkdir trim

#trim salmon
cat  SalmonSeqIds | \
parallel -j 30 \
  'fastp -i salmon_reads/{}_R1.fastq.gz \
  -I salmon_reads/{}_R2.fastq.gz \
  -o trim/{}_R1.trimmed.fastq.gz \
  -O trim/{}_R2.trimmed.fastq.gz \
  --thread 4 '

#trim charr
cat  CharrSeqIds | \
parallel -j 30 \
  'fastp -i reads_charr/{}_R1.fastq.gz \
  -I reads_charr/{}_R2.fastq.gz \
  -o trim/{}_R1.trimmed.fastq.gz \
  -O trim/{}_R2.trimmed.fastq.gz \
  --thread 4 '

conda deactivate
```

Prepare alignment environment and activate

```
conda create -n align samtools htslib bwa bwa-mem2 -c bioconda
conda activate align
```

Index reference genome

```
bwa index reference/CIGENE-ICSASG_v2.fa
samtools faidx reference/CIGENE-ICSASG_v2.fa
```

Align

```
mkdir align

while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
  SMID=$(echo $ind | sed 's/NS.*i5.//') ;
  LBID=$(echo $ind | sed 's/.UDP.*//');
  bwa mem \
  -t 100 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
  reference/CIGENE-ICSASG_v2.fa \
  trim/$ind\_R1.trimmed.fastq.gz  trim/$ind\_R2.trimmed.fastq.gz | \
  samtools sort -o align/$ind\.sorted.bam -T $ind -@ 100 -m 20G ;
  done < SalmonSeqIds


while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
  SMID=$(echo $ind | sed 's/NS.*i5.//') ;
  LBID=$(echo $ind | sed 's/.UDP.*//');
  bwa mem \
  -t 100 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
  reference/CIGENE-ICSASG_v2.fa \
  trim/$ind\_R1.trimmed.fastq.gz  trim/$ind\_R2.trimmed.fastq.gz | \
  samtools sort -o align/$ind\.sorted.bam -T $ind -@ 100 -m 20G ;
  done < CharrSeqIds 

conda deactivate

```
