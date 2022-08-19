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

Deduplicate in parallel

```
cat CharrSeqIds SalmonSeqIds | \
  parallel --jobs 100 'gatk --java-options "-Xmx4G" \
  MarkDuplicates \
  I=align/{}.sorted.bam \
  O=align/{}.deDup.bam M=align/{}_deDupMetrics.txt \
  REMOVE_DUPLICATES=true'

``` 

Index after deduplication
```
while read  ind; 
  do samtools index align/$ind.deDup.bam ; 
done < SalmonSeqIds
```

Depth analysis

```
conda create -n mosdepth mosdepth -c bioconda
conda activate mosdepth

cat SalmonSeqIds |
parallel --jobs 50 'mosdepth -n --fast-mode \
--by 5000 depths/Salmon{} align/{}.deDup.bam '


cat CharrSeqIds |
parallel --jobs 50 'mosdepth -n --fast-mode \
--by 5000 depths/Charr{} align/{}.deDup.bam '

```
Extract per window depths

```
 ls *bed.gz | \
  sed s'/.bed.gz//' | \
  parallel 'gzip -dc {}.bed.gz | cut -f4 > {}.5k.depths ' 

```
Check to ensure same number of lines occur per file, extract chrom and window information, make species-level depth matrix files

```
wc -l *depths

gzip -dc CharrNS.1779.001.UDP0289_i7---UDP0289_i5.PG001_2013_1M.regions.bed.gz | cut -f1,2,3 > locs 

paste locs Charr*depths > Charr.5kdepths
paste locs Salmo*depths > Salmo.5kdepths

```
 Now, go to salmon_charr_depthcompare.R[] to generate Charr_salm_top10_deltacov.bed
