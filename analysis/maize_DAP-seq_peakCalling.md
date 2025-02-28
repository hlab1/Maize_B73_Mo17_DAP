# Mapping maize DAP-seq reads to the B73v5 genome  

#### analysis is typically done on an HPC  
#### example is shown for BHLH57  

## Trim reads using Trimmomatic  
```
java -jar /home/mjg328/software/trimmomatic-0.39.jar PE -threads 16 -phred33 BHLH57_B73_R1.fastq.gz BHLH57_B73_R2.fastq.gz BHLH57_B73_R1_paired.fastq.gz BHLH57_B73_R1_unpaired.fastq.gz BHLH57_B73_R2_paired.fastq.gz BHLH57_B73_R2_unpaired.fastq.gz ILLUMINACLIP:/home/mjg328/software/Trimmomatic/Trimmomatic-main/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```  

## Map trimmed reads to B73v5 genome using bowtie2 with default parameters (you will need to have a bowtie2 index already generated for B73v5)  

```
bowtie2 -x B73_NAM_v5_chr/B73ref_NAM_v5 -1 BHLH57_B73_R1_paired.fastq.gz -2 BHLH57_B73_R2_paired.fastq.gz -U BHLH57_B73_R1_unpaired.fastq.gz,BHLH57_B73_R2_unpaired.fastq.gz -S BHLH57_B73_B73v5_pairAndUnpair.sam -p 16 >& BHLH57_B73_B73v5_bowtie2_pairAndUnpair.log
```
#### convert .sam file to .bam file  
```
samtools view -bS BHLH57_B73_B73v5_pairAndUnpair.sam | samtools sort -o BHLH57_B73v5.bam
```
#### filter the mapped reads to keep only those with a MAPQ score greater than 30  
```
samtools view -b -q 30 BHLH57_B73v5.bam > BHLH57_B73v5_Q30.bam
```
```
samtools index BHLH57_B73v5_Q30.bam
```

## Call peaks using GEM3 and a threshold of 0.00001 (--q 5)  
```
java -Xmx10G -jar /apps/gem/3.3/gem.jar --d /apps/gem/3.3/Read_Distribution_default.txt --g B73_NAM_v5_chr/B73v5_chr_sizes.txt --genome B73_NAM_v5_chr/chr_files_GEM --expt BHLH57_B73v5_Q30.bam --ctrl HALO-GST_NAMv5_accepted_hits_Q30.bam --ex B73v5_blocklist6_3-1-2022.txt --f SAM --out BHLH57_B73v5_Q30_qvale5_finalBl --t 16 --k_min 5 --k_max 14 --outNP --outMEME --q 5 --k_neg_dinu_shuffle --print_bound_seqs --print_aligned_seqs
```

## Create a .bigwig file that can viewed in a genome browser such as JBrowse or IGV. Exlude any reads falling within regions assigned to the B73v5 blacklist  
```
bamCoverage --bam BHLH57_B73v5_Q30.bam --binSize 10 -p max --normalizeUsing RPKM --bl B73v5_blocklist6_3-1-2022_nooverlaps2.bed -o BHLH57_B73v5_Q30_bl.bigwig
```


# Mapping maize DAP-seq reads to the Mo17 CAU1.0 genome  

#### example is shown for SBP30  
## Trim reads using Trimmomatic  
```
java -jar /home/mjg328/software/trimmomatic-0.39.jar PE -threads 16 -phred33 SBP30_Mo17_R1.fastq.gz SBP30_Mo17_R2.fastq.gz SBP30_Mo17_R1_paired.fastq.gz SBP30_Mo17_R1_unpaired.fastq.gz SBP30_Mo17_R2_paired.fastq.gz SBP30_Mo17_R2_unpaired.fastq.gz ILLUMINACLIP:/home/mjg328/software/Trimmomatic/Trimmomatic-main/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```

## Map trimmed reads to Mo17 genome using bowtie2 with default parameters (you will need to have a bowtie2 index already generated for Mo17 CAU1.0)  

```
bowtie2 -x ZmMo17_CAU1.0 -1 SBP30_Mo17_R1_paired.fastq.gz -2 SBP30_Mo17_R2_paired.fastq.gz -U SBP30_Mo17_R1_unpaired.fastq.gz,SBP30_Mo17_R2_unpaired.fastq.gz -S SBP30_Mo17_pairAndUnpair.sam -p 16 >& SBP30_Mo17_bowtie2_pairAndUnpair.log
```
#### convert .sam file to .bam file  
```
samtools view -bS SBP30_Mo17_pairAndUnpair.sam | samtools sort -o SBP30_Mo17.bam
```
#### filter the mapped reads to keep only those with a MAPQ score greater than 30  
```
samtools view -b -q 30 SBP30_Mo17.bam > SBP30_Mo17_Q30.bam
```
```
samtools index SBP30_Mo17_Q30.bam
```

## Call peaks using GEM3 and a threshold of 0.00001 (--q 5)  
```
java -Xmx10G -jar /apps/gem/3.3/gem.jar --d /apps/gem/3.3/Read_Distribution_default.txt --g Zm-Mo17-REFERENCE-CAU-1.0_chr_sizes.txt --genome Zm-Mo17-REFERENCE-CAU-1.0.id_chr_nuc --expt *Q30.bam --ctrl HALO_GST_Mo17CAU_accepted_hits_Q30.bam --ex Mo17_blocklist6_3-7-22.txt --f SAM --out SBP30_Mo17_Q30_qvale5_finalBl --t 16 --q 5 --k_min 5 --k_max 14 --outNP --outMEME --k_neg_dinu_shuffle --print_bound_seqs --print_aligned_seqs
```

## Create a .bigwig file that can viewed in a genome browser such as JBrowse or IGV. Exlude any reads falling within regions assigned to the Mo17 blacklist  
```
bamCoverage --bam SBP30_Mo17_Q30.bam --binSize 10 -p max --normalizeUsing RPKM --bl Mo17_blocklist6_3-7-22_nooverlaps.bed -o SBP30_Mo17_Q30.bam_bl.bigwig
```
