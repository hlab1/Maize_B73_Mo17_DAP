## Calculate overlap of shared and specific DAP-seq peaks with SV, SNPs, indels etc.

#### Use SyRI to compare structural features between B73v5 and Mo17 CAU1.0  
> https://github.com/schneebergerlab/syri

#### convert anchorwave .maf to .sam for analysis with syri to detect structural variation, snps, and indels  
```
python2 maf-convert sam anchorwave_Mo17toB73v5.maf | sed 's/[0-9]\+H//g' > anchorwave_Mo17toB73v5.sam  
```
#### run syri to create .vcf file  
```
syri -c anchorwave_analysis_B73v5/anchorwave_Mo17toB73v5.sam -F S --prefix anchorwave_Mo17toB73v5_sam_ --cigar -f --log DEBUG -r Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa -q Zm-Mo17-REFERENCE-CAU-1.0.id_chr_nuc_mg.fa
```
#### Parse .vcf file to identify SNPs, INDELs (<50bp), SVs (INDELs greater than 50bp), and DUPs (duplication regions)  
```
grep -E '#|INS' anchorwave_Mo17toB73v5_sam_syri.vcf > anchorwave_Mo17toB73v5_sam_syri_INSs.vcf
```
```
grep -E '#|DEL' anchorwave_Mo17toB73v5_sam_syri.vcf > anchorwave_Mo17toB73v5_sam_syri_DELs.vcf
```

#create file with Insertions less than 50bp or more than 50bp; column 5 shows actual nucleotide seq of insertion so can count how long using 'length'  
```
awk '{ if (length($5) < 50) print }' anchorwave_Mo17toB73v5_sam_syri_INSs.vcf > anchorwave_Mo17toB73v5_sam_syri_INSs_less50bp.vcf
```
```
awk '{ if (length($5) >= 50) print }' anchorwave_Mo17toB73v5_sam_syri_INSs.vcf > anchorwave_Mo17toB73v5_sam_syri_INSs_more50bp.vcf
```
#create file with Deletions less than 50bp or more than 50bp; column 4 shows actual nucleotide seq of deletion  
```
grep -E '#|DEL' anchorwave_Mo17toB73v5_sam_syri_INSDELs_no0s.vcf > anchorwave_Mo17toB73v5_sam_syri_DELs_no0s.vcf
```
```
grep -E '#|INS' anchorwave_Mo17toB73v5_sam_syri_INSDELs_no0s.vcf > anchorwave_Mo17toB73v5_sam_syri_INSs_no0s.vcf
```
```
awk '{ if (length($4) >= 50) print }' anchorwave_Mo17toB73v5_sam_syri_DELs_no0s.vcf > anchorwave_Mo17toB73v5_sam_syri_DELs_no0s_more50bp.vcf #this file does not have entries with a 0 start position which causes a bedtools error
```
```
awk '{ if (length($5) >= 50) print }' anchorwave_Mo17toB73v5_sam_syri_INSs_no0s.vcf > anchorwave_Mo17toB73v5_sam_syri_INSs_no0s_more50bp.vcf #this file does not have entries with a 0 start position which causes a bedtools error
```
#### Combine files for Insertions and Deletions  
##### small indels
```
cat anchorwave_Mo17toB73v5_sam_syri_DELs_less50bp.vcf anchorwave_Mo17toB73v5_sam_syri_INSs_less50bp_noHeader.vcf > anchorwave_Mo17toB73v5_sam_syri_INSDELs_less50bp.vcf
```
##### Structural Variants (SVs)
```
cat vcf_header.txt anchorwave_Mo17toB73v5_sam_syri_DELs_no0s_more50bp.vcf anchorwave_Mo17toB73v5_sam_syri_INSs_no0s_more50bp.vcf > anchorwave_Mo17toB73v5_sam_syri_INSDELs_no0s_more50bp.vcf
```

#### Create file with SNPs  
```
grep -E '#|SNP' anchorwave_Mo17toB73v5_sam_syri.vcf > anchorwave_Mo17toB73v5_sam_syri_SNPs.vcf
```

#### Create file with Duplications (DUPs) and Inverted duplications (INVDP) regions  
```
grep -E '<DUP>|<INVDP>' anchorwave_Mo17toB73v5_sam_syri.vcf > anchorwave_Mo17toB73v5_sam_syri_DUP_INVDP_B73coords.bed
```
#### modify DUP_INVDP_B73coords.bed file to create four column .bed file containing  
`chr  start  end (from END=)  name`  

#### Create file with Not Aligned Regions (NOTAL)  
```
grep 'NOTAL' anchorwave_Mo17toB73v5_sam_syri.vcf > anchorwave_Mo17toB73v5_sam_syri_NOTAL_B73coords.bed
```
#### modify NOTAL_B73coords.bed file to create four column .bed file containing  
`chr  start  end (from END=)  name`  

## Check overlap of peaks with SNPs, INDELs, and SVs  
#### example shown for ARF16 shared and B73-specific peaks (top 20% of peaks)    

#### check overlap with SNPs
```
bedtools intersect -a ARF16_B73v5-specific.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_SNPs.vcf -wa -c | awk '$16 !=0' | awk '$16 !=201' | wc -l
```
```
bedtools intersect -a ARF16_B73v5_shared_nodups.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_SNPs.vcf -wa -c | awk '$16 !=0' | awk '$16 !=201' | wc -l
```

#### check overlap with small indels  
```
bedtools intersect -a ARF16_B73v5-specific.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_INSDELs_less50bp.vcf -wa -c | awk '$16 !=0' | wc -l
```
```
bedtools intersect -a ARF16_B73v5_shared_nodups.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_INSDELs_less50bp.vcf -wa -c | awk '$16 !=0' | wc -l
```

#### check overlap with SVs  
```
bedtools intersect -a ARF16_B73v5-specific.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_INSDELs_no0s_more50bp.vcf -wa -c | awk '$16 !=0' | wc -l
```
```
bedtools intersect -a ARF16_B73v5_shared_nodups.narrowPeak -b anchorwave_Mo17toB73v5_sam_syri_INSDELs_no0s_more50bp.vcf -wa -c | awk '$16 !=0' | wc -l
```

#### analysis done for all 200 TFs  
#### MG *_shared and *_specific files use awk $16; SCH *_shared and *_specific files use awk $11
