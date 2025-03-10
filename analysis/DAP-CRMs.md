# Annotation of DAP-CRMs and analysis

#### B73v5 DAP-CRMs were generated using the 66 TF diversity panel TFs.
#### Peak summits within 300bp of each other were merged and those with three or more TF binding events were selected to form DAP-CRM as follows:

#### label all B73v5 .narrowPeak files from panel with name of TF  
#example shown for ABI19  
```
awk '{OFS="\t"; $6 = "ABI19"; print $0}' ABI19_B73v5_Q30_default_finalBl.GEM_events.narrowPeak > TF_labeled_files/ABI19_B73v5_default.narrowPeak
```

#### create .narrowPeak file with only summit; repeat for all files in TF diversity panel  
```
awk '{OFS="\t";size=$3-$2;half=size/2;summit=half+$2;summit_int=int(summit);print $1,summit_int,summit_int+1,$4,$5,$6,$7,$8,$9,$10}' TF_labeled_files/ABI19_B73v5_default.narrowPeak > summit_TF_labeled_files/ABI19_TFlabeled_summit.narrowPeak
```

#### create file with DAP-CRMs consisting of three or more TFs
```
cat *summit.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 300 -c 2,6 -o count,collapse | awk -F "\t" '{ if($4 > 2) { print $0 } }' > TF_diversityPanel66_summit_merged300bp_3ormoreTFs.txt
```

#### calculate size of DAP-CRMs  
```
awk '{ sum += $6 } END { if (NR > 0) print sum / NR }' TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMsize.txt
```
`343.663`  
#### range of sizes of DAP-CRMs
```
cut -f 6 TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMsize.txt | sort -k1,1n | uniq -c
```

#### total bp of DAP-CRMS:  
`77404941`  

#### TOTAL % of genome  
`77404941/sum chrs`  
`77404941/2131846805 = 3.6%`  

#### calculate avg number of TFs per DAP-CRM
```
awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMsize.txt
```
`5.27965 TFs per DAP-CRM`  

#### counts of min and max number of TFs per DAP-CRM
```
cut -f 4 TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMsize.txt | sort -k1,1n | uniq -c
```      
#### run Chipseeker on TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed using standard setttings to determine overlap with gene features

#### label each DAP-CRM with unique ID
```
cat -n TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMsize.txt | awk '{OFS="\t"; $1 = "CRM3TF_"$1; print $2,$3,$4,$1,$5,$6,$7}' - > TF_diversityPanel66_summit_merged300bp_3ormoreTFs_CRMannod_CRMsize.txt 
```

### overlap with orthogonal datasets and comparison to randomly shuffled DAP-CRMs  
#### create randomly shuffled DAP-CRM dataset
```
bedtools shuffle -excl B73v5_blocklist6_3-1-2022_nooverlaps2.bed -i TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -g B73v5_chr_sizes.txt > shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed
```

#### count overlap of DAP-CRMs with published ATAC-seq
```
bedtools intersect -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_ATAC_Ear_rep1_and_leafNAM_merged.bed -u | wc -l
```
`38521`  
#### count shuffled
```
bedtools intersect -a shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_ATAC_Ear_rep1_and_leafNAM_merged.bed -u | wc -l
```
`6076`  

#### count overlap of DAP-CRMs with MOA-seq
```
bedtools intersect -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b MOA-seq_coverage_peaks_chr1-10_ear_mg.bed -u | wc -l
```
`24420`  
#### count shuffled  
```
bedtools intersect -a shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed -b MOA-seq_coverage_peaks_chr1-10_ear_mg.bed -u | wc -l
```
`3836`  

#### count overlap of DAP-CRMs with UMRs
```
bedtools intersect -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_UMRs_MGDB.bed -u | wc -l
```
`95104`  
#### count shuffled
```
bedtools intersect -a shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_UMRs_MGDB.bed -u | wc -l
```
`23130`  

#### count overlap of DAP-CRMs with CNS
```
bedtools intersect -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b sorghum_CNS_B73v4tov5_ensembl_chr.bed -u | wc -l
```
`58314`  
#### count shuffled
```
bedtools intersect -a shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed -b ../../sorghum_CNS_B73v4tov5_ensembl_chr.bed -u | wc -l
```
`14227`  

#GWAS hits data from Wallace et al., 2014 downloaded from maizeGDB B73v5
#### count overlap of DAP-CRMs with GWAS hits
```
bedtools intersect -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_GWAS_SNPs_chr1-10_mbdg.bed -u | wc -l
```
`5109`  
#### count shuffled
```
bedtools intersect -a shuffled_TFdiversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_GWAS_SNPs_chr1-10_mbdg.bed -u | wc -l
```
`1470`  

#### check for significance of overlap with each set of orthogonal data (i.e. ATAC)  
#### Must use bedtools2.30 version or greater
```
bedtools fisher -a TF_diversityPanel66_summit_merged300bp_3ormoreTFs.bed -b B73v5_ATAC_Ear_rep1_and_leafNAM_merged.bed -g B73v5_chr_sizes.txt
```



