# Comparative-analysis-of-TF-binding-in-B73-and-Mo17

# General workflow  

## Normalize mapped reads for B73 and Mo17 DAP-seq files  
#### Check non-normalized B73v5 and Mo17 bam files for number of mapped reads using samtools flagstat (check for each single TF)  
```
samtools flagstat ARF16_B73v5_Q30.bam
```  
```
samtools flagstat ARF16_Mo17_Q30.bam
```  

#### If read numbers are more than >90% similar, use as is (no downsampling needed)

#### If datasets differ by more than 90% mapped reads, downsample reads with samtools -s so that higher dataset has same number of reads as lower dataset (sometimes B73 needs to be downsampled, sometimes Mo17)
```
samtools view -b -s 0.75 ARF16_B73v5_Q30.bam > ARF16_B73v5_Q30_downsampled.bam
```

#### run GEM3 peak caller on downsampled dataset using a threshold similar to non-downsampled dataset with the goal being to get a similar number of peaks for the matching B73 and Mo17 TFs.  
```
java -Xmx10G -jar /apps/gem/3.3/gem.jar --d /apps/gem/3.3/Read_Distribution_default.txt --g B73v5_chr_sizes.txt --genome B73_NAM_v5_chr/chr_files_GEM --expt ARF16_B73v5_Q30_bamDownsampled.bam --ctrl /scratch/mjg328/DAP_all_B73v5_mapping/HALO-GST_NAMv5_accepted_hits_Q30.bam --ex /scratch/mjg328/DAP_all_B73v5_mapping/B73v5_blocklist6_3-1-2022.txt --f SAM --out ARF16_B73v5_Q30_bamDownsampled_qvale5_finalBl --q 5 --t 16 --k_min 5 --k_max 14 --outNP --outMEME --k_neg_dinu_shuffle --print_bound_seqs --print_aligned_seqs
```

## Run crossmap on normalized narrowPeak files to convert B73 peaks to Mo17 coordinates and Mo17 peaks to B73v5 coordinates  
> First need to create whole genome alignment of B73v5 and Mo17 and generate .chain file
> See file 'Anchorwave_B73_Mo17_alignment'

#### Crossmap analysis  
> https://crossmap.readthedocs.io/en/latest/  
> See file 'Crossmap_analysis_B73_Mo17'

## Identify shared and genotype-specific peaks

#### Using top 20% of B73 peaks (or Mo17), run `bedtools intersect` with B73peaks and Mo17peaksWithB73coordinates to identify peaks that are not present in Mo17 (or vice versa)  
#### Identify peaks that are unique and shared  

