# Identification of shared and genotype-specific peaks

#### To compare peaks from B73 and Mo17, need to have the peaks coordinates correspond to the same genome. This can be done by converting Mo17 peak coordinates to B73v5 (and vice versa). 

#### Use Crossmap to liftover coordinates  
>https://crossmap.readthedocs.io/en/latest/

#### Example using ARF16
###### using bam normalized and peak threshold-adjusted peaks, confirm that number of peaks for each genome is roughly similar  
39784 ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events.narrowPeak  
32750 ARF16_Mo17_Q30_qvale5_finalBl.GEM_events.narrowPeak

#### prepare files for running Crossmap to convert coordinates between B73v5 and Mo17 and vice versa  
#### create four column bed files from narrowPeak file  
```
awk -v OFS='\t' '{print $1, $2, $3, $4}' ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events.narrowPeak > ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events_4col.bed
```  

```
awk -v OFS='\t' '{print $1, $2, $3, $4}' ARF16_Mo17_Q30_qvale5_finalBl.GEM_events.narrowPeak > ARF16_Mo17_Q30_qvale5_finalBl.GEM_events_4col.bed
```  

#### run crossMap coordinate conversion (B73v5 to Mo17; chain file name is misleading but is correct file for this orientation)  
```
CrossMap.py bed --chromid a --unmap-file ARF16_B73v5_Q30_downsampled_qvale5_to_Mo17coords.crossmap.unmap.bed anchorwave_Mo17toB73v5.chain ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events_4col.bed ARF16_B73v5_Q30_downsampled_qvale5_to_Mo17coords.crossmap.bed
```

#### run crossMap coordinate conversion (Mo17 to B73v5; chain file name is misleading but is correct file for this orientation)  
```
CrossMap.py bed --chromid a --unmap-file ARF16_Mo17_Q30_qvale5_to_B73v5coords.crossmap.unmap.bed anchorwave_Mo17toB73v5_swap.chain ARF16_Mo17_Q30_qvale5_finalBl.GEM_events_4col.bed ARF16_Mo17_Q30_qvale5_to_B73v5coords.crossmap.bed
```

#### Identify B73-specific peaks (i.e. not in Mo17) and B73-Mo17-shared peaks (i.e. present in both B73 and Mo17)  

#calculate 20% of total B73 peak dataset (this ensures only highest quality peaks are used in analysis)  
```
python -c 'print(0.2*39784)'
```  
`7956`

#### Create file with B73v5-specific peaks (top 20% peaks):  
```
head -7956 ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events.narrowPeak | bedtools intersect -a stdin -b ARF16_Mo17_Q30_qvale5_to_B73v5coords.crossmap.bed -wao | awk '$15 == 0' > ARF16_B73v5-specific.narrowPeak
```

#### Create file with B73v5-Mo17-shared peaks (top 20% peaks):  
```
head -7956 ARF16_B73v5_Q30_downsampled_qvale5_finalBl.GEM_events.narrowPeak | bedtools intersect -a stdin -b ARF16_Mo17_Q30_qvale5_to_B73v5coords.crossmap.bed -wao | awk '$15 != 0' | awk '!A[$4]++' > ARF16_B73v5-shared.narrowPeak
```

## Identify Mo17-specific and Mo17-B73-shared peaks

#### calculate how many Mo17 peaks are unique (i.e. not in B73) and Mo17-B73-shared peaks (i.e. present in both Mo17 and B73)

#calculate 20% of total Mo17 peak dataset (this ensures only highest quality peaks are used in analysis)  
```
python -c 'print(0.2*32750)'
```  
`6550`

#### Create file with Mo17-specific peaks (top 20% peaks):  
```
head -6550 ARF16_Mo17_Q30_qvale5_finalBl.GEM_events.narrowPeak | bedtools intersect -a stdin -b ARF16_B73v5_Q30_downsampled_qvale5_to_Mo17coords.crossmap.bed -wao | awk '$15 == 0' > ARF16_Mo17-specific.narrowPeak
```

#### Create file with Mo17-B73-shared peaks (top 20% peaks):  
```
head -6550 ARF16_Mo17_Q30_qvale5_finalBl.GEM_events.narrowPeak | bedtools intersect -a stdin -b ARF16_B73v5_Q30_downsampled_qvale5_to_Mo17coords.crossmap.bed -wao | awk '$15 != 0' | awk '!A[$4]++' > ARF16_Mo17_shared.narrowPeak
```




