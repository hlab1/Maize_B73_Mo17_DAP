# Positional Variation analysis  
#### Analysis describes identification of B73 and Mo17 conserved peaks whose relative position changes in relation to the transcription start site (TSS)  

#### Target gene file preparation: start by calculating B73-Mo17 shared peaks (see also Crossmap pipeline)  
##### example shown for ARF16  

```
bedtools intersect -a ARF16_B73v5_Q30_downsampled.narrowPeak -b ARF16_Mo17_Q30_to_B73v5coords.crossmap.bed -wao | awk '$15 != 0' | awk '!A[$4]++' | head
```
`confirm that column 15 contains a '1' to indicate that the B73 peak is shared i.e. overlaps with a Mo17 peak converted to B73v5 coordinates`  

#### save output to file
```
bedtools intersect -a ARF16_B73v5_Q30_downsampled.narrowPeak -b ARF16_Mo17_Q30_to_B73v5coords.crossmap.bed -wao | awk '$15 != 0' | awk '!A[$4]++' > ARF16_B73v5_allShared_bamNorm.narrowPeak
```

#### output file is a 15 column .narrowPeak file but Chipseeker will only take a 10 column .narrowPeak; also need to use file with peak summit not total peak range
#### also will need Mo17 peak coordinate info from `ARF16_B73v5_allShared_bamNorm.narrowPeak` file (Step 2 below); cut columns $11-$15 containing Mo17 peak info and create new file (Mo17_5col.bed) for use later  

```
cut -f11-15 ARF16_B73v5_allShared_bamNorm.narrowPeak > allShared_Mo17peakinfo_5col.bed
```

#### make ARF16_B73v5_allShared_bamNorm.narrowPeak file that has only 10 columns for use with Chipseeker and has coordinates for B73v5 peak summit
```
cut -f1-10 ARF16_B73v5_allShared_bamNorm.narrowPeak > ARF16_allShared_bamNorm_10col.narrowPeak
```
#### use above file to create new file with peak summit and 10 columns
```
awk '{OFS="\t";size=$3-$2;half=size/2;summit=half+$2;summit_int=int(summit);print $1,summit_int,summit_int+1,$4,$5,$6,$7,$8,$9,$10}' ARF16_B73v5_allShared_bamNorm_10col.narrowPeak | head -4
```
> chr6	179136084	179136085	6:179136084	1000	.	2027.5	999.00	999.00	100  
> chr2	168390595	168390596	2:168390595	887	.	1773.0	999.00	999.00	100  
> chr3	196929778	196929779	3:196929778	816	.	1614.0	999.00	999.00	100  
> chr5	169357513	169357514	5:169357513	782	.	1538.5	999.00	999.00	100  

```
awk '{OFS="\t";size=$3-$2;half=size/2;summit=half+$2;summit_int=int(summit);print $1,summit_int,summit_int+1,$4,$5,$6,$7,$8,$9,$10}' $file > ARF16_allShared_summit_10col.narrowPeak
```

#### file will be used to call target genes using Chipseeker  

## Step1. Identify shared B73v5 target genes and determine distance of peak from TSS  

#### in R, run chipseeker on ARF16_allShared_summit_10col.narrowPeak

```
library(ChIPseeker)  
library(ChIPpeakAnno) #use for loadDb function.  
library(GenomicRanges)  
library(GenomicFeatures)  

#B73v5 db  
setwd("/Users/marygalli/Desktop/")  
datadirGFFv5 <- "/Users/marygalli/Desktop/"  
datafileGFFv5 <- paste(datadirGFFv5, "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", sep="/")  
txdb.B73v5_gff <- makeTxDbFromGFF(file=datafileGFFv5, format="auto")  

#put all *_allShared_summit_10col.narrowPeak files (i.e. ARF16_allShared_summit_10col.narrowPeak) in same directory  
setwd("/Users/marygalli/Desktop/PosV_2-25-25_allSharedPeaks_analysis/B73v5_allShared_SUMMIT_10col_narrowPeaks/")  
filename <- system("ls /Users/marygalli/Desktop/PosV_2-25-25_allSharedPeaks_analysis/B73v5_allShared_SUMMIT_10col_narrowPeaks/",intern=TRUE)  

for(i in 1:length(filename)){
  
  B73file <- paste(filename[i]) ## 
  peakfile_B73 <- toGRanges(B73file, format="narrowPeak")
  peakAnno_B73 <- annotatePeak(peakfile_B73, tssRegion=c(-10000, 1),TxDb=txdb.B73v5_gff)
  TF_NP_B73 <- as.data.frame(peakAnno_B73)
  write.table(as.data.frame(TF_NP_B73), file=paste("/Users/marygalli/Desktop/PosV_2-25-25_allSharedPeaks_analysis/B73v5_shared_allPeaks_SUMMIT_Chipseeker_files/",paste(filename[i],"Chipseeker_10kb1bp.txt",sep="_"), sep="/"), sep = "\t", row.names = FALSE)
}
```

## Step 2. Identify shared Mo17 target genes and determine distance of peak to TSS

#### Mo17 peak summit is in Mo17_5col.bed files; note that other coords are B73  

`head ARF16_B73v5_allShared_Mo17peakinfo_5col.bed`  
chr6	179135985	179136186	6:166771515	200  
chr2	168390482	168390683	2:169198039	188  
chr3	196929673	196929874	3:196911817	196  
chr5	169357410	169357611	5:168232278	198  
chr1	104943842	104943997	1:105644298	155  
chr3	222458012	222458213	3:225808429	199  

#### create three col .bed file using  column 4 of 5col.bed file and run Chipseeker bed with it
#### use column 4 of 5col.bed (original Mo17 coordinate peak summit) to create 3 column bed file of Mo17 peak summit file of B73-Mo17 shared peaks (with actual Mo17 coordinates).

```
cut -f4 ARF16_B73v5_allShared_Mo17peakinfo_5col.bed | sed -e 's/:/\t/g' | awk '{OFS = "\t"; $1="chr"$1;$3=$2+1; print $1,$2,$3}' > ARF16_B73v5_allShared_MatchedMo17peakSummit_3col.bed
```

#### Run Chipseeker with this file and B73v5 gene models lifted to Mo17 coordinates

#### Run ChIPseeker with Mo17 liftoff in R
library(ChIPseeker)  
library(ChIPpeakAnno) #use for loadDb function  
library(GenomicRanges)  
library(GenomicFeatures)  

### db of B73v5 gene annotations lifted to Mo17 with liftoff (Shumate, A. & Salzberg, S.L. Liftoff: accurate mapping of gene annotations. Bioinformatics 37, 1639-1643 (2021))  

```
datadirGFFMo17 <- "/Users/marygalli/Desktop/B73_Mo17_positional_variation_analysis/"  
datafileGFFMo17_liftoff <- paste(datadirGFFMo17, "B73v5genesToMo17coords_liftoff.gff3", sep="/")  
txdb.Mo17_gff3_liftoff <- makeTxDbFromGFF(file=datafileGFFMo17_liftoff, format="auto")  
```

#### loop thru files in folder with all 3_col bed files 
```
filename <- system("ls /B73v5_allShared_narrowPeaks/B73v5_allShared_MatchedMo17peakSummit_3col_beds/",intern=TRUE)

for(i in 1:length(filename)){
  
  Mo17file <- paste(filename[i]) 
  peakfile_Mo17 <- toGRanges(Mo17file, format="BED")
  peakAnno_Mo17 <- annotatePeak(peakfile_Mo17, tssRegion=c(-10000, 1),TxDb=txdb.Mo17_gff3_liftoff) #promoter defined as -10kb and +1bp
  TF_NP_Mo17 <- as.data.frame(peakAnno_Mo17)
  write.table(as.data.frame(TF_NP_Mo17), file=paste("/Mo17_Matched_B73v5lifted_allShared_Chipseeker_files/",paste("Chipseeker",filename[i],sep="_"), sep="/"), sep = "\t", row.names = FALSE)
}
```

#### ChIPseeker puts "" around certain fields; remove with sed loop  
```
sed -i -e 's/"//g' Chipseeker_ARF16_B73v5_allShared_MatchedMo17peakSummit_3col.bed
```

#### ARF16 Chipseeker output example:  
seqnames	start	end	width	strand	annotation	geneChr	geneStart	geneEnd	geneLength	geneStrand	geneId	transcriptId	distanceToTSS  
chr6	166771516	166771516	1	*	Promoter (5-6kb)	6	166762728	166766442	3715	2	Zm00001eb297400	Zm00001eb297400_T001	-5074  
chr2	169198040	169198040	1	*	5' UTR	2	169195908	169198156	2249	2	Zm00001eb095710	Zm00001eb095710_T001	116  
chr3	196911818	196911818	1	*	Promoter (3-4kb)	3	196915251	196919401	4151	1	Zm00001eb151670	Zm00001eb151670_T001	-3433  
chr5	168232279	168232279	1	*	Distal Intergenic	5	168187073	168191586	4514	2	Zm00001eb242060	Zm00001eb242060_T001	-40693  

#### add Mo17 chr:summit column info so can track peaks (add as last column)  
```
awk -F $'\t' '{OFS = "\t"; $15=$1":"$2; print $0}' Chipseeker_ARF16_B73v5_allShared_MatchedMo17peakSummit_3col.bed > Chipseeker_ARF16_B73v5_allShared_Mo17peakCoordsLiftoffmodels.txt
```

#### final output is list of allshared peaks with Mo17 coordinates, closest lifted B73v5 geneID and distance to TSS

#### ARF16 example of output:  
seqnames	start	end	width	strand	annotation	geneChr	geneStart	geneEnd	geneLength	geneStrand	geneId	transcriptId	distanceToTSS	seqnames:start  
chr6	166771516	166771516	1	*	Promoter (5-6kb)	6	166762728	166766442	3715	2	Zm00001eb297400	Zm00001eb297400_T001	-5074	chr6:166771516  
chr2	169198040	169198040	1	*	5' UTR	2	169195908	169198156	2249	2	Zm00001eb095710	Zm00001eb095710_T001	116	chr2:169198040  
chr3	196911818	196911818	1	*	Promoter (3-4kb)	3	196915251	196919401	4151	1	Zm00001eb151670	Zm00001eb151670_T001	-3433	chr3:196911818  
chr5	168232279	168232279	1	*	Distal Intergenic	5	168187073	168191586	4514	2	Zm00001eb242060	Zm00001eb242060_T001	-40693	chr5:168232279  

## Step 3. Combine B73 and Mo17-lifted target genes and calculate relative distances from TSS

#### paste B73 and Mo17_lifted Chipseeker files together  
```
paste -d'\t' Chipseeker_ARF16_B73v5shared_top20_summit_B73v5peakCoords.txt Chipseeker_ARF16_Mo17_B73v5shared_top20_summitMo17peakCoordsLiftoffmodels.txt > B73v5_and_Mo17_combined_Chipseeker_files/ChIPseeker_ARF16_B73andMo17_B73v5shared_top20percent_compareDistToTSS.txt
```

#### subtract distance between B73 and Mo17, select geneID matches, select peaks with Promoter and 5'UTR annotation, and then record peaks with relative distance to TSS greater than or less than 500 bp  
#### greater than 500bp  
```
awk -F $'\t' '{OFS = "\t"; diff=$19-$34; print $0,diff}' ARF16_ChIPseeker_B73andMo17_B73v5shared_all_SUMMIT_compareDistToTSS.txt | awk -F $'\t' '{OFS = "\t"; if ($17 == $32) print $0,"geneMatch",$17; else print $0,"noMatch"}' - | grep "geneMatch" - | grep -E "5' UTR|Promoter" - | awk -F $'\t' '{OFS = "\t"; if ($36 > 500 || $36 < -500) print $0}' - > ARF16_all_shared_summit_posV_greaterThan500.txt
```
#### less than 500bp  
```
awk -F $'\t' '{OFS = "\t"; diff=$19-$34; print $0,diff}' ARF16_ChIPseeker_B73andMo17_B73v5shared_all_SUMMIT_compareDistToTSS.txt | awk -F $'\t' '{OFS = "\t"; if ($17 == $32) print $0,"geneMatch",$17; else print $0,"noMatch"}' - | grep "geneMatch" - | grep -E "5' UTR|Promoter" - | awk -F $'\t' '{OFS = "\t"; if ($36 <= 500 && $36 >= -500) print $0}' - > ARF16_all_shared_summit_posV_lessThan500.txt
```
#### example of output file; last column contains gene ID to use for differential gene expression analysis  
`head ARF16_all_shared_summit_posV_greaterThan500.txt`  
chr7	2363846	2363847	2	*	648	1235.6	999	999	100	Promoter (5-6kb)	7	2369632	2371252	1621	1	Zm00001eb298790	Zm00001eb298790_T001	-5785	chr7:2363846	chr7	2332170	2332170	1	*	Promoter (3-4kb)	7	2335782	2337449	1668	1	Zm00001eb298790	Zm00001eb298790_T001	-3612	chr7:2332170	-2173	geneMatch	Zm00001eb298790  
chr1	38798171	38798172	2	*	452	794	230	226.93	100	Promoter (3-4kb)	1	38792743	38794651	1909	2	Zm00001eb011940	Zm00001eb011940_T001	-3520	chr1:38798171	chr1	40199810	40199810	1	*	Promoter (2-3kb)	1	40195736	40197653	1918	2	Zm00001eb011940	Zm00001eb011940_T001	-2157	chr1:40199810	-1363	geneMatch	Zm00001eb011940  
chr1	301597492	301597493	2	*	464	821	225.96	222.91	100	Promoter (5-6kb)	1	301603401	301605367	1967	1	Zm00001eb063070	Zm00001eb063070_T001	-5908	chr1:301597492	chr1	299397156	299397156	1	*	Distal Intergenic	1	299418232	299420205	1974	1	Zm00001eb063070	Zm00001eb063070_T001	-21076	chr1:299397156	15168	geneMatch	Zm00001eb063070  
chr1	255205937	255205938	2	*	444	776	214.23	211.3	100	Distal Intergenic	1	255223470	255228015	4546	1	Zm00001eb050010	Zm00001eb050010_T001	-17532	chr1:255205937	chr1	253193966	253193966	1	*	Promoter (6-7kb)	1	253200116	253204549	4434	1	Zm00001eb050010	Zm00001eb050010_T001	-6150	chr1:253193966	-11382	geneMatch	Zm00001eb050010  
chr10	11527261	11527262	2	*	425	733	206.39	203.51	100	Distal Intergenic	10	11508098	11513323	5226	2	Zm00001eb408430	Zm00001eb408430_T001	-13938	chr10:11527261	chr10	11713328	11713328	1	*	Promoter (2-3kb)	10	11708590	11711240	2651	2	Zm00001eb408430	Zm00001eb408430_T001	-2088	chr10:11713328	-11850	geneMatch	Zm00001eb408430  


## Step 4. Use B73 and Mo17 DEG lists to determine enrichment for positional varaition effect.
> Considering only the DAP-seq TF targets that were positionally variant between inbreds, two-by-two continency tables (a gene is DEG or not DEG vs. a target gene is associated with PosV greater than 500-bp or less than 500-bp) were created for each of the 65 TFs and seven tissue types, which were used in one-sided Fisher’s exact tests by the fisher.test function in R. The reported P values were then adjusted for multiple hypothesis testing by the Benjamini–Hochberg procedure.
> 




