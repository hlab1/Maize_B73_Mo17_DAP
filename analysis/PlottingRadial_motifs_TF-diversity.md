## Plotting motifs for TF diversity panel TFs (Fig. 2a)  

#### use meme to extract motif from top 1000 peaks in each dataset  
#### sort data by adjusted p-value  
```
sort -k9,9nr -k7,7nr TF_B73NAMv5_Q30_default_finalBl.GEM_events.narrowPeak | head -1000 > meme_analysis/TF_GEM_events.top600.narrowPeak
```
#### get fasta sequence for top 1000 peaks  
```
bedtools getfasta -fi /scratch/mjg328/Genomes/B73v5_NAM_1-2022/B73v5_NAM_1-2022/Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa -bed TF_B73NAMv5_GEM_events.top1000.narrowPeak -fo TF_B73NAMv5_GEM_events.top1000.fasta
```

#### Run meme analysis (can use webserver or cluster; on amarel cluster need 'module load openmpi/2.1.1')  
```
meme-chip -oc . -time 240 -ccut 100 -fdesc description -dna -order 2 -minw 6 -maxw 10 -db db/motif_databases/ARABD/ArabidopsisDAPv1.meme -meme-mod zoops -meme-nmotifs 3 -meme-searchsize 100000 -streme-pvt 0.05 -streme-align center -streme-totallength 4000000 -centrimo-score 5.0 -centrimo-ethresh 10.0 TF_B73v5_GEM_events_top1000.fasta
```
#### access count matrix data from meme and save as TF_counts.txt
#### convert TF_counts.txt to the following format:
> position	A	C	G	T  
pos_1	380	217	243	160  
pos_2	251	352	272	125  
pos_3	92	722	161	25  
pos_4	418	539	43	0  
pos_5	9	0	991	0  
pos_6	991	0	9	0  
pos_7	0	1000	0	0  
pos_8	980	20	0	0  
pos_9	346	154	259	241  
pos_10	162	281	511	46
> 

#### motifStack for maize diversity panel TFs
`library(motifStack)`  
`library(reshape2)`

#### convert meme PWM to proper format  
`setwd("/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/")`  
`filename <- system("ls /Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/",intern=TRUE)`  

```
for(i in 1:length(filename)){  
  
  file <- read.table(filename[i],header=TRUE)  
  #mean <- mean(file$V4)  
  melt <- melt(file, id.vars = c("position"),  
               variable.name = "new_base", value.name = "prop")  
  dcast1 <- dcast(melt, new_base ~ position, value.var = "prop")  
  dcast1 <- dcast1[,2:ncol(dcast1)]  
  rownames(dcast1) <- c("A","C","G","T")  
  write.table(dcast1,file=paste("/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts",paste("reshaped_2",filename[i],sep="."),sep="/"), sep = "\t")  
}  
```

#### copy output files to a new directory called reshaped_motifs  
#### change file names to *.pcm and use with readPCM function  

#### view motifs as stack
```
pcm_TFs<-readPCM(file.path("/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/reshaped_motifs/"),"pcm$")  
motifTFs<-lapply(pcm_TFs,pcm2pfm)  
pdf(file = "/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/diversityTFs_stack.pdf", width = 6, height = 60)  
motifStack(motifTFs, layout="stack", ncex=0.5)  
dev.off()  
```
```
pdf(file = "/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/diversityTFs_stack_1.pdf", width = 6, height = 60)  
motifStack(motifTFs, layout="stack", ncex=0.01)# to mask labels  
dev.off()
```
```
motifTFs_rc<-lapply(pcm_TFs,matrixReverseComplement)
```

`library(RColorBrewer)`  

#### plot logo stack with radial style  
```
pcm_TFs<-readPCM(file.path("/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/reshaped_motifs/"),"pcm$")  
motifTFs<-lapply(pcm_TFs,pcm2pfm)  
pfms <- lapply(pcm_TFs,pcm2pfm)  
color <- brewer.pal(12, "Set3")  
pdf(file = "/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/meme_diversity_panel/meme_motif_counts/diversityTFs_radial_fixed.pdf", width = 10, height = 12)  
motifStack(motifTFs, layout="radialPhylog")  
dev.off()
```
