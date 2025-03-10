## Pearson correlation of B73v5 DAP-seq TF datasets 

#### use deeptools2 to create multiBigwigSummary to examine Pearson correlation of TF binding events
> https://deeptools.readthedocs.io/en/latest/
> 
#### copy all TF.bigwig files to single directory

`cd /projects/f_waksman_1/Gallavotti_Lab/marygalli/DAP-seq_data/B73v5_bigwig_all_final/`

#### create shell script for sbatch submission

`multiBigwigSumm_B73v5all.sh`

```
#!/bin/bash
#SBATCH --partition=p_waksman_1     # Partition (job queue)

#SBATCH --requeue                   # Return job to the queue if preempted

#SBATCH --job-name=bigwigsum        # Assign a short name to your job

#SBATCH --nodes=1                   # Number of nodes you require

#SBATCH --ntasks=1                  # Total # of tasks across all nodes

#SBATCH --cpus-per-task=16          # Cores per task (>1 if multithread tasks)

#SBATCH --mem=80000                 # Real memory (RAM) required (MB)

#SBATCH --time=72:00:00             # Total run time limit (HH:MM:SS)

#SBATCH --output=slurm.%N.%j.out    # STDOUT output file

#SBATCH --error=slurm.%N.%j.err     # STDERR output file (optional)

srun multiBigwigSummary bins --bwfiles /projects/f_waksman_1/Gallavotti_Lab/marygalli/DAP-seq_data/B73v5_bigwig_all_final/*bigwig -bs 100 --region chr10 --blackListFileName B73v5_blocklist6_3-1-2022_nooverlaps2.bed -p max -o DAP_B73v5_TFsall_finalbl_bs100.npz  
#use 100bp bin size and use only chr 10 reads  
```
#### use deeptools2 plotCorrelation feature to plot Pearson correlation of TF binding events  
```
plotCorrelation -in DAP_B73v5_TFsall_finalbl_bs100.npz -c pearson -p heatmap -o DAP_B73V5_TFsall_finalbl_bs100_skipZeros.pdf --labels ABI19 ABI34 ABI8 ARR5 BAF1 BHLH106_ZmBIM2 BHLH125_ZmSPT1 BHLH149_MYC2A BHLH168_ZmSPT2 BHLH17 BHLH57_ZmJAM1 BHLH67_BA1 BHLH85_ZmHEC3 BHLH91_MYC2B BHLH99_ZmJAM2 BZIP105 BZIP106 BZIP111 BZIP112 BZIP115 BZIP124 BZIP128_OHP1 BZIP12 BZIP13 BZIP160 BZIP1_O2 BZIP21 BZIP23 BZIP25 BZIP27 BZIP2_DLF1 BZIP34 BZIP39 BZIP40 BZIP41 BZIP48 BZIP52 BZIP54 BZIP58 BZIP61 BZIP62 BIZP68_ZmABF2 BIZP70 BZIP9 BZIP84 BZIP86 BZIP87 BZIP90_LRS1 BZIP91 BZIP92 BZIP94_OHP2 BZIP99 BZR1 CA2P13:CA3P1:CA5P7 CA2P13:CA3P4:CA5P7 CCHH2_ZmIDD15 DOF1_PBF1 DOF24 EREB126 EREB145 EREB149 EREB183_ZmPUCHI EREB201 EREB2 EREB36 EREB46_ZmSHINE EREB49 EREB4 EREB53_BBM EREB87 GATA20_TSH1 GBP10 GBP11 GBP4 GLK17 GLK20 GLK21 GLK44 GLK47 GLK48 GLK53 GLK56 GRF2 GRF3 GRF6 GRF7 GRF9 HB102 HB117_WOX9C HB122_WUS2 HB125_NS1 HB15_OCL4 HB22 HB54 HB67_WUS1 HB68_GT1 HSF10 HSF24 HSF25 HSF7 MADS12_BDE1 MADS43 MADS68 MADS69 MADS72_TU1 MADS73 MYB131 MYB138 MYB40 MYB50 MYB83 MYBR101 MYBR111 MYBR17 MYBR52 MYBR58 NAC115 NAC116 NAC12 NAC20 NAC25 NAC27 NAC28 NAC33 NAC3 NAC42 NAC49 NAC71 NAC74 NAC86 NFYA1:LEC1:CA5P16 NLP14 SBP15_LG1 SBP1_TGA1 SBP2_THS4 SBP30_UB3 SBP32_NOT1 TCP11 TCP15 TCP2 TCP30 TCP4 TCP5 TCP8 THX16 THX20 THX26 WRKY11 WRKY125 WRKY53 WRKY57 WRKY64 WRKY82 YAB2 ZHD15 ZHD1 ZHD21 ZHD4 ZHD5 ZIM2 --colorMap Blues --skipZeros --outFileCorMatrix DAP_B73v5_TFsall_finalbl_bs100_skipZeros.txt
```
#### plot data in R using NMF aheatmap
```
pearson_b <- read.table("/Users/marygalli/Desktop/DAP_B73v5_TFsall_finalbl_bs100_skipZeros.tab", header = TRUE, sep = "\t", row.names = 1)`
```

#### annotate each TF by family type
```
familyTF <- read.table("/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/BigwigSum_B73v5_annCol.txt", header = TRUE) #make list that in same order as in big dataframe, the header will be used as the label  
```

```
Var1 = c("#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#a6cee3","#33a02c","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f","#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef","#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221")  
ann_colors = list(Var1)
```

#### plot
`library(NMF)`
```
pdf(file= "/Users/marygalli/Desktop/DAP_B73_Mo17_manuscript_analysis_files/BigwigSum_B73v5_bs100_chr10_Blues_annC.pdf", width = 13, height = 12)
aheatmap(pearson_b, color = "Blues:100", revC=TRUE, annColors = ann_colors, annCol = familyTF, annRow = familyTF, Colv = FALSE, Rowv = FALSE, fontsize = 10)
dev.off()
```
