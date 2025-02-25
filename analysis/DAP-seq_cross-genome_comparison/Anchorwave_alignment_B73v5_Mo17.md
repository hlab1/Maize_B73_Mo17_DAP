# Whole genome alignment with anchorwave (Song et al., 2022)

## using anchorwave to create whole genome alignement of B73v5 (ref) and Mo17_CAU1.0 (query): this pipeline converts B73v5 coordinates to Mo17 coordinates.

B73v5 genome file:  
Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa  
Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_nucOnly_mg.gff3  

Mo17 genome file:  
Zm-Mo17-REFERENCE-CAU-1.0.id_chr_nuc_mg.fa  
#both use chr1, chr2, etc as chr names (anchorwave will look for matching chr IDs to pair alignments)  

```
anchorwave gff2seq -i Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_nucOnly_mg.gff3 -r Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa -o B73v5_cds.fa
```  

```
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-Mo17-REFERENCE-CAU-1.0.id_chr_nuc_mg.fa B73v5_cds.fa > B73v5_cds.sam
```  

```
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa B73v5_cds.fa > B73v5_ref.sam
```  

##### create script file <genoAli_1thread.sh> with script below
```
#!/bin/bash

#SBATCH --partition=main            # Partition (job queue)

#SBATCH --requeue                   # Return job to the queue if preempted

#SBATCH --job-name=genomeAlign         # Assign a short name to your job

#SBATCH --nodes=1                   # Number of nodes you require

#SBATCH --ntasks=1                 # Total # of tasks across all nodes

#SBATCH --cpus-per-task=1           # Cores per task (>1 if multithread tasks)

#SBATCH --mem=120000                # Real memory (RAM) required (MB)

#SBATCH --time=72:00:00             # Total run time limit (HH:MM:SS)

#SBATCH --output=slurm.%N.%j.out    # STDOUT output file

#SBATCH --error=slurm.%N.%j.err     # STDERR output file (optional)

cd anchorwave_analysis_B73v5/

#module purge
module load intel/19.0.3
module use /projects/community/modulefiles
module load gcc/9.2.0-gc563
module load cmake

srun anchorwave genoAli -i Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_nucOnly_mg.gff3 -as B73v5_cds.fa -r Zm-B73-REFERENCE-NAM-5.0.id_chrs_mg.fa -a B73v5_cds.sam -ar B73v5_ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.id_chr_nuc_mg.fa -n anchors -o anchorwave_Mo17toB73v5.maf -f anchorwave_B73v5.f.maf -IV
```

##### use last to create chain file that is needed to run Crossmap and for viewing synteny alignment in JBrowse2  
> https://gitlab.com/mcfrith/last

```
maf-convert chain anchorwave_Mo17toB73v5.maf > anchorwave_Mo17toB73v5.chain
```  
###### chain file check  
`CrossMap.py viewchain anchorwave_Mo17toB73v5.chain > chainCheck_1.tab`  
`head chainCheck_1.tab`  
chr1	0	1	+	chr1	0	1	+  
chr1	1	5	+	chr1	40421	40425	+  
chr1	16618	16621	+	chr1	40425	40428	+  
chr1	16623	18723	+	chr1	40428	42528	+  
chr1	18723	21690	+	chr1	42529	45496	+  
chr1	21690	23929	+	chr1	54627	56866	+  
chr1	23929	25155	+	chr1	65849	67075	+  
chr1	25156	25366	+	chr1	67075	67285	+  
chr1	25366	28249	+	chr1	67286	70169	+  
chr1	28249	30176	+	chr1	70170	72097	+  

##### Name of chain file for converting B73v5 to Mo17 coordinates (note this file converts B73v5 to Mo17 coordinates; file name is somewhat misleading)  
`anchorwave_Mo17toB73v5.chain`  

#### To convert Mo17 coordinates to B73v5 coordinates (i.e. swap the reference and query sequences -- no need to re-run AnchorWave)

#### Run the script anchorwave-map-swap (part of anchorwave scripts)

`cd anchorwave_analysis_B73v5`  
`module use /projects/community/modulefiles`  
`module load gcc/9.2.0-gc563`  
`module load cmake`  
`module load python/2.7.12`  
```
cat anchorwave_Mo17toB73v5.maf | python2 anchorwave/scripts/anchorwave-maf-swap.py > anchorwave_Mo17toB73v5_swap.maf
```  
#### run maf-convert on swapped file with last

```
maf-convert chain anchorwave_Mo17toB73v5_swap.maf > anchorwave_Mo17toB73v5_swap.chain
```

#### .chain files can be viewed directly in JBrowse2 as synteny files or can be used in CrossMap.py for coordinate conversion

##### use `anchorwave_Mo17toB73v5_swap.chain` to convert Mo17read files to B73 coordinates (note the file name is misleading).






****************************************************************************************************************************
