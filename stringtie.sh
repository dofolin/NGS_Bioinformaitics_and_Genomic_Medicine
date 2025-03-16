#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J stringtie         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_stringtie.log          # Path to the standard output file 
#SBATCH -e err_stringtie.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END


# Please enter the R1 & R2 file name and your username
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Stringtie/stringtie_v2.2.1/bin:$PATH

samples=(
SRR13310285
SRR13310286
SRR13310291
SRR13310326
SRR13310331
SRR13310342
SRR13310354
SRR13310356
SRR13310358
SRR13310368
SRR13310370
SRR13310388
SRR13310390
SRR13310434
SRR13310448
SRR13310462
SRR13310482
SRR13310293
SRR13310328
SRR13310355
SRR13310387
SRR13310423
SRR13310436
SRR13310437
SRR13310463
)

sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

user=dofolin77

cd /work/${user}/final

## Set up the environment for running fastqc
#module load biology/FastQC

## Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
#stringtie SRR13147325.bam -o tra_SRR13147325.gtf -p 8 -G ref_SRR13147325.gtf
stringtie /work/dofolin77/new/${sample}.bam -e -B -o /work/dofolin77/new/${sample}.gtf -p 8 -G merged26.gtf
