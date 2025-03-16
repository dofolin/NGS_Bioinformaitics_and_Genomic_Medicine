#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J fastqc         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_fastqc.log          # Path to the standard output file 
#SBATCH -e err_fastqc.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
R1=SRR13147308_1
R2=SRR13147308_2
user=dofolin77

cd /work/${user}/final

## Set up the environment for running fastqc
module load biology/FastQC

## Analyzing your sample's sequence QC by fastqc
fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
