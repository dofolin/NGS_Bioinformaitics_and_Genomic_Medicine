#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J sam2bam        # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_sam2bam.log          # Path to the standard output file 
#SBATCH -e err_sam2bam.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# 
user=dofolin77

cd /work/${user}/final

## Set up the environment for running samtool
module load biology/SAMTOOLS/1.18

## Converting your sample's SAM to BAM by samtool
samtools view -bS SRR13147325.sam | samtools sort -o SRR13147325.bam
