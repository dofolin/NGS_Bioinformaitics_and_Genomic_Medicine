#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J hisat2        # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_hisat2.log          # Path to the standard output file 
#SBATCH -e err_hisat2.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
sample=SRR1331046
R1=${sample}_1
R2=${sample}_2
user=dofolin77


cd /work/${user}/final

## Set up the environment for running hisat2
module load old-module
module load biology/HISAT2/2.2.1
module load biology/SAMTOOLS/1.18
export PATH=/opt/ohpc/Taiwania3/pkg/biology/Stringtie/stringtie_v2.2.1/bin:$PATH

## Mapping your sample's sequence reads by hisat2
hisat2 -x genome -1 /work/dofolin77/new/${R1}.fastq.gz -2 /work/dofolin77/new/${R2}.fastq.gz -S /work/dofolin77/new/${sample}.sam

samtools view -bS /work/dofolin77/new/${sample}.sam | samtools sort -o /work/dofolin77/new/${sample}.bam

stringtie /work/dofolin77/new/${sample}.bam -o /work/dofolin77/new/${sample}.gtf -p 8 -G Homo_sapiens.GRCh38.112.gtf
