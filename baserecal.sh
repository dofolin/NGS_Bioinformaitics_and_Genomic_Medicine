#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J baserecal         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_baserecal.log          # Path to the standard output file 
#SBATCH -e err_baserecal.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
#R1=SRR13147308_1
#R2=SRR13147308_2
user=dofolin77

cd /work/${user}/final

## Set up the environment for running fastqc
module load old-module
module load biology/GATK/4.2.3.0
module load miniconda3

## Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
gatk BaseRecalibrator -I dedup_SRR13147305.bam -R resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --known-sites Homo_sapiens_assembly38.dbsnp138.vcf 	--known-sites 1000G_phase1.snps.high_confidence.hg38.vcf --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf	-O SRR13147305_recal.table
