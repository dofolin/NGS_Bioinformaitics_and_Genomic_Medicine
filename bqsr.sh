#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J bqsr         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_bqsr.log          # Path to the standard output file 
#SBATCH -e err_bqsr.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
fasta=Homo_sapiens_assembly38.fasta
user=dofolin77
sample_name=SRR13147308
known_indel=1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_mills=Homo_sapiens_assembly38.dbsnp138.vcf
known_dbsnp=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

cd /work/${user}/final

## Set up the environment for running fastqc
module load old-module
module load biology/Python/3.12.2
module load biology/GATK/4.2.3.0
module load biology/SAMTOOLS/1.18

## Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
gatk BaseRecalibrator \
	    -R ${fasta} \
	        -I ${sample_name}.markdup.bam \
		    --known-sites ${known_indel} \
		        --known-sites ${known_mills} \
			    --known-sites ${known_dbsnp} \
			    --known-sites 1000G_omni2.5.hg38.vcf.gz \
			    --known-sites hapmap_3.3.hg38.vcf.gz \
			        -O ${sample_name}.recal_data.table
gatk ApplyBQSR \
	    --bqsr-recal-file ${sample_name}.recal_data.table \
	        -R ${fasta} \
		    -I ${sample_name}.markdup.bam \
		        -O ${sample_name}.recaled.bam
samtools index ${sample_name}.recaled.bam
