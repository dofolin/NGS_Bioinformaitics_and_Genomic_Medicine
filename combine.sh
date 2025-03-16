#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J markdup         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_markdup.log          # Path to the standard output file 
#SBATCH -e err_markdup.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
#fasta=Homo_sapiens_assembly38.fasta
#fastq_1=SRR13147308_1.fastq
#fastq_2=SRR13147308_2.fastq
user=dofolin77
#group=S
#sample=7S1
#platform=illumina
#nt=2
#output=SRR13147308

cd /work/${user}/final

## Set up the environment for running fastqc
module load old-module
module load biology/Python/3.12.2
module load biology/GATK/4.2.3.0
#module load biology/SAMTOOLS/1.18

## Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
gatk CombineGVCFs \
	   -R Homo_sapiens_assembly38.fasta \
	      -V SRR13147325.g.vcf \
	         -V SRR13147319.g.vcf \
		    -V SRR13147305.g.vcf \
		                  -V SRR13147318.g.vcf \
				           -V SRR13147312.g.vcf \
				                       -V SRR13147316.g.vcf \
						                     -V SRR13147321.g.vcf \
							                      -V SRR13147308.g.vcf \
		       -O combine.g.vcf
