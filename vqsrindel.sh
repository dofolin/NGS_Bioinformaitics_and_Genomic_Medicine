#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J vqsrindel        # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_vqsr.log          # Path to the standard output file 
#SBATCH -e err_vqsr.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

#
user=dofolin77

cd /work/${user}/final

##
module load old-module
module load biology/Python/3.12.2
module load biology/GATK/4.2.3.0
###
gatk VariantRecalibrator \
   -R Homo_sapiens_assembly38.fasta \
      -V snp.vcf \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
       -an QD -an FS -an ReadPosRankSum -an MQRankSum \
	        -mode INDEL \
		   -O indel.recal \
      --tranches-file indel.tranches \
       --rscript-file indel.plots.R
gatk ApplyVQSR \
	    -R Homo_sapiens_assembly38.fasta \
	        -V snp.vcf \
		    --ts-filter-level 99.0 \
		        --tranches-file indel.tranches \
:			    --recal-file indel.recal \
			        --mode INDEL \
				    -O recalibrated.vcf
