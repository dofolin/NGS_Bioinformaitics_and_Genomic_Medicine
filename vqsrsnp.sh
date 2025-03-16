#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J vqsrsnp        # Job name
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
      -V final_variants.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
	     -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	        -mode SNP \
		   -O snp.recal \
      --tranches-file snp.tranches \
       --rscript-file snp.plots.R
gatk ApplyVQSR \
	    -R Homo_sapiens_assembly38.fasta \
	        -V final_variants.vcf \
		    --ts-filter-level 99.0 \
		        --tranches-file snp.tranches \
			    --recal-file snp.recal \
			        --mode SNP \
				    -O snp.vcf
