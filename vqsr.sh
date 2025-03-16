#!/bin/bash
#SBATCH -A ACD113027       # Account name/project number
#SBATCH -J gvcf      # Job name
#SBATCH -p ngscourse           # Partition Name
#SBATCH -c 2               # core preserved
#SBATCH --mem=13G           # memory used
#SBATCH -o out_vqsr.log          # Path to the standard output file 
#SBATCH -e err_vqsr.log          # Path to the standard error output file
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=FAIL,END

# Please enter the sample name (e.g., SRR13076392)
sample=pCR
# Please enter your username
user=dofolin77

## Change the path below if your HW1 folder isn't like this
HW1dir=/work/dofolin77/HW3/pipelineB


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
## Path of forward and backward reads
#sampleR1=${HW1dir}/${sample}_1.fastq
#sampleR2=${HW1dir}/${sample}_2.fastq

mkdir -p ./pipelineB
cd pipelineB

#echo "$(date '+%Y-%m-%d %H:%M:%S') Job started"

# Create the environment for alignment, mark duplicate and variant calling
module load old-module
module load biology/Python/3.12.2
module load biology/SAMTOOLS/1.18
module load biology/BWA/0.7.17
module load biology/GATK/4.2.0.0
PICARD=/work/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.27.4/share/picard-2.27.4-0/picard.jar

# Paths for BQSR
#known_indel=/work/r12455009/GATK_BQSR/Homo_sapiens_assembly38.known_indels.vcf.gz
#known_mills=/work/r12455009/GATK_BQSR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#known_dbsnp=/work/r12455009/GATK_BQSR/dbsnp_146.hg38.vcf.gz


##############################
# Mapping reads with BWA-MEM #
##############################
#bwa mem -M -R "@RG\tID:GP_${sample}\tSM:SM_${sample}\tPL:ILLUMINA" \
#    -t 40 -K 1000000 ${ref} ${sampleR1} ${sampleR2} > ${sample}.sam

#echo "Mapping Reads: Finished"

######################
# Sorting & Indexing #
######################
#samtools view -@ 2 -S -b ${sample}.sam > ${sample}.bam
#samtools sort -@ 2 ${sample}.bam -o ${sample}.sorted.bam
#samtools index -@ 20 ${sample}.sorted.bam

#rm ${sample}.sam
#echo "Sorting and Indexing: Finished"

###################
# Mark duplicates #
###################
#java -jar ${PICARD} MarkDuplicates \
#      -I ${sample}.sorted.bam \
#      -O ${sample}.sorted.markdup.bam \
#      -M ${sample}_markdup_metrics.txt \
#      --CREATE_INDEX true

#echo "Mark duplicates: Finished"

######################
# Base recalibration #
######################

gatk VariantRecalibrator \
   -R ${ref} \
   -V ${sample}.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --max-gaussians 4 \
   -O ${sample}.SNP.recal \
   --tranches-file ${sample}.SNP.tranches \
   --rscript-file ${sample}.SNP.plots.R

gatk ApplyVQSR \
   -R ${ref} \
   -V ${sample}.vcf.gz \
   --ts-filter-level 99.0 \
   --tranches-file ${sample}.SNP.tranches \
   --recal-file ${sample}.SNP.recal \
   --mode SNP \
   -O ${sample}.SNP.vcf.gz

gatk VariantRecalibrator \
   -R ${ref} \
   -V ${sample}.SNP.vcf.gz \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an FS -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   --max-gaussians 4 \
   -O ${sample}.INDEL.recal \
   --tranches-file ${sample}.INDEL.tranches \
   --rscript-file ${sample}.INDEL.plots.R

gatk ApplyVQSR \
   -R ${ref} \
   -V ${sample}.SNP.vcf.gz \
   --ts-filter-level 99.0 \
   --tranches-file ${sample}.INDEL.tranches \
   --recal-file ${sample}.INDEL.recal \
   --mode INDEL \
   -O ${sample}.recalibrated.vcf.gz


#echo "Base recalibration: Finished"

############################################
# Calling variants by GATK HaplotypeCaller #
############################################
#gatk HaplotypeCaller \
#    -R ${ref} \
#    -I ${sample}.recaled.bam \
#    -O ${sample}.g.vcf \
#    -ERC GVCF

#echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"
