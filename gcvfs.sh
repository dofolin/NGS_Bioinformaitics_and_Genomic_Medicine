#!/bin/bash
#SBATCH -A ACD113027       # Account name/project number
#SBATCH -J variant_calling      # Job name
#SBATCH -p ngscourse           # Partition Name
#SBATCH -c 2               # core preserved
#SBATCH --mem=13G           # memory used
#SBATCH -o out_vc_pipelineB.log          # Path to the standard output file 
#SBATCH -e err_vc_pipelineB.log          # Path to the standard error output file
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=FAIL,END

# Please enter the sample name (e.g., SRR13076392)
sample=SRR13147308
# Please enter your username
user=dofolin77

## Change the path below if your HW1 folder isn't like this
HW1dir=/work/dofolin77/final


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
## Path of forward and backward reads
sampleR1=${HW1dir}/${sample}_1.fastq
sampleR2=${HW1dir}/${sample}_2.fastq

mkdir -p ./pipelineB
cd pipelineB

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started"

# Create the environment for alignment, mark duplicate and variant calling
module load biology/SAMTOOLS/1.18
module load biology/BWA/0.7.17
module load biology/GATK/4.2.0.0
PICARD=/work/opt/ohpc/Taiwania3/pkg/biology/Picard/picard_v2.27.4/share/picard-2.27.4-0/picard.jar

# Paths for BQSR
known_indel=/work/r12455009/GATK_BQSR/Homo_sapiens_assembly38.known_indels.vcf.gz
known_mills=/work/r12455009/GATK_BQSR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_dbsnp=/work/r12455009/GATK_BQSR/dbsnp_146.hg38.vcf.gz

gatk CombineGVCFs \
	   -R ${ref} \
	   -v SRR13147325.g.vcf \
	   -v SRR13147319.g.vcf \
	   -v SRR13147305.g.vcf \
	   -v SRR13147318.g.vcf \
	   -v SRR13147312.g.vcf \
	   -v SRR13147316.g.vcf \
	   -v SRR13147321.g.vcf \
	   -v SRR13147308.g.vcf \
	   -O adjusted.g.vcf

gatk GenotypeGVCFs \
	   -R ${ref} \
	   -v adjusted.g.vcf \
	   -O adjusted.vcf.gz
##############################
# Mapping reads with BWA-MEM #
##############################
#bwa mem -M -R "@RG\tID:GP_${sample}\tSM:SM_${sample}\tPL:ILLUMINA" \
#    -t 40 -K 1000000 ${ref} ${sampleR1} ${sampleR2} > ${sample}.sam

echo "Mapping Reads: Finished"

######################
# Sorting & Indexing #
######################
#samtools view -@ 2 -S -b ${sample}.sam > ${sample}.bam
#samtools sort -@ 2 ${sample}.bam -o ${sample}.sorted.bam
#samtools index -@ 20 ${sample}.sorted.bam

#rm ${sample}.sam
echo "Sorting and Indexing: Finished"

###################
# Mark duplicates #
###################
#java -jar ${PICARD} MarkDuplicates \
#      -I ${sample}.sorted.bam \
#      -O ${sample}.sorted.markdup.bam \
#      -M ${sample}_markdup_metrics.txt \
#      --CREATE_INDEX true

echo "Mark duplicates: Finished"

######################
# Base recalibration #
######################
#gatk BaseRecalibrator \
#   -I ${sample}.sorted.markdup.bam \
#   -R ${ref} \
#   --known-sites ${known_indel} \
#   --known-sites ${known_mills} \
#   --known-sites ${known_dbsnp} \
#   -O ${sample}.recal_data.table

#gatk ApplyBQSR \
#   -R ${ref} \
#   -I ${sample}.sorted.markdup.bam \
#   --bqsr-recal-file ${sample}.recal_data.table \
#   -O ${sample}.recaled.bam

echo "Base recalibration: Finished"

############################################
# Calling variants by GATK HaplotypeCaller #
############################################
#gatk HaplotypeCaller \
#    -R ${ref} \
#    -I ${sample}.recaled.bam \
#    -O ${sample}.g.vcf \
#    -ERC GVCF

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"
