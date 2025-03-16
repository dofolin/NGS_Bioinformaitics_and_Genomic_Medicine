#!/bin/bash
#SBATCH -A ACD113027       # Account name/project number
#SBATCH -J variant_calling      # Job name
#SBATCH -p ngscourse           # Partition Name
#SBATCH -c 2               # core preserved
#SBATCH --mem=13G           # memory used
#SBATCH -o out_vc_pipelineA.log          # Path to the standard output file 
#SBATCH -e err_vc_pipelineA.log          # Path to the standard error output file
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=FAIL,END

# Please enter the sample name (e.g., SRR13076392)
sample=SRR13076392
# Please enter your username
user=dofolin77

## Change the path below if your HW1 folder isn't like this
HW1dir=/home/dofolin77/HW1


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
## Path of forward and backward reads
sampleR1=${HW1dir}/${sample}*R1_001.fastq.gz
sampleR2=${HW1dir}/${sample}*R2_001.fastq.gz

mkdir -p ./pipelineA
cd pipelineA

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started"

# Create the environment for alignment and variant calling
module load biology/SAMTOOLS/1.18
module load biology/BWA/0.7.17
module load biology/GATK/4.2.0.0
set -euo pipefail

##############################
# Mapping reads with BWA-MEM #
##############################
bwa mem -M -R "@RG\tID:GP_${sample}\tSM:SM_${sample}\tPL:ILLUMINA" \
    -t 40 -K 1000000 ${ref} ${sampleR1} ${sampleR2} > ${sample}.sam

echo "Mapping Reads: Finished"

######################################################
# Preparing for variant calling (sorting & indexing) #
######################################################
samtools view -@ 2 -S -b ${sample}.sam > ${sample}.bam
samtools sort -@ 2 ${sample}.bam -o ${sample}.sorted.bam
samtools index -@ 20 ${sample}.sorted.bam

rm ${sample}.sam
echo "Sorting and Indexing: Finished"

############################################
# Calling variants by GATK HaplotypeCaller #
############################################
gatk HaplotypeCaller \
    -R ${ref} \
    -I ${sample}.sorted.bam \
    -O ${sample}.HC.vcf.gz

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"
