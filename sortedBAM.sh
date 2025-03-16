#!/bin/bash
#SBATCH -A ACD113027       # Account name/project number
#SBATCH -J variant_calling      # Job name
#SBATCH -p ngscourse           # Partition Name
#SBATCH -c 2               # core preserved
#SBATCH --mem=13G           # memory used
#SBATCH -o out_vc_pipelineA.log          # Path to the standard output file 
#SBATCH -e err_vc_pipelineA.log          # Path to the standard error output file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL,END

# Please enter the sample name (e.g., SRR13076392)
sample=sample1
# Please enter your username
user=username

## Change the path below if your HW1 folder isn't like this
HW1dir=/home/${user}/final


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Reference: Homo_sapiens_assembly38.fasta
ref=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
## Path of the existing BAM file
input_bam=${HW1dir}/${sample}.bam
input_vcf=${HW1dir}/${sample}.vcf.gz

mkdir -p ./pipelineA
## Copy BAM and VCF files to pipelineA directory
cp ${input_bam} ./pipelineA
cp ${input_vcf} ./pipelineA

cd pipelineA

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started"

# Create the environment for alignment and variant calling
module load biology/SAMTOOLS/1.18
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
samtools sort -@ 2 ${input_bam} -o ${sample}.sorted.bam
samtools index -@ 20 ${sample}.sorted.bam

echo "Sorting and Indexing: Finished"

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished"