#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J vep        # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=FAIL,END

# Please enter the sample name (e.g., SRR13076392)
sample=sensitive.recalibrated
# Please enter which pipeline you are going to run (pipelineA or pipelineB)
pipeline=pipelineB


# ------------------------------------ #
# Please don't change the script below #
# ------------------------------------ #
## Set up the environment and path for running VEP
VEP_PATH=/opt/ohpc/Taiwania3/pkg/biology/Ensembl-VEP/ensembl-vep/vep
VEP_CACHE_DIR=/opt/ohpc/Taiwania3/pkg/biology/DATABASE/VEP/Cache
VEP_FASTA=/opt/ohpc/Taiwania3/pkg/biology/reference/Homo_sapiens/GATK/hg38/Homo_sapiens_assembly38.fasta
BCFTOOLS=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools

module load old-module
module load biology/Perl/5.28.1
module load pkg/Anaconda3
export PATH=${PATH}:/opt/ohpc/Taiwania3/pkg/biology/HTSLIB/htslib_v1.13/bin:/opt/ohpc/Taiwania3/pkg/biology/SAMTOOLS/samtools_v1.15.1/bin
set -euo pipefail

## TWB SNV/indel custom annotation file
twb_snv=/work/r12455009/hg38_TWB_official_AF/TWB_official_snv_indel_AF.vcf.gz
custom_snv=${twb_snv},TWB_official_SNV_indel,vcf,exact,0,AF

## Input data path
wkdir=$(realpath ${pipeline})
cd ${wkdir}

## Log file settings
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_vep_${sample}.log

## Redirect standard output and error to the log file
exec > "$logfile" 2>&1

printf "#############################################################################\n"
printf "###                Work started:   $(date +%Y-%m-%d:%H:%M:%S)             ###\n"
printf "#############################################################################\n"

# split multiallelic
${BCFTOOLS} norm -m -any ${sample}.vcf.gz \
    -Oz \
    -o ${sample}.HC_normed.vcf.gz
${BCFTOOLS} index -t -f ${sample}.HC_normed.vcf.gz

echo "Split multiallelic: Finished"

INPUT_VCF=${sample}.HC_normed.vcf.gz
SAMPLE_ID=${sample}.HC.VEP

#############################
# Variant annotation by VEP #
#############################
${VEP_PATH} --cache --offline \
    --cache_version 108 \
    --dir_cache ${VEP_CACHE_DIR} \
    --assembly GRCh38 \
    --fasta ${VEP_FASTA} \
    --fork 4 \
    -i ${INPUT_VCF} \
	--custom ${custom_snv} \
    --check_existing \
    --af_gnomade \
    --af_gnomadg \
    --vcf \
    -o ${SAMPLE_ID}.vcf \
    --force_overwrite

echo "VEP annotation: Finished"

# Generate tsv
echo -e "CHROM\tPOS\tREF\tALT\tDP\t$(${BCFTOOLS} +split-vep -l ${SAMPLE_ID}.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > ${SAMPLE_ID}.tsv
${BCFTOOLS} +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%CSQ\n' -A tab ${SAMPLE_ID}.vcf >> ${SAMPLE_ID}.tsv

awk 'NR==1 || ($1 ~ /^chr[1-9]$|^chr10$/) {
    gsub(/,.*$/, "", $6)
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $9 "\t" $10 "\t" $34 "\t" $44 "\t" $50 "\t" $54
}' ${SAMPLE_ID}.tsv > ${SAMPLE_ID}_filtered.tsv

echo "Format changing & column filtering: Finished"

printf "#############################################################################\n"
printf "###               Work completed: $(date +%Y-%m-%d:%H:%M:%S)              ###\n"
printf "#############################################################################\n"

rm -f ../slurm*.out
