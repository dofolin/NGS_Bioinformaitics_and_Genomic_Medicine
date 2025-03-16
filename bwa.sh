#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J BWA         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_BWA.log          # Path to the standard output file 
#SBATCH -e err_BWA.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
fasta=Homo_sapiens_assembly38.fasta
fastq_1=SRR13147308_1.fastq
fastq_2=SRR13147308_2.fastq
user=dofolin77
group=S
sample=7S
platform=illumina
nt=2
sample_name=SRR13147308

cd /work/${user}/final

## Set up the environment for running fastqc
module load old-module
module load biology/BWA
module load biology/SAMTOOLS/1.18

## Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq ${R2}.fastq -o /work/${user}/final
bwa mem -M \
	-R "@RG\tID:${group}\tSM:${sample}\tPL:${platform}" \
	-t ${nt} -K 1000000 ${fasta} ${fastq_1} ${fastq_2} | \
	samtools sort -o ${sample_name}.sorted.bam -

