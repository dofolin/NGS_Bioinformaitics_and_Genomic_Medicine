#!/usr/bin/sh
#SBATCH -A ACD113027        # Account name/project number
#SBATCH -J sradownload         # Job name
#SBATCH -p ngscourse           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 2               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=13g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out_sra.log          # Path to the standard output file 
#SBATCH -e err_sra.log 
#SBATCH --mail-user=dofolin@gmail.com
#SBATCH --mail-type=END

# Please enter the R1 & R2 file name and your username
#R1=SRR13076392_S14_L002_R1_001
#R2=SRR13076392_S14_L002_R2_001
user=dofolin77

cd /work/${user}/new

## Set up the environment for running fastqc
module load old-module
module load biology/SRAToolkit/2.11.1

##
SRR_LIST=(
    SRR13310284
    SRR13310285
    SRR13310286
    SRR13310291
    SRR13310326
    SRR13310331
    SRR13310342
    SRR13310354
    SRR13310356
    SRR13310358
    SRR13310368
    SRR13310370
    SRR13310388
    SRR13310390
    SRR13310434
    SRR13310448
    SRR13310462
    SRR13310482
    SRR13310293
    SRR13310328
    SRR13310355
    SRR13310387
    SRR13310423
    SRR13310436
    SRR13310437
    SRR13310463
		)

### Analyzing your sample's sequence QC by fastqc
#fastqc ${R1}.fastq.gz ${R2}.fastq.gz -o /home/${user}/HW1
for SRR in "${SRR_LIST[@]}"; do
    fasterq-dump --split-files ${SRR}

    gzip ${SRR}_1.fastq
    gzip ${SRR}_2.fastq 
done
