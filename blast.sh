#!/usr/bin/bash
#SBATCH -A MST112120     
#SBATCH -J blast
#SBATCH -p ct56           
#SBATCH -c 32       
#SBATCH --mem=64g        
#SBATCH -o /home/dofolin77/blast/slurm_io/bed.std.log
#SBATCH -e /home/dofolin77/blast/slurm_io/bed.err.log

cd /home/dofolin77/blast/
module load pkg/Anaconda3
conda activate /home/dofolin77/blast/conda

which conda
conda list

date
makeblastdb -dbtype nucl -parse_seqids -in data/CYP2D6.haplotypes.fasta -out dataBase/Haplotypes
blastn -query HG002.1.CYP2D6.fa.fai -db dataBase/Haplotypes -task blastn -evalue 1e-6 -outfmt 0 -out blast2.tsv
date
