#!/usr/bin/bash
#SBATCH -A MST112120     
#SBATCH -J bed
#SBATCH -p ct56           
#SBATCH -c 10          
#SBATCH --mem=64g        
#SBATCH -o /home/dofolin77/blast/slurm_io/bed.std.log
#SBATCH -e /home/dofolin77/blast/slurm_io/bed.err.log

cd /home/dofolin77/blast/
module load pkg/Anaconda3
conda activate /home/dofolin77/blast/conda

which conda
conda list

date
makeblastdb -dbtype nucl -parse_seqids -in data/HG002.1.fa -out dataBase/HG002.1
date
