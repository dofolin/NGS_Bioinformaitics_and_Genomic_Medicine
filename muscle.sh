#!/usr/bin/bash
#SBATCH -A MST112120     
#SBATCH -J muscle
#SBATCH -p ct56           
#SBATCH -c 32       
#SBATCH --mem=64g        
#SBATCH -o /home/dofolin77/blast/slurm_io/muscle.std.log
#SBATCH -e /home/dofolin77/blast/slurm_io/muscle.err.log

cd /home/dofolin77/blast/
module load pkg/Anaconda3
conda activate /home/dofolin77/blast/conda

which conda
conda list

date
muscle -in data/CYP2D6.haplotypes_summer.modified.fasta -out cyp2d6.summer.msa.afa
muscle -maketree -in cyp2d6.summer.msa.afa -out cyp2d6.summer.msa.phy
date
