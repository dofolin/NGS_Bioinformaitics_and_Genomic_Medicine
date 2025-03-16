#!/usr/bin/bash
#SBATCH -A MST112120
#SBATCH -J blast_job
#SBATCH -p ct56
#SBATCH -c 32
#SBATCH --mem=64g
#SBATCH -o /home/dofolin77/blast/slurm_io/blast_job.std.log
#SBATCH -e /home/dofolin77/blast/slurm_io/blast_job.err.log

cd /home/dofolin77/blast/
module load pkg/Anaconda3 
conda activate ./conda

which conda
conda list

date
{cmd}
date
