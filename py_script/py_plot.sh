#!/bin/bash
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=1  
#SBATCH --time=06:00:00    
#SBATCH --error=plot.err
#SBATCH --output=plot.out
#SBATCH --job-name=plot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=1059291645@qq.com
module load anaconda/3.7 
python obs.py