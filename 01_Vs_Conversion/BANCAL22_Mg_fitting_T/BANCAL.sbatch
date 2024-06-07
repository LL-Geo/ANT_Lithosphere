#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --mem=100G
source activate simpeg
python main.py >> ausrem_run.txt