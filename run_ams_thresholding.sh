#!/bin/bash
#SBATCH --job-name=thresholding_add_module_score
#SBATCH --partition=compms-cpu-big
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --time=15:00:00
#SBATCH --output=/nfs/home/students/i.kaciran/SysBioMed-PLAs/slurm_logs/thresholding_ams_%j.out
#SBATCH --error=/nfs/home/students/i.kaciran/SysBioMed-PLAs/slurm_logs/thresholding_ams_%j.err

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Running on node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"

Rscript src/thresholding/tresholding_addmodulescore.R