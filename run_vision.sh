#!/bin/bash
#SBATCH --job-name=VISION_Bcell
#SBATCH --partition=compms-cpu-big
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/home/students/i.kaciran/SysBioMed-PLAs/slurm_logs/vision_%j.out
#SBATCH --error=/nfs/home/students/i.kaciran/SysBioMed-PLAs/slurm_logs/vision_%j.err

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Running on node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"

Rscript src/Vision.R
