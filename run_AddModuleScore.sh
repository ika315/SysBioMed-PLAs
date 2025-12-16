#!/bin/bash
#SBATCH --job-name=AddModuleScore
#SBATCH --output=slurm_logs/addmodulescore_%j.out
#SBATCH --error=slurm_logs/addmodulescore_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# Conda/Mamba Umgebung aktivieren
source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Starte R-Skript: AddModuleScore.R"
echo "Verwende Umgebung: $CONDA_DEFAULT_ENV"

# R-Skript ausführen
$(which Rscript) src/AddModuleScore.R

echo "Job AddModuleScore beendet."
