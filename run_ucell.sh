#!/bin/bash
#SBATCH --job-name=UCell_Score
#SBATCH --output=slurm_logs/ucell_%j.out
#SBATCH --error=slurm_logs/ucell_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8       # Anzahl der Kerne für BiocParallel
#SBATCH --mem=64G               # Speicher in GB
#SBATCH --time=04:00:00         # Maximale Laufzeit

# Conda/Mamba Umgebung aktivieren
source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Starte R-Skript: UCell.R"
echo "Verwende Umgebung: $CONDA_DEFAULT_ENV"

# R-Skript ausführen
$(which Rscript) src/UCell.R

echo "Job UCell beendet."
