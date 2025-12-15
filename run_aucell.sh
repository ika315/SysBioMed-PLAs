#!/bin/bash
#SBATCH --job-name=AUCell_Score
#SBATCH --output=slurm_logs/aucell_%j.out
#SBATCH --error=slurm_logs/aucell_%j.err
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8      # WICHTIG: Anzahl der Kerne für BiocParallel
#SBATCH --mem=64G              # WICHTIG: Speicher in GB (64 GB sind oft gut für große Seurat-Objekte)
#SBATCH --time=04:00:00        # Maximale Laufzeit (4 Stunden)

# Conda/Mamba Umgebung aktivieren
source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Starte R-Skript: AUCell.R"
echo "Verwende Umgebung: $CONDA_DEFAULT_ENV"

# Ausführen des R-Skripts mit dem korrekten Rscript-Pfad aus dem Conda-Environment
$(which Rscript) src/AUCell.R

echo "Job AUCell beendet."
