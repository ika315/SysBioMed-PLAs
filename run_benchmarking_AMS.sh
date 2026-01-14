#!/bin/bash
#SBATCH --job-name=BenchmarkingASM_Bcells
#SBATCH --output=slurm_logs/benchmark_%j.out
#SBATCH --error=slurm_logs/benchmark_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10      # Hochgesetzt für die Scoring-Methoden
#SBATCH --mem=80G               # Angemessener RAM für große Seurat-Objekte
#SBATCH --time=06:00:00         # Maximale Laufzeit

source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

echo "Starte R-Skript: Benchmarking_AddModuleScore.R"
echo "Verwende Umgebung: $CONDA_DEFAULT_ENV"

# Ausführen des R-Skripts mit dem korrekten Rscript-Pfad aus dem Conda-Environment
$(which Rscript) src/Benchmarking_AddModuleScore.R

echo "Job Benchmakring beendet."
