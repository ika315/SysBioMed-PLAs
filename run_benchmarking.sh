#!/bin/bash
#SBATCH --job-name=PLA_Benchmarking
#SBATCH --output=slurm_logs/benchmark_%j.out
#SBATCH --error=slurm_logs/benchmark_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=06:00:00

# Conda Umgebung aktivieren
source /opt/mambaforge/etc/profile.d/conda.sh
conda activate pla

METHOD=$1
SIG_NAME=$2
SIG_FILE=$3
USE_EXT=$4
THRESH=$5

echo "Starte R-Skript: Benchmarking_Master_v2.R"
echo "Methode: $METHOD | Signature: $SIG_NAME | Extension: $USE_EXT"
echo "Job-ID: $SLURM_JOB_ID"

# R-Skript ausführen
$(which Rscript) src/Benchmarking_Master_v2.R "$METHOD" "$SIG_NAME" "$SIG_FILE" "$USE_EXT" "$THRESH"

echo "Job für $METHOD beendet."
