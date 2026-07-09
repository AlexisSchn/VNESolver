#!/bin/bash
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=0-2
#SBATCH --exclusive
#SBATCH --job-name=VNE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --partition=COMPUTE-SHORT

PROJ_DIR="$HOME/vne/VNESolver/experiments"

sizes=(uni three five)
SIZE="${sizes[$SLURM_ARRAY_TASK_ID]}"

ALGO="FLOWFORMULATION"

# To read the note - if no note is given, default to "No note provided."
XP_NOTE="${1:-No note provided.}"

# create folder for results
TIMESTAMP=$(date +"%Y-%m-%d")
XP_DIR="$PROJ_DIR/results/xp_${TIMESTAMP}_${SLURM_ARRAY_JOB_ID}_$ALGO"
mkdir -p "$XP_DIR"

# Note into the folder. 
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "=========================================" > "$XP_DIR/note.txt"
    echo "Experiment Note:" >> "$XP_DIR/note.txt"
    echo "$XP_NOTE" >> "$XP_DIR/note.txt"
    echo "Timestamp: $(date)" >> "$XP_DIR/note.txt"
    echo "Job ID: $SLURM_ARRAY_JOB_ID" >> "$XP_DIR/note.txt"
    echo "=========================================" >> "$XP_DIR/note.txt"
fi

# tiny sleep to ensure the folder and note are completely settled
sleep 2

srun julia --project=$PROJ_DIR scripts/flow_formulation/test_flow.jl "$SIZE" "$XP_DIR" \
    > "$XP_DIR/output-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt" \
    2> "$XP_DIR/error-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"