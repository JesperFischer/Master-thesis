#!/bin/bash
#SBATCH --job-name=Hierarchical_power
#SBATCH --cpus-per-task=10
#SBATCH --mem=5000

# Define your R script and its arguments
R_SCRIPT="Make_datasets_slurm.R"
SUBJECTS="0"          # Set your default value for subjects
TRIALS="0"           # Set your default value for trials
EFFECT_SIZE_ALPHA="0"      # Set your default value for effect size
EFFECT_SIZE_BETA="0"      # Set your default value for effect size


# If command line arguments are provided, use them
if [ "$#" -eq 4 ]; then
  SUBJECTS=$1
  TRIALS=$2
  EFFECT_SIZE_ALPHA=$3
  EFFECT_SIZE_BETA=$4
fi

# Run R script with the specified arguments
Rscript $R_SCRIPT $SUBJECTS $TRIALS $EFFECT_SIZE_ALPHA $EFFECT_SIZE_BETA